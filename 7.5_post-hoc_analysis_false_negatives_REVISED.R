library(data.table)
library(dplyr)
library(reshape2)
library(tidyr)
library(purrr)
library(stringr)
library(Matrix)
library(progress)
library(parallel)
library(caret)

site <- "Zambezia"

###########################################
# Input full training genomic data and metadata
PAIRS_METADATA <- readRDS(paste0("PAIRS_METADATA_", site, ".RDS"))
PAIRS_GENOMIC <- readRDS(paste0("PAIRS_GENOMIC_", site, ".RDS"))

PAIRS_GENOMIC <- as.data.table(PAIRS_GENOMIC)
PAIRS_METADATA <- as.data.table(PAIRS_METADATA)

PAIRS_GENOMIC <- left_join(PAIRS_GENOMIC, PAIRS_METADATA, by = "PairsID") 

# Input features testing data for subsetting
test_data <- readRDS(paste0(site, "_test_data.RDS"))
TEST_META <- test_data %>%
  mutate(PairsID = as.character(PairsID),
         labels = factor(labels, levels = c("NI", "R"))) %>%
  select(PairsID, labels, pair_type = eCOI_pairs) %>%  # Renamed to pair_type for consistency
  as.data.table()

# Subset metadata using test data
# DO NOT subset PAIRS_GENOMIC (it removes alleles)
# Instead, subset only the metadata to get test pairs
PAIRS_METADATA <- PAIRS_METADATA[PAIRS_METADATA$PairsID %in% TEST_META$PairsID,]

###########################################
# Apply dropout function to PAIRS_GENOMIC
# Remove alleles from each timepoint of each PairsID, but only when there's more than 1 per loci

dropout_rate = 0.75

# Function to apply dropout to genomic data
apply_dropout <- function(data, dropout_rate = 0) {
  # Group by PairsID, time_point, and locus
  data_grouped <- data %>%
    group_by(PairsID, time_point, locus) %>%
    mutate(
      allele_count = n(),
      # Only apply dropout when more than 1 allele exists at a locus
      dropout = ifelse(allele_count > 1, 
                       sample(c(TRUE, FALSE), n(), replace = TRUE, 
                              prob = c(dropout_rate, 1 - dropout_rate)), 
                       FALSE)
    ) %>%
    filter(!dropout) %>%
    select(-allele_count, -dropout) %>%
    ungroup()
  
  return(data_grouped)
}

# Apply dropout to PAIRS_GENOMIC
PAIRS_GENOMIC <- apply_dropout(PAIRS_GENOMIC, dropout_rate = dropout_rate)

###########################################
# Input model and optimal thresholds for each pair_type
model <- readRDS(paste0("LogReg_model_", site, ".RDS"))
thresholds <- read.csv(paste0("training_Results_LogReg_", site, ".csv"))
thresholds <- thresholds %>% rename(pair_type = eCOI_pairs)

###########################################
# Calculate features
calculate_features_optimized <- function(sample1, sample2) {
  # 1) Unique alleles & loci
  alleles1 <- unique(sample1$allele)
  alleles2 <- unique(sample2$allele)
  all_alleles <- union(alleles1, alleles2)
  
  # 2) Allele‐level intersection & union via set operations
  inter_cnt   <- length(intersect(alleles1, alleles2))
  union_cnt   <- length(union(alleles1, alleles2))
  jaccard     <- if (union_cnt > 0) inter_cnt/union_cnt else 0
  retention   <- if (length(alleles1) > 0) inter_cnt/length(alleles1) else 0
  allele_gain <- if (union_cnt > 0) length(setdiff(alleles2, alleles1))/union_cnt else 0
  allele_loss <- if (union_cnt > 0) length(setdiff(alleles1, alleles2))/union_cnt else 0
  
  # 3) Transition asymmetry
  trans_asym <- if ((allele_gain + allele_loss) > 0)
    (allele_gain - allele_loss)/(allele_gain + allele_loss) else 0
  
  # 4) Prepare locus‐grouped allele lists once
  split1 <- split(sample1$allele, sample1$locus)
  split2 <- split(sample2$allele, sample2$locus)
  loci   <- union(names(split1), names(split2))
  n_loci <- length(loci)
  
  # 5) Compute discordant loci count
  discordant_loci <- sum(vapply(
    loci,
    function(l) { length(intersect(split1[[l]] %||% character(0),
                                   split2[[l]] %||% character(0))) == 0 },
    logical(1)
  ))
  
  # 6) Compute replacement pattern sum
  replacement_pattern_sum <- sum(vapply(
    loci,
    function(l) {
      a1 <- split1[[l]] %||% character(0)
      a2 <- split2[[l]] %||% character(0)
      if      (length(a2) == 0)         1
      else if (all(a2 %in% a1))       0
      else                             length(setdiff(a2, a1)) / length(a2)
    },
    numeric(1)
  ))
  
  # 7) Final locus‐level rates
  locus_discordance_rate    <- discordant_loci / n_loci
  replacement_pattern_score <- replacement_pattern_sum / n_loci
  
  # 8) Return all features
  c(
    jaccard_similarity          = jaccard,
    allele_retention_rate       = retention,
    allele_gain                 = allele_gain,
    locus_discordance_rate      = locus_discordance_rate,
    allele_transition_asymmetry = trans_asym,
    replacement_pattern_score   = replacement_pattern_score
  )
}

# Extract unique pairs once
unique_pairs <- unique(PAIRS_METADATA$PairsID)

# Create a named vector for mapping PairsID to D0 and Dx for faster access
pairs_to_d0 <- PAIRS_METADATA$NIDA1[match(unique_pairs, PAIRS_METADATA$PairsID)]
pairs_to_dx <- PAIRS_METADATA$NIDA2[match(unique_pairs, PAIRS_METADATA$PairsID)]

# Preprocess merged_dfs_filtered into a list for fast access
merged_dfs_filtered <- split(PAIRS_GENOMIC, PAIRS_GENOMIC$NIDA)
merged_dfs_filtered <- lapply(merged_dfs_filtered, function(df) { # reduces df size significantly
  df %>% select(-PairsID, -time_point) %>% distinct()
})

# Parallel processing setup
num_cores <- detectCores() - 1  # Using all cores except one for system stability
chunk_size <- 1000  # Adjust based on available memory

# Initialize the progress bar
pb <- progress_bar$new(
  format = "Processing [:bar] :percent :elapsedfull",
  total = ceiling(length(unique_pairs) / chunk_size),
  width = 100
)

# Output file for writing results
output_file <- paste0("lalala", site, "_", ".csv")

# Function to process each pair
process_pair <- function(pair_id) {
  sample1 <- merged_dfs_filtered[[pairs_to_d0[pair_id]]] 
  sample2 <- merged_dfs_filtered[[pairs_to_dx[pair_id]]]
  
  # Calculate features
  feats <- as.data.frame(t(calculate_features_optimized(sample1, sample2)))
  feats$PairsID <- unique_pairs[pair_id]
  
  return(feats)
}

# Function to write results in chunks to avoid memory bloat
write_chunk_to_csv <- function(data, file_path, append = FALSE) {
  if (append) {
    fwrite(data, file_path, append = TRUE)
  } else {
    fwrite(data, file_path)
  }
}

# Clear memory
gc()

# Process unique_pairs in chunks
for (i in seq(1, length(unique_pairs), by = chunk_size)) {
  # Define chunk range
  chunk_pair_ids <- i:min(i + chunk_size - 1, length(unique_pairs))
  
  print(paste0("Processing PairsID ", unique_pairs[chunk_pair_ids[1]], 
               " to ", unique_pairs[chunk_pair_ids[length(chunk_pair_ids)]], "..."))
  
  # Parallel processing of the current chunk
  results <- mclapply(chunk_pair_ids, process_pair, mc.cores = num_cores)
  
  # Filter out NULL results
  metrics_list_chunk <- Filter(Negate(is.null), results)
  
  # Check if there are valid results to write
  if (length(metrics_list_chunk) > 0) {
    # Combine results into a data table
    delta_features_df_chunk <- rbindlist(metrics_list_chunk, fill = TRUE)
    
    # Incrementally write results to CSV
    write_chunk_to_csv(delta_features_df_chunk, output_file, append = (i > 1))
  }
  
  pb$tick()
}

# Read and finalize the results
delta_metrics_df_final <- read.csv(paste0("lalala", site, "_", ".csv"))
delta_metrics_df_final <- delta_metrics_df_final %>%
  select(PairsID, everything())
write.csv(delta_metrics_df_final, paste0("lalala", site, "_", ".csv"), row.names = FALSE)

#####################################################
# CHECK if test_data vs delta_metrics_df_final match when DROPOUT = 0
# The dataframes should be identical

if (dropout_rate == 0){
  
  cols <- colnames(delta_metrics_df_final)
  test_data_sub <- test_data[cols]
  
  # Check if all elements are equal
  comparison <- delta_metrics_df_final == test_data_sub
  all_equal <- all(comparison, na.rm = TRUE)
  
  if (all_equal) {
    print("Validation successful: All data matches between test_data and delta_metrics_df_final")
  } else {
    # Find which columns have differences
    col_equality <- colSums(comparison, na.rm = TRUE) == nrow(comparison)
    mismatch_cols <- names(col_equality[!col_equality])
    print(paste("Validation failed: Differences found in columns:", paste(mismatch_cols, collapse = ", ")))
  }
  
}



#####################################################
# Predict probabilities on the test (holdout) set
TEST_META_threshs <- left_join(TEST_META, thresholds[c("pair_type", "decision_threshold")], by = "pair_type")
TEST_META_threshs$preds_prob <- predict(model, newdata = delta_metrics_df_final, type = "prob")[, "R"]
TEST_META_threshs$preds <- ifelse(TEST_META_threshs$preds_prob >= TEST_META_threshs$decision_threshold, "R", "NI")

# Initialize an empty results dataframe
results <- data.frame(pair_type = character(),
                      sensitivity = numeric(),
                      specificity = numeric(),
                      R_pairs = numeric(),
                      NI_pairs = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each unique pair_type combination
for (strain_comb in unique(TEST_META_threshs$pair_type)) {
  
  # Subset the TEST and TEST_labels based on the current combination
  subset_indices <- TEST_META_threshs$pair_type == strain_comb
  subset_TEST_labels <- TEST_META_threshs$labels[subset_indices]
  
  # Count R and NI pairs
  r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
  ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
  
  # Get the corresponding predictions for the current subset
  subset_preds <- TEST_META_threshs$preds[subset_indices]
  
  # Evaluate the confusion matrix for the current subset
  cm <- confusionMatrix(as.factor(subset_preds), as.factor(subset_TEST_labels), positive = "R")
  
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  # Append results
  results <- rbind(results, data.frame(pair_type = strain_comb,
                                       sensitivity = sens,
                                       specificity = spec,
                                       R_pairs = r,
                                       NI_pairs = ni,
                                       stringsAsFactors = FALSE))
}

# Join with thresholds
results <- left_join(results, thresholds[c("pair_type", "decision_threshold")], by = "pair_type")

results

# # Write final results
# write.csv(results, paste0("results_", site, ".csv"), row.names = FALSE)