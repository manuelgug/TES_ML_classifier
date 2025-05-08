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
library(ggplot2)

# Set site
site <- "Tete"

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
# Function to apply dropout to genomic data
apply_dropout <- function(data, dropout_rate = 0) {
  if (dropout_rate == 0) return(data)
  data_copy <- copy(data)
  
  # flag/drop…
  data_flagged <- data_copy %>%
    group_by(PairsID, time_point, locus) %>%
    mutate(
      allele_count = n(),
      dropout = ifelse(allele_count > 1,
                       sample(c(TRUE, FALSE), n(), replace = TRUE,
                              prob = c(dropout_rate, 1 - dropout_rate)),
                       FALSE)
    ) %>%
    ungroup()
  
  # # **Here** is your post-dropout summary:
  # cat("=== post-dropout allele counts ===\n")
  # data_flagged %>%
  #   filter(!dropout) %>%
  #   group_by(PairsID, time_point, locus) %>%
  #   summarize(count_after = n(), .groups = "drop") %>%
  #   ungroup() %>%
  #   count(count_after) %>%
  #   print()
  
  # now drop and return
  data_flagged %>%
    filter(!dropout) %>%
    select(-allele_count, -dropout)
}


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

# Function to calculate performance metrics
calculate_metrics <- function(TEST_META_threshs) {
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
  return(results)
}

# Function to run a single Monte Carlo iteration
run_montecarlo_iteration <- function(dropout_rate, iteration) {
  # Apply dropout to PAIRS_GENOMIC
  PAIRS_GENOMIC_dropped <- apply_dropout(PAIRS_GENOMIC, dropout_rate = dropout_rate)
  
  # Initialize the progress bar for this iteration
  pb <- progress_bar$new(
    format = paste0("Dropout Rate: ", dropout_rate, ", Iteration: ", iteration, 
                    " [:bar] :percent :elapsedfull"),
    total = ceiling(length(unique_pairs) / chunk_size),
    width = 100
  )
  
  # Output file for writing results
  output_file <- paste0("temp_", site, "_dropout_", dropout_rate, "_iter_", iteration, ".csv")
  
  # Update the merged_dfs_filtered with the new dropped data
  merged_dfs_filtered_iter <- split(PAIRS_GENOMIC_dropped, PAIRS_GENOMIC_dropped$NIDA)
  merged_dfs_filtered_iter <- lapply(merged_dfs_filtered_iter, function(df) {
    df %>% select(-PairsID, -time_point) %>% distinct()
  })
  
  # Process unique_pairs in chunks
  for (i in seq(1, length(unique_pairs), by = chunk_size)) {
    # Define chunk range
    chunk_pair_ids <- i:min(i + chunk_size - 1, length(unique_pairs))
    
    # Parallel processing of the current chunk
    results <- mclapply(chunk_pair_ids, function(pair_id) {
      # Use merged_dfs_filtered_iter instead of global merged_dfs_filtered
      sample1 <- merged_dfs_filtered_iter[[pairs_to_d0[pair_id]]] 
      sample2 <- merged_dfs_filtered_iter[[pairs_to_dx[pair_id]]]
      
      # Calculate features
      feats <- as.data.frame(t(calculate_features_optimized(sample1, sample2)))
      feats$PairsID <- unique_pairs[pair_id]
      
      return(feats)
    }, mc.cores = num_cores)
    
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
  delta_metrics_df_final <- fread(output_file)
  delta_metrics_df_final <- delta_metrics_df_final %>%
    select(PairsID, everything())
  
  # Predict probabilities on the test (holdout) set
  TEST_META_threshs <- left_join(TEST_META, thresholds[c("pair_type", "decision_threshold")], by = "pair_type")
  TEST_META_threshs$preds_prob <- predict(model, newdata = delta_metrics_df_final, type = "prob")[, "R"]
  TEST_META_threshs$preds <- ifelse(TEST_META_threshs$preds_prob >= TEST_META_threshs$decision_threshold, "R", "NI")
  
  # Calculate metrics
  results <- calculate_metrics(TEST_META_threshs)
  
  # Add dropout rate and iteration information
  results$dropout_rate <- dropout_rate
  results$iteration <- iteration
  
  # Clean up temporary file
  file.remove(output_file)
  
  return(results)
}

###########################################
# Run Monte Carlo simulation for different dropout rates
###########################################

# Define dropout rates to test
dropout_rates <- seq(0, 0.5, by = 0.25)
mc_iterations <- 2  # 10 iterations for each dropout rate

# Initialize empty dataframe to store all Monte Carlo results
all_results <- data.frame()

# Run simulations for each dropout rate
for (dropout_rate in dropout_rates) {
  if (dropout_rate == 0) {
    # For dropout rate 0, run once without Monte Carlo
    cat("Processing dropout rate 0 (no Monte Carlo needed)...\n")
    
    # Apply dropout to PAIRS_GENOMIC (which does nothing at rate 0)
    PAIRS_GENOMIC_dropped <- apply_dropout(PAIRS_GENOMIC, dropout_rate = 0)
    
    # Initialize the progress bar
    pb <- progress_bar$new(
      format = "Processing dropout rate 0 [:bar] :percent :elapsedfull",
      total = ceiling(length(unique_pairs) / chunk_size),
      width = 100
    )
    
    # Output file for writing results
    output_file <- paste0("temp_", site, "_dropout_0.csv")
    
    # Update merged_dfs_filtered with dropout applied data (no change at rate 0)
    merged_dfs_filtered_rate0 <- split(PAIRS_GENOMIC_dropped, PAIRS_GENOMIC_dropped$NIDA)
    merged_dfs_filtered_rate0 <- lapply(merged_dfs_filtered_rate0, function(df) {
      df %>% select(-PairsID, -time_point) %>% distinct()
    })
    
    # Process unique_pairs in chunks
    for (i in seq(1, length(unique_pairs), by = chunk_size)) {
      # Define chunk range
      chunk_pair_ids <- i:min(i + chunk_size - 1, length(unique_pairs))
      
      # Parallel processing of the current chunk
      results <- mclapply(chunk_pair_ids, function(pair_id) {
        # Use the local merged_dfs_filtered_rate0
        sample1 <- merged_dfs_filtered_rate0[[pairs_to_d0[pair_id]]] 
        sample2 <- merged_dfs_filtered_rate0[[pairs_to_dx[pair_id]]]
        
        # Calculate features
        feats <- as.data.frame(t(calculate_features_optimized(sample1, sample2)))
        feats$PairsID <- unique_pairs[pair_id]
        
        return(feats)
      }, mc.cores = num_cores)
      
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
    delta_metrics_df_final <- fread(output_file)
    delta_metrics_df_final <- delta_metrics_df_final %>%
      select(PairsID, everything())
    
    # CHECK if test_data vs delta_metrics_df_final match when DROPOUT = 0
    # The dataframes should be identical
    cols <- colnames(delta_metrics_df_final)
    test_data_sub <- test_data[cols]
    
    # Convert both to data.frame for comparison
    delta_metrics_df_final_df <- as.data.frame(delta_metrics_df_final)
    test_data_sub_df <- as.data.frame(test_data_sub)
    
    # Check if all elements are equal
    all_equal <- TRUE
    mismatch_cols <- character(0)
    
    # First check if dimensions match
    if (!identical(dim(delta_metrics_df_final_df), dim(test_data_sub_df))) {
      all_equal <- FALSE
      cat("Validation failed: Dimensions don't match\n")
    } 
    # Then check if column names match
    else if (!identical(colnames(delta_metrics_df_final_df), colnames(test_data_sub_df))) {
      all_equal <- FALSE
      cat("Validation failed: Column names don't match\n")
    } 
    # Then check column by column
    else {
      for (col in colnames(delta_metrics_df_final_df)) {
        col1 <- delta_metrics_df_final_df[[col]]
        col2 <- test_data_sub_df[[col]]
        
        # Handle different column types
        if (is.numeric(col1) && is.numeric(col2)) {
          # For numeric columns, use tolerance
          if (!all(abs(col1 - col2) < 1e-10, na.rm = TRUE)) {
            all_equal <- FALSE
            mismatch_cols <- c(mismatch_cols, col)
          }
        } 
        else if (is.character(col1) && is.character(col2) ||
                 is.factor(col1) && is.factor(col2) ||
                 is.logical(col1) && is.logical(col2)) {
          # For character, factor, or logical columns
          if (!all(col1 == col2, na.rm = TRUE)) {
            all_equal <- FALSE
            mismatch_cols <- c(mismatch_cols, col)
          }
        }
        else {
          # If types don't match, convert to character
          if (!all(as.character(col1) == as.character(col2), na.rm = TRUE)) {
            all_equal <- FALSE
            mismatch_cols <- c(mismatch_cols, col)
          }
        }
      }
    }
    
    if (all_equal) {
      cat("Validation successful: All data matches between test_data and delta_metrics_df_final\n")
    } else {
      if (length(mismatch_cols) > 0) {
        cat(paste("Validation failed: Differences found in columns:", paste(mismatch_cols, collapse = ", "), "\n"))
      }
    }
    
    # Predict probabilities on the test (holdout) set
    TEST_META_threshs <- left_join(TEST_META, thresholds[c("pair_type", "decision_threshold")], by = "pair_type")
    TEST_META_threshs$preds_prob <- predict(model, newdata = delta_metrics_df_final, type = "prob")[, "R"]
    TEST_META_threshs$preds <- ifelse(TEST_META_threshs$preds_prob >= TEST_META_threshs$decision_threshold, "R", "NI")
    
    # Calculate metrics
    results <- calculate_metrics(TEST_META_threshs)
    
    # Add dropout rate and iteration information
    results$dropout_rate <- 0
    results$iteration <- 0  # No iteration for dropout rate 0
    
    # Add to all results
    all_results <- rbind(all_results, results)
    
    # Clean up temporary file
    file.remove(output_file)
    
  } else {
    # For other dropout rates, run Monte Carlo iterations
    cat(paste0("Processing dropout rate ", dropout_rate, " with ", mc_iterations, " Monte Carlo iterations...\n"))
    
    for (iter in 1:mc_iterations) {
      cat(paste0("  - Running iteration ", iter, " of ", mc_iterations, "\n"))
      
      # Set seed for reproducibility but different for each iteration
      set.seed(1000 * dropout_rate + iter)
      
      # Run Monte Carlo iteration
      iter_results <- run_montecarlo_iteration(dropout_rate, iter)
      
      # Add to all results
      all_results <- rbind(all_results, iter_results)
    }
  }
}

# Save all Monte Carlo results
write.csv(all_results, paste0(site, "_montecarlo_results.csv"), row.names = FALSE)

###########################################
# Plot Results
###########################################

# Calculate mean and standard deviation of sensitivity and specificity for each dropout rate and pair_type
summary_results <- all_results %>%
  group_by(dropout_rate, pair_type) %>%
  summarize(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    mean_specificity = mean(specificity, na.rm = TRUE),
    sd_specificity = sd(specificity, na.rm = TRUE),
    n_iterations = n(),
    .groups = "drop"
  )

# Plot sensitivity by dropout rate for each pair_type
sensitivity_plot <- ggplot(summary_results, aes(x = dropout_rate, y = mean_sensitivity, color = pair_type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_sensitivity - sd_sensitivity, 
                    ymax = mean_sensitivity + sd_sensitivity), 
                width = 0.01) +
  labs(title = paste0(site, " - Sensitivity by Dropout Rate"),
       x = "Dropout Rate",
       y = "Sensitivity",
       color = "Pair Type") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    legend.background = element_rect(fill = "white")
  )+
  ylim(0,1)

# Plot specificity by dropout rate for each pair_type
specificity_plot <- ggplot(summary_results, aes(x = dropout_rate, y = mean_specificity, color = pair_type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_specificity - sd_specificity, 
                    ymax = mean_specificity + sd_specificity), 
                width = 0.01) +
  labs(title = paste0(site, " - Specificity by Dropout Rate"),
       x = "Dropout Rate",
       y = "Specificity",
       color = "Pair Type") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    legend.background = element_rect(fill = "white")
  )+
  ylim(0,1)

# Save plots
ggsave(paste0(site, "_sensitivity_by_dropout.png"), sensitivity_plot, width = 10, height = 6, bg = "white")
ggsave(paste0(site, "_specificity_by_dropout.png"), specificity_plot, width = 10, height = 6, bg = "white")

# Create a combined plot for all pair types
# Long format for easier faceting
plot_data <- summary_results %>%
  pivot_longer(
    cols = c(mean_sensitivity, mean_specificity),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = case_when(
      metric == "mean_sensitivity" ~ "Sensitivity",
      metric == "mean_specificity" ~ "Specificity"
    ),
    sd = case_when(
      metric == "Sensitivity" ~ sd_sensitivity,
      metric == "Specificity" ~ sd_specificity
    )
  )

# Combined plot
combined_plot <- ggplot(plot_data, aes(x = dropout_rate, y = value, color = pair_type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.01, alpha = 0.7) +
  facet_grid(metric ~ ., scales = "free_y") +
  labs(title = paste0(site, " - Performance Metrics by Dropout Rate"),
       x = "Dropout Rate",
       y = "Value",
       color = "Pair Type") +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    legend.background = element_rect(fill = "white")
  )+
  ylim(0,1)

# Save combined plot
ggsave(paste0(site, "_combined_metrics_by_dropout.png"), combined_plot, width = 12, height = 8, bg = "white")

cat("\nMonte Carlo simulation complete. Results saved to:\n")
cat(paste0("  - ", site, "_montecarlo_results.csv\n"))
cat(paste0("  - ", site, "_sensitivity_by_dropout.png\n"))
cat(paste0("  - ", site, "_specificity_by_dropout.png\n"))
cat(paste0("  - ", site, "_combined_metrics_by_dropout.png\n"))







######################################################################################################

# Example for a single pair in the “2__1” group:
pair_id <- as.numeric(TEST_META$PairsID[TEST_META$pair_type == "2__1"][11])

# Extract the raw genomic tables for that pair:
sample1_orig <- merged_dfs_filtered[[ pairs_to_d0[pair_id] ]]
sample2_orig <- merged_dfs_filtered[[ pairs_to_dx[pair_id] ]]

# Do two independent dropouts:
drop1 <- apply_dropout(PAIRS_GENOMIC[PAIRS_GENOMIC$NIDA %in% c(pairs_to_d0[pair_id], pairs_to_dx[pair_id])], 0.25)
drop2 <- apply_dropout(PAIRS_GENOMIC[PAIRS_GENOMIC$NIDA %in% c(pairs_to_d0[pair_id], pairs_to_dx[pair_id])], 0.25)

# Split them again by NIDA:
s1a <- drop1[ drop1$NIDA == pairs_to_d0[pair_id], ]
s2a <- drop1[ drop1$NIDA == pairs_to_dx[pair_id], ]
s1b <- drop2[ drop2$NIDA == pairs_to_d0[pair_id], ]
s2b <- drop2[ drop2$NIDA == pairs_to_dx[pair_id], ]

# Calculate features:
feat1 <- calculate_features_optimized(s1a, s2a)
feat2 <- calculate_features_optimized(s1b, s2b)
print(feat1)
print(feat2)



apply_dropout <- function(data, dropout_rate = 0) {
  if (dropout_rate == 0) return(data)
  # 1) keep only one row per allele-locus
  alleles_only <- data %>%
    distinct(PairsID, time_point, locus, allele)
  # 2) flag dropout on those unique alleles
  alleles_flagged <- alleles_only %>%
    group_by(PairsID, time_point, locus) %>%
    mutate(drop = sample(c(TRUE, FALSE), n(), replace = TRUE,
                         prob = c(dropout_rate, 1 - dropout_rate))) %>%
    ungroup()
  # 3) keep only non-dropped alleles
  kept_alleles <- alleles_flagged %>% filter(!drop) %>% select(-drop)
  return(kept_alleles)
}

drop1 <- apply_dropout(PAIRS_GENOMIC[PAIRS_GENOMIC$NIDA %in% c(pairs_to_d0[pair_id], pairs_to_dx[pair_id])], 0.25)
drop2 <- apply_dropout(PAIRS_GENOMIC[PAIRS_GENOMIC$NIDA %in% c(pairs_to_d0[pair_id], pairs_to_dx[pair_id])], 0.25)
s1a <- drop1[ drop1$NIDA == pairs_to_d0[pair_id], ]
s2a <- drop1[ drop1$NIDA == pairs_to_dx[pair_id], ]
s1b <- drop2[ drop2$NIDA == pairs_to_d0[pair_id], ]
s2b <- drop2[ drop2$NIDA == pairs_to_dx[pair_id], ]
feat1 <- calculate_features_optimized(s1a, s2a)
feat2 <- calculate_features_optimized(s1b, s2b)
