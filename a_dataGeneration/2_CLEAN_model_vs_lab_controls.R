
library(dplyr)
library(ggplot2)
library(purrr)

# createggplot2# create model with TRUTH

#####################################################################################
############################# TRUTH ##############################
#####################################################################################

TRUTH_CLEAN <- read.csv("TRUTH_CLEAN.csv")

############# SCRIPT 4 MODED FOR TRUTH STRAINS @@@@@@@@@@@@

# --- 3. Create All Mixes ----
min_coi <- 1
max_coi <- length(unique(TRUTH_CLEAN$Strain))

nidas <- unique(TRUTH_CLEAN$Strain)


# Max number of combinations to keep per COI level (to avoid memory explosion)
MAX_COMBOS <- 10000

create_combinations_df_safe <- function(vec, k, max_combos = MAX_COMBOS) {
  total_combos <- choose(length(vec), k)
  
  if (total_combos > max_combos) {
    # Sample random combos without generating full combn matrix
    sampled_combos <- replicate(max_combos, sort(sample(vec, k)), simplify = FALSE)
    comb_df <- as.data.frame(do.call(rbind, sampled_combos), stringsAsFactors = FALSE)
  } else {
    comb_df <- as.data.frame(t(combn(vec, k)), stringsAsFactors = FALSE)
  }
  
  setNames(comb_df, paste0("strain_", 1:k))
}


# Build mixes with sampling protection
coi_values <- seq(min_coi+1, max_coi, 1)
strain_mixes <- map(coi_values, ~ create_combinations_df_safe(nidas, .x, MAX_COMBOS))

# Add mix_1 (monoclonals)
strain_mixes <- setNames(c(list(mix_1 = data.frame(strain_1 = nidas)),
                           strain_mixes),
                         c("mix_1", paste0("mix_", coi_values)))

# subset given the lab_controls data!
strain_mixes <- strain_mixes[c("mix_1", "mix_2", "mix_3", "mix_5")]


############################################

#just to not crash the script
TRUTH_CLEAN$reads <- 1

initial_sample_size <- 10000



# --- 4. Subsample Mixes ----
strain_mixes_subsampled <- lapply(strain_mixes, \(mix_df) {
  sample_n(mix_df, size = min(nrow(mix_df), initial_sample_size))
})

# --- 5. Create Metadata & Genomic Mixes ----
max_coi_len <- max_coi
metadata_list <- list()
genomic_list <- list()

meta_counter <- 1
geno_counter <- 1

for (mix_num in seq(min_coi, max_coi)) {
  mix_name <- paste0("mix_", mix_num)
  if (!mix_name %in% names(strain_mixes_subsampled)) next
  
  current_mix <- strain_mixes_subsampled[[mix_name]]
  
  for (i in seq_len(nrow(current_mix))) {
    strains <- unlist(current_mix[i, ], use.names = FALSE)
    strains <- strains[strains != ""]  # Remove blanks
    
    selected <- TRUTH_CLEAN %>%
      filter(Strain %in% strains)
    
    total_reads <- sum(selected$reads)
    
    prop_reads <- selected %>%
      group_by(Strain) %>%
      summarise(prop_reads = sum(reads) / total_reads, .groups = 'drop')
    
    mixID <- paste0("mix", length(strains), "_ID", i)
    
    # Build metadata row
    row_metadata <- c(
      mixID,
      strains,
      rep(NA, max_coi_len - length(strains)),
      prop_reads$prop_reads,
      rep(NA, max_coi_len - length(strains))
    )
    
    metadata_list[[meta_counter]] <- row_metadata
    meta_counter <- meta_counter + 1
    
    # Genomic info
    rows_genomic <- selected %>%
      group_by(locus, allele) %>%
      summarise(read_counts = sum(reads), .groups = "drop") %>%
      mutate(mixID = mixID) %>%
      group_by(mixID, locus) %>%
      mutate(
        norm.reads.locus = read_counts / sum(read_counts),
        n.alleles = n_distinct(allele)
      ) %>%
      ungroup()
    
    genomic_list[[geno_counter]] <- rows_genomic
    geno_counter <- geno_counter + 1
  }
}

# Finalize Metadata
strain_cols <- paste0("strain", seq_len(max_coi_len))
prop_cols <- paste0("strain_prop", seq_len(max_coi_len))

MIXES_METADATA <- do.call(rbind, metadata_list) %>%
  as.data.frame(stringsAsFactors = FALSE)
colnames(MIXES_METADATA) <- c("NIDA", strain_cols, prop_cols)

MIXES_METADATA[strain_cols] <- lapply(MIXES_METADATA[strain_cols], as.character)
MIXES_METADATA[prop_cols] <- lapply(MIXES_METADATA[prop_cols], as.numeric)

# Finalize Genomic
MIXES_GENOMIC <- bind_rows(genomic_list)

# --- 6. Validate ----
stopifnot(nrow(MIXES_METADATA) == length(unique(MIXES_GENOMIC$mixID)))
stopifnot(all(MIXES_METADATA$NIDA %in% MIXES_GENOMIC$mixID))

# check
length(unique(MIXES_GENOMIC$mixID))
nrow(MIXES_METADATA)


# --- 7. Save ----
saveRDS(MIXES_METADATA, paste0("MIXES_METADATA_TRUTH_CLEAN.RDS"))
saveRDS(MIXES_GENOMIC, paste0("MIXES_GENOMIC_TRUTH_CLEAN.RDS"))

# --- 8. EDA Plot ----
alleles_per_mix <- MIXES_GENOMIC %>%
  group_by(mixID) %>%
  summarise(alleles_per_mix = n_distinct(allele), .groups = "drop") %>%
  separate(mixID, into = c("mix_type", "ID"), sep = "_", remove = FALSE)

allele_plot <- ggplot(alleles_per_mix, aes(x = alleles_per_mix, fill = mix_type)) +
  geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
  theme_minimal() +
  labs(x = "Number of Unique Alleles", y = "Count", fill = "Mix Type") +
  xlim(0, NA)

allele_plot

ggsave(paste0("mixes_EDA_TRUTH_CLEAN.png"), allele_plot, width = 9, height = 6, dpi = 300, bg = "white")




############# SCRIPT 5 MODED FOR TRUTH STRAINS @@@@@@@@@@@@

#------------------------------------------------
# 2. Create pairs
#------------------------------------------------

MIXES_METADATA <- readRDS(paste0("MIXES_METADATA_TRUTH_CLEAN.RDS"))
MIXES_GENOMIC <- readRDS(paste0("MIXES_GENOMIC_TRUTH_CLEAN.RDS"))

# Generate all combinations
unique_combos <- expand.grid(
  offset_naive_coi_D0 = 1:7,
  offset_naive_coi_Dx = 1:7
)

nidas_all <- MIXES_METADATA$NIDA

pairs_df <- expand.grid(NIDA1 = nidas_all, NIDA2 = nidas_all, stringsAsFactors = FALSE)

pairs_df <- pairs_df[rowSums(sapply(1:nrow(unique_combos), function(i) {
  grepl(unique_combos$offset_naive_coi_D0[i], pairs_df$NIDA1) &
    grepl(unique_combos$offset_naive_coi_Dx[i], pairs_df$NIDA2)
})) > 0, ]

pairs_df <- pairs_df %>%
  mutate(PairsID = row_number()) %>%
  pivot_longer(cols = c(NIDA1, NIDA2), names_to = "time_point", values_to = "NIDA") %>%
  mutate(time_point = ifelse(time_point == "NIDA1", "D0", "Dx")) %>%
  select(PairsID, NIDA, time_point)

colnames(MIXES_GENOMIC)[colnames(MIXES_GENOMIC) == "mixID"] <- "NIDA"
merged_dfs <- left_join(pairs_df, MIXES_GENOMIC, by = "NIDA")

length(unique(merged_dfs$PairsID))

gc()
saveRDS(merged_dfs, paste0("PAIRS_GENOMIC_TRUTH_CLEAN.RDS"))

#------------------------------------------------
# 3. Compare allele content of pairs
#------------------------------------------------

alleles <- merged_dfs %>%
  group_by(PairsID, NIDA, time_point) %>%
  summarize(alleles = list(allele), .groups = "drop") %>%
  arrange(PairsID, time_point)

gc()

alleles_shared_prop <- alleles %>%
  group_by(PairsID) %>%
  summarize(
    NIDA1 = NIDA[1],
    NIDA2 = NIDA[2],
    shared_count = length(intersect(alleles[[1]], alleles[[2]])),
    union_count = length(union(alleles[[1]], alleles[[2]])),
    shared_prop = shared_count / union_count,
    .groups = "drop"
  ) %>%
  mutate(
    NIDA1_trimmed = sub("_.*", "", NIDA1),
    NIDA2_trimmed = sub("_.*", "", NIDA2),
    pair_type = paste0(NIDA1_trimmed, "_", NIDA2_trimmed)
  ) %>%
  select(-NIDA1_trimmed, -NIDA2_trimmed)

#------------------------------------------------
# 4. Label data
#------------------------------------------------

MIXES_METADATA <- MIXES_METADATA %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  rowwise() %>%
  mutate(
    STRAINS = list(na.omit(c_across(matches("^strain\\d?$"))))
  ) %>%
  ungroup()

PAIRS <- merge(pairs_df, MIXES_METADATA[c("NIDA", "STRAINS")], by = "NIDA") %>%
  arrange(PairsID, time_point)

labels <- PAIRS %>%
  group_by(PairsID) %>%
  summarise(
    labels = ifelse(
      any(unlist(STRAINS[time_point == "Dx"]) %in% unlist(STRAINS[time_point == "D0"])),
      "R", "NI"
    )
  )

PAIRS <- PAIRS %>%
  rowwise() %>%
  mutate(
    nstrains = length(STRAINS),
    sample = paste(STRAINS, collapse = "_")
  )

PAIRS_wide <- PAIRS %>%
  pivot_wider(
    names_from = time_point,
    values_from = c(sample, nstrains),
    names_glue = "{time_point}_{.value}"
  )

PAIRS_metadata <- PAIRS_wide %>%
  group_by(PairsID) %>%
  select(-NIDA, -STRAINS) %>%
  reframe(
    D0_sample = first(na.omit(D0_sample)),
    Dx_sample = first(na.omit(Dx_sample)),
    D0_nstrains = first(na.omit(D0_nstrains)),
    Dx_nstrains = first(na.omit(Dx_nstrains))
  )

PAIRS_metadata <- inner_join(PAIRS_metadata, labels, by = "PairsID")
PAIRS_metadata <- inner_join(PAIRS_metadata, alleles_shared_prop, by = "PairsID")
PAIRS_metadata$pair_type <- paste0(PAIRS_metadata$D0_nstrains, "__", PAIRS_metadata$Dx_nstrains)

saveRDS(PAIRS_metadata, paste0("PAIRS_METADATA_TRUTH_CLEAN.RDS"))

#------------------------------------------------
# 5. Pair metadata mini EDA
#------------------------------------------------

PAIRS_summary <- PAIRS_metadata %>%
  group_by(pair_type) %>%
  summarise(
    n_pairs = n(),
    NI_prop = mean(labels == "NI"),
    R_prop = mean(labels == "R"),
    # median_shared_prop = median(shared_prop, na.rm = TRUE),
    NI_size = NI_prop * n_pairs,
    R_size = R_prop * n_pairs
  )

PAIRS_summary

write.csv(PAIRS_summary, paste0("PAIRS_SUMMARY_TRUTH_CLEAN.csv"), row.names = FALSE)




############# SCRIPT 6 MODED FOR TRUTH STRAINS @@@@@@@@@@@@

library(data.table)
library(dplyr)
library(dcifer)
library(reshape2)
library(tidyr)
library(purrr)
library(stringr)
library(Matrix)
library(progress)
library(parallel)


#select data type betweem "TRAINING_DATA" or "REAL_DATA"
DATA_TYPE = "TRAINING_DATA"


if (DATA_TYPE == "TRAINING_DATA") {
  
  # for the training data:
  PAIRS_METADATA <- readRDS(paste0("PAIRS_METADATA_TRUTH_CLEAN.RDS"))
  PAIRS_GENOMIC <- readRDS(paste0("PAIRS_GENOMIC_TRUTH_CLEAN.RDS"))
  
  ### REMOVE THE IMPISSIBLE PAIRS (IN THIS CASE, SINCE THERE ARE ONLY 8 STRAINS, EVERYTHING ABOVE 9 WUEH ADDING UP d0 AND dX IS NOT USEFUL)
  max_coi_TRUTH <- max(PAIRS_METADATA$D0_nstrains, PAIRS_METADATA$Dx_nstrains)
  remove_pairs <- PAIRS_METADATA[PAIRS_METADATA$D0_nstrains + PAIRS_METADATA$Dx_nstrains > max_coi, ]$PairsID
  
  PAIRS_METADATA <- PAIRS_METADATA[!PAIRS_METADATA$PairsID %in% remove_pairs,]
  PAIRS_GENOMIC <- PAIRS_GENOMIC[!PAIRS_GENOMIC$PairsID %in% remove_pairs,]
  
  
  PAIRS_GENOMIC <- as.data.table(PAIRS_GENOMIC)
  PAIRS_METADATA <- as.data.table(PAIRS_METADATA)
  
  PAIRS_GENOMIC <- left_join(PAIRS_GENOMIC, PAIRS_METADATA, by = "PairsID") 
  
} else if (DATA_TYPE == "REAL_DATA") {
  
  #the actual data
  PAIRS_GENOMIC <- read.csv(paste0("genomic_updated_",site,".csv"))
  PAIRS_METADATA <- read.csv(paste0("metadata_updated_",site,".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
  
  PAIRS_METADATA <- PAIRS_METADATA[!is.na(PAIRS_METADATA$time_point),]
  
  PAIRS_GENOMIC <- as.data.table(PAIRS_GENOMIC)
  PAIRS_METADATA <- as.data.table(PAIRS_METADATA)
  
  PAIRS_GENOMIC <- PAIRS_GENOMIC %>% rename(read_counts = reads)  
  PAIRS_GENOMIC <- PAIRS_GENOMIC %>% rename(NIDA = sampleID) # format genomic file
  suppressWarnings(PAIRS_GENOMIC <- PAIRS_GENOMIC %>%
                     separate(NIDA, into = c("NIDA", "run"), sep = "__", remove = TRUE))
  PAIRS_GENOMIC <- inner_join(PAIRS_METADATA, PAIRS_GENOMIC, by = "NIDA")
  
  # Step 1: Extract offset_naive_coi at D0 and Dx per PairsID
  coi_summary <- PAIRS_GENOMIC[time_point %in% c("D0", "Dx"), 
                               .(offset_naive_coi_D0 = unique(offset_naive_coi[time_point == "D0"]),
                                 offset_naive_coi_Dx = unique(offset_naive_coi[time_point == "Dx"])),
                               by = PairsID]
  
  # Step 2: Create pair_type as "D0__Dx"
  coi_summary[, pair_type := paste0(offset_naive_coi_D0, "__", offset_naive_coi_Dx)]
  
  # Step 3: Merge back into the full dataset
  PAIRS_GENOMIC <- merge(PAIRS_GENOMIC, coi_summary[, .(PairsID, pair_type)], by = "PairsID", all.x = TRUE)
  
  PAIRS_METADATA <- PAIRS_METADATA %>% select(PairsID, NIDA, time_point) # format metadata file
  
  PAIRS_METADATA <- PAIRS_METADATA %>%
    pivot_wider(names_from = time_point, values_from = NIDA, names_prefix = "NIDA") %>%
    rename(NIDA1 = NIDAD0, NIDA2 = NIDADx)
  
  PAIRS_METADATA <- left_join(PAIRS_METADATA, coi_summary, by = "PairsID")
  
  
} else {
  
  print("Incorrect data type. Options are 'TRAINING_DATA' and 'REAL_DATA'.")
  
}

###

calculate_features_optimized <- function(sample1, sample2) {
  # 1) Unique alleles & loci
  alleles1 <- unique(sample1$allele)
  alleles2 <- unique(sample2$allele)
  all_alleles <- union(alleles1, alleles2)
  
  # 2) Allele-level intersection & union via set operations
  inter_cnt   <- length(intersect(alleles1, alleles2))
  union_cnt   <- length(all_alleles)
  jaccard     <- if (union_cnt > 0) inter_cnt / union_cnt else 0
  retention   <- if (length(alleles1) > 0) inter_cnt / length(alleles1) else 0
  allele_gain <- if (union_cnt > 0) length(setdiff(alleles2, alleles1)) / union_cnt else 0
  allele_loss <- if (union_cnt > 0) length(setdiff(alleles1, alleles2)) / union_cnt else 0
  
  # 3) Transition asymmetry
  trans_asym <- if ((allele_gain + allele_loss) > 0)
    (allele_gain - allele_loss) / (allele_gain + allele_loss) else 0
  
  # 4) Prepare locus-grouped allele lists
  split1 <- split(sample1$allele, sample1$locus)
  split2 <- split(sample2$allele, sample2$locus)
  loci   <- union(names(split1), names(split2))
  n_loci <- length(loci)
  
  # 5) Compute discordant loci count
  discordant_loci <- sum(vapply(
    loci,
    function(l) {
      length(intersect(split1[[l]] %||% character(0),
                       split2[[l]] %||% character(0))) == 0
    },
    logical(1)
  ))
  
  # 6) Compute replacement pattern sum
  replacement_pattern_sum <- sum(vapply(
    loci,
    function(l) {
      a1 <- split1[[l]] %||% character(0)
      a2 <- split2[[l]] %||% character(0)
      if      (length(a2) == 0)         1
      else if (all(a2 %in% a1))         0
      else                              length(setdiff(a2, a1)) / length(a2)
    },
    numeric(1)
  ))
  
  # 7) Final locus-level rates
  locus_discordance_rate    <- discordant_loci / n_loci
  replacement_pattern_score <- replacement_pattern_sum / n_loci
  
  # 8) Return all features
  c(
    jaccard_similarity          = jaccard,
    allele_retention_rate       = retention,
    allele_gain                 = allele_gain,
    allele_loss                 = allele_loss,
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

gc()

# Parallel processing setup
num_cores <- detectCores() - 0
chunk_size <- 1000  # Adjust based on available memory

# Initialize the progress bar
pb <- progress_bar$new(
  format = "Processing [:bar] :percent :elapsedfull",
  total = ceiling(length(unique_pairs) / chunk_size),
  width = 100
)

# Temporary file for writing output
output_file <- paste0("delta_features_TRUTH_CLEAN_", DATA_TYPE, ".csv")

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

delta_metrics_df_final <- read.csv(paste0("delta_features_TRUTH_CLEAN_", DATA_TYPE, ".csv"))

delta_metrics_df_final <- delta_metrics_df_final %>%
  select(PairsID, everything())

write.csv(delta_metrics_df_final, paste0("delta_features_TRUTH_CLEAN_", DATA_TYPE, ".csv"), row.names = F)



##### MERGE FEATURES AND OUTPUT ------

# FEATURES <- left_join(dres0_long_final_summarized, delta_metrics_df_final, by = "PairsID")

FEATURES <- left_join(delta_metrics_df_final, PAIRS_METADATA, by = "PairsID")

#coi change as feature
if (DATA_TYPE == "REAL_DATA"){
  
  FEATURES$coi_change <- FEATURES$offset_naive_coi_Dx - FEATURES$offset_naive_coi_D0
  
} else {
  
  FEATURES$coi_change <- FEATURES$Dx_nstrains - FEATURES$D0_nstrains
  
}



write.csv(FEATURES, paste0("TRUTH_CLEAN_", DATA_TYPE,".csv"), row.names = F)



############# SCRIPT 7 MODED FOR TRUTH STRAINS @@@@@@@@@@@@

library(caret)    
library(dplyr)    
library(tidyr)    
library(ggplot2)  
library(broom)
library(glmnet)


site <- "TRUTH_CLEAN"


### 1) IMPORT TRAINING AND REAL DATA ----------

TRAINING_DATA <- read.csv(paste0(site, "_TRAINING_DATA.csv"), row.names = 1) 
TRAINING_DATA$eCOI_pairs <- paste0(TRAINING_DATA$D0_nstrains, "__", TRAINING_DATA$Dx_nstrains) # Create a new variable by combining 'D0nstrains' and 'Dxnstrains'
LABELS <- data.frame(labels = TRAINING_DATA$labels)
LABELS$labels <- as.factor(LABELS$labels)

#REAL_DATA <- read.csv(paste0(site, "_REAL_DATA.csv"), stringsAsFactors = FALSE, colClasses = c(NIDA1 = "character", NIDA2= "character")) 

features_to_use <- colnames(TRAINING_DATA)[!colnames(TRAINING_DATA) %in% c("PairsID", "NIDA1", "NIDA2", "pair_type", "IBD_estimate",
                                                                           "offset_naive_coi_D0", "offset_naive_coi_Dx", 
                                                                           "replacement_pattern_score", "locus_discordance_rate", 
                                                                           "D0_sample", "Dx_sample","D0_nstrains", "Dx_nstrains", "eCOI_pairs", "labels", "shared_count","union_count", "shared_prop")]

corrplot::corrplot(cor(TRAINING_DATA %>% select(features_to_use), use = "complete.obs"), "pie")


### 2) SPLIT DATA ------------

set.seed(420)
# Step 1: Stratified Sampling by `pair_type`
train_indices <- createDataPartition(TRAINING_DATA$eCOI_pairs, p = 0.7, list = FALSE)
train_data <- TRAINING_DATA[train_indices, ]
test_data <- TRAINING_DATA[-train_indices, ]


## output TEST for later use in false negative tests
test_data$PairsID <- rownames(test_data)
test_data %>% select(PairsID, everything(), -eCOI_pairs)
saveRDS(test_data, paste0(site, "_test_data.RDS"))

# Step 2: Extract Metadata and Labels
TRAIN_META <- train_data %>% select(-all_of(features_to_use))
TRAIN <- train_data %>% select(all_of(features_to_use))
TRAIN_labels <- LABELS[train_indices, ]

TEST_META <- test_data %>% select(-all_of(features_to_use))
TEST <- test_data %>% select(all_of(features_to_use))
TEST_labels <- LABELS[-train_indices, ]

# Step 3: Verify Stratification
prop.table(table(TRAIN_META$eCOI_pairs))
prop.table(table(TRAIN_labels))
table(TRAIN_labels)

prop.table(table(TEST_META$eCOI_pairs))
prop.table(table(TEST_labels))
table(TEST_labels)



##### 4) TEST MODEL USING IBD ONLY (BEST FEATURE)--------

# Create training and test data frames for the current feature
df_train_IBD <- data.frame(TRAIN[features_to_use], label = as.factor(TRAIN_labels))
df_test_IBD <- data.frame(TEST[features_to_use], label = as.factor(TEST_labels))

# # Set up 10-fold cross-validation
# ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
# 
# # Re-level factor so that "R" is the positive class
# df_train_IBD$label <- relevel(df_train_IBD$label, ref = "R")
# 
# # Train logistic regression model with cross-validation
# fit_IBD <- train(label ~ ., 
#                  data = df_train_IBD, 
#                  method = "glm", 
#                  family = "binomial", 
#                  trControl = ctrl, 
#                  metric = "ROC")
# 
# print(fit_IBD)
# 
# fit_IBD$finalModel
# 
# ### feature importance
# # Extract coefficients from the final model
# coefs <- summary(fit_IBD$finalModel)$coefficients
# coefs_df <- as.data.frame(coefs)
# coefs_df$Variable <- rownames(coefs_df)
# 
# # Remove intercept for feature importance plot
# coefs_df <- coefs_df[coefs_df$Variable != "(Intercept)", ]
# 
# # Calculate absolute coefficient values to rank by importance
# coefs_df$AbsEstimate <- log(abs(coefs_df$Estimate))
# 
# # Sort by absolute coefficient value
# coefs_df <- coefs_df[order(coefs_df$AbsEstimate, decreasing = TRUE), ]
# 
# # Create a color vector (positive coefficients in blue, negative in red)
# coefs_df$Color <- ifelse(coefs_df$Estimate > 0, "positive", "negative")
# 
# 
# # Add significance stars
# coefs_df$Significance <- ifelse(coefs_df$`Pr(>|z|)` < 0.001, "***",
#                                 ifelse(coefs_df$`Pr(>|z|)` < 0.01, "**",
#                                        ifelse(coefs_df$`Pr(>|z|)` < 0.05, "*", "")))
# 
# # Plot with significance indicators
# importance <- ggplot(coefs_df, aes(x = reorder(Variable, AbsEstimate), y = AbsEstimate, fill = Color)) +
#   geom_bar(stat = "identity") +
#   geom_text(aes(label = Significance, hjust = ifelse(AbsEstimate < 0, 1.2, -0.2))) +
#   coord_flip() +
#   scale_fill_manual(values = c("positive" = "steelblue", "negative" = "firebrick")) +
#   theme_minimal() +
#   labs(#title = "Feature Importance in Logistic Regression Model",
#     #subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
#     x = "Features",
#     y = "log(Absolute Coefficient Value)",
#     fill = "Coefficient Direction") +
#   theme(legend.position = "bottom")
# 
# importance
# 
# ggsave(paste0("feat_importance_", site, ".png"), importance, bg = "white", dpi = 300, height = 5, width = 8)
# 
# 
# 
# 
# # Predict probabilities on the test (holdout) set
# preds_prob <- predict(fit_IBD, newdata = df_test_IBD, type = "prob")[, "R"]
# 
# # Define the range of decision_thresholds
# decision_thresholds <- seq(0, 1, by = 0.05)
# 
# # Initialize an empty results dataframe
# results <- data.frame(eCOI_pairs = character(),
#                       decision_threshold = numeric(),
#                       sensitivity = numeric(),
#                       specificity = numeric(),
#                       R_pairs = numeric(),
#                       NI_pairs = numeric(),
#                       stringsAsFactors = FALSE)
# 
# #decision_thresholds <- 0.5 # if wanting to use only 0.5 decision threshold for all pair types...
# 
# # Loop through each decision_threshold
# for (thresh in decision_thresholds) {
#   
#   # Convert probabilities to binary predictions at the current decision_threshold
#   preds <- ifelse(preds_prob >= thresh, "R", "NI")
#   
#   # Loop through each unique eCOI_pairs combination
#   for (strain_comb in unique(TEST_META$eCOI_pairs)) {
#     
#     # Subset the TEST and TEST_labels based on the current combination
#     subset_indices <- TEST_META$eCOI_pairs == strain_comb
#     subset_TEST_labels <- TEST_labels[subset_indices]
#     
#     # Count R and NI pairs
#     r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
#     ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
#     
#     # Get the corresponding predictions for the current subset
#     subset_preds <- preds[subset_indices]
#     
#     # Evaluate the confusion matrix for the current subset
#     cm <- confusionMatrix(as.factor(subset_preds), as.factor(subset_TEST_labels), positive = "R")
#     
#     sens <- cm$byClass["Sensitivity"]
#     spec <- cm$byClass["Specificity"]
#     
#     # Append results
#     results <- rbind(results, data.frame(eCOI_pairs = strain_comb,
#                                          decision_threshold = thresh,
#                                          sensitivity = sens,
#                                          specificity = spec,
#                                          R_pairs = r,
#                                          NI_pairs = ni,
#                                          stringsAsFactors = FALSE))
#   }
# }
# 
# 
# # Select the best decision_threshold per eCOI_pair based on balance between sensitivity and specificity
# best_decision_thresholds <- results %>%
#   mutate(youden_j = sensitivity + specificity - 1) %>%  # Compute Youden’s J
#   group_by(eCOI_pairs) %>%
#   slice_max(youden_j) %>%  # Select rows with the largest Youden's J
#   #slice_min(decision_threshold, with_ties = FALSE) %>%  # If ties, pick the lowest decision_threshold
#   slice_min(abs(decision_threshold - 0.5), with_ties = FALSE) %>%  # Pick decision_threshold closest to 0.5
#   ungroup() 
# 
# 
# print(best_decision_thresholds)
# 
# # remove those that sum 8 (total number of strains)
# best_decision_thresholds <- best_decision_thresholds[!best_decision_thresholds$eCOI_pairs %in% c("3__5", "5__3"), ]
# 
# 
# # Reshape data to long format for easy plotting
# best_decision_thresholds_long <- best_decision_thresholds %>%
#   select(eCOI_pairs, sensitivity, specificity) %>%
#   pivot_longer(cols = c(sensitivity, specificity), 
#                names_to = "Metric", 
#                values_to = "Value")
# 
# # Create bar plot
# metrics <- ggplot(best_decision_thresholds_long, aes(x = eCOI_pairs, y = Value, fill = Metric)) +
#   geom_bar(stat = "identity", position = "dodge") +  # Dodge separates bars for clarity
#   labs(title = "",
#        x = "Pair Type",
#        y = "Value") +
#   theme_minimal() +
#   scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
#   geom_hline(yintercept = 0.9, linetype = "solid", color = "black")+
#   geom_hline(yintercept = 0.80, linetype = "dashed", color = "black")
# 
# metrics
# 
# 
# # save model and model results
# saveRDS(fit_IBD, paste0("LogReg_model_", site, ".RDS"))
# write.csv(best_decision_thresholds, paste0("training_Results_LogReg_", site, ".csv"), row.names = F)
# ggsave(paste0(site, "_model_results_LogReg.png"), metrics, bg = "white", dpi = 300, height = 7, width = 14)
# 
# 
# 
# # Reshape data into long format for plotting
# results_long <- results %>%
#   pivot_longer(cols = c(sensitivity, specificity), 
#                names_to = "Metric", 
#                values_to = "Value")
# 
# # Merge best_decision_thresholds to get the specific decision_threshold for each eCOI_pair
# best_decision_thresholds_long <- best_decision_thresholds %>%
#   select(eCOI_pairs, decision_threshold, sensitivity, specificity) %>%
#   pivot_longer(cols = c(sensitivity, specificity), 
#                names_to = "Metric", 
#                values_to = "Value")
# 
# # Plot both Sensitivity and Specificity in the same graph
# sens_spec_plot <- ggplot(results_long, aes(x = decision_threshold, y = Value, linetype = Metric)) +
#   geom_line() +
#   #geom_point(data = best_decision_thresholds_long, aes(x = decision_threshold, y = log(Value)), shape = 19, size = 3, stroke = 1.5) + 
#   geom_vline(data = best_decision_thresholds, aes(xintercept = decision_threshold), color = "red", linetype = "solid") +
#   facet_wrap(~eCOI_pairs) +
#   labs(title = "",
#        x = "Decision Threshold",
#        y = "Value") +
#   theme_minimal()
# 
# sens_spec_plot


# --------------------------
# 1) SETUP
# --------------------------
ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

df_train_IBD$label <- relevel(df_train_IBD$label, ref = "R")

# --------------------------
# 2) TRAINING WITH GLMNET
# --------------------------
fit_IBD <- train(
  label ~ .,
  data = df_train_IBD,
  method = "glmnet",
  family = "binomial",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 10
)

print(fit_IBD)

# --------------------------
# 3) FEATURE IMPORTANCE (GLMNET)
# --------------------------
best_lambda <- fit_IBD$bestTune$lambda
coefs <- coef(fit_IBD$finalModel, s = best_lambda)
coefs_df <- as.data.frame(as.matrix(coefs))
colnames(coefs_df) <- "Estimate"
coefs_df$Variable <- rownames(coefs_df)

coefs_df <- coefs_df[coefs_df$Variable != "(Intercept)", ]
coefs_df$AbsEstimate <- log(abs(coefs_df$Estimate) + 1e-8)  # avoid log(0)
coefs_df$Color <- ifelse(coefs_df$Estimate > 0, "positive", "negative")

importance <- ggplot(coefs_df, aes(x = reorder(Variable, Estimate), y = Estimate, fill = Color)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("positive" = "steelblue", "negative" = "firebrick")) +
  theme_minimal() +
  labs(x = "Features", y = "Coefficient Value", fill = "Coefficient Direction") +
  theme(legend.position = "bottom")

ggsave(paste0("feat_importance_", site, ".png"), importance, bg = "white", dpi = 300, height = 5, width = 8)

# --------------------------
# 4) PREDICT ON TEST SET
# --------------------------
preds_prob <- predict(fit_IBD, newdata = df_test_IBD, type = "prob")[, "R"]

decision_thresholds <- seq(0, 1, by = 0.05)
results <- data.frame(eCOI_pairs = character(),
                      decision_threshold = numeric(),
                      sensitivity = numeric(),
                      specificity = numeric(),
                      R_pairs = numeric(),
                      NI_pairs = numeric(),
                      stringsAsFactors = FALSE)

for (thresh in decision_thresholds) {
  preds <- ifelse(preds_prob >= thresh, "R", "NI")
  
  for (strain_comb in unique(TEST_META$eCOI_pairs)) {
    subset_indices <- TEST_META$eCOI_pairs == strain_comb
    subset_TEST_labels <- TEST_labels[subset_indices]
    
    r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
    ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
    
    subset_preds <- preds[subset_indices]
    cm <- confusionMatrix(as.factor(subset_preds), as.factor(subset_TEST_labels), positive = "R")
    
    sens <- cm$byClass["Sensitivity"]
    spec <- cm$byClass["Specificity"]
    
    results <- rbind(results, data.frame(eCOI_pairs = strain_comb,
                                         decision_threshold = thresh,
                                         sensitivity = sens,
                                         specificity = spec,
                                         R_pairs = r,
                                         NI_pairs = ni,
                                         stringsAsFactors = FALSE))
  }
}

# --------------------------
# 5) BEST THRESHOLDS (YOUDEN'S J)
# --------------------------
best_decision_thresholds <- results %>%
  mutate(youden_j = sensitivity + specificity - 1) %>%
  group_by(eCOI_pairs) %>%
  slice_max(youden_j) %>%
  slice_min(abs(decision_threshold - 0.5), with_ties = FALSE) %>%
  ungroup()

# remove those that sum 8 (total number of strains)
best_decision_thresholds <- best_decision_thresholds[!best_decision_thresholds$eCOI_pairs %in% c("3__5", "5__3"), ]

best_decision_thresholds_long <- best_decision_thresholds %>%
  select(eCOI_pairs, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity),
               names_to = "Metric",
               values_to = "Value")

metrics <- ggplot(best_decision_thresholds_long, aes(x = eCOI_pairs, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Pair Type", y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0.9, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "black")

# --------------------------
# 6) SAVE RESULTS
# --------------------------
saveRDS(fit_IBD, paste0("GLMNET_model_", site, ".RDS"))
write.csv(best_decision_thresholds, paste0("training_Results_GLMNET_", site, ".csv"), row.names = F)
ggsave(paste0(site, "_model_results_GLMNET.png"), metrics, bg = "white", dpi = 300, height = 5, width = 7)



# Reshape data into long format for plotting
results_long <- results %>%
  pivot_longer(cols = c(sensitivity, specificity), 
               names_to = "Metric", 
               values_to = "Value")

# Merge best_decision_thresholds to get the specific decision_threshold for each eCOI_pair
best_decision_thresholds_long <- best_decision_thresholds %>%
  select(eCOI_pairs, decision_threshold, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity), 
               names_to = "Metric", 
               values_to = "Value")

# Plot both Sensitivity and Specificity in the same graph
sens_spec_plot <- ggplot(results_long, aes(x = decision_threshold, y = Value, linetype = Metric)) +
  geom_line() +
  geom_vline(data = best_decision_thresholds, aes(xintercept = decision_threshold), color = "red", linetype = "solid") +
  facet_wrap(~eCOI_pairs) +
  labs(title = "",
       x = "Decision Threshold",
       y = "Value") +
  theme_minimal()



###### 5) BENCHMAKR AGAINST A DUMMY RANDOM CLASSIFIER ----------------

set.seed(420) 

# Create a dummy classifier that randomly assigns labels based on class distribution
class_probs <- prop.table(table(TRAIN_labels))
TEST$dummy_random <- sample(names(class_probs), size = nrow(TEST), replace = TRUE, prob = class_probs)

# Evaluate performance of both dummy classifiers per eCOI_pair
dummy_results <- data.frame(eCOI_pairs = character(),
                            model = character(),
                            sensitivity = numeric(),
                            specificity = numeric())

for (strain_comb in unique(TEST_META$eCOI_pairs)) {
  
  subset_indices <- TEST_META$eCOI_pairs == strain_comb
  subset_TEST_labels <- TEST_labels[subset_indices]
  
  preds <- TEST$dummy_random[subset_indices]
  
  cm <- confusionMatrix(as.factor(preds), as.factor(subset_TEST_labels), positive = "R")
  
  dummy_results <- rbind(dummy_results, data.frame(
    eCOI_pairs = strain_comb,
    model = "dummy",
    sensitivity = cm$byClass["Sensitivity"],
    specificity = cm$byClass["Specificity"]
  ))
}

# Reshape data to long format for easy plotting
dummy_results_long <- dummy_results %>%
  select(eCOI_pairs, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity), 
               names_to = "Metric", 
               values_to = "Value")

dummy_comparison <- ggplot(dummy_results_long, aes(x = eCOI_pairs, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +  # Dodge separates bars for clarity
  labs(title = "",
       x = "Pair Type",
       y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept = 0.9, linetype = "solid", color = "black")+
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "black")+
  ylim(0,1)

dummy_comparison

ggsave(paste0(site, "_sensitivity_dummy_model_comparison_LR.png"), dummy_comparison, bg = "white", dpi = 300, height = 7, width = 14)



#####################################################################################
############################# LAB_CONTROLS ##############################
#####################################################################################


CONTROLS_ALL <- read.csv("LAB_CONTROLS_CLEAN.csv")


#####################################################################################3
### 2. CREATE CONTROL METADATA (PAIRS) ----

sampleIDs <- unique(CONTROLS_ALL$sampleID)

# Create all possible pairs of sampleIDs
pairs_df <- expand.grid(NIDA1 = sampleIDs, NIDA2 = sampleIDs, stringsAsFactors = FALSE)

# #just remove pairs of the same sample
# pairs_df <- pairs_df %>%
#   filter(NIDA1 != NIDA2)


# Create a unique PairsID for each pair
pairs_df <- pairs_df %>%
  mutate(PairsID = row_number()) %>%
  pivot_longer(cols = c(NIDA1, NIDA2), names_to = "time_point", values_to = "NIDA")

# Assign time_point based on the column names
pairs_df <- pairs_df %>%
  mutate(time_point = ifelse(time_point == "NIDA1", "D0", "Dx"))

# Rearrange columns to match the desired format
control_pairs <- pairs_df %>%
  select(PairsID, NIDA, time_point)

# add metadata to the same df

# Step 1: Extract 'run' from NIDA column
control_pairs <- control_pairs %>%
  mutate(run = str_extract(NIDA, "(?<=___).*"))



#####################################################################################3
####### 3. FORMATTING #######------------------

#paired_data import and formatting
# control_pairs$NIDA <- gsub("\\.", "_", control_pairs$NIDA) #change "." for "_" on NIDA column
# control_pairs$NIDA <- ifelse(startsWith(control_pairs$NIDA, "N"), 
#                              control_pairs$NIDA, 
#                              paste0("N", control_pairs$NIDA)) #add initial N to all NIDAs
# 
# 
# Remove "_S*" from sampleID only if "_S" is not present in control_pairs$NIDA
# contains_S <- any(grepl("_S", control_pairs$NIDA))
# if (!contains_S) {
#   CONTROLS_ALL$sampleID <- sub("_S.*$", "", CONTROLS_ALL$sampleID)
# }
colnames(CONTROLS_ALL)[1] <- "NIDA"
# CONTROLS_ALL$NIDA <- gsub("\\.", "_", CONTROLS_ALL$NIDA) #change "." for "_" on NIDA column
# CONTROLS_ALL$NIDA <- ifelse(startsWith(CONTROLS_ALL$NIDA, "N"), 
#                             CONTROLS_ALL$NIDA, 
#                             paste0("N", CONTROLS_ALL$NIDA)) #add initial N to all NIDAs

#are all the NIDAs in paired_data present in CONTROLS_ALL?
control_pairs_NIDA <- unique(control_pairs$NIDA)
CONTROLS_ALL_NIDA <- unique(CONTROLS_ALL$NIDA)
missing_NIDAs <- control_pairs_NIDA[!control_pairs_NIDA %in% CONTROLS_ALL_NIDA]

if (length(missing_NIDAs) > 0) {
  print("These NIDAs from control_pairs are not present in CONTROLS_ALL:")
  print(missing_NIDAs)
  
  #remove paired samples with missing NIDAs
  missing_NIDAs_PairsID<-subset(control_pairs, NIDA %in% missing_NIDAs)[2]
  control_pairs <- control_pairs[!(control_pairs$PairsID %in% missing_NIDAs_PairsID$PairsID), ]
  
  
} else {
  print("All NIDAs from control_pairs are present in CONTROLS_ALL.")
}

# ### remove asv, pseudo_cigar and run columns because ugh
# if ("run" %in% colnames(CONTROLS_ALL)) {
#   # If present, select all columns except 'asv', 'pseudo_cigar', and 'run'
#   CONTROLS_ALL <- CONTROLS_ALL %>%
#     select(-pseudo_cigar, -run, -n.alleles)
# } else {
#   # If not present, keep only the 'run' column (if it exists)
#   CONTROLS_ALL <- CONTROLS_ALL %>%
#     select(-pseudo_cigar, -n.alleles)
# }
# 
# if ("run" %in% colnames(control_pairs)) {
#   # If present, remove 'run' column
#   control_pairs <- control_pairs %>%
#     select(-run)
# }

CONTROLS_ALL <- CONTROLS_ALL %>% select(-run)

#merge control_pairs and CONTROLS_ALL by NIDA
merged_dfs <-  inner_join(control_pairs, CONTROLS_ALL, by = "NIDA")

### remove asv column because ugh
merged_dfs <- merged_dfs %>%
  arrange(PairsID, time_point)


#check
all(unique(control_pairs$NIDA) %in% unique(merged_dfs$NIDA))
all(unique(control_pairs$PairsID) %in% unique(merged_dfs$PairsID))

all(unique(merged_dfs$NIDA) %in% unique(control_pairs$NIDA))
all(unique(merged_dfs$PairsID) %in% unique(control_pairs$PairsID))


# #### RAM DUMP
# rm(controls_ucsf, controls_ucsf_kathryn, data_all, data)
# gc()

merged_dfs <- merged_dfs[merged_dfs$PairsID %in% control_pairs$PairsID,]

#check
all(unique(control_pairs$NIDA) %in% unique(merged_dfs$NIDA))
all(unique(control_pairs$PairsID) %in% unique(merged_dfs$PairsID))

all(unique(merged_dfs$NIDA) %in% unique(control_pairs$NIDA))
all(unique(merged_dfs$PairsID) %in% unique(control_pairs$PairsID))





## format metadata. add labels. add nstrains for each time

# 1) import control metadata
meta <- read.csv("control_metadata.csv")

#optional: remove low density controls
# above_100_paras <- na.omit(meta[meta$parasitemia > 100,]$NIDA)
# 
# meta <- meta[meta$NIDA %in% above_100_paras,]
# 
# couts <- table(control_pairs[control_pairs$NIDA %in% above_100_paras, ]$PairsID)
# pairs_above_100_paras <- names(couts[couts == 2])
# 
# merged_dfs <- merged_dfs[merged_dfs$PairsID %in% pairs_above_100_paras,]
# control_pairs <- control_pairs[control_pairs$PairsID %in% pairs_above_100_paras,]

# 2) count strains for each mix
library(stringr)

# Step 1: Convert strain columns to character
meta <- meta %>%
  mutate(across(matches("^strain[1-9]$"), ~ as.character(.)))

# Step 2: Replace any NA with ""
meta <- meta %>%
  mutate(across(matches("^strain[1-9]$"), ~ replace_na(., "")))

# Step 3: Recompute n_strains
meta <- meta %>%
  rowwise() %>%
  mutate(
    n_strains = sum(c_across(matches("^strain[1-9]$")) != "")
  ) %>%
  ungroup()

# 3) merge with control_pairs to get the n_strains
control_pairs <- left_join(control_pairs, meta[c("n_strains", "NIDA")], by = "NIDA")

# 4) reformat
library(tidyr)

# Pivot the data wider by time_point (D0 vs Dx)
formatted_pairs <- control_pairs %>%
  mutate(row = row_number()) %>%
  pivot_wider(
    id_cols = PairsID,
    names_from = time_point,
    values_from = c(NIDA, n_strains),
    names_sep = "_"
  ) %>%
  # Rename columns as requested
  rename(
    D0_nstrains = n_strains_D0,
    Dx_nstrains = n_strains_Dx,
    D0_sample = NIDA_D0,
    Dx_sample = NIDA_Dx) %>%
  select(PairsID,  D0_sample, Dx_sample, D0_nstrains, Dx_nstrains)

formatted_pairs$NIDA1 <- formatted_pairs$D0_sample
formatted_pairs$NIDA2 <- formatted_pairs$Dx_sample


## 4) assign labels
library(purrr)

# Step 1: Get all unique NIDA values
nida_list <- unique(meta$NIDA)

# Step 2: Create all pairwise combinations (including both orders if needed)
nida_pairs <- expand.grid(NIDA1 = nida_list, NIDA2 = nida_list) %>%
  filter(NIDA1 != NIDA2)

# # Optional: Remove duplicate reversed pairs
# nida_pairs <- nida_pairs %>%
#   rowwise() %>%
#   mutate(pair_id = paste(sort(c(NIDA1, NIDA2)), collapse = "_")) %>%
#   distinct(pair_id, .keep_all = TRUE) %>%
#   select(-pair_id) %>%
#   ungroup()

# Step 3: Create a helper function to get strains from a given NIDA
get_strains <- function(nida_name) {
  meta %>%
    filter(NIDA == nida_name) %>%
    select(starts_with("strain")) %>%
    unlist(use.names = FALSE) %>%
    na.omit() %>%
    trimws() %>%
    .[. != ""]
}

# Step 4: Compare strain sets for each pair
compare_strains <- nida_pairs %>%
  rowwise() %>%
  mutate(
    shared = length(intersect(get_strains(NIDA1), get_strains(NIDA2))) > 0,
    labels = ifelse(shared, "R", "NI")
  ) %>%
  ungroup() %>%
  select(NIDA1, NIDA2, labels)

# View result
print(compare_strains)

### 5) merge all to create final metadata
final_metadata <- formatted_pairs %>%
  left_join(compare_strains, by = c("NIDA1", "NIDA2"))

final_metadata$labels <- ifelse(final_metadata$NIDA1 == final_metadata$NIDA2, "R", final_metadata$labels)

# 6) final touches
# reoving some pairs that don't appear in the metadata. this is becauase metadata was originally done manually and i added more runs later. all is good. no need to add new runs anyways
final_metadata <- final_metadata[!is.na(final_metadata$labels),]

# subset genomic
merged_dfs <- merged_dfs[merged_dfs$PairsID %in% final_metadata$PairsID,]

# 1. Check PairsID matches both ways
all(merged_dfs$PairsID %in% final_metadata$PairsID)     # A in B
all(final_metadata$PairsID %in% merged_dfs$PairsID)     # B in A

# 2. Check NIDA matches both ways (against NIDA1 and NIDA2)
all(merged_dfs$NIDA %in% c(final_metadata$NIDA1, final_metadata$NIDA2))  # A in B
all(c(final_metadata$NIDA1, final_metadata$NIDA2) %in% merged_dfs$NIDA)  # B in A


## check:
nrow(final_metadata)
length(unique(merged_dfs$PairsID))

################################## export

saveRDS(merged_dfs,  "Lab_Controls_genomic.RDS")
saveRDS(final_metadata,  "Lab_Controls_metadata.RDS")


#################################################################################### AQUÍ VOY!!! 11/JUN/2025


library(data.table)
library(dplyr)
library(dcifer)
library(reshape2)
library(tidyr)
library(purrr)
library(stringr)
library(Matrix)
library(progress)
library(parallel)


# for the training data:
PAIRS_METADATA <- readRDS(paste0("Lab_Controls_metadata.RDS"))
PAIRS_GENOMIC <- readRDS(paste0("Lab_Controls_genomic.RDS"))

PAIRS_GENOMIC <- as.data.table(PAIRS_GENOMIC)
PAIRS_METADATA <- as.data.table(PAIRS_METADATA)

PAIRS_GENOMIC <- left_join(PAIRS_GENOMIC, PAIRS_METADATA, by = "PairsID") 



calculate_features_optimized <- function(sample1, sample2) {
  # 1) Unique alleles & loci
  alleles1 <- unique(sample1$allele)
  alleles2 <- unique(sample2$allele)
  all_alleles <- union(alleles1, alleles2)
  
  # 2) Allele-level intersection & union via set operations
  inter_cnt   <- length(intersect(alleles1, alleles2))
  union_cnt   <- length(all_alleles)
  jaccard     <- if (union_cnt > 0) inter_cnt / union_cnt else 0
  retention   <- if (length(alleles1) > 0) inter_cnt / length(alleles1) else 0
  allele_gain <- if (union_cnt > 0) length(setdiff(alleles2, alleles1)) / union_cnt else 0
  allele_loss <- if (union_cnt > 0) length(setdiff(alleles1, alleles2)) / union_cnt else 0
  
  # 3) Transition asymmetry
  trans_asym <- if ((allele_gain + allele_loss) > 0)
    (allele_gain - allele_loss) / (allele_gain + allele_loss) else 0
  
  # 4) Prepare locus-grouped allele lists
  split1 <- split(sample1$allele, sample1$locus)
  split2 <- split(sample2$allele, sample2$locus)
  loci   <- union(names(split1), names(split2))
  n_loci <- length(loci)
  
  # 5) Compute discordant loci count
  discordant_loci <- sum(vapply(
    loci,
    function(l) {
      length(intersect(split1[[l]] %||% character(0),
                       split2[[l]] %||% character(0))) == 0
    },
    logical(1)
  ))
  
  # 6) Compute replacement pattern sum
  replacement_pattern_sum <- sum(vapply(
    loci,
    function(l) {
      a1 <- split1[[l]] %||% character(0)
      a2 <- split2[[l]] %||% character(0)
      if      (length(a2) == 0)         1
      else if (all(a2 %in% a1))         0
      else                              length(setdiff(a2, a1)) / length(a2)
    },
    numeric(1)
  ))
  
  # 7) Final locus-level rates
  locus_discordance_rate    <- discordant_loci / n_loci
  replacement_pattern_score <- replacement_pattern_sum / n_loci
  
  # 8) Return all features
  c(
    jaccard_similarity          = jaccard,
    allele_retention_rate       = retention,
    allele_gain                 = allele_gain,
    allele_loss                 = allele_loss,
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

gc()

# Parallel processing setup
num_cores <- detectCores() - 0
chunk_size <- 1000  # Adjust based on available memory

# Initialize the progress bar
pb <- progress_bar$new(
  format = "Processing [:bar] :percent :elapsedfull",
  total = ceiling(length(unique_pairs) / chunk_size),
  width = 100
)

# Temporary file for writing output
output_file <- paste0("delta_features_Lab_Controls_TRAINING_DATA.csv")

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

delta_metrics_df_final <- read.csv(paste0("delta_features_Lab_Controls_TRAINING_DATA.csv"))

delta_metrics_df_final <- delta_metrics_df_final %>%
  select(PairsID, everything())

write.csv(delta_metrics_df_final, paste0("delta_features_Lab_Controls_TRAINING_DATA.csv"), row.names = F)



##### MERGE FEATURES AND OUTPUT ------

# FEATURES <- left_join(dres0_long_final_summarized, delta_metrics_df_final, by = "PairsID")

FEATURES <- left_join(delta_metrics_df_final, PAIRS_METADATA, by = "PairsID")

# ECOI CHANGE
FEATURES$coi_change <- FEATURES$Dx_nstrains - FEATURES$D0_nstrains

# some NAs for whatever reason...
FEATURES <- FEATURES[!is.na(FEATURES$coi_change),]


write.csv(FEATURES, paste0("Lab_Controls_TRAINING_DATA.csv"), row.names = F)



####################### EDA @@@@@@@@@@@@@@@@@@

library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(uwot)
library(ggpubr)
library(rstatix)
library(cluster)

#----------------------------------------------------------
# Data Import & Preparation
#----------------------------------------------------------
library(corrplot)

TRAINING_DATA <- read.csv("Lab_Controls_TRAINING_DATA.csv", row.names = 1)
TRAINING_DATA$eCOI_pairs <- paste0(TRAINING_DATA$D0_nstrains, "__", TRAINING_DATA$Dx_nstrains)
LABELS <- data.frame(labels = TRAINING_DATA$labels)
LABELS$labels <- as.factor(LABELS$labels)

#REAL_DATA <- read.csv(paste0(site, "_REAL_DATA.csv"), stringsAsFactors = FALSE, colClasses = c(NIDA1 = "character", NIDA2= "character")) 

features_to_use <- colnames(TRAINING_DATA)[!colnames(TRAINING_DATA) %in% c("PairsID", "NIDA1", "NIDA2", "pair_type", "IBD_estimate",
                                                                           "offset_naive_coi_D0", "offset_naive_coi_Dx", 
                                                                           "replacement_pattern_score", "locus_discordance_rate", 
                                                                           "D0_sample", "Dx_sample","D0_nstrains", "Dx_nstrains", "eCOI_pairs", "labels", "shared_count","union_count", "shared_prop")]

#----------------------------------------------------------
# 1) FEATURE CORRELATIONS
#----------------------------------------------------------
corr_data <- TRAINING_DATA %>% select(all_of(features_to_use))
corrplot(cor(corr_data, use = "complete.obs"), method = "pie", tl.cex = 0.7)

#----------------------------------------------------------
# 2) FEATURE DISTRIBUTIONS: GLOBAL & STRATIFIED BY eCOI_pairs
#----------------------------------------------------------
# Convert to long format
training_long <- TRAINING_DATA %>%
  select(all_of(features_to_use), eCOI_pairs, labels) %>%
  pivot_longer(cols = all_of(features_to_use), names_to = "feature", values_to = "value")

# Global distributions (density plots)
ggplot(training_long, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  facet_wrap(~ feature, scales = "free") +
  theme_minimal() +
  labs(title = "", x = "Value", y = "Density")

# Stratified by eCOI_pairs
distros_strat <- ggplot(training_long, aes(x = value, fill = labels, color = labels)) +
  geom_density(alpha = 0.5) +
  facet_grid(feature~ eCOI_pairs, scales = "free") +
  theme_minimal() +
  labs(title = "", x = "Value", y = "Density")+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8) )

ggsave(paste0("distributions_stratified_Lab_Controls.png"), distros_strat, height = 12, width = 14, dpi = 300, bg= "white")

#----------------------------------------------------------
# 3) BOXPLOTS OF FEATURES: GLOBAL & STRATIFIED BY eCOI_pairs WITH WILCOX TEST P-VALUES
#----------------------------------------------------------
# Boxplots by label (global)
p1 <- ggplot(training_long, aes(x = labels, y = value, fill = labels)) +
  geom_boxplot(alpha = 1, color = "black", outlier.alpha = 0.05) +  # hide outliers; outlier.shape = NA
  facet_grid(feature ~ eCOI_pairs, scales = "free") +
  # stat_compare_means(method = "wilcox.test", label = "p.signif", label.y.npc = "top") +
  #scale_fill_brewer(palette = "Set3") +
  theme_minimal(base_size = 14) +
  labs(title = "", 
       x = "Label", 
       y = "Value") +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )

ggsave(paste0("boxplots_stratified_Lab_Controls.png"), p1, height = 12, width = 17, dpi = 300, bg= "white")


# p2 <- ggplot(training_long, aes(x = eCOI_pairs, y = value, fill = eCOI_pairs)) +
#   geom_boxplot(alpha = 0.8, color = "black", outlier.alpha = 0.05) +  # hide outliers; outlier.shape = NA
#   facet_wrap(~ feature, scales = "free") +
#   #stat_compare_means(method = "wilcox.test", label = "p.signif", label.y.npc = "top") +
#   scale_fill_brewer(palette = "Set3") +
#   theme_minimal(base_size = 14) +
#   labs(title = "", 
#        x = "Pair Type", 
#        y = "Value") +
#   theme(legend.position = "none",
#         strip.background = element_rect(fill = "lightgrey", color = "black"),
#         strip.text = element_text(face = "bold"),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.border = element_rect(color = "black", fill = NA, size = 0.8) )
# 
# ggsave(paste0("boxplots_global_", site, ".png"), p2, height = 8, width = 13, dpi = 300, bg= "white")



#----------------------------------------------------------
# 4) UMAP OF FEATURES: COLORED BY LABEL, SHAPED BY eCOI_pairs
#----------------------------------------------------------
# Run UMAP using uwot
umap_res <- uwot::umap(TRAINING_DATA %>% select(all_of(features_to_use)),
                       n_neighbors = 15,
                       seed = 420, 
                       n_threads = 20, 
                       n_components = 3,
                       verbose = T)

# Convert the result into a data frame and rename columns
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2", "UMAP3")

# Add label and eCOI_pairs information
umap_df$label <- LABELS$labels
umap_df$eCOI_pairs <- TRAINING_DATA$eCOI_pairs

# Plot UMAP embedding colored by label and shaped by eCOI_pairs
umap_all <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, size = UMAP3, color = label)) +
  geom_point(size = 3, alpha = 0.3) +
  theme_minimal() +
  labs(title = "", x = "UMAP1", y = "UMAP2", 
       color = "Label", shape = "eCOI_pairs")

#umap_all

ggsave(paste0("umap_global_Lab_Controls.png"), umap_all, height = 7, width = 8, dpi = 300, bg= "white")




### PREDICT WITH THE MODEL CREATED FROM TRUTH STRAINS!

model <- readRDS("GLMNET_model_TRUTH_CLEAN.RDS")

thresholds <- read.csv("training_Results_GLMNET_TRUTH_CLEAN.csv")
thresholds <- thresholds %>% rename(pair_type = eCOI_pairs)


# metatada
PAIRS_METADATA <- readRDS("Lab_Controls_metadata.RDS")

# Input features testing data for subsetting
test_data <- read.csv("Lab_Controls_TRAINING_DATA.csv")
test_data$pair_type <- paste0(test_data$D0_nstrains, "__", test_data$Dx_nstrains)
TEST_META <- test_data %>%
  mutate(PairsID = as.character(PairsID),
         labels = factor(labels, levels = c("NI", "R"))) %>%
  select(PairsID, labels, pair_type) %>%  # Renamed to pair_type for consistency
  as.data.table()

# subset metadata to only account for test data
PAIRS_METADATA <- PAIRS_METADATA[PAIRS_METADATA$PairsID %in% TEST_META$PairsID,]

# # add parasitemia
# left_join(PAIRS_METADATA, meta, by = "")

TEST_META_threshs <- left_join(TEST_META, thresholds[c("pair_type", "decision_threshold")], by = "pair_type")
TEST_META_threshs$preds_prob <- predict(model, newdata = test_data, type = "prob")[, "R"]
TEST_META_threshs$preds <- ifelse(TEST_META_threshs$preds_prob >= TEST_META_threshs$decision_threshold, "R", "NI")


# Function to calculate performance metrics
calculate_metrics <- function(TEST_META_threshs) {
  results <- data.frame(pair_type = character(),
                        sensitivity = numeric(),
                        specificity = numeric(),
                        R_pairs = numeric(),
                        NI_pairs = numeric(),
                        stringsAsFactors = FALSE)
  
  for (strain_combo in unique(TEST_META_threshs$pair_type)) {
    subset_indices <- TEST_META_threshs$pair_type == strain_combo
    subset_TEST_labels <- TEST_META_threshs$labels[subset_indices]
    subset_preds <- TEST_META_threshs$preds[subset_indices]
    
    r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
    ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
    
    # Only calculate if both prediction and label vectors are non-empty
    if (length(subset_preds) > 0 && length(subset_TEST_labels) > 0) {
      cm <- caret::confusionMatrix(
        factor(subset_preds, levels = c("R", "NI")),
        factor(subset_TEST_labels, levels = c("R", "NI")),
        positive = "R"
      )
      
      sens <- cm$byClass["Sensitivity"]
      spec <- cm$byClass["Specificity"]
    } else {
      sens <- NA
      spec <- NA
    }
    
    results <- rbind(results, data.frame(pair_type = strain_combo,
                                         sensitivity = sens,
                                         specificity = spec,
                                         R_pairs = r,
                                         NI_pairs = ni,
                                         stringsAsFactors = FALSE))
  }
  
  results <- left_join(results, thresholds[c("pair_type", "decision_threshold")], by = "pair_type")
  return(results)
}


# Calculate metrics
results <- calculate_metrics(TEST_META_threshs)

results <- results[complete.cases(results),]


# Reshape data to long format for easy plotting
best_decision_thresholds_long <- results %>%
  select(pair_type, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity), 
               names_to = "Metric", 
               values_to = "Value")

# Create bar plot
metrics <- ggplot(best_decision_thresholds_long, aes(x = pair_type, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +  # Dodge separates bars for clarity
  labs(title = "",
       x = "Pair Type",
       y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept = 0.9, linetype = "solid", color = "black")+
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "black")

metrics

ggsave(paste0("TRUTH_model_vs_Lab_Controls_GLMNET.png"), metrics, bg = "white", dpi = 300, height = 5, width = 7)


# TO DO:
## TRIPLE CHECK BOTH THE LAB CONTROLS AND THE TRUTH DATA. MAYBE THERE ARE SOME SAMPLES/MIXES IN THE LAB CONTROLS THAT MAY BE BETTER OFF NOT INCLUDING DUE TO UNCERTAIN COMPOSITION
## IS THE AMOUNT OF PAIRS TOO LOW FOR THE TRUTH DATA?
## MAYBE DON'T SPLIT THE TRUTH DATA THIS TIME?? IS NOT NEEDED GIVEN THAT THE TSTING WILL BE DONE WOTH THE LAB CONTROLS...
# REMOVE MIXES WITH DD2k? AND W2? APPARENTLY BOTH ARE THE SAME AS DD2. THEY ARE TAGGED AS DD2, BUT MAYBE THERE COULD BE SOME TAGGING ISSUES IN SOME LABS...
# MAYBE RESTRICT TO >0.98 MIXES??? AT THE MOMENT IS THOUGHT AS IF IT WERE >0.99



### model using only lab data, no truth.


############# SCRIPT 7 MODED FOR TRUTH STRAINS @@@@@@@@@@@@

library(caret)    
library(dplyr)    
library(tidyr)    
library(ggplot2)  
library(broom)


site <- "Lab_Controls"


### 1) IMPORT TRAINING AND REAL DATA ----------

TRAINING_DATA <- read.csv(paste0(site, "_TRAINING_DATA.csv"), row.names = 1) 
TRAINING_DATA$eCOI_pairs <- paste0(TRAINING_DATA$D0_nstrains, "__", TRAINING_DATA$Dx_nstrains) # Create a new variable by combining 'D0nstrains' and 'Dxnstrains'
LABELS <- data.frame(labels = TRAINING_DATA$labels)
LABELS$labels <- as.factor(LABELS$labels)

#REAL_DATA <- read.csv(paste0(site, "_REAL_DATA.csv"), stringsAsFactors = FALSE, colClasses = c(NIDA1 = "character", NIDA2= "character")) 

features_to_use <- colnames(TRAINING_DATA)[!colnames(TRAINING_DATA) %in% c("PairsID", "NIDA1", "NIDA2", "pair_type", "IBD_estimate",
                                                                           "offset_naive_coi_D0", "offset_naive_coi_Dx", 
                                                                           "replacement_pattern_score", "locus_discordance_rate", 
                                                                           "D0_sample", "Dx_sample","D0_nstrains", "Dx_nstrains", "eCOI_pairs", "labels", "shared_count","union_count", "shared_prop")]

corrplot::corrplot(cor(TRAINING_DATA %>% select(features_to_use), use = "complete.obs"), "pie")


### 2) SPLIT DATA ------------

set.seed(420)
# Step 1: Stratified Sampling by `pair_type`
train_indices <- createDataPartition(TRAINING_DATA$eCOI_pairs, p = 0.7, list = FALSE)
train_data <- TRAINING_DATA[train_indices, ]
test_data <- TRAINING_DATA[-train_indices, ]


## output TEST for later use in false negative tests
test_data$PairsID <- rownames(test_data)
test_data %>% select(PairsID, everything(), -eCOI_pairs)
saveRDS(test_data, paste0(site, "_test_data.RDS"))

# Step 2: Extract Metadata and Labels
TRAIN_META <- train_data %>% select(-all_of(features_to_use))
TRAIN <- train_data %>% select(all_of(features_to_use))
TRAIN_labels <- LABELS[train_indices, ]

TEST_META <- test_data %>% select(-all_of(features_to_use))
TEST <- test_data %>% select(all_of(features_to_use))
TEST_labels <- LABELS[-train_indices, ]

# Step 3: Verify Stratification
prop.table(table(TRAIN_META$eCOI_pairs))
prop.table(table(TRAIN_labels))
table(TRAIN_labels)

prop.table(table(TEST_META$eCOI_pairs))
prop.table(table(TEST_labels))
table(TEST_labels)



##### 4) TEST MODEL USING IBD ONLY (BEST FEATURE)--------

# Create training and test data frames for the current feature
df_train_IBD <- data.frame(TRAIN[features_to_use], label = as.factor(TRAIN_labels))
df_test_IBD <- data.frame(TEST[features_to_use], label = as.factor(TEST_labels))

# # Set up 10-fold cross-validation
# ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
# 
# # Re-level factor so that "R" is the positive class
# df_train_IBD$label <- relevel(df_train_IBD$label, ref = "R")
# 
# # Train logistic regression model with cross-validation
# fit_IBD <- train(label ~ ., 
#                  data = df_train_IBD, 
#                  method = "glm", 
#                  family = "binomial", 
#                  trControl = ctrl, 
#                  metric = "ROC")
# 
# print(fit_IBD)
# 
# fit_IBD$finalModel
# 
# ### feature importance
# # Extract coefficients from the final model
# coefs <- summary(fit_IBD$finalModel)$coefficients
# coefs_df <- as.data.frame(coefs)
# coefs_df$Variable <- rownames(coefs_df)
# 
# # Remove intercept for feature importance plot
# coefs_df <- coefs_df[coefs_df$Variable != "(Intercept)", ]
# 
# # Calculate absolute coefficient values to rank by importance
# coefs_df$AbsEstimate <- log(abs(coefs_df$Estimate))
# 
# # Sort by absolute coefficient value
# coefs_df <- coefs_df[order(coefs_df$AbsEstimate, decreasing = TRUE), ]
# 
# # Create a color vector (positive coefficients in blue, negative in red)
# coefs_df$Color <- ifelse(coefs_df$Estimate > 0, "positive", "negative")
# 
# 
# # Add significance stars
# coefs_df$Significance <- ifelse(coefs_df$`Pr(>|z|)` < 0.001, "***",
#                                 ifelse(coefs_df$`Pr(>|z|)` < 0.01, "**",
#                                        ifelse(coefs_df$`Pr(>|z|)` < 0.05, "*", "")))
# 
# # Plot with significance indicators
# importance <- ggplot(coefs_df, aes(x = reorder(Variable, AbsEstimate), y = AbsEstimate, fill = Color)) +
#   geom_bar(stat = "identity") +
#   geom_text(aes(label = Significance, hjust = ifelse(AbsEstimate < 0, 1.2, -0.2))) +
#   coord_flip() +
#   scale_fill_manual(values = c("positive" = "steelblue", "negative" = "firebrick")) +
#   theme_minimal() +
#   labs(#title = "Feature Importance in Logistic Regression Model",
#     #subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
#     x = "Features",
#     y = "log(Absolute Coefficient Value)",
#     fill = "Coefficient Direction") +
#   theme(legend.position = "bottom")
# 
# importance
# 
# ggsave(paste0("feat_importance_", site, ".png"), importance, bg = "white", dpi = 300, height = 5, width = 8)
# 
# 
# 
# 
# # Predict probabilities on the test (holdout) set
# preds_prob <- predict(fit_IBD, newdata = df_test_IBD, type = "prob")[, "R"]
# 
# # Define the range of decision_thresholds
# decision_thresholds <- seq(0, 1, by = 0.05)
# 
# # Initialize an empty results dataframe
# results <- data.frame(eCOI_pairs = character(),
#                       decision_threshold = numeric(),
#                       sensitivity = numeric(),
#                       specificity = numeric(),
#                       R_pairs = numeric(),
#                       NI_pairs = numeric(),
#                       stringsAsFactors = FALSE)
# 
# #decision_thresholds <- 0.5 # if wanting to use only 0.5 decision threshold for all pair types...
# 
# # Loop through each decision_threshold
# for (thresh in decision_thresholds) {
#   
#   # Convert probabilities to binary predictions at the current decision_threshold
#   preds <- ifelse(preds_prob >= thresh, "R", "NI")
#   
#   # Loop through each unique eCOI_pairs combination
#   for (strain_comb in unique(TEST_META$eCOI_pairs)) {
#     
#     # Subset the TEST and TEST_labels based on the current combination
#     subset_indices <- TEST_META$eCOI_pairs == strain_comb
#     subset_TEST_labels <- TEST_labels[subset_indices]
#     
#     # Count R and NI pairs
#     r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
#     ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
#     
#     # Get the corresponding predictions for the current subset
#     subset_preds <- preds[subset_indices]
#     
#     # Evaluate the confusion matrix for the current subset
#     cm <- confusionMatrix(as.factor(subset_preds), as.factor(subset_TEST_labels), positive = "R")
#     
#     sens <- cm$byClass["Sensitivity"]
#     spec <- cm$byClass["Specificity"]
#     
#     # Append results
#     results <- rbind(results, data.frame(eCOI_pairs = strain_comb,
#                                          decision_threshold = thresh,
#                                          sensitivity = sens,
#                                          specificity = spec,
#                                          R_pairs = r,
#                                          NI_pairs = ni,
#                                          stringsAsFactors = FALSE))
#   }
# }
# 
# results <- results[complete.cases(results),]
# 
# # Select the best decision_threshold per eCOI_pair based on balance between sensitivity and specificity
# best_decision_thresholds <- results %>%
#   mutate(youden_j = sensitivity + specificity - 1) %>%  # Compute Youden’s J
#   group_by(eCOI_pairs) %>%
#   slice_max(youden_j) %>%  # Select rows with the largest Youden's J
#   #slice_min(decision_threshold, with_ties = FALSE) %>%  # If ties, pick the lowest decision_threshold
#   slice_min(abs(decision_threshold - 0.5), with_ties = FALSE) %>%  # Pick decision_threshold closest to 0.5
#   ungroup() 
# 
# 
# print(best_decision_thresholds)
# 
# 
# # Reshape data to long format for easy plotting
# best_decision_thresholds_long <- best_decision_thresholds %>%
#   select(eCOI_pairs, sensitivity, specificity) %>%
#   pivot_longer(cols = c(sensitivity, specificity), 
#                names_to = "Metric", 
#                values_to = "Value")
# 
# # Create bar plot
# metrics <- ggplot(best_decision_thresholds_long, aes(x = eCOI_pairs, y = Value, fill = Metric)) +
#   geom_bar(stat = "identity", position = "dodge") +  # Dodge separates bars for clarity
#   labs(title = "",
#        x = "Pair Type",
#        y = "Value") +
#   theme_minimal() +
#   scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
#   geom_hline(yintercept = 0.9, linetype = "solid", color = "black")+
#   geom_hline(yintercept = 0.80, linetype = "dashed", color = "black")
# 
# metrics
# 
# ggsave(paste0("Lab_controls_model_vs_Lab_Controls_LogReg.png"), metrics, bg = "white", dpi = 300, height = 7, width = 14)
# 


# --------------------------
# 1) SETUP
# --------------------------
ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

df_train_IBD$label <- relevel(df_train_IBD$label, ref = "R")

# --------------------------
# 2) TRAINING WITH GLMNET
# --------------------------
fit_IBD <- train(
  label ~ .,
  data = df_train_IBD,
  method = "glmnet",
  family = "binomial",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 10
)

print(fit_IBD)

# --------------------------
# 3) FEATURE IMPORTANCE (GLMNET)
# --------------------------
best_lambda <- fit_IBD$bestTune$lambda
coefs <- coef(fit_IBD$finalModel, s = best_lambda)
coefs_df <- as.data.frame(as.matrix(coefs))
colnames(coefs_df) <- "Estimate"
coefs_df$Variable <- rownames(coefs_df)

coefs_df <- coefs_df[coefs_df$Variable != "(Intercept)", ]
coefs_df$AbsEstimate <- log(abs(coefs_df$Estimate) + 1e-8)  # avoid log(0)
coefs_df$Color <- ifelse(coefs_df$Estimate > 0, "positive", "negative")

importance <- ggplot(coefs_df, aes(x = reorder(Variable, Estimate), y = Estimate, fill = Color)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("positive" = "steelblue", "negative" = "firebrick")) +
  theme_minimal() +
  labs(x = "Features", y = "Coefficient Value", fill = "Coefficient Direction") +
  theme(legend.position = "bottom")

ggsave(paste0("feat_importance_", site, ".png"), importance, bg = "white", dpi = 300, height = 5, width = 8)

# --------------------------
# 4) PREDICT ON TEST SET
# --------------------------
preds_prob <- predict(fit_IBD, newdata = df_test_IBD, type = "prob")[, "R"]

decision_thresholds <- seq(0, 1, by = 0.05)
results <- data.frame(eCOI_pairs = character(),
                      decision_threshold = numeric(),
                      sensitivity = numeric(),
                      specificity = numeric(),
                      R_pairs = numeric(),
                      NI_pairs = numeric(),
                      stringsAsFactors = FALSE)

for (thresh in decision_thresholds) {
  preds <- ifelse(preds_prob >= thresh, "R", "NI")
  
  for (strain_comb in unique(TEST_META$eCOI_pairs)) {
    subset_indices <- TEST_META$eCOI_pairs == strain_comb
    subset_TEST_labels <- TEST_labels[subset_indices]
    
    r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
    ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
    
    subset_preds <- preds[subset_indices]
    cm <- confusionMatrix(as.factor(subset_preds), as.factor(subset_TEST_labels), positive = "R")
    
    sens <- cm$byClass["Sensitivity"]
    spec <- cm$byClass["Specificity"]
    
    results <- rbind(results, data.frame(eCOI_pairs = strain_comb,
                                         decision_threshold = thresh,
                                         sensitivity = sens,
                                         specificity = spec,
                                         R_pairs = r,
                                         NI_pairs = ni,
                                         stringsAsFactors = FALSE))
  }
}

# --------------------------
# 5) BEST THRESHOLDS (YOUDEN'S J)
# --------------------------
best_decision_thresholds <- results %>%
  mutate(youden_j = sensitivity + specificity - 1) %>%
  group_by(eCOI_pairs) %>%
  slice_max(youden_j) %>%
  slice_min(abs(decision_threshold - 0.5), with_ties = FALSE) %>%
  ungroup()

best_decision_thresholds_long <- best_decision_thresholds %>%
  select(eCOI_pairs, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity),
               names_to = "Metric",
               values_to = "Value")

metrics <- ggplot(best_decision_thresholds_long, aes(x = eCOI_pairs, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Pair Type", y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0.9, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "black")

metrics
