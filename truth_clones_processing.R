

library(dplyr)
library(ggdendro)
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(viridis)


TRUTH <- read.csv("pool1A_truth.csv")


## subset pfPHAST amplicons
madhito_amps <- read.csv("../../madhito_20amps.csv")
TRUTH <- TRUTH[TRUTH$locus %in% madhito_amps$locus,]


# Ignore masking, turn it into ref (.)
TRUTH$pseudo_cigar <- gsub("\\d+\\+[^N]*N", "", TRUTH$pseudo_cigar) # Remove masking
TRUTH$pseudo_cigar <- ifelse(TRUTH$pseudo_cigar == "" | is.na(TRUTH$pseudo_cigar), ".", TRUTH$pseudo_cigar) # If empty, add "." since it was reference


# Aggregate unmasked sequences
TRUTH <- unique(TRUTH)


# Create allele column
TRUTH$allele <- paste0(TRUTH$locus, "__", TRUTH$pseudo_cigar)


# Remove indels
TRUTH <- TRUTH[!grepl("I=", TRUTH$allele), ] # Remove alleles with I (insertion)
TRUTH <- TRUTH[!grepl("D=", TRUTH$allele), ] # Remove alleles with D (deletion)



# loci check
TRUTH %>% group_by(Strain) %>% summarise(length(unique(locus)))

# allele check
TRUTH %>% group_by(Strain) %>% summarise(length(allele))


##################################################################

# SHARED ALLELES COMPARISON

# Step 1: Create the presence/absence matrix
allele_matrix <- TRUTH %>%
  distinct(Strain, allele) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Strain, values_from = present, values_fill = 0)

# Step 2: Extract matrix (remove allele column)
allele_data <- as.data.frame(allele_matrix)
allele_only <- as.matrix(allele_data[ , -1])  # remove 'allele' column

# Step 3: Get strain names from column names
strain_names <- colnames(allele_only)

# Step 4: Compute pairwise shared alleles matrix
shared_matrix <- t(allele_only) %*% allele_only  # dot product

# Step 5: Assign row/column names (based on matrix dimensions)
rownames(shared_matrix) <- strain_names
colnames(shared_matrix) <- strain_names

# Sort strain names by total shared alleles (row sums)
strain_order <- shared_matrix %>%
  rowSums() %>%
  sort(decreasing = TRUE) %>%
  names()

# Convert to long format
shared_df <- as.data.frame(as.table(shared_matrix))
colnames(shared_df) <- c("Strain1", "Strain2", "SharedAlleles")

# Convert strain factors to sorted order
shared_df <- shared_df %>%
  mutate(Strain1 = factor(Strain1, levels = strain_order),
         Strain2 = factor(Strain2, levels = strain_order))

# Normalize shared alleles to similarity [0,1]
similarity_matrix <- shared_matrix / max(shared_matrix)

# Convert similarity to distance
distance_matrix <- 1 - similarity_matrix

# Convert to 'dist' object for clustering
dist_obj <- as.dist(distance_matrix)


# Define a viridis color ramp like ggplot2
viridis_colors <- colorRamp2(
  seq(0, max(shared_matrix), length.out = 100),
  viridis(100)
)

Heatmap(shared_matrix,
        name = "Shared Alleles",
        col = viridis_colors,  # equivalent to scale_fill_viridis_c()
        clustering_distance_rows = dist_obj,
        clustering_distance_columns = dist_obj,
        clustering_method_rows = "average",
        clustering_method_columns = "average",
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_names_side = "left",
        column_names_rot = 45)



### REMOVE GENETICALLY IDENTICAL CLONES

remove_clones <- c("W2", "DD2k")
TRUTH_CLEAN <- TRUTH[!TRUTH$Strain %in% remove_clones,]


##  EXPORT CLEAN DATASET

write.csv(TRUTH_CLEAN, "TRUTH_CLEAN.csv", row.names = F)




## proceed to scripts 4 (mixes)  to 7 (model)



############# SCRIPT 4 MODED FOR TRUTH STRAINS @@@@@@@@@@@@

# --- 3. Create All Mixes ----
min_coi <- 1
max_coi <- 7

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


############################################

#just to not crash the script
TRUTH_CLEAN$reads <- 1

initial_sample_size <- 200



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
train_indices <- createDataPartition(TRAINING_DATA$eCOI_pairs, p = 0.9, list = FALSE)
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

# Set up 10-fold cross-validation
ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)

# Re-level factor so that "R" is the positive class
df_train_IBD$label <- relevel(df_train_IBD$label, ref = "R")

# Train logistic regression model with cross-validation
fit_IBD <- train(label ~ ., 
                 data = df_train_IBD, 
                 method = "glm", 
                 family = "binomial", 
                 trControl = ctrl, 
                 metric = "ROC")

print(fit_IBD)

fit_IBD$finalModel

### feature importance
# Extract coefficients from the final model
coefs <- summary(fit_IBD$finalModel)$coefficients
coefs_df <- as.data.frame(coefs)
coefs_df$Variable <- rownames(coefs_df)

# Remove intercept for feature importance plot
coefs_df <- coefs_df[coefs_df$Variable != "(Intercept)", ]

# Calculate absolute coefficient values to rank by importance
coefs_df$AbsEstimate <- log(abs(coefs_df$Estimate))

# Sort by absolute coefficient value
coefs_df <- coefs_df[order(coefs_df$AbsEstimate, decreasing = TRUE), ]

# Create a color vector (positive coefficients in blue, negative in red)
coefs_df$Color <- ifelse(coefs_df$Estimate > 0, "positive", "negative")


# Add significance stars
coefs_df$Significance <- ifelse(coefs_df$`Pr(>|z|)` < 0.001, "***",
                                ifelse(coefs_df$`Pr(>|z|)` < 0.01, "**",
                                       ifelse(coefs_df$`Pr(>|z|)` < 0.05, "*", "")))

# Plot with significance indicators
importance <- ggplot(coefs_df, aes(x = reorder(Variable, AbsEstimate), y = AbsEstimate, fill = Color)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Significance, hjust = ifelse(AbsEstimate < 0, 1.2, -0.2))) +
  coord_flip() +
  scale_fill_manual(values = c("positive" = "steelblue", "negative" = "firebrick")) +
  theme_minimal() +
  labs(#title = "Feature Importance in Logistic Regression Model",
    #subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
    x = "Features",
    y = "log(Absolute Coefficient Value)",
    fill = "Coefficient Direction") +
  theme(legend.position = "bottom")

importance

ggsave(paste0("feat_importance_", site, ".png"), importance, bg = "white", dpi = 300, height = 5, width = 8)




# Predict probabilities on the test (holdout) set
preds_prob <- predict(fit_IBD, newdata = df_test_IBD, type = "prob")[, "R"]

# Define the range of decision_thresholds
decision_thresholds <- seq(0, 1, by = 0.05)

# Initialize an empty results dataframe
results <- data.frame(eCOI_pairs = character(),
                      decision_threshold = numeric(),
                      sensitivity = numeric(),
                      specificity = numeric(),
                      R_pairs = numeric(),
                      NI_pairs = numeric(),
                      stringsAsFactors = FALSE)

#decision_thresholds <- 0.5 # if wanting to use only 0.5 decision threshold for all pair types...

# Loop through each decision_threshold
for (thresh in decision_thresholds) {
  
  # Convert probabilities to binary predictions at the current decision_threshold
  preds <- ifelse(preds_prob >= thresh, "R", "NI")
  
  # Loop through each unique eCOI_pairs combination
  for (strain_comb in unique(TEST_META$eCOI_pairs)) {
    
    # Subset the TEST and TEST_labels based on the current combination
    subset_indices <- TEST_META$eCOI_pairs == strain_comb
    subset_TEST_labels <- TEST_labels[subset_indices]
    
    # Count R and NI pairs
    r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
    ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
    
    # Get the corresponding predictions for the current subset
    subset_preds <- preds[subset_indices]
    
    # Evaluate the confusion matrix for the current subset
    cm <- confusionMatrix(as.factor(subset_preds), as.factor(subset_TEST_labels), positive = "R")
    
    sens <- cm$byClass["Sensitivity"]
    spec <- cm$byClass["Specificity"]
    
    # Append results
    results <- rbind(results, data.frame(eCOI_pairs = strain_comb,
                                         decision_threshold = thresh,
                                         sensitivity = sens,
                                         specificity = spec,
                                         R_pairs = r,
                                         NI_pairs = ni,
                                         stringsAsFactors = FALSE))
  }
}


# Select the best decision_threshold per eCOI_pair based on balance between sensitivity and specificity
best_decision_thresholds <- results %>%
  mutate(youden_j = sensitivity + specificity - 1) %>%  # Compute Youden’s J
  group_by(eCOI_pairs) %>%
  slice_max(youden_j) %>%  # Select rows with the largest Youden's J
  #slice_min(decision_threshold, with_ties = FALSE) %>%  # If ties, pick the lowest decision_threshold
  slice_min(abs(decision_threshold - 0.5), with_ties = FALSE) %>%  # Pick decision_threshold closest to 0.5
  ungroup() 


print(best_decision_thresholds)


# Reshape data to long format for easy plotting
best_decision_thresholds_long <- best_decision_thresholds %>%
  select(eCOI_pairs, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity), 
               names_to = "Metric", 
               values_to = "Value")

# Create bar plot
metrics <- ggplot(best_decision_thresholds_long, aes(x = eCOI_pairs, y = Value, fill = Metric)) +
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


# save model and model results
saveRDS(fit_IBD, paste0("LogReg_model_", site, ".RDS"))
write.csv(best_decision_thresholds, paste0("training_Results_LogReg_", site, ".csv"), row.names = F)
ggsave(paste0(site, "_model_results_LogReg.png"), metrics, bg = "white", dpi = 300, height = 5, width = 7)



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
  #geom_point(data = best_decision_thresholds_long, aes(x = decision_threshold, y = log(Value)), shape = 19, size = 3, stroke = 1.5) + 
  geom_vline(data = best_decision_thresholds, aes(xintercept = decision_threshold), color = "red", linetype = "solid") +
  facet_wrap(~eCOI_pairs) +
  labs(title = "",
       x = "Decision Threshold",
       y = "Value") +
  theme_minimal()

sens_spec_plot



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

ggsave(paste0(site, "_sensitivity_dummy_model_comparison_LR.png"), dummy_comparison, bg = "white", dpi = 300, height = 5, width = 7)



#### AQUÍ VOY!!! jun 11, 2025
