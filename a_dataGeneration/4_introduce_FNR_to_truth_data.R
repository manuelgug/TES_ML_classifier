
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(tibble)

#false negative rates
FNRs <- read.csv("FNRs.csv")

FNRs <- FNRs %>%
  rowwise() %>%
  mutate(FNR = list({
    ranks <- c_across(starts_with("rank_"))
    fnrs <- c()
    fnrs <- c(fnrs, rep(0.14, sum(ranks >= 0.02 & ranks < 0.03, na.rm = TRUE))) ### FNR of 0.14 for strains below 0.03 (experimentally checked with the lab controls dataset!) == DROP 2 ALLELES
    fnrs <- c(fnrs, rep(0.15, sum(ranks < 0.02, na.rm = TRUE))) ### FNR of 0.14 for strains below 0.03 (experimentally checked with the lab controls dataset!) == DROP 3 ALLELES
    fnrs
  })) %>%
  ungroup()

# final values: number of elements in FNR are the strains that should have given FNR for that mix because there were many strains below 0.05 in the real data)
FNRs <- FNRs %>%
  distinct(COI, FNR)

#delete empty vectors (NO FNR APPLIED)
FNRs <- FNRs %>% filter(lengths(FNR) > 0)


# original TRUTH data for the model
PAIRS_METADATA <- readRDS("PAIRS_METADATA_TRUTH_CLEAN.RDS")
PAIRS_GENOMIC <- readRDS("PAIRS_GENOMIC_TRUTH_CLEAN.RDS")

## EXTRACT PAIRS THAT WILL BE FNRed
FNR_pairs <- PAIRS_METADATA[PAIRS_METADATA$D0_nstrains %in% FNRs$COI | PAIRS_METADATA$Dx_nstrains %in% FNRs$COI ,]$PairsID
FNR_pairs_GENOMIC <- PAIRS_GENOMIC[PAIRS_GENOMIC$PairsID %in% FNR_pairs,]
FNR_pairs_METADATA <- PAIRS_METADATA[PAIRS_METADATA$PairsID %in% FNR_pairs,]

FNR_pairs_METADATA <- FNR_pairs_METADATA %>%
  select(PairsID, D0_sample, Dx_sample, D0_nstrains, Dx_nstrains) %>%
  pivot_longer(
    cols = c(D0_sample, Dx_sample, D0_nstrains, Dx_nstrains),
    names_to = c("time_point", ".value"),
    names_sep = "_"
  ) %>%
  rename(sample = sample, n_strains = nstrains)

FNR_pairs_METADATA <- FNR_pairs_METADATA %>%
  mutate(sample = str_split(sample, "_"))

# # add COI
# FNR_pairs_GENOMIC <- FNR_pairs_GENOMIC %>%
#   mutate(COI = as.integer(stringr::str_extract(NIDA, "(?<=mix)\\d+(?=_)")))


## CREATE THE NEW TRUTH + FNR DATASETS
TRUTH <- read.csv("TRUTH_CLEAN.csv")

# Get unique FNR values
unique_fnr_values <- unique(unlist(FNRs$FNR))

set.seed(420)

# Apply FNR to TRUTH
TRUTH_filtered_list <- map(unique_fnr_values, function(fnr) {
  TRUTH %>%
    group_by(Strain) %>%
    group_modify(~ {
      n_to_keep <- ceiling((1 - fnr) * nrow(.x))
      .x[sample(nrow(.x), n_to_keep), , drop = FALSE]
    }) %>%
    ungroup()
})

# Name each element by its FNR
names(TRUTH_filtered_list) <- as.character(unique_fnr_values)

# check
map_dfr(
  names(TRUTH_filtered_list),
  function(fnr_name) {
    TRUTH_filtered_list[[fnr_name]] %>%
      group_by(Strain) %>%
      summarise(n_alleles = n(), .groups = "drop") %>%
      mutate(FNR = fnr_name)
  })


## remake the pairs with the FNR clones/strains

# fill FNR vecgor with TRUTH for easier handling
FNRs <- FNRs %>%
  mutate(
    FNR = map2(FNR, COI, ~ {
      len <- length(.x)
      if (len < .y) {
        c(.x, rep("TRUTH", .y - len))
      } else {
        .x
      }
    })
  )

# FNR_pairs_GENOMIC
# FNR_pairs_METADATA
# 
# TRUTH_filtered_list
# TRUTH

# ADD TRUTH DF TO THE FILTEREDLIST as element 1 FOR EASIER MANAGING
TRUTH_filtered_list <- c(list(as.tibble(TRUTH)), TRUTH_filtered_list)
names(TRUTH_filtered_list)[1] <- "TRUTH"


### TO REMAKE THE PAIRS I NEED TO:
## 1) i'm iterating over FNRs combos
## 2) grab combo 1 
## 3) subset pairids that have the coi
## 4) extract the strains from TRUTH that doesn't have the coi
## 5) extract strains (could be whichever) from the correspondingTRUTH_filtered_list (0.08 or 0.14, depending on the FNR(s))
## 6) put verything together into a df which it's going into a list with the name paste0( FNRs$COI, FNRs$FNR)
## 7) proceed to next FNR combo.
### NOTE!! manage instances with where FNR needs to be applied to many strains have (eg: COI = 5, FNR= c(0.08, 0.14, 0.14) [one strain needs to have 0.08 and two 0.14])
  

# Initialize storage objects
final_sample_df <- tibble()
TRUTH_samples_by_row <- list()

# Iterate over each row in FNRs
for (numba in seq(1, nrow(FNRs), 1)) {
  
  # Get iteration data
  iteration <- FNRs[numba, ]
  
  print(paste0("PROCESSING: COI = ", iteration$COI, "; FNR = ", paste(unlist(iteration$FNR), collapse = "_")))
  
  
  # Subset FNR_pairs_METADATA based on COI
  coi_subset <- FNR_pairs_METADATA %>%
    group_by(PairsID) %>%
    filter(any(n_strains == iteration$COI)) %>%
    ungroup() %>%
    rename(COI = n_strains) %>%
    full_join(iteration, by = "COI")
  
  # Add unique ID to PairsID
  coi_subset$PairsID <- paste0(coi_subset$PairsID, "_FNR", numba)
  
  # Function to process one row
  process_row_optimized <- function(pid, tp, strains, fnrs) {
    # Handle NULL FNR case
    if (is.null(fnrs)) {
      fnrs <- rep("TRUTH", length(strains))
    }
    
    # Create strain-fnr pairs for vectorized processing
    strain_fnr_df <- data.frame(
      strain = strains,
      fnr = fnrs,
      stringsAsFactors = FALSE
    )
    
    # Process all strain-fnr combinations at once
    sampled_data <- strain_fnr_df %>%
      pmap_dfr(function(strain, fnr) {
        TRUTH_filtered_list[[fnr]] %>%
          filter(Strain == strain)
      })
    
    # Create FNR suffix
    fnr_suffix <- paste(fnrs, collapse = "_")
    
    # Return results
    list(
      list_name = paste0(pid, "_", tp, "__", fnr_suffix),
      list_data = sampled_data,
      df_data = tibble(
        PairsID = paste0(pid, "__", fnr_suffix),
        time_point = tp,
        locus = sampled_data$locus,
        pseudo_cigar = sampled_data$pseudo_cigar,
        allele = sampled_data$allele
      )
    )
  }
  
  # Apply processing function to all rows in this iteration
  all_results <- pmap(coi_subset, function(PairsID, time_point, sample, COI, FNR) {
    process_row_optimized(PairsID, time_point, sample, FNR)
  })
  
  # Append TRUTH sample list
  TRUTH_samples_by_row <- c(
    TRUTH_samples_by_row,
    set_names(
      map(all_results, "list_data"),
      map_chr(all_results, "list_name")
    )
  )
  
  # Append to final_sample_df
  final_sample_df <- bind_rows(
    final_sample_df,
    map_dfr(all_results, "df_data")
  )
}

# Deduplicate and enrich final_sample_df
final_sample_df <- final_sample_df %>%
  distinct() %>%
  mutate(NIDA = paste0(PairsID, "__", time_point)) %>%
  rename(FNR_scheme = PairsID) %>%
  mutate(PairsID = sub("__.*", "", FNR_scheme)) %>%
  relocate(PairsID, .before = FNR_scheme)

# Create final_sample_df_METADATA
final_sample_df_METADATA <- final_sample_df %>%
  select(PairsID, time_point, NIDA) %>%
  distinct() %>%
  pivot_wider(
    names_from = time_point,
    values_from = NIDA,
    names_prefix = "NIDA"
  ) %>%
  rename(
    NIDA1 = NIDAD0,
    NIDA2 = NIDADx
  )

  


#########3 CALCULATE FEATURES OF FNR DATA
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

PAIRS_GENOMIC_FNR <- as.data.table(final_sample_df)
PAIRS_METADATA_FNR <- as.data.table(final_sample_df_METADATA)

PAIRS_GENOMIC_FNR <- left_join(PAIRS_GENOMIC_FNR, PAIRS_METADATA_FNR, by = "PairsID") 


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
unique_pairs <- unique(PAIRS_METADATA_FNR$PairsID)

# Create a named vector for mapping PairsID to D0 and Dx for faster access
pairs_to_d0 <- PAIRS_METADATA_FNR$NIDA1[match(unique_pairs, PAIRS_METADATA_FNR$PairsID)]
pairs_to_dx <- PAIRS_METADATA_FNR$NIDA2[match(unique_pairs, PAIRS_METADATA_FNR$PairsID)]

# Preprocess merged_dfs_filtered into a list for fast access
merged_dfs_filtered <- split(PAIRS_GENOMIC_FNR, PAIRS_GENOMIC_FNR$NIDA)
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
output_file <- paste0("TRUTH_FNR_TRAINING_DATA.csv")

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

delta_metrics_df_final <- read.csv(paste0("TRUTH_FNR_TRAINING_DATA.csv"))

delta_metrics_df_final <- delta_metrics_df_final %>%
  select(PairsID, everything())

write.csv(delta_metrics_df_final, paste0("TRUTH_FNR_TRAINING_DATA.csv"), row.names = F)



##### MERGE FEATURES AND OUTPUT ------

# FEATURES <- left_join(dres0_long_final_summarized, delta_metrics_df_final, by = "PairsID")

FEATURES <- left_join(delta_metrics_df_final, PAIRS_METADATA_FNR, by = "PairsID")


### AQUÍ VOY 30 jun 25:
### MERGEAR CON ORIGINAL PAIRS_METADATA BY PAIRSID (split pairsid IN FEATURES!)

FEATURES <- FEATURES %>%
  mutate(PairsID_clean = sub("_.*", "", PairsID))

FEATURES$PairsID_clean <- as.numeric(FEATURES$PairsID_clean)


### MERGE WITH CLEAN DATA'S METADATA
TRUTH_CLEAN<- read.csv("TRUTH_CLEAN_TRAINING_DATA.csv")

# Perform the join using the cleaned ID
merged_df <- FEATURES %>%
  left_join(TRUTH_CLEAN[c("D0_sample", "Dx_sample", "D0_nstrains", "Dx_nstrains", "PairsID", "coi_change", "labels")], by = c("PairsID_clean" = "PairsID")) %>%
  select(-PairsID_clean)

merged_df$PairsID <- as.character(merged_df$PairsID)


#output final features with TRUTH + TRUTH_FNR data
write.csv(merged_df, paste0("TRUTH_FNR_TRAINING_DATA.csv"), row.names = F)


### CONCAT WITH CLEAN DATA
# Find common columns between the two dataframes
common_cols <- intersect(names(merged_df), names(TRUTH_CLEAN))

TRUTH_CLEAN$PairsID <- as.character(TRUTH_CLEAN$PairsID)

# Bind the data using only the common columns
CLEAN_and_FNR_final <- bind_rows(
  merged_df %>% select(all_of(common_cols)),
  TRUTH_CLEAN %>% select(all_of(common_cols))
)

# delete columns with NA
CLEAN_and_FNR_final <- CLEAN_and_FNR_final %>% 
  filter(if_all(everything(), ~ !is.na(.)))


write.csv(CLEAN_and_FNR_final, "TRUTH_CLEAN_and_FNR_FINAL_TRAINING_DATA.csv", row.names = F)

  



### EDA

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

TRAINING_DATA <- read.csv("TRUTH_CLEAN_and_FNR_FINAL_TRAINING_DATA.csv")
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

ggsave(paste0("distributions_stratified_TRUTH_CLEAN_and_FNR.png"), distros_strat, height = 12, width = 14, dpi = 300, bg= "white")

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

ggsave(paste0("boxplots_stratified_TRUTH_CLEAN_and_FNR.png"), p1, height = 12, width = 17, dpi = 300, bg= "white")


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

ggsave(paste0("umap_global_TRUTH_CLEAN_and_FNR.png"), umap_all, height = 7, width = 8, dpi = 300, bg= "white")








### MODEL with FNR
############# SCRIPT 7 MODED FOR TRUTH STRAINS @@@@@@@@@@@@

library(caret)    
library(dplyr)    
library(tidyr)    
library(ggplot2)  
library(broom)


site <- "TRUTH_CLEAN_and_FNR_FINAL"


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
# train_indices <- createDataPartition(TRAINING_DATA$eCOI_pairs, p = 0.7, list = FALSE)

# Create a combined stratification factor
strat_group <- paste0(TRAINING_DATA$eCOI_pairs, "_", ifelse(grepl("FNR", TRAINING_DATA$PairsID), "FNR", "TRUTH"))
train_indices <- createDataPartition(strat_group, p = 0.7, list = FALSE)

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
ggsave(paste0(site, "_model_results_LogReg.png"), metrics, bg = "white", dpi = 300, height = 7, width = 14)



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

ggsave(paste0(site, "_sensitivity_dummy_model_comparison_LR.png"), dummy_comparison, bg = "white", dpi = 300, height = 7, width = 14)





  
  


### PREDICT WITH THE MODEL CREATED FROM TRUTH STRAINS!

model <- readRDS("LogReg_model_TRUTH_CLEAN_and_FNR_FINAL.RDS")

thresholds <- read.csv("training_Results_LogReg_TRUTH_CLEAN_and_FNR_FINAL.csv")
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

ggsave(paste0("TRUTH_CLEAN_and_FNR_model_vs_Lab_Controls_LogReg.png"), metrics, bg = "white", dpi = 300, height = 7, width = 14)


  
  
  
