

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


site <- "Tete"


#select data type betweem "TRAINING_DATA" or "REAL_DATA"
DATA_TYPE = "TRAINING_DATA"


if (DATA_TYPE == "TRAINING_DATA") {
  
  # for the training data:
  PAIRS_METADATA <- readRDS(paste0("PAIRS_METADATA_",site,".RDS"))
  PAIRS_GENOMIC <- readRDS(paste0("PAIRS_GENOMIC_",site,".RDS"))
  
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
  
  # Step 1: Extract naive_coi at D0 and Dx per PairsID
  coi_summary <- PAIRS_GENOMIC[time_point %in% c("D0", "Dx"), 
                               .(naive_coi_D0 = unique(naive_coi[time_point == "D0"]),
                                 naive_coi_Dx = unique(naive_coi[time_point == "Dx"])),
                               by = PairsID]
  
  # Step 2: Create pair_type as "D0__Dx"
  coi_summary[, pair_type := paste0(naive_coi_D0, "__", naive_coi_Dx)]
  
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



# # ####### DCIFER'S IBD #######------------------  
# 
# dsmp <- formatDat(PAIRS_GENOMIC, svar = "NIDA", lvar = "locus", avar = "allele")
# 
# lrank <- 2
# coi   <- getCOI(dsmp, lrank = lrank)
# 
# afreq <- calcAfreq(dsmp, coi, tol = 1e-5)
# 
# dres0 <- ibdDat(dsmp, coi, afreq, pval = TRUE, confint = TRUE, rnull = 0,
#                 alpha = 0.05, nr = 1e3)
# 
# gc()
# 
# suppressWarnings({
#   dres0_long <- melt(dres0)
# })
# dres0_long$value <- ifelse(dres0_long$Var1 == dres0_long$Var2, 1, dres0_long$value) # put 1 if the sample is compared with itself
# 
# #need extra foramtting for real data after passing through dcifer because the nidas...
# if (DATA_TYPE == "REAL_DATA"){
# 
#   dres0_long$Var1 <- as.character(dres0_long$Var1)
#   dres0_long$Var2 <- as.character(dres0_long$Var2)
#   dres0_long <- dres0_long %>%
#     mutate(Var1 = if_else(!str_detect(Var1, "\\."), paste0(Var1, ".0"), Var1),
#            Var2 = if_else(!str_detect(Var2, "\\."), paste0(Var2, ".0"), Var2))  # If no ".", add ".0"
# }
# 
# dres0_long <- dres0_long[dres0_long$Var3 == "estimate" & !is.na(dres0_long$value),]
# dres0_long <- dres0_long %>% select(-Var3)
# colnames(dres0_long) <- c("infection1", "infection2", "IBD_estimate")
# 
# 
# # First match: D0 with infection1 and Dx with infection2
# match1 <- PAIRS_METADATA %>%
#   left_join(dres0_long, by = c("NIDA1" = "infection1", "NIDA2" = "infection2")) %>%
#   select(PairsID, NIDA1, NIDA2, IBD_estimate)
# 
# match1 <- match1[!is.na(match1$IBD_estimate),]
# 
# 
# # Second match: Dx with infection1 and D0 with infection2
# match2 <- PAIRS_METADATA %>%
#   left_join(dres0_long, by = c("NIDA2" = "infection1", "NIDA1" = "infection2")) %>%
#   select(PairsID, NIDA1, NIDA2, IBD_estimate)
# 
# match2 <- match2[!is.na(match2$IBD_estimate),]
# 
# dres0_long_final<- rbind(match1, match2)
# dres0_long_final <- distinct(dres0_long_final)
# 
# dres0_long_final_summarized <- dres0_long_final %>%
#   select(PairsID, IBD_estimate) %>%
#   arrange(PairsID)
# 
# dres0_long_final_summarized <- dres0_long_final_summarized[complete.cases(dres0_long_final_summarized),]
# 
# # MERGE WITH METADATA
# dres0_long_final_summarized <- merge(dres0_long_final_summarized, PAIRS_METADATA, by = "PairsID")

# if (DATA_TYPE == "REAL_DATA"){
# 
#   metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
# 
#   metadata_updated <- metadata_updated[!is.na(metadata_updated$time_point),]
#   
#   metadata_updated$naive_coi <- round(metadata_updated$naive_coi)
# 
#   metadata_updated_wide <- metadata_updated %>%
#     pivot_wider(
#       id_cols = PairsID,
#       names_from = time_point,
#       values_from = c(NIDA, naive_coi),
#       names_glue = "{.value}_{time_point}"
#     )
# 
#   metadata_updated_wide$pair_type <- paste0(metadata_updated_wide$naive_coi_D0, "__", metadata_updated_wide$naive_coi_Dx)
# 
#   dres0_long_final_summarized <- merge(dres0_long_final_summarized, metadata_updated_wide[c("PairsID", "pair_type")], by = "PairsID")
# 
# }


####### DIVERSITY/DELTA FEATURES #######------------------

# # MORE COMPLEX, 6-FEATURE FUNCTION (CURRENT)
calculate_features_optimized <- function(sample1, sample2) {
  # 1) Unique alleles & loci
  alleles1 <- unique(sample1$allele)
  alleles2 <- unique(sample2$allele)
  all_alleles <- union(alleles1, alleles2)
  
  # 2) Allele‐level intersection & union via set operations
  inter_cnt   <- length(intersect(alleles1, alleles2))
  union_cnt   <- length(union(alleles1, alleles2))
  jaccard     <- if (union_cnt>0) inter_cnt/union_cnt else 0
  retention   <- if (length(alleles1)>0) inter_cnt/length(alleles1) else 0
  allele_gain <- if (union_cnt>0) length(setdiff(alleles2, alleles1))/union_cnt else 0
  allele_loss <- if (union_cnt>0) length(setdiff(alleles1, alleles2))/union_cnt else 0
  
  # 3) Transition asymmetry
  trans_asym <- if ((allele_gain+allele_loss)>0)
    (allele_gain - allele_loss)/(allele_gain+allele_loss) else 0
  
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
      if      (length(a2)==0)         1
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


# ### POR SI ACASO... NO BORAR, PUEDE SER ÚTIL CUANDO INTRODUZCA ERRORES:
#
# calculate_features_optimized <- function(sample1, sample2) {
#   # 1) Unique alleles & loci
#   alleles1   <- unique(sample1$allele)
#   alleles2   <- unique(sample2$allele)
#   all_alleles<- union(alleles1, alleles2)
#   
#   # 2) Allele‐level intersection & union
#   inter_cnt  <- length(intersect(alleles1, alleles2))
#   union_cnt  <- length(all_alleles)
#   jaccard    <- if (union_cnt>0) inter_cnt/union_cnt else 0
#   retention  <- if (length(alleles1)>0) inter_cnt/length(alleles1) else 0
#   allele_gain<- if (union_cnt>0) length(setdiff(alleles2, alleles1))/union_cnt else 0
#   allele_loss<- if (union_cnt>0) length(setdiff(alleles1, alleles2))/union_cnt else 0
#   
#   # 3) Transition asymmetry
#   trans_asym <- if ((allele_gain+allele_loss)>0)
#     (allele_gain - allele_loss)/(allele_gain+allele_loss) else 0
#   
#   # 4) Locus‐grouped allele lists
#   split1 <- split(sample1$allele, sample1$locus)
#   split2 <- split(sample2$allele, sample2$locus)
#   loci   <- union(names(split1), names(split2))
#   n_loci <- length(loci)
#   
#   # 5) Per locus: shared count and replacement score
#   shared_counts <- integer(n_loci)
#   share_ratios  <- numeric(n_loci)
#   zero_bits     <- integer(n_loci)
#   replacement_scores <- numeric(n_loci)
#   
#   for (i in seq_along(loci)) {
#     l <- loci[i]
#     a1 <- unique(split1[[l]] %||% character(0))
#     a2 <- unique(split2[[l]] %||% character(0))
#     
#     S  <- length(intersect(a1, a2))
#     U  <- length(union(a1, a2))
#     shared_counts[i] <- S
#     share_ratios[i]  <- if (U>0) S/U else 0
#     zero_bits[i]     <- as.integer(S == 0)
#     
#     # replacement score
#     replacement_scores[i] <- if      (length(a2)==0)         1
#     else if (all(a2 %in% a1))       0
#     else                             sum(!a2 %in% a1)/length(a2)
#   }
#   
#   # 6) Locus‐level metrics
#   discordant_loci        <- sum(zero_bits)
#   replacement_pattern_sum<- sum(replacement_scores)
#   locus_discordance_rate <- discordant_loci / n_loci
#   replacement_pattern_score <- replacement_pattern_sum / n_loci
#   
#   # 7) Original allele‐level features
#   # (jaccard, retention, allele_gain already computed above)
#   
#   # 8) New content‐based features
#   avg_locus_share <- mean(share_ratios)
#   min_locus_share <- min(share_ratios)
#   
#   # core_prop: alleles shared at *every* locus
#   core_set <- Reduce(intersect, lapply(loci, function(l) {
#     intersect(split1[[l]] %||% character(0),
#               split2[[l]] %||% character(0))
#   }))
#   core_prop <- if (n_loci>0) length(core_set)/n_loci else 0
#   
#   richness_diff <- abs(length(alleles1) - length(alleles2))
#   coi_ratio     <- if (length(alleles1)>0) length(alleles2)/length(alleles1) else 0
#   
#   # zero_run_max: longest run of unshared loci
#   r <- rle(zero_bits)
#   zero_run_max <- if (any(r$values==1)) max(r$lengths[r$values==1]) else 0
#   
#   # 9) Return all features
#   c(
#     # original
#     jaccard_similarity          = jaccard,
#     allele_retention_rate       = retention,
#     allele_gain                 = allele_gain,
#     locus_discordance_rate      = locus_discordance_rate,
#     allele_transition_asymmetry = trans_asym,
#     replacement_pattern_score   = replacement_pattern_score,
#     
#     # new
#     avg_locus_share             = avg_locus_share,
#     min_locus_share             = min_locus_share,
#     core_prop                   = core_prop,
#     richness_diff               = richness_diff,
#     coi_ratio                   = coi_ratio,
#     zero_run_max                = zero_run_max
#   )
# }


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
output_file <- paste0("delta_features_",site,"_", DATA_TYPE, ".csv")

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

delta_metrics_df_final <- read.csv(paste0("delta_features_",site,"_", DATA_TYPE, ".csv"))

delta_metrics_df_final <- delta_metrics_df_final %>%
  select(PairsID, everything())

write.csv(delta_metrics_df_final, paste0("delta_features_",site,"_", DATA_TYPE, ".csv"), row.names = F)



##### MERGE FEATURES AND OUTPUT ------

# FEATURES <- left_join(dres0_long_final_summarized, delta_metrics_df_final, by = "PairsID")

FEATURES <- left_join(delta_metrics_df_final, PAIRS_METADATA, by = "PairsID")


write.csv(FEATURES, paste0(site,"_", DATA_TYPE,".csv"), row.names = F)
