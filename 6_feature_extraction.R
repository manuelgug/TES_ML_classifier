

library(data.table)
library(dplyr)
library(dcifer)
library(reshape2)
library(tidyr)
library(purrr)
library(stringr)
library(igraph)
library(Matrix)


site <- "Cabo_Delgado"


#select data type betweem "TRAINING_DATA" or "REAL_DATA"
DATA_TYPE = "TRAINING_DATA"


if (DATA_TYPE == "TRAINING_DATA") {
  
  # for the training data:
  PAIRS_METADATA <- readRDS(paste0("PAIRS_METADATA_",site,".RDS"))
  PAIRS_GENOMIC <- readRDS(paste0("PAIRS_GENOMIC_",site,".RDS"))
  
  PAIRS_GENOMIC <- as.data.table(PAIRS_GENOMIC)
  PAIRS_METADATA <- as.data.table(PAIRS_METADATA)
  
} else if (DATA_TYPE == "REAL_DATA") {
  
  #the actual data
  PAIRS_GENOMIC <- read.csv(paste0("genomic_updated_",site,".csv"))
  PAIRS_METADATA <- read.csv(paste0("metadata_updated_",site,".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
  
  PAIRS_GENOMIC <- as.data.table(PAIRS_GENOMIC)
  PAIRS_METADATA <- as.data.table(PAIRS_METADATA)

  PAIRS_GENOMIC <- PAIRS_GENOMIC %>% rename(read_counts = reads)  
  PAIRS_GENOMIC <- PAIRS_GENOMIC %>% rename(NIDA = sampleID) # format genomic file
  suppressWarnings(PAIRS_GENOMIC <- PAIRS_GENOMIC %>%
    separate(NIDA, into = c("NIDA", "run"), sep = "__", remove = TRUE))
  PAIRS_GENOMIC <- inner_join(PAIRS_METADATA, PAIRS_GENOMIC, by = "NIDA")
  
  PAIRS_METADATA <- PAIRS_METADATA %>% select(PairsID, NIDA, time_point) # format metadata file
  PAIRS_METADATA <- PAIRS_METADATA %>%
    pivot_wider(names_from = time_point, values_from = NIDA, names_prefix = "NIDA") %>%
    rename(NIDA1 = NIDAD0, NIDA2 = NIDADx)
  
} else {
  
  print("Incorrect data type. Options are 'TRAINING_DATA' and 'REAL_DATA'.")
  
}




####### DCIFER'S IBD #######------------------  

dsmp <- formatDat(PAIRS_GENOMIC, svar = "NIDA", lvar = "locus", avar = "allele")

#use already calculated coi instead?
lrank <- 2
coi   <- getCOI(dsmp, lrank = lrank)

afreq <- calcAfreq(dsmp, coi, tol = 1e-5) 

dres0 <- ibdDat(dsmp, coi, afreq, pval = TRUE, confint = TRUE, rnull = 0, 
                alpha = 0.05, nr = 1e3)   

gc()

suppressWarnings({
  dres0_long <- melt(dres0)
})
dres0_long$value <- ifelse(dres0_long$Var1 == dres0_long$Var2, 1, dres0_long$value) # put 1 if the sample is compared with itself

#need extra foramtting for real data after passing through dcifer because the nidas...
if (DATA_TYPE == "REAL_DATA"){
  
  dres0_long$Var1 <- as.character(dres0_long$Var1)
  dres0_long$Var2 <- as.character(dres0_long$Var2)
  dres0_long <- dres0_long %>%
    mutate(Var1 = if_else(!str_detect(Var1, "\\."), paste0(Var1, ".0"), Var1),
           Var2 = if_else(!str_detect(Var2, "\\."), paste0(Var2, ".0"), Var2))  # If no ".", add ".0"
}

dres0_long <- dres0_long[dres0_long$Var3 == "estimate" & !is.na(dres0_long$value),]
dres0_long <- dres0_long %>% select(-Var3)
colnames(dres0_long) <- c("infection1", "infection2", "IBD_estimate")


# First match: D0 with infection1 and Dx with infection2
match1 <- PAIRS_METADATA %>%
  left_join(dres0_long, by = c("NIDA1" = "infection1", "NIDA2" = "infection2")) %>%
  select(PairsID, NIDA1, NIDA2, IBD_estimate)

match1 <- match1[!is.na(match1$IBD_estimate),]


# Second match: Dx with infection1 and D0 with infection2
match2 <- PAIRS_METADATA %>%
  left_join(dres0_long, by = c("NIDA2" = "infection1", "NIDA1" = "infection2")) %>%
  select(PairsID, NIDA1, NIDA2, IBD_estimate)

match2 <- match2[!is.na(match2$IBD_estimate),]

dres0_long_final<- rbind(match1, match2)
dres0_long_final <- distinct(dres0_long_final)

dres0_long_final_summarized <- dres0_long_final %>%
  select(PairsID, IBD_estimate) %>%
  arrange(PairsID)

dres0_long_final_summarized <- dres0_long_final_summarized[complete.cases(dres0_long_final_summarized),]

# MERGE WITH METADATA
dres0_long_final_summarized <- merge(dres0_long_final_summarized, PAIRS_METADATA, by = "PairsID")

if (DATA_TYPE == "REAL_DATA"){
  
  metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
  
  metadata_updated$naive_coi <- round(metadata_updated$naive_coi)
  
  metadata_updated_wide <- metadata_updated %>%
    pivot_wider(
      id_cols = PairsID, 
      names_from = time_point, 
      values_from = c(NIDA, naive_coi), 
      names_glue = "{.value}_{time_point}"
    )
 
  metadata_updated_wide$eCOI_pairs <- paste0(metadata_updated_wide$naive_coi_D0, "__", metadata_updated_wide$naive_coi_Dx)
   
  dres0_long_final_summarized <- merge(dres0_long_final_summarized, metadata_updated_wide[c("PairsID", "eCOI_pairs")], by = "PairsID")
  
}


####### DIVERSITY/DELTA FEATURES #######------------------


calculate_features_optimized <- function(sample1, sample2) {
  
  # Extract alleles and loci
  alleles1 <- unique(sample1$allele)
  loci1 <- unique(sample1$locus)
  alleles2 <- unique(sample2$allele)
  loci2 <- unique(sample2$locus)
  
  # Global presence/absence vectors
  all_alleles <- unique(c(alleles1, alleles2))
  m1 <- as.integer(all_alleles %in% alleles1)
  m2 <- as.integer(all_alleles %in% alleles2)
  
  # Global similarity metrics
  intersection <- sum(m1 & m2)
  union_alleles <- sum(m1 | m2)
  num_alleles1 <- length(alleles1)
  num_alleles2 <- length(alleles2)
  
  # Allele turnover features
  retention_rate <- intersection / num_alleles1
  allele_gain <- sum(m2 & !m1) / union_alleles
  
  # Locus-based calculations
  all_loci <- unique(c(loci1, loci2))
  locus_concordance <- 0
  
  # 4. Locus Concordance Rate & 5. Max Locus Similarity
  for(locus in all_loci) {
    a1 <- unique(sample1$allele[sample1$locus == locus])
    a2 <- unique(sample2$allele[sample2$locus == locus])
    
    # Locus concordance
    if(length(intersect(a1, a2)) > 0) locus_concordance <- locus_concordance + 1
  
  return(c(
    jaccard_similarity = intersection / union_alleles,
    allele_retention_rate = retention_rate,
    allele_gain = allele_gain,
    locus_concordance_rate = locus_concordance / length(all_loci)
  ))
  }
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

library(progress)
library(parallel)

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

FEATURES <- left_join(dres0_long_final_summarized, delta_metrics_df_final, by = "PairsID")

write.csv(FEATURES, paste0(site,"_", DATA_TYPE,".csv"), row.names = F)

