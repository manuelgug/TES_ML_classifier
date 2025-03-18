

library(data.table)
library(dplyr)
library(progress)
library(parallel)
library(dcifer)
library(igraph)
library(Matrix)
library(reshape2)
library(tidyr)
library(purrr)


#select data type betweem "training_data" or "real_data"
DATA_TYPE = "REAL_DATA"


if (DATA_TYPE == "TRAINING_DATA") {
  
  # for the training data:
  PAIRS_METADATA <- readRDS("PAIRS_METADATA.RDS")
  PAIRS_GENOMIC <- readRDS("PAIRS_GENOMIC.RDS")
  
  PAIRS_GENOMIC <- as.data.table(PAIRS_GENOMIC)
  PAIRS_METADATA <- as.data.table(PAIRS_METADATA)
  
} else if (DATA_TYPE == "REAL_DATA") {
  
  #the actual data
  PAIRS_GENOMIC <- read.csv("genomic_updated.csv")
  PAIRS_METADATA <- read.csv("metadata_updated.csv", stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
  
  PAIRS_GENOMIC <- as.data.table(PAIRS_GENOMIC)
  PAIRS_METADATA <- as.data.table(PAIRS_METADATA)
  
  
  PAIRS_GENOMIC <- PAIRS_GENOMIC %>% rename(NIDA = sampleID) # format genomic file
  PAIRS_GENOMIC <- PAIRS_GENOMIC %>%
    separate(NIDA, into = c("NIDA", "run"), sep = "__", remove = TRUE)
  PAIRS_GENOMIC <- inner_join(PAIRS_METADATA, PAIRS_GENOMIC, by = "NIDA")
  
  PAIRS_METADATA <- PAIRS_METADATA %>% select(PairsID, NIDA, time_point) # format metadata file
  PAIRS_METADATA <- PAIRS_METADATA %>%
    pivot_wider(names_from = time_point, values_from = NIDA, names_prefix = "NIDA") %>%
    rename(NIDA1 = NIDAD0, NIDA2 = NIDADx)
  
} else {
  
  print("Incorrect data type. Options are 'TRAINING_DATA' and 'REAL_DATA'.")
  
}





####### DIVERSITY/DELTA FEATURES #######------------------

calculate_features_optimized <- function(sample1, sample2) {
  
  # Extract alleles and read counts
  alleles1 <- sample1$allele
  alleles2 <- sample2$allele
  read_counts1 <- sample1$read_counts
  read_counts2 <- sample2$read_counts
  
  # Get unique alleles across both samples
  all_alleles <- unique(c(alleles1, alleles2))
  
  # Create presence/absence vectors
  m1 <- as.integer(all_alleles %in% alleles1)
  m2 <- as.integer(all_alleles %in% alleles2)
  
  # # Create frequency vectors
  # f1 <- numeric(length(all_alleles))
  # f2 <- numeric(length(all_alleles))
  # 
  # match_alleles1 <- match(alleles1, all_alleles)
  # match_alleles2 <- match(alleles2, all_alleles)
  # 
  # f1[match_alleles1] <- read_counts1
  # f2[match_alleles2] <- read_counts2
  # 
  # # Normalize read counts (relative abundance per sample)
  # f1_norm <- f1 / sum(f1, na.rm = TRUE)
  # f2_norm <- f2 / sum(f2, na.rm = TRUE)
  
  # ---- Genetic Similarity Measures ---- #
  
  # Jaccard similarity
  shared_alleles_mask <- m1 & m2
  jaccard <- sum(shared_alleles_mask) / sum(m1 | m2)
  
  # # Manhattan distance
  # manhattan_dist <- sum(abs(f1_norm - f2_norm), na.rm = TRUE)
  
  # # Cosine similarity
  # cosine_sim <- sum(f1_norm * f2_norm, na.rm = TRUE) / (sqrt(sum(f1_norm^2, na.rm = TRUE)) * sqrt(sum(f2_norm^2, na.rm = TRUE)))

  # Bray-Curtis dissimilarity
  # bray_curtis <- 1 - (2 * sum(pmin(f1_norm, f2_norm), na.rm = TRUE) / (sum(f1_norm, na.rm = TRUE) + sum(f2_norm, na.rm = TRUE)))
  
  # ---- Allele Dynamics Features ---- #
  
  # Number of alleles per sample
  num_alleles1 <- length(unique(alleles1))
  num_alleles2 <- length(unique(alleles2))
  
  # Allele retention rate
  retention_rate <- sum(shared_alleles_mask) / num_alleles1
  
  # Allele gain and loss
  allele_gain <- sum(m2 & !m1) / sum(m1 | m2)  # Fraction of new alleles
  # allele_loss <- sum(m1 & !m2) / sum(m1 | m2)  # Fraction of lost alleles
  
  # ---- Diversity Change Features ---- #
  
  # # Shannon diversity index
  # shannon_div1 <- -sum(f1_norm * log(f1_norm + 1e-10), na.rm = TRUE)
  # shannon_div2 <- -sum(f2_norm * log(f2_norm + 1e-10), na.rm = TRUE)
  # shannon_diversity_change <- shannon_div2 - shannon_div1   # Change in SHannon div (ΔShannon)
  # 
  # # Expected heterozygosity (He) per locus
  # heterozygosity1 <- 1 - mean(f1_norm^2, na.rm = TRUE)  # Mean He across loci for sample 1
  # heterozygosity2 <- 1 - mean(f2_norm^2, na.rm = TRUE)  # Mean He across loci for sample 2
  # heterozygosity_change <- heterozygosity2 - heterozygosity1   # Change in heterozygosity (ΔHe)
  
  # ---- Return All Features ---- #
  
  return(c(
    jaccard_similarity = jaccard,
    # manhattan_distance = manhattan_dist,
    # cosine_similarity = cosine_sim,
    #bray_curtis_dissimilarity = bray_curtis,
    allele_retention_rate = retention_rate,
    allele_gain = allele_gain
    #allele_loss = allele_loss,
    #shannon_diversity_change = shannon_diversity_change,
    #heterozygosity_change = heterozygosity_change
  ))
}


# Extract unique pairs once
unique_pairs <- unique(PAIRS_METADATA$PairsID)

# Create a named vector for mapping PairsID to D0 and Dx for faster access
pairs_to_d0 <- PAIRS_METADATA$NIDA1[match(unique_pairs, PAIRS_METADATA$PairsID)]
pairs_to_dx <- PAIRS_METADATA$NIDA2[match(unique_pairs, PAIRS_METADATA$PairsID)]

# Preprocess merged_dfs_filtered into a list for fast access
merged_dfs_filtered <- split(PAIRS_GENOMIC, PAIRS_GENOMIC$NIDA)

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
output_file <- paste0("delta_features_", DATA_TYPE, ".csv")

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

delta_metrics_df_final <- read.csv(paste0("delta_features_", DATA_TYPE, ".csv"))

delta_metrics_df_final <- delta_metrics_df_final %>%
  select(PairsID, everything())

write.csv(delta_metrics_df_final, paste0("delta_features_", DATA_TYPE, ".csv"), row.names = F)

# free RAM
rm(pairs_to_d0, pairs_to_dx, delta_metrics_df_final)
gc()


####### NETWORK FEATURES #######------------------

# Pre-set the key for faster subsetting (using the data.table index)
setkey(PAIRS_GENOMIC, PairsID)

# Function to process each pair
process_pair <- function(pair_id, PAIRS_GENOMIC) {
  # Subset the data efficiently by PairsID (using data.table indexing)
  merged_dfs_filtered_count_sub <- PAIRS_GENOMIC[.(pair_id)]
  
  # Use table function to create the co-occurrence matrix
  w <- table(merged_dfs_filtered_count_sub$allele, merged_dfs_filtered_count_sub$time_point)
  
  return(w)
}

# Function to process a chunk of pairs
process_chunk <- function(pairs_chunk, PAIRS_GENOMIC) {
  result_chunk <- lapply(pairs_chunk, process_pair, PAIRS_GENOMIC)
  names(result_chunk) <- pairs_chunk  # Assign names within the chunk
  return(result_chunk)
}

chunk_size <- 2000  # Define the chunk size
pair_chunks <- split(unique(PAIRS_GENOMIC$PairsID), ceiling(seq_along(unique(PAIRS_GENOMIC$PairsID)) / chunk_size))
num_cores <- detectCores() - 0  # Use all but one core

# Process each chunk in parallel
dcasted_list <- list()

for (chunk in pair_chunks) {
  # Parallel processing of the current chunk
  print(paste0("Processing pairs ",  min(chunk), " to ", max(chunk)))
  
  chunk_result <- mclapply(chunk, process_pair, PAIRS_GENOMIC, mc.cores = num_cores)
  
  # Convert the result into a named list and combine with the main result
  names(chunk_result) <- chunk
  dcasted_list <- c(dcasted_list, chunk_result)
  
}

cat("Calculated co-occurrence matrices for", length(dcasted_list), "pairs")


#checkpoint
saveRDS(dcasted_list, paste0("coocurrence_matrices_", DATA_TYPE, ".RDS"))
# dcasted_list <- readRDS(paste0(outprefix, "_coocurrence_matrices.RDS"))

pair_ids <- names(dcasted_list)

# Function to calculate network metrics
network_metrics <- function(graph) {
  list(
    density = edge_density(graph),
    clustering_coefficient = transitivity(graph, type = "average"),
    # eigenvector_centrality = mean(eigen_centrality(graph)$vector),
    avg_path_length = mean_distance(graph, directed = FALSE),  # Faster version
    # assortativity = assortativity_degree(graph),
    betweenness_centrality = mean(betweenness(graph, normalized = TRUE)),
    # degree_centrality = mean(degree(graph)),
    modularity_value = modularity(cluster_fast_greedy(graph)),
    #avg_node_strength = mean(strength(graph)),
    edge_betweenness = mean(edge_betweenness(graph)),
    # motifs_triangles = count_motifs(graph, size = 3),  # Counting triangles
    # global_efficiency = 1 / mean_distance(graph, directed = FALSE),
    constraint = mean(constraint(graph))
  )
}


# Optimized function to process each PairsID with error handling
process_pair_id <- function(pair_id) {
  tryCatch({
    x <- dcasted_list[[as.character(pair_id)]]
    x[is.na(x)] <- 0
    x <- (x > 0) + 0  # Convert to binary
    
    # Co-occurrence matrix as sparse matrix
    v <- as(Matrix(x, sparse = TRUE) %*% t(Matrix(x, sparse = TRUE)), "sparseMatrix")
    diag(v) <- 0
    
    # Create graph and compute metrics
    g <- graph_from_adjacency_matrix(v, mode = "undirected", weighted = TRUE, diag = FALSE)
    g <- simplify(g)
    metrics <- network_metrics(g)
    metrics$PairsID <- pair_id
    
    return(metrics)
  }, error = function(e) {
    message(paste("Error processing PairsID:", pair_id, "-", e$message))
    return(NULL)
  })
}

# Parallel processing setup
num_cores <- detectCores() - 0
chunk_size <- 1000  # Adjust based on available memory

# Initialize the progress bar
pb <- progress_bar$new(
  format = "Processing [:bar] :percent :elapsedfull",
  total = ceiling(length(pair_ids) / chunk_size),
  width = 100
)

# Batch write results to disk in chunks to avoid memory bloat
write_chunk_to_csv <- function(data, file_path, append = FALSE) {
  if (append) {
    fwrite(data, file_path, append = TRUE)
  } else {
    fwrite(data, file_path)
  }
}

gc()

# Temporary file for writing
output_file <- paste0("network_features_", DATA_TYPE, ".csv")

# Process pair_ids in chunks
for (i in seq(1, length(pair_ids), by = chunk_size)) {
  chunk_pair_ids <- pair_ids[i:min(i + chunk_size - 1, length(pair_ids))]
  
  print(paste0("Processing PairsID ", first(chunk_pair_ids), " to ", last(chunk_pair_ids), "..."))
  
  # Parallel processing of chunk
  results <- mclapply(chunk_pair_ids, process_pair_id, mc.cores = num_cores)
  
  # Filter out NULL results
  metrics_list_chunk <- Filter(Negate(is.null), results)
  
  if (length(metrics_list_chunk) > 0) {
    # Combine results into a data table
    network_metrics_df_chunk <- rbindlist(metrics_list_chunk, fill = TRUE)
    
    # Incrementally write results to CSV
    write_chunk_to_csv(network_metrics_df_chunk, output_file, append = (i > 1))
  }
  
  pb$tick()
}

network_metrics_df_final <- read.csv(paste0("network_features_", DATA_TYPE, ".csv"))

network_metrics_df_final <- network_metrics_df_final %>%
  select(PairsID, everything())

write.csv(network_metrics_df_final, paste0("network_features_", DATA_TYPE, ".csv"), row.names = F)

gc()



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

dres0_long_final_summarized[complete.cases(dres0_long_final_summarized),]

write.csv(dres0_long_final_summarized, paste0("IBD_features_", DATA_TYPE, ".csv"), row.names = F)



###### MERGE FEATURES ###### ----------------

network_features <- read.csv(paste0("network_features_", DATA_TYPE, ".csv"))
delta_features <- read.csv(paste0("delta_features_", DATA_TYPE, ".csv"))
ibd_features <- read.csv(paste0("IBD_features_", DATA_TYPE, ".csv"))

# List of data frames to be merged
all_Feats <- list(delta_features, network_features, ibd_features)
all_Feats <- map(all_Feats, ~ .x %>% mutate(PairsID = as.character(PairsID)))

# Merge all data frames by PairsID
all_Feats_merged_df <- reduce(all_Feats, full_join, by = "PairsID")

#extra formatting if training data in order to add metadata (labels basically)
if(DATA_TYPE == "TRAINING_DATA"){
  
  # Merge features with metadata
  PAIRS_METADATA$PairsID <- as.character(PAIRS_METADATA$PairsID)
  all_Feats_merged_df <- inner_join(all_Feats_merged_df, PAIRS_METADATA, by = "PairsID")
  
}

## OUTPUT
write.csv(all_Feats_merged_df, paste0(DATA_TYPE,".csv"), row.names = F)
