
library(tibble)
library(dplyr)
library(purrr)
library(data.table)
library(tidyr)
library(stringr)
library(parallel)
library(progress)


# Set site
site <- "Zambezia"


#### INPUTS

#clones genomic data
CLONES <- read.csv(paste0("clones_genomic_data_", site, ".csv"))

CLONES <- CLONES %>% select(sampleID, locus, allele)


# metatada
PAIRS_METADATA <- readRDS(paste0("PAIRS_METADATA_", site, ".RDS"))

# Input features testing data for subsetting
test_data <- readRDS(paste0(site, "_test_data.RDS"))
TEST_META <- test_data %>%
  mutate(PairsID = as.character(PairsID),
         labels = factor(labels, levels = c("NI", "R"))) %>%
  select(PairsID, labels, pair_type = eCOI_pairs) %>%  # Renamed to pair_type for consistency
  as.data.table()

# subset metadata to only account for test data
PAIRS_METADATA <- PAIRS_METADATA[PAIRS_METADATA$PairsID %in% TEST_META$PairsID,]


########
generate_evenness_table <- function(max_clones = 7) {
  
  # Normalized Shannon index
  shannon_index <- function(p) {
    if (length(p) == 1) return(0)
    -sum(p * log(p)) / log(length(p))
  }
  
  # Generator for HIGH (equal proportions)
  generate_high <- function(n) rep(1/n, n)
  
  # Generator for skewed (medium, low, ultra low), with unique proportions
  generate_skewed <- function(n, level) {
    if (n == 1) return(c(1))
    
    dominant <- switch(level,
                       "medium" = ifelse(n == 2, 0.6, 0.5),
                       "low" = 0.8,
                       "ultra_low" = 0.9,
                       stop("Invalid level"))
    
    rest_weights <- seq(from = n - 1, to = 1)
    rest_props <- rest_weights / sum(rest_weights) * (1 - dominant)
    
    p <- c(dominant, rest_props)
    sort(p, decreasing = TRUE)
  }
  
  tibble::tibble(
    n_clones = rep(1:max_clones, each = 4),
    evenness_level = rep(c("high", "medium", "low", "ultra_low"), times = max_clones)
  ) %>%
    dplyr::mutate(
      clone_proportions = purrr::map2(n_clones, evenness_level, 
                                      ~ round(
                                        if (.y == "high") generate_high(.x) else generate_skewed(.x, .y), 
                                        3)
      ),
      has_duplicates = purrr::map_lgl(clone_proportions, ~ any(duplicated(.x))),
      shannon_evenness = purrr::map_dbl(clone_proportions, shannon_index),
      proportions_key = purrr::map_chr(clone_proportions, ~ paste(.x, collapse = "-"))
    ) %>%
    filter(evenness_level == "high" | !has_duplicates) %>%
    distinct(n_clones, proportions_key, .keep_all = TRUE) %>%
    select(n_clones, evenness_level, clone_proportions, shannon_evenness)
}

evenness_tbl <- generate_evenness_table(max_clones = 7)

print(evenness_tbl)

# Load progress package
library(progress)

# Pre-compute constants and lookups
eveness_levels <- unique(evenness_tbl$evenness_level)
pairsID       <- PAIRS_METADATA$PairsID

# Pre-compute evenness lookup for faster access
evenness_lookup <- split(evenness_tbl, evenness_tbl$evenness_level)

# Pre-compute PAIRS_METADATA lookups to avoid repeated filtering
pairs_lookup <- setNames(
  lapply(seq_along(pairsID), function(idx) {
    pair <- pairsID[idx]
    pair_data <- PAIRS_METADATA[PAIRS_METADATA$PairsID == pair, ]
    
    # Process D0 samples
    D0_split <- gsub("(?<=\\d)_(?=\\d)", ",", pair_data$D0_sample, perl = TRUE)
    D0_vec   <- unlist(strsplit(D0_split, ","))
    D0_named <- setNames(D0_vec, rep(pair_data$NIDA1, length(D0_vec)))
    
    # Process Dx samples
    Dx_split <- gsub("(?<=\\d)_(?=\\d)", ",", pair_data$Dx_sample, perl = TRUE)
    Dx_vec   <- unlist(strsplit(Dx_split, ","))
    Dx_named <- setNames(Dx_vec, rep(pair_data$NIDA2, length(Dx_vec)))
    
    list(D0 = D0_named, Dx = Dx_named)
  }),
  pairsID
)

# Initialize result list for better memory management
FNR_list <- vector("list", length(eveness_levels))
names(FNR_list) <- eveness_levels

# Optimized FNR function (unchanged)
simulate_detected_alleles <- function(df, proportions = NULL, alpha = 1, beta = 1) {
  clones <- unique(df$sampleID)
  N      <- length(clones)
  
  if (is.null(proportions)) {
    proportions <- rep(1, N)
  } else if (length(proportions) != N) {
    stop("Length of proportions does not match number of clones.")
  }
  
  detection_probs <- 1 - (((1 - proportions) * (N - 1) / N) ^ alpha) * beta
  detection_lookup <- setNames(detection_probs, clones)
  
  df$detection_prob <- detection_lookup[df$sampleID]
  set.seed(42)
  df$detected <- runif(nrow(df)) < df$detection_prob
  
  detected_df <- df[df$detected, c("sampleID", "locus", "allele", "time_point",
                                   grep("NIDA", names(df), value = TRUE))]
  detected_df <- detected_df[!duplicated(detected_df), ]
  
  detected_loci <- unique(detected_df$locus)
  all_loci      <- unique(df$locus)
  missing_loci  <- setdiff(all_loci, detected_loci)
  
  if (length(missing_loci) > 0) {
    missing_df <- df[df$locus %in% missing_loci, ]
    additional_alleles <- do.call(rbind, lapply(missing_loci, function(loc) {
      loc_data <- missing_df[missing_df$locus == loc, ]
      max_idx  <- which.max(loc_data$detection_prob)
      loc_data[max_idx, c("sampleID", "locus", "allele", "time_point",
                          grep("NIDA", names(loc_data), value = TRUE))]
    }))
    additional_alleles <- additional_alleles[!duplicated(additional_alleles), ]
    detected_df       <- rbind(detected_df, additional_alleles)
  }
  
  detected_df
}

# Main processing loop
for (level_idx in seq_along(eveness_levels)) {
  level_select <- eveness_levels[level_idx]
  cat("Processing evenness level:", level_select, "\n")
  
  # Progress bar
  pairs_pb <- progress_bar$new(
    format = paste0("  [", level_select, "] Processing pairs [:bar] :percent :current/:total ETA: :eta"),
    total = length(pairsID), clear = FALSE, width = 80
  )
  
  freq_scheme   <- evenness_lookup[[level_select]]
  level_results <- vector("list", length(pairsID))
  
  for (pair_idx in seq_along(pairsID)) {
    pair <- pairsID[pair_idx]
    pairs_pb$tick()
    
    # Fetch precomputed maps
    pair_data <- pairs_lookup[[pair_idx]]
    D0         <- pair_data$D0
    Dx         <- pair_data$Dx
    
    # Subset CLONES once each
    CLONES_D0 <- CLONES[CLONES$sampleID %in% D0, ]
    CLONES_D0$time_point <- "D0"
    CLONES_D0$NIDA       <- names(D0)[1]
    
    CLONES_Dx <- CLONES[CLONES$sampleID %in% Dx, ]
    CLONES_Dx$time_point <- "Dx"
    CLONES_Dx$NIDA       <- names(Dx)[1]
    
    # Counts
    n0 <- length(unique(CLONES_D0$sampleID))
    nx <- length(unique(CLONES_Dx$sampleID))
    
    # **Safe** proportion lookup with fallback to 1
    idx0 <- which(freq_scheme$n_clones == n0)
    idxx <- which(freq_scheme$n_clones == nx)
    p0   <- if (length(idx0)) freq_scheme$clone_proportions[[idx0]] else 1
    px   <- if (length(idxx)) freq_scheme$clone_proportions[[idxx]] else 1
    
    # Simulate
    D0_fnr <- simulate_detected_alleles(CLONES_D0, proportions = p0, alpha = 3, beta = 1.5)
    Dx_fnr <- simulate_detected_alleles(CLONES_Dx, proportions = px, alpha = 3, beta = 1.5)
    
    # Dedup
    D0_fnr <- D0_fnr[!duplicated(D0_fnr[c("locus", "allele")]), ]
    Dx_fnr <- Dx_fnr[!duplicated(Dx_fnr[c("locus", "allele")]), ]
    
    fnr <- rbind(D0_fnr, Dx_fnr)
    fnr$PairsID        <- pair
    fnr$evenness_level <- level_select
    
    level_results[[pair_idx]] <- fnr
  }
  
  pairs_pb$terminate()
  FNR_list[[level_idx]] <- do.call(rbind, level_results)
  cat("Completed evenness level:", level_select, "\n\n")
}

# Final assembly
FNR_ALL <- do.call(rbind, FNR_list)
cat("Processing complete. Total rows in FNR_ALL:", nrow(FNR_ALL), "\n")

length(unique(FNR_ALL$PairsID))


## PUT EVERYTHING INTO FNR_ALL

############################ LOOP END!!!!!!!! ###############################



###########################################
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


# Extract unique pairs once
unique_pairs <- unique(PAIRS_METADATA$PairsID)

# Create a named vector for mapping PairsID to D0 and Dx for faster access
pairs_to_d0 <- PAIRS_METADATA$NIDA1[match(unique_pairs, PAIRS_METADATA$PairsID)]
pairs_to_dx <- PAIRS_METADATA$NIDA2[match(unique_pairs, PAIRS_METADATA$PairsID)]


#############################

for (evenness_level in eveness_levels){
  
  print(paste0("Processing: ", evenness_level))
  
  ### SUBSET LEVEL
  FNR <- FNR_ALL[FNR_ALL$evenness_level == evenness_level,]
  ####################33
  
  
  # Preprocess merged_dfs_filtered into a list for fast access
  # PAIRS_GENOMIC <- PAIRS_GENOMIC[PAIRS_GENOMIC$PairsID == 22,]
  # merged_dfs_filtered <- split(PAIRS_GENOMIC, PAIRS_GENOMIC$NIDA)
  merged_dfs_filtered <- split(FNR, FNR$NIDA)
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
  
  # Temporary file for writing output
  output_file <- paste0("FNR_features_", evenness_level, "_evenness_" ,site, ".csv")
  
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
  
}





###################33
#PREDICTIONS
#####################



