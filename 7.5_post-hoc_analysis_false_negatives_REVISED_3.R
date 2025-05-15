
library(tibble)
library(dplyr)
library(purrr)



# Set site
site <- "Tete"


#### INPUTS

#clones genomic data
CLONES <- read.csv("clones_genomic_data_Tete.csv")

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
                       "medium" = ifelse(n == 2, 0.6, 0.5),  # exception for 2 clones
                       "low" = 0.8,
                       "ultra low" = 0.9,
                       stop("Invalid level"))
    
    rest_weights <- seq(n - 1, 1)  # descending to create unique values
    rest_props <- rest_weights / sum(rest_weights) * (1 - dominant)
    
    p <- c(dominant, rest_props)
    sort(p, decreasing = TRUE)
  }
  
  # Create the full table with 4 levels now
  tibble::tibble(
    n_clones = rep(1:max_clones, each = 4),
    evenness_level = rep(c("high", "medium", "low", "ultra low"), times = max_clones)
  ) %>%
    mutate(
      clone_proportions = purrr::map2(n_clones, evenness_level, 
                                      ~ if (.y == "high") generate_high(.x) else generate_skewed(.x, .y)),
      has_duplicates = purrr::map_lgl(clone_proportions, ~ any(duplicated(round(.x, 5)))),
      shannon_evenness = purrr::map_dbl(clone_proportions, shannon_index),
      proportions_key = purrr::map_chr(clone_proportions, ~ paste(round(.x, 5), collapse = "-"))
    ) %>%
    filter(evenness_level == "high" | !has_duplicates) %>%
    distinct(n_clones, proportions_key, .keep_all = TRUE) %>%
    select(n_clones, evenness_level, clone_proportions, shannon_evenness)
}

evenness_tbl <- generate_evenness_table(max_clones = 7)
print(evenness_tbl)
########


#### TEST PAIR 5

pairsID <- PAIRS_METADATA$PairsID

pair <- pairsID[5]

D0 <- PAIRS_METADATA %>%
  filter(PairsID == pair) %>%
  mutate(D0_sample_split = str_replace_all(D0_sample, "(?<=\\d)_(?=\\d)", ",")) %>%
  mutate(D0_vec = str_split(D0_sample_split, ",")) %>%
  select(NIDA1, D0_vec) %>%
  unnest(D0_vec) %>%
  pull(D0_vec) %>%
  setNames(rep(PAIRS_METADATA %>% filter(PairsID == pair) %>% pull(NIDA1), 
               each = length(.)))
Dx <- PAIRS_METADATA %>%
  filter(PairsID == pair) %>%
  mutate(Dx_sample_split = str_replace_all(Dx_sample, "(?<=\\d)_(?=\\d)", ",")) %>%
  mutate(Dx_vec = str_split(Dx_sample_split, ",")) %>%
  select(NIDA1, Dx_vec) %>%
  unnest(Dx_vec) %>%
  pull(Dx_vec) %>%
  setNames(rep(PAIRS_METADATA %>% filter(PairsID == pair) %>% pull(NIDA2), 
               each = length(.)))

CLONES_D0<- CLONES[CLONES$sampleID %in% D0,]
CLONES_D0$time_point <- "D0"
CLONES_D0$NIDA <- unique(names(D0))

CLONES_Dx<- CLONES[CLONES$sampleID %in% Dx,]
CLONES_Dx$time_point <- "Dx"
CLONES_Dx$NIDA <- unique(names(Dx))

# number of clones, useful for selecting a wsaf scheme
n_clones_D0 <- length(unique(CLONES_D0$sampleID))
n_clones_Dx <- length(unique(CLONES_Dx$sampleID))

# FNR function
simulate_detected_alleles <- function(df, proportions = NULL, alpha = 1, beta = 1) {
  
  # 1. Get unique clones
  clones <- unique(df$sampleID)
  N <- length(clones)
  
  # 2. Assign proportions (if not provided, assume no dropout aka complete detection of alleles)
  if (is.null(proportions)) {
    proportions <- rep(1, N)
  } else if (length(proportions) != N) {
    stop("Length of proportions does not match number of clones.")
  }
  
  # 3. Create clone detection probabilities with alpha and beta
  detection_probs <- 1 - (( (1 - proportions) * (N - 1) / N ) ^ alpha) * beta
  clone_df <- data.frame(sampleID = clones,
                         proportion = proportions,
                         detection_prob = detection_probs,
                         stringsAsFactors = FALSE)
  
  print(detection_probs)
  
  # 4. Merge detection probabilities back into data
  df <- df %>%
    left_join(clone_df, by = "sampleID")
  
  # 5. Simulate allele detection
  set.seed(42)  # For reproducibility
  df$detected <- runif(nrow(df)) < df$detection_prob
  
  # 6. Ensure that at least one allele per locus is kept
  # For each locus, if no alleles are detected, keep one from any clone
  detected_df <- df %>%
    filter(detected) %>%
    select(sampleID, locus, allele, time_point, matches("NIDA")) %>%
    distinct()
  
  # 7. Identify loci that need to be filled (no detected alleles)
  missing_loci <- setdiff(unique(df$locus), unique(detected_df$locus))
  
  if (length(missing_loci) > 0) {
    # Keep one allele from each locus where no alleles were detected, from the clone with the highest proportion
    additional_alleles <- df %>%
      filter(locus %in% missing_loci) %>%
      group_by(locus) %>%
      arrange(desc(proportion)) %>%  # Sort by proportion, highest first
      slice(1) %>%  # Keep the allele from the clone with the highest proportion
      select(sampleID, locus, allele, time_point, matches("NIDA")) %>%
      distinct()
    
    # Append the additional alleles to the detected_df
    detected_df <- bind_rows(detected_df, additional_alleles)
  }
  
  return(detected_df)
}

freq_schemes <-111111 ###########################

# Apply the function
# the lower the alpha, the more penalized are the clones in low proportion and less penalized the clones in high proportion (steepness)
# beta is a scaling factor; overall severity: lower betas increase the prob of detection for all clones (less severe overall)
D0_fnr <- simulate_detected_alleles(CLONES_D0, proportions = c(0.9, 0.1), alpha = 3, beta = 1.5) 
Dx_fnr <- simulate_detected_alleles(CLONES_Dx, proportions = c(1), alpha = 3, beta = 1.5) 

# TURN INTO A MIX
D0_fnr <- D0_fnr %>%
  distinct(locus, allele, .keep_all = TRUE)

Dx_fnr <- Dx_fnr %>%
  distinct(locus, allele, .keep_all = TRUE)

## concat + add pairsid column
fnr <- rbind(D0_fnr, Dx_fnr)
fnr$PairsID <- pair

fnr <- as_tibble(FNR)


### MEGRE ALL fnrs INTO capitlized FNR
FNR <- bind_rows(fnr)


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
output_file <- paste0("LALALA_features_",site, ".csv")

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
