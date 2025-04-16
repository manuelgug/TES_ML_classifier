library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)

site <- "Zambezia"

# Load data
metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), 
                             stringsAsFactors = FALSE, 
                             colClasses = c(NIDA = "character"))

data <- read.csv(paste0("genomic_updated_", site, ".csv"), 
                 stringsAsFactors = FALSE, 
                 colClasses = c(sampleID = "character")) %>%
  filter(data_type == "tes") %>%
  mutate(sampleID = gsub("__.*", "", sampleID))

coi_stats <- read.csv(paste0("coi_stats_", site, ".csv"), 
                      stringsAsFactors = FALSE, 
                      colClasses = c(NIDA = "character")) %>%
  mutate(NIDA = gsub("__.*", "", NIDA))

# Join time_point info
data <- left_join(data, metadata_updated[c("NIDA", "time_point")], 
                  by = c("sampleID" = "NIDA"))

# coi_stats <- left_join(coi_stats, metadata_updated[c("NIDA", "time_point")], 
#                        by = "NIDA") %>%
#   filter(time_point == "D0")

# Allele frequencies
AF <- moire::summarize_allele_freqs(readRDS(paste0("coi_mcmc_", site, ".RDS"))) %>%
  select(locus, allele, post_allele_freqs_mean)


### 1) Separate monoclonal infections ----
clones <- coi_stats %>%
  filter(naive_coi < 1.1) %>%
  pull(NIDA)

clones_genomic <- data %>%
  filter(sampleID %in% clones)


### 2) Create artificial clones from polyclonal infections ----
N_CLONES <- 1000 - length(clones)
data_polyclonal <- data %>% filter(!sampleID %in% clones)

polyclonal_samples <- unique(data_polyclonal$sampleID)
iterations <- ceiling(N_CLONES / length(polyclonal_samples))

random_draw <- function(locus_data) {
  slice_sample(locus_data, n = 1, weight_by = post_allele_freqs_mean) %>%
    pull(allele)
}

library(purrr)
library(data.table)  # for faster data manipulation

# Convert AF to data.table for speed
AF_dt <- as.data.table(AF)
setkey(AF_dt, locus, allele)

# Pre-allocate list to store clone data
sampled_monoclonals_list <- vector("list", length(polyclonal_samples))

set.seed(42069)  # Set seed once, outside the loop

for (idx in seq_along(polyclonal_samples)) {
  sample_id <- polyclonal_samples[idx]
  samp <- filter(data_polyclonal, sampleID == sample_id)
  
  allele_combinations <- samp %>%
    select(locus, allele) %>%
    distinct() %>%
    left_join(AF, by = c("locus", "allele"))
  
  # Split by locus for faster access
  split_locus <- split(allele_combinations, allele_combinations$locus)
  
  clones_per_sample <- map_dfr(seq_len(iterations), function(i) {
    clone_id <- paste0(sample_id, "__clone_", i)
    
    # Sample one allele per locus
    sampled_alleles <- map_dfr(split_locus, function(df) {
      df %>%
        slice_sample(n = 1, weight_by = post_allele_freqs_mean) %>%
        select(locus, allele)
    })
    
    clone_result <- sampled_alleles %>%
      mutate(sampleID = clone_id) %>%
      left_join(samp[c("allele", "reads", "norm.reads.locus", "data_type", "time_point")], 
                by = "allele") %>%
      select(sampleID, locus, allele, reads, norm.reads.locus, data_type, time_point)
    
    return(clone_result)
  })
  
  sampled_monoclonals_list[[idx]] <- clones_per_sample
}

sampled_monoclonals <- bind_rows(sampled_monoclonals_list)
all_clones <- bind_rows(sampled_monoclonals, clones_genomic)

cat("Unique clones retained:", length(unique(all_clones$sampleID)), "\n")


### 3) Calculate pairwise proportion of shared alleles ----
alleles <- all_clones %>%
  group_by(sampleID) %>%
  summarize(alleles = list(allele), .groups = "drop")

n <- nrow(alleles)
comparison_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(alleles$sampleID, alleles$sampleID))

# Precompute all unique index combinations (i < j)
pair_indices <- combn(n, 2)

# Function to compute Jaccard similarity (or any overlap metric)
compute_shared_prop <- function(i, j) {
  a_i <- alleles$alleles[[i]]
  a_j <- alleles$alleles[[j]]
  shared <- length(intersect(a_i, a_j))
  total <- length(unique(c(a_i, a_j)))
  shared / total
}

# Apply function across all combinations
shared_values <- apply(pair_indices, 2, function(x) compute_shared_prop(x[1], x[2]))

# Fill the upper and lower triangle
for (k in seq_along(shared_values)) {
  i <- pair_indices[1, k]
  j <- pair_indices[2, k]
  comparison_matrix[i, j] <- shared_values[k]
  comparison_matrix[j, i] <- shared_values[k]
}

# Fill diagonal with 1s (since each sample fully overlaps with itself)
diag(comparison_matrix) <- 1

comparison_df <- as.data.frame(comparison_matrix) %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
  filter(Var1 < Var2) %>%
  arrange(value)

# Plot histogram
hist <- ggplot(comparison_df, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.6, color = "black") +
  labs(title = "", x = "Proportion of Shared Alleles", y = "Density of Monoclonal Comparisons") +
  theme_minimal() +
  xlim(0, 1)

ggsave(paste0("hist_shared_alleles_", site, ".png"), hist, 
       height = 5, width = 8, bg = "white", dpi = 300)


### 4) Clean monoclonal data: remove redundant clone groups ----
same_clones <- comparison_df %>% filter(value > 0.99)

groups <- list()

for (i in seq_len(nrow(same_clones))) {
  v1 <- same_clones$Var1[i]
  v2 <- same_clones$Var2[i]
  matched <- FALSE
  
  for (j in seq_along(groups)) {
    if (v1 %in% groups[[j]] || v2 %in% groups[[j]]) {
      groups[[j]] <- unique(c(groups[[j]], v1, v2))
      matched <- TRUE
      break
    }
  }
  
  if (!matched) {
    groups <- append(groups, list(c(v1, v2)))
  }
}

# Keep only one representative per group
clones_to_remove <- unlist(lapply(groups, function(g) g[-1]))
all_clones <- all_clones %>% filter(!sampleID %in% clones_to_remove)

cat("Final unique clones after filtering:", length(unique(all_clones$sampleID)), "\n")


### 5) Export cleaned clone data ----
write.csv(all_clones, paste0("clones_genomic_data_", site, ".csv"), row.names = FALSE)
