library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)


site <- "Inhambane"
clone_cap <- 10
initial_sample_size <- 200
set.seed(69420)  # set once globally

# --- 1. Load Data ----
clones_genomic <- read.csv(paste0("clones_genomic_data_", site, ".csv"),
                           stringsAsFactors = FALSE,
                           colClasses = c(sampleID = "character"))

metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"),
                             stringsAsFactors = FALSE,
                             colClasses = c(NIDA = "character"))

# --- 2. Subsample Clones ----
all_clones <- unique(clones_genomic$sampleID)
tesclones <- all_clones[!grepl("_clone_", all_clones)]
n_synthetic_clones <- clone_cap - length(tesclones)
N_CLONES <- length(all_clones)

#minimum amount of clones possible is the n of tes clones. 
subsampled_clones <- sample(all_clones, ifelse(N_CLONES > clone_cap, 
                                               ifelse(n_synthetic_clones < 0, 0, n_synthetic_clones)
                                               , N_CLONES))

keep_clones <- union(tesclones, subsampled_clones)


clones_genomic <- filter(clones_genomic, sampleID %in% keep_clones)

# Check: all clones are monoallelic per locus
stopifnot(all(
  clones_genomic %>%
    group_by(sampleID, locus) %>%
    summarise(n = n_distinct(allele), .groups = "drop") %>%
    pull(n) == 1
))

# --- 3. Create All Mixes ----
min_coi <- round(min(metadata_updated$naive_coi))
max_coi <- round(max(metadata_updated$naive_coi))

nidas <- unique(clones_genomic$sampleID)

# Create all mixes
create_combinations_df <- function(vec, k) {
  as.data.frame(t(combn(vec, k)), stringsAsFactors = FALSE) |>
    setNames(paste0("strain_", 1:k))
}

#########################################3

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

# Unique COI values > 1
coi_values <- sort(unique(round(metadata_updated$naive_coi)))
coi_values <- coi_values[coi_values > 1]


# Build mixes with sampling protection
strain_mixes <- map(coi_values, ~ create_combinations_df_safe(nidas, .x, MAX_COMBOS))

# Add mix_1 (monoclonals)
strain_mixes <- setNames(c(list(mix_1 = data.frame(strain_1 = nidas)),
                           strain_mixes),
                         c("mix_1", paste0("mix_", coi_values)))

############################################

# # Define mix sizes based on unique naive_coi values > 1
# coi_values <- sort(unique(metadata_updated$naive_coi[metadata_updated$naive_coi > 1]))
# 
# strain_mixes <- lapply(coi_values, \(k) create_combinations_df(nidas, k))
# 
# # Add mix_1 manually
# strain_mixes <- setNames(c(list(mix_1 = data.frame(strain_1 = nidas)),
#                            strain_mixes),
#                          c("mix_1", paste0("mix_", coi_values)))

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
    
    selected <- clones_genomic %>%
      filter(sampleID %in% strains)
    
    total_reads <- sum(selected$reads)
    
    prop_reads <- selected %>%
      group_by(sampleID) %>%
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
saveRDS(MIXES_METADATA, paste0("MIXES_METADATA_", site, ".RDS"))
saveRDS(MIXES_GENOMIC, paste0("MIXES_GENOMIC_", site, ".RDS"))

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

ggsave(paste0("mixes_EDA_", site, "_", N_CLONES, "_clones.png"), allele_plot, width = 9, height = 6, dpi = 300, bg = "white")
