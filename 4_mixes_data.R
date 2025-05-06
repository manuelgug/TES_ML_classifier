library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(vegan)


site <- "Zambezia"
initial_sample_size <- 200
set.seed(69420) 

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
synthetic_clones <- all_clones[grepl("_clone_", all_clones)]

#extract natural clones
clones_genomic_TES <- filter(clones_genomic, sampleID %in% tesclones)

paste0(length(tesclones), " natural clones found in the data.")


if (length(tesclones) >= 3){ ## IF THERE ARE AT LEAST 3 CLONES IN THE DATA (can do rarefaction):
  print("Performing allele accumulationc curbve.")
  
  # Create the matrix
  allele_matrix <- clones_genomic_TES %>%
    select(sampleID, allele) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = allele, values_from = present, values_fill = 0)
  
  # Save sampleIDs and convert to matrix
  sample_ids <- allele_matrix$sampleID
  allele_matrix <- allele_matrix %>% select(-sampleID)
  allele_matrix <- as.matrix(allele_matrix)
  rownames(allele_matrix) <- sample_ids
  
  # make curve
  spec_accum <- specaccum(allele_matrix, method = "random", )
  
  # 1. Build your accumulation data.frame
  df_accum <- data.frame(
    sites = spec_accum$sites,
    rich  = spec_accum$rich,
    sd    = spec_accum$sd
  )
  
  # 2. Fit the SSasymp model
  fit_asymp <- nls(rich ~ SSasymp(sites, Asym, R0, lrc), data = df_accum)
  params    <- coef(fit_asymp)
  Asym      <- params["Asym"]
  R0        <- params["R0"]
  lrc       <- params["lrc"]
  
  # 3. Curve completeness threshold
  threshold <- 0.95
  target_richness <- threshold * Asym
  
  # 4. Solve for the required # of clones (sites) to hit that target:
  s_needed <- -log((target_richness - Asym)/(R0 - Asym)) / exp(lrc)
  additional <- s_needed - max(df_accum$sites)
  additional_clones_needed <- round(additional, 0)
  additional_clones_needed <- ifelse(additional_clones_needed < 0, 0, additional_clones_needed) # if curve is complete, no need for more clones
  
  # 5. Generate model predictions up to s_needed
  new_sites <- seq(0, s_needed, length.out = 200)
  df_model  <- data.frame(
    sites = new_sites,
    rich  = predict(fit_asymp, newdata = data.frame(sites = new_sites))
  )
  
  # 6. Plot everything
  CURVE <- ggplot(df_accum, aes(x = sites, y = rich)) +
    # ±1 SD ribbon
    geom_ribbon(aes(ymin = rich - sd, ymax = rich + sd),
                fill = "steelblue", alpha = 0.3) +
    # observed curve
    geom_line(color = "steelblue", size = 1) +
    # fitted model
    geom_line(data = df_model, aes(x = sites, y = rich),
              color = "red", size = 1) +
    # vertical line at s_needed
    geom_vline(xintercept = s_needed, linetype = "dashed") +
    # annotate how many more clones
    annotate("text",
             x = s_needed, 
             y = min(df_accum$rich),
             label = paste0(additional_clones_needed, " more clones"),
             angle = 90, vjust = 1.2) +
    labs(
      x     = "Number of clones sampled",
      y     = "Cumulative unique alleles",
      title = paste0("Rarefaction + Asymptote (", threshold*100, "%)"),
      subtitle = paste0("Asymptote ≈ ", round(Asym,1),
                        " | Observed final ≈ ", tail(df_accum$rich,1))
    ) +
    theme_minimal()
  
  CURVE
  
  ggsave(paste0(site, "_allele_curve.png"), CURVE, dpi = 300, height = 7, width = 8, bg = "white")
  
  
  #### extract samples synthetic clones needed to complete the curve
  additional_synthetic_clones <- sample(synthetic_clones, additional_clones_needed)
  clones_genomic_synthetic <- filter(clones_genomic, sampleID %in% additional_synthetic_clones)
  
  ### put everything together
  clones_genomic <- rbind(clones_genomic_TES,clones_genomic_synthetic)
  
  
} else { # IF LESS THAN 2 TES CLONES (can't do ratefaction), INCLUDING NONE:
  
  clone_cap <- 30 # sample up to 30 synthetic clones by default
  
  print(paste0("No clones in the data. Subsampling ", clone_cap, " synthetic clones."))
  
  n_clones <- clone_cap - length(tesclones)
  
  additional_synthetic_clones <- sample(synthetic_clones, n_clones)
  clones_genomic_synthetic <- filter(clones_genomic, sampleID %in% additional_synthetic_clones)
  
  
  ### put everything together
  clones_genomic <- rbind(clones_genomic_TES,clones_genomic_synthetic)
}

N_CLONES <- length(unique(clones_genomic$sampleID))
paste0("FINAL CLONES: ", N_CLONES)

write.csv(clones_genomic, paste0(site, "_clones_genomic_SELECTED_FOR_TRAINING.csv"), row.names = F)


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
metadata_updated_tes <- metadata_updated[!is.na(metadata_updated$PairsID),] # only tes data
coi_values <- sort(unique(round(metadata_updated_tes$naive_coi)))
coi_values <- coi_values[coi_values > 1]


# Build mixes with sampling protection
strain_mixes <- map(coi_values, ~ create_combinations_df_safe(nidas, .x, MAX_COMBOS))

# Add mix_1 (monoclonals)
strain_mixes <- setNames(c(list(mix_1 = data.frame(strain_1 = nidas)),
                           strain_mixes),
                         c("mix_1", paste0("mix_", coi_values)))

############################################

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
