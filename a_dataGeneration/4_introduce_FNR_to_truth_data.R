
library(dplyr)
library(stringr)
library(ggplot2)


######################################################################################################
##### DATA CLEANING
######################################################################################################


# -------------------- 0) PARAMETERS --------------------

site <- "Tete"
cum_curve_threshold <- 0.999
main_dir <- "."
metadata_file <- paste0("metadata_tes_", site, ".csv")
maf_filter <- 0.01
min_allele_read_count <- 10
OPTIM_AMPSET <- FALSE # FALSE uses pfPHAST pool 1A's 20 amps; TRUE uses all 165 amps from pool 1A

# -------------------- 1) IMPORT DATA --------------------

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
filtered_dirs <- list.dirs(main_dir, recursive = FALSE, full.names = TRUE)

list_of_dfs <- list()
for (dir in filtered_dirs) {
  folder_name <- gsub("_RESULTS.*", "", basename(dir))
  csv_file <- file.path(dir, "allele_data_global_max_0_filtered.csv")
  if (file.exists(csv_file)) {
    print(paste0("Importing: ", dir))
    df <- read.csv(csv_file)
    df$run <- folder_name
    list_of_dfs[[folder_name]] <- df
  }
}
merged_dfs <- bind_rows(list_of_dfs)

# -------------------- 2) FORMAT DATA --------------------

merged_dfs <- merged_dfs %>%
  mutate(
    sampleID = gsub("N|_S.*$", "", sampleID),
    sampleID = if_else(str_detect(sampleID, "_"), sampleID, paste0(sampleID, ".0")),
    sampleID = gsub("_", ".", sampleID)
  )

metadata <- metadata %>%
  mutate(
    NIDA = gsub("N|_S.*$", "", NIDA),
    NIDA = if_else(!str_detect(NIDA, "\\."), paste0(NIDA, ".0"), NIDA)
  )

merged_dfs <- merged_dfs[merged_dfs$sampleID %in% metadata$NIDA,]
merged_dfs$sampleID <- paste0(merged_dfs$sampleID, "__", merged_dfs$run)
merged_dfs <- merged_dfs %>% select(sampleID, locus, pseudo_cigar, reads, norm.reads.locus, run)

# -------------------- 3) CLEAN TES DATA --------------------

clean_data <- function(df, maf_filter, min_read_count) {
  good_sampleID <- df %>%
    group_by(sampleID, locus) %>%
    summarise(reads = sum(reads), .groups = "drop") %>%
    filter(reads >= 100) %>%
    group_by(sampleID) %>%
    summarise(n_loci = n_distinct(locus), .groups = "drop") %>%
    filter(n_loci >= 50) %>%
    pull(sampleID)
  
  df <- df[df$sampleID %in% good_sampleID, ]
  
  high_read_samples <- df %>%
    group_by(sampleID) %>%
    summarise(total_reads = sum(reads), .groups = "drop") %>%
    filter(total_reads > 10000) %>%
    pull(sampleID)
  
  df <- df[df$sampleID %in% high_read_samples, ]
  df <- df[grepl("-1A$", df$locus), ]
  df$pseudo_cigar <- gsub("\\d+\\+[^N]*N", "", df$pseudo_cigar)
  df$pseudo_cigar[df$pseudo_cigar == "" | is.na(df$pseudo_cigar)] <- "."
  
  df <- df %>%
    group_by(sampleID, locus, pseudo_cigar) %>%
    summarise(reads = sum(reads), norm.reads.locus = sum(norm.reads.locus), .groups = "drop") %>%
    mutate(allele = paste0(locus, "__", pseudo_cigar)) %>%
    filter(!grepl("I=|D=", allele)) %>%
    filter(norm.reads.locus > maf_filter & reads > min_read_count) %>%
    select(-pseudo_cigar)
  
  return(df)
}

merged_dfs_agg <- clean_data(merged_dfs, maf_filter, min_allele_read_count)

# -------------------- 4) KEEP SAMPLES WITH PAIRS --------------------

ids <- data.frame(NIDA = unique(gsub("__.*$", "", merged_dfs_agg$sampleID)))
ids_merged <- merge(ids, metadata, by = "NIDA")

pairs_counts <- ids_merged %>% group_by(PairsID) %>% summarise(n_samples = n(), .groups = "drop") %>%
  filter(n_samples == 2)

metadata_updated <- metadata[metadata$PairsID %in% pairs_counts$PairsID, ]
merged_dfs_agg <- merged_dfs_agg[gsub("__.*$", "", merged_dfs_agg$sampleID) %in% metadata_updated$NIDA, ]

# -------------------- 5) CLEAN SITE DATA --------------------

genomic_site <- read.csv(paste0("genomic_site_", site, ".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA2 = "character")) %>%
  rename(sampleID = NIDA2)

genomic_site_agg <- clean_data(genomic_site, maf_filter, min_allele_read_count)

genomic_site_agg <- genomic_site_agg[!sapply(genomic_site_agg$sampleID, function(x) any(grepl(x, merged_dfs_agg$sampleID, fixed = TRUE))), ]

# -------------------- 6) MERGE DATA --------------------

merged_dfs_agg$data_type <- "tes"
genomic_site_agg$data_type <- "site"
data_all <- rbind(merged_dfs_agg, genomic_site_agg)

# -------------------- 7) SHARED AMPLICONS --------------------

amps_data <- data_all %>%
  group_by(locus) %>%
  summarise(in_n_samples = n_distinct(sampleID), .groups = "drop") %>%
  arrange(desc(in_n_samples))

selected_amps <- amps_data$in_n_samples == max(amps_data$in_n_samples)
data_all <- data_all[data_all$locus %in% amps_data[selected_amps, ]$locus, ]

# -------------------- 8) SELECT AMPS --------------------

heterozygosity_per_sample <- data_all %>%
  group_by(sampleID, locus) %>%
  mutate(freq = norm.reads.locus / sum(norm.reads.locus, na.rm = TRUE)) %>%
  summarise(heterozygosity = 1 - sum(freq^2, na.rm = TRUE), .groups = "drop")

variability_per_locus <- heterozygosity_per_sample %>%
  group_by(locus) %>%
  summarise(mean_He = mean(heterozygosity, na.rm = TRUE), sd_He = sd(heterozygosity, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_He))

if (OPTIM_AMPSET) {
  variability_per_locus <- variability_per_locus %>%
    mutate(prob_identity = 1 - mean_He)
  
  n_loci <- nrow(variability_per_locus)
  cumulative_He <- sapply(1:n_loci, function(i) 1 - prod(variability_per_locus$prob_identity[1:i]))
  cum_curve <- data.frame(loci_included = 1:n_loci, multilocus_He = cumulative_He)
  amps_var_data <- cbind(variability_per_locus, cum_curve)
  
  #write.csv(amps_var_data, paste0("amps_variation_", site, ".csv"), row.names = FALSE)
  
  top_n_amps <- length(cum_curve$multilocus_He[cum_curve$multilocus_He < cum_curve_threshold])
  top_var_amps <- variability_per_locus$locus[1:top_n_amps]
  data_all <- data_all[data_all$locus %in% top_var_amps, ]
  
} else {
  
  madhito_amps <- read.csv("madhito_20amps.csv")
  data_all <- data_all[data_all$locus %in% madhito_amps$locus, ]
  
  # used_amps <- unique(data_all$locus)
  # variability_per_locus2 <- left_join(variability_per_locus, madhito_amps, by = "locus") %>%
  #   mutate(used = locus %in% used_amps, color = "selected_madhito_amps")
  # 
  # selected_madhito_amps <- variability_per_locus2[!is.na(variability_per_locus2$color), ] %>%
  #   arrange(desc(mean_He)) %>%
  #   mutate(prob_identity = 1 - mean_He)
  # 
  # n_loci <- nrow(selected_madhito_amps)
  # cumulative_He <- sapply(1:n_loci, function(i) 1 - prod(selected_madhito_amps$prob_identity[1:i]))
  # cum_curve <- data.frame(loci_included = 1:n_loci, multilocus_He = cumulative_He)
  # amps_var_data <- cbind(selected_madhito_amps, cum_curve) %>% select(locus, mean_He, multilocus_He)
  # 
  #write.csv(amps_var_data, paste0("amps_variation_", site, ".csv"), row.names = FALSE)
}

# -------------------- 9) CHECKS --------------------

# samples (failure pairs + site)
length(unique(data_all$sampleID))
# shared locus
length(unique(data_all$locus))
# failure pairs
length(unique(data_all[data_all$data_type == "tes", ]$sampleID)) / 2

# -------------------- 10) OUTPUTS --------------------

write.csv(data_all, paste0("genomic_updated_", site, ".csv"), row.names = FALSE)
write.csv(metadata_updated, paste0("metadata_updated_", site, ".csv"), row.names = FALSE)




######################################################################################################
##### OFFSET NAIVE COI CALCULATION
######################################################################################################

# calculate offset naive coi directly, without moire
coi_stats <- data_all %>%
  group_by(sampleID, locus) %>%
  summarise(n_alleles = n_distinct(allele), .groups = "drop") %>%
  group_by(sampleID) %>%
  summarise(
    offset_naive_coi = {
      vals <- sort(n_alleles, decreasing = T)
      if (length(vals) >= 2) vals[2] else vals[1] # OFFSET = THE 2ND VALUE
    }
  ) %>% rename(NIDA = sampleID)


# update metadata
coi_stats$NIDA <- gsub("__.*", "", coi_stats$NIDA)
metadata_updated2 <- merge(metadata_updated, coi_stats, by = c("NIDA"))
metadata_updated2 <- metadata_updated2 %>% arrange(SampleID)


write.csv(metadata_updated2, paste0("metadata_updated_", site, ".csv"), row.names = F)




######################################################################################################
##### FNR CALCULATION
######################################################################################################

# Load data
metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), 
                             stringsAsFactors = FALSE, 
                             colClasses = c(NIDA = "character"))

metadata_updated$time_point <- ifelse(is.na(metadata_updated$time_point), "D0", metadata_updated$time_point) # avoid issues with NA in the site samples

data <- read.csv(paste0("genomic_updated_", site, ".csv"),
                 stringsAsFactors = FALSE,
                 colClasses = c(sampleID = "character")) %>%
  #filter(data_type == "tes") %>%
  mutate(sampleID = gsub("__.*", "", sampleID))


# 2) PROPORTIONS
# for proportions, coi should be the numver of strains, so i'm picking all loci that have n_alleles == coi for each mix and then averaging out the norm.reads.locus
allele_count_per_locus <- data %>% 
  group_by(sampleID, locus) %>%
  summarise(n_alleles = length(unique(allele)))

# Merge the two data frames by NIDA
merged_df <- left_join(allele_count_per_locus, metadata_updated[c("NIDA", "offset_naive_coi")], by = c("sampleID" = "NIDA"))
#merged_df <- merged_df[merged_df$offset_naive_coi > 1,] # no COI = 1

# Filter rows where COI matches n_alleles: here are the loci that will be used for the proportions calculation for each mix (NIDA)
matched_df <- merged_df %>%
  filter(n_alleles == offset_naive_coi)

# Subset labcontrols_genomic by rowwise matching of NIDA and locus from matched_df
subset_labcontrols_genomic <- semi_join(data, matched_df, by = c("sampleID", "locus"))

# keep uniques
subset_labcontrols_genomic <- unique(subset_labcontrols_genomic)

# Step 1: Rank norm.reads.locus descendingly within each NIDA and locus
ranked_df <- subset_labcontrols_genomic %>%
  group_by(sampleID, locus) %>%
  arrange(desc(norm.reads.locus), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Step 2: Compute the average of norm.reads.locus by rank per NIDA
average_by_rank <- ranked_df %>%
  group_by(sampleID, rank) %>%
  summarise(avg_norm_reads = mean(norm.reads.locus, na.rm = TRUE), .groups = "drop") %>%
  arrange(sampleID, rank)

library(tidyr)

average_by_rank_wide <- average_by_rank %>%
  pivot_wider(
    names_from = rank,
    values_from = avg_norm_reads,
    names_prefix = "rank_"
  )


FNRs <- average_by_rank_wide %>%
  rowwise() %>%
  mutate(FNR = list({
    ranks <- c_across(starts_with("rank_"))
    fnrs <- c()
    fnrs <- c(fnrs, rep(0.11, sum(ranks >= 0.02 & ranks < 0.03, na.rm = TRUE))) ### FNR of 0.14 for strains below 0.03 (experimentally checked with the lab controls dataset!) == DROP 2 ALLELES
    fnrs <- c(fnrs, rep(0.16, sum(ranks < 0.02, na.rm = TRUE))) ### FNR of 0.15 for strains below 0.03 (experimentally checked with the lab controls dataset!) == DROP 3 ALLELES
    fnrs
  })) %>%
  ungroup()

FNRs <- left_join(FNRs, metadata_updated[c("NIDA", "offset_naive_coi")], by = c("sampleID" = "NIDA"))

# final values: number of elements in FNR are the strains that should have given FNR for that mix because there were many strains below 0.05 in the real data)
FNRs <- FNRs %>%
  distinct(offset_naive_coi, FNR)

#delete empty vectors or coi = 1 (NO FNRs APPLIED for those)
FNRs <- FNRs %>% filter(lengths(FNR) > 0)
FNRs <- FNRs[FNRs$offset_naive_coi > 1,]

saveRDS(FNRs, paste0(site, "_FNRs.RDS"))



######################################################################################################
##### GENERATE CLONE DATA
######################################################################################################


library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(purrr)
library(data.table)

# Load data
metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), 
                             stringsAsFactors = FALSE, 
                             colClasses = c(NIDA = "character"))

metadata_updated$time_point <- ifelse(is.na(metadata_updated$time_point), "D0", metadata_updated$time_point) # avoid issues with NA in the site samples

data <- read.csv(paste0("genomic_updated_", site, ".csv"), 
                 stringsAsFactors = FALSE, 
                 colClasses = c(sampleID = "character")) %>%
  #filter(data_type == "tes") %>%
  mutate(sampleID = gsub("__.*", "", sampleID))

coi_stats <- metadata_updated %>% select(SampleID, PairsID, NIDA, offset_naive_coi, time_point)


# Join time_point info
data <- left_join(data, metadata_updated[c("NIDA", "time_point", "offset_naive_coi")], 
                  by = c("sampleID" = "NIDA"))


# samples with only 1 allele for each loci, regardless of offset naive coi, which sometimes considers samples with more than 1 allele because of the offset thing
clones <- data %>%
  group_by(sampleID, locus) %>%
  summarise(n = n_distinct(allele), .groups = "drop") %>%
  group_by(sampleID) %>%
  summarise(all_single = all(n == 1), .groups = "drop") %>%
  filter(all_single) %>%
  pull(sampleID)

clones_genomic <- data %>%
  filter(sampleID %in% clones)

clones_genomic$time_point <- ifelse(is.na(clones_genomic$time_point), "D0", clones_genomic$time_point) # avoid issues with NA in the site samples

clones_genomic <- clones_genomic[clones_genomic$time_point == "D0",] # only D0

length(unique(clones_genomic$sampleID))


### 2) Create artificial clones from polyclonal infections ----
N_CLONES <- 1000 - length(clones)
data_polyclonal <- data %>% filter(!sampleID %in% clones & time_point == "D0") # only D0

polyclonal_samples <- unique(data_polyclonal$sampleID)
iterations <- ceiling(N_CLONES / length(polyclonal_samples))

# testing new synthetic clone sampling;

clone_result <- data.frame()  # Initialize empty result data frame

# Loop through each unique polyclonal sample
for (sid in polyclonal_samples) {
  sample_data <- data_polyclonal[data_polyclonal$sampleID == sid, ]
  
  for (i in 1:iterations) {
    # Sample one allele per locus using norm.reads.locus as weights
    sampled <- sample_data %>%
      group_by(locus) %>%
      slice_sample(n = 1, weight_by = norm.reads.locus, replace = TRUE) %>%
      ungroup()
    
    # Add unique clone ID
    sampled$sampleID <- paste0(sid, "__clone_", i)
    
    # Append to result
    clone_result <- bind_rows(clone_result, sampled)
  }
}

all_clones <- bind_rows(clone_result, clones_genomic)

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




######################################################################################################
##### GENERATE MIX DATA
######################################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(vegan)


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
  print("Performing allele accumulationc curve.")
  
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
min_coi <- round(min(metadata_updated$offset_naive_coi))
max_coi <- round(max(metadata_updated$offset_naive_coi))

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
coi_values <- sort(unique(round(metadata_updated_tes$offset_naive_coi)))
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

ggsave(paste0("mixes_EDA_", site, "_", N_CLONES, "_clones.png"), allele_plot, width = 9, height = 6, dpi = 300, bg = "white")




######################################################################################################
##### GENERATE PAIRS DATA
######################################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
metadata_updated <- metadata_updated[!is.na(metadata_updated$PairsID),] # only tes data
#------------------------------------------------
# 1. Determine pairs of mixes based on COIs
#------------------------------------------------

metadata_updated$offset_naive_coi <- round(metadata_updated$offset_naive_coi)

metadata_updated_wide <- metadata_updated %>%
  pivot_wider(
    id_cols = PairsID,
    names_from = time_point,
    values_from = c(NIDA, offset_naive_coi),
    names_glue = "{.value}_{time_point}"
  )

unique_combos <- metadata_updated_wide %>%
  distinct(offset_naive_coi_D0, offset_naive_coi_Dx) %>%
  mutate(
    offset_naive_coi_D0 = paste0("mix", offset_naive_coi_D0),
    offset_naive_coi_Dx = paste0("mix", offset_naive_coi_Dx)
  ) %>%
  arrange(offset_naive_coi_D0, offset_naive_coi_Dx)

#------------------------------------------------
# 2. Create pairs
#------------------------------------------------

MIXES_METADATA <- readRDS(paste0("MIXES_METADATA_", site, ".RDS"))
MIXES_GENOMIC <- readRDS(paste0("MIXES_GENOMIC_", site, ".RDS"))

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

saveRDS(merged_dfs, paste0("PAIRS_GENOMIC_", site, ".RDS"))

#------------------------------------------------
# 3. Compare allele content of pairs
#------------------------------------------------

alleles <- merged_dfs %>%
  group_by(PairsID, NIDA, time_point) %>%
  summarize(alleles = list(allele), .groups = "drop") %>%
  arrange(PairsID, time_point)

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

saveRDS(PAIRS_metadata, paste0("PAIRS_METADATA_", site, ".RDS"))

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

write.csv(PAIRS_summary, paste0("PAIRS_SUMMARY_", site, ".csv"), row.names = FALSE)



######################################################################################################
##### CALCULATE FEATURES
######################################################################################################

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

options(datatable.warn = FALSE)

chunk_size      <- 1000
num_cores       <- detectCores()
output_suffixes <- c("TRAINING_DATA", "REAL_DATA")

# Helper: load & standardize both PAIRS_GENOMIC & PAIRS_METADATA
load_pairs <- function(type) {
  if (type == "TRAINING_DATA") {
    meta <- readRDS(paste0("PAIRS_METADATA_", site, ".RDS"))
    geno <- readRDS(paste0("PAIRS_GENOMIC_" , site, ".RDS"))
    
    setDT(meta); meta <- copy(meta)
    setDT(geno); geno <- copy(geno)
    geno <- merge(geno, meta, by = "PairsID", all.x = TRUE)
    
  } else if (type == "REAL_DATA") {
    geno <- fread(paste0("genomic_updated_", site, ".csv"))
    meta <- fread(paste0("metadata_updated_", site, ".csv"),
                  colClasses = c(NIDA = "character")) %>%
      filter(!is.na(time_point))
    
    setDT(geno); geno <- copy(geno)
    setDT(meta); meta <- copy(meta)
    
    # rename + split
    setnames(geno, c("reads", "sampleID"), c("read_counts", "NIDA"))
    geno[, c("NIDA", "run") := tstrsplit(NIDA, "__", fixed = TRUE)]
    
    # merge + derive pair_type
    geno <- merge(meta, geno, by = "NIDA", all = FALSE)
    coi_sum <- geno[
      time_point %in% c("D0","Dx"),
      .(
        offset_naive_coi_D0 = unique(offset_naive_coi[time_point=="D0"]),
        offset_naive_coi_Dx = unique(offset_naive_coi[time_point=="Dx"])
      ),
      by = PairsID
    ]
    coi_sum[, pair_type := paste0(offset_naive_coi_D0, "__", offset_naive_coi_Dx)]
    geno <- merge(geno, coi_sum[, .(PairsID, pair_type)], by = "PairsID", all.x = TRUE)
    
    # wide metadata, then convert to data.table and copy to avoid selfref warnings
    meta <- meta %>%
      select(PairsID, NIDA, time_point) %>%
      pivot_wider(
        names_from  = time_point,
        values_from = NIDA,
        names_prefix = "NIDA"
      ) %>%
      rename(NIDA1 = NIDAD0, NIDA2 = NIDADx) %>%
      left_join(coi_sum, by = "PairsID")
    setDT(meta); meta <- copy(meta)
    
  } else {
    stop("Unknown type: ", type)
  }
  list(META = meta, GENO = geno)
}

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

# Process both types
# … your library loads, parameter setup, load_pairs, calculate_features_optimized, etc …

# Process both types, silencing self‑ref warnings
suppressWarnings({
  for (TYPE in output_suffixes) {
    cat("=== Processing", TYPE, "===\n")
    
    parts           <- load_pairs(TYPE)
    PAIRS_METADATA  <- parts$META
    PAIRS_GENOMIC   <- parts$GENO
    
    # Pre-split genomic by NIDA
    merged_list <- split(PAIRS_GENOMIC, by = "NIDA")
    merged_list <- lapply(merged_list, function(dt) {
      dt[, c("PairsID","time_point") := NULL]
      unique(dt)
    })
    gc()
    
    unique_pairs <- unique(PAIRS_METADATA$PairsID)
    pairs_to_d0  <- PAIRS_METADATA$NIDA1[match(unique_pairs, PAIRS_METADATA$PairsID)]
    pairs_to_dx  <- PAIRS_METADATA$NIDA2[match(unique_pairs, PAIRS_METADATA$PairsID)]
    
    out_file <- paste0("delta_features_", site, "_", TYPE, ".csv")
    
    # progress bar
    pb <- progress_bar$new(
      format = "  [:bar] :percent ETA: :eta",
      total  = ceiling(length(unique_pairs) / chunk_size),
      width  = 60
    )
    
    first_chunk <- TRUE
    for (i in seq(1, length(unique_pairs), by = chunk_size)) {
      idx <- i:min(i + chunk_size - 1, length(unique_pairs))
      cat(" chunk", i, "–", idx[length(idx)], "\n")
      
      results <- mclapply(idx, function(k) {
        s1 <- merged_list[[ pairs_to_d0[k] ]]
        s2 <- merged_list[[ pairs_to_dx[k] ]]
        feats <- calculate_features_optimized(s1, s2)
        data.frame(t(feats), PairsID = unique_pairs[k], check.names = FALSE)
      }, mc.cores = num_cores)
      
      dt_chunk <- rbindlist(Filter(Negate(is.null), results), fill = TRUE)
      fwrite(dt_chunk, out_file, append = !first_chunk)
      first_chunk <- FALSE
      gc()
      pb$tick()
    }
    
    # Final merge back
    delta_df <- fread(out_file) %>%
      select(PairsID, everything()) %>%
      left_join(PAIRS_METADATA, by = "PairsID")
    
    if (TYPE == "REAL_DATA") {
      delta_df[, coi_change := offset_naive_coi_Dx - offset_naive_coi_D0]
    } else {
      delta_df[, coi_change := Dx_nstrains - D0_nstrains]
    }
    
    write.csv(delta_df,
              paste0(site, "_", TYPE, ".csv"),
              row.names = FALSE)
    
    cat("=> done:", out_file, "\n\n")
  }
})


###


######################################################################################################
##### TARGET AUGMENTATION
######################################################################################################

FNRs <- readRDS(paste0(site, "_FNRs.RDS"))

if (nrow(FNRs)>0){
  
  
  PAIRS_METADATA <- readRDS(paste0("PAIRS_METADATA_", site, ".RDS"))
  PAIRS_GENOMIC <- readRDS(paste0("PAIRS_GENOMIC_" , site, ".RDS"))
  
  ## EXTRACT PAIRS THAT WILL BE FNRed
  FNR_pairs <- PAIRS_METADATA[PAIRS_METADATA$D0_nstrains %in% FNRs$offset_naive_coi | PAIRS_METADATA$Dx_nstrains %in% FNRs$offset_naive_coi ,]$PairsID
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
    mutate(sample = str_split(sample, "(?<=\\d)_(?=\\d)"))
  
  # clone data
  TRUTH <- read.csv(paste0("clones_genomic_data_", site, ".csv"))
  
  # # remove synthetic clones? NO, they also make the pairs!
  # TRUTH <- TRUTH[!grepl("__clone", TRUTH$sampleID),]
  
  
  # Get unique FNR values
  unique_fnr_values <- unique(unlist(FNRs$FNR))
  
  set.seed(420)
  
  library(purrr)
  
  # Apply FNR to TRUTH
  TRUTH_filtered_list <- map(unique_fnr_values, function(fnr) {
    TRUTH %>%
      group_by(sampleID) %>%
      group_modify(~ {
        n_to_keep <- ceiling((1 - fnr) * nrow(.x))
        .x[sample(nrow(.x), n_to_keep), , drop = FALSE]
      }) %>%
      ungroup()
  })
  
  # Name each element by its FNR
  names(TRUTH_filtered_list) <- as.character(unique_fnr_values)
  
  library(tibble)
  
  # check
  map_dfr(
    names(TRUTH_filtered_list),
    function(fnr_name) {
      TRUTH_filtered_list[[fnr_name]] %>%
        group_by(sampleID) %>%
        summarise(n_alleles = n(), .groups = "drop") %>%
        mutate(FNR = fnr_name)
    })
  
  # fill FNR vecgor with TRUTH for easier handling
  FNRs <- FNRs %>%
    mutate(
      FNR = map2(FNR, offset_naive_coi, ~ {
        len <- length(.x)
        if (len < .y) {
          c(.x, rep("TRUTH", .y - len))
        } else {
          .x
        }
      })
    )
  
  # ADD TRUTH DF TO THE FILTEREDLIST as element 1 FOR EASIER MANAGING
  TRUTH_filtered_list <- c(list(as.tibble(TRUTH)), TRUTH_filtered_list)
  names(TRUTH_filtered_list)[1] <- "TRUTH"
  
  
  # Initialize storage objects
  final_sample_df <- tibble()
  TRUTH_samples_by_row <- list()
  
  # Iterate over each row in FNRs
  for (numba in seq(1, nrow(FNRs), 1)) {
    
    # Get iteration data
    iteration <- FNRs[numba, ]
    
    print(paste0("PROCESSING: offset_naive_coi = ", iteration$offset_naive_coi, "; FNR = ", paste(unlist(iteration$FNR), collapse = "_")))
    
    
    # Subset FNR_pairs_METADATA based on offset_naive_coi
    coi_subset <- FNR_pairs_METADATA %>%
      group_by(PairsID) %>%
      filter(any(n_strains == iteration$offset_naive_coi)) %>%
      ungroup() %>%
      rename(offset_naive_coi = n_strains) %>%
      full_join(iteration, by = "offset_naive_coi")
    
    
    #####  #####  #####  #####  #####  #####
    
    # randomly sample 10% of the pairs for FNR scheme
    n_size <- round(length(unique(coi_subset$PairsID))*0.1,1)
    
    set.seed(42069)
    sample_pairsid <- sample(unique(coi_subset$PairsID), n_size)
    
    coi_subset <- coi_subset[coi_subset$PairsID %in% sample_pairsid,]
    #####  #####  #####  #####  #####  #####
    
    
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
            filter(sampleID == strain)
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
    all_results <- pmap(coi_subset, function(PairsID, time_point, sample, offset_naive_coi, FNR) {
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
  output_file <- paste0(site, "_TRUTH_FNR_TRAINING_DATA.csv")
  
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
  
  delta_metrics_df_final_FNR <- read.csv(paste0(site, "_TRUTH_FNR_TRAINING_DATA.csv"))
  
  delta_metrics_df_final_FNR <- delta_metrics_df_final_FNR %>%
    select(PairsID, everything())
  
  write.csv(delta_metrics_df_final_FNR, paste0(site, "_TRUTH_FNR_TRAINING_DATA.csv"), row.names = F)
  
  
  
  ### APPEND FNR DATA TO ORIGINAL DATA
  
  FEATURES <- left_join(delta_metrics_df_final_FNR, PAIRS_METADATA_FNR, by = "PairsID")
  
  FEATURES <- FEATURES %>%
    mutate(PairsID_clean = sub("_.*", "", PairsID))
  
  FEATURES$PairsID_clean <- as.numeric(FEATURES$PairsID_clean)
  
  ### MERGE WITH CLEAN DATA'S METADATA
  TRUTH_CLEAN<- read.csv(paste0(site, "_TRAINING_DATA.csv"))
  
  merged_df <- FEATURES %>%
    left_join(TRUTH_CLEAN[c("D0_sample", "Dx_sample", "D0_nstrains", "Dx_nstrains", "PairsID", "coi_change", "labels")], by = c("PairsID_clean" = "PairsID")) %>%
    select(-PairsID_clean)
  
  merged_df$PairsID <- as.character(merged_df$PairsID)
  
  #output final features with TRUTH + TRUTH_FNR data
  write.csv(merged_df, paste0(site, "_FNR_TRAINING_DATA.csv"), row.names = F)
  
  
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
  
  CLEAN_and_FNR_final <- unique(CLEAN_and_FNR_final)     # why duplicates at this point??????
  
  # OVERWRITE ORIHINAL TRAINING DATA WITH THE APPENDED ORIGINAL + FNR DATA.
  write.csv(CLEAN_and_FNR_final, paste0(site, "_TRAINING_DATA.csv"), row.names = F)
  
  
} else {
  
  print("ALl strain proportions from D0 failure samples are > 0.03. No target augmentation needed.")
  
}


######################################################################################################
##### ML MODEL
######################################################################################################

library(caret)    
library(dplyr)    
library(tidyr)    
library(ggplot2)  
library(broom)
library(glmnet)


### 1) IMPORT TRAINING AND REAL DATA ----------
TRAINING_DATA <- read.csv(paste0(site, "_TRAINING_DATA.csv"), row.names = 1) 
TRAINING_DATA$eCOI_pairs <- paste0(TRAINING_DATA$D0_nstrains, "__", TRAINING_DATA$Dx_nstrains) # Create a new variable by combining 'D0nstrains' and 'Dxnstrains'
LABELS <- data.frame(labels = TRAINING_DATA$labels)
LABELS$labels <- as.factor(LABELS$labels)

REAL_DATA <- read.csv(paste0(site, "_REAL_DATA.csv"), stringsAsFactors = FALSE, colClasses = c(NIDA1 = "character", NIDA2= "character")) 

features_to_use <- colnames(REAL_DATA)[!colnames(REAL_DATA) %in% c("PairsID", "NIDA1", "NIDA2", "pair_type", "IBD_estimate",
                                                                   "offset_naive_coi_D0", "offset_naive_coi_Dx", 
                                                                   "replacement_pattern_score", "locus_discordance_rate")]

corrplot::corrplot(cor(TRAINING_DATA %>% select(features_to_use), use = "complete.obs"), "pie")


### 2) SPLIT DATA ------------

set.seed(420)
# Step 1: Stratified Sampling by `pair_type`
train_indices <- createDataPartition(TRAINING_DATA$eCOI_pairs, p = 0.7, list = FALSE)
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

# # Set up 10-fold cross-validation
# ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
# 
# # Re-level factor so that "R" is the positive class
# df_train_IBD$label <- relevel(df_train_IBD$label, ref = "R")
# 
# # Train logistic regression model with cross-validation
# fit_IBD <- train(label ~ ., 
#                  data = df_train_IBD, 
#                  method = "glm", 
#                  family = "binomial", 
#                  trControl = ctrl, 
#                  metric = "ROC")
# 
# # print(fit_IBD)
# # 
# # fit_IBD$finalModel
# 
# ### feature importance
# coefs <- summary(fit_IBD$finalModel)$coefficients
# coefs_df <- as.data.frame(coefs)
# coefs_df$Variable <- rownames(coefs_df)
# 
# coefs_df <- coefs_df[coefs_df$Variable != "(Intercept)", ]
# 
# coefs_df$AbsEstimate <- log(abs(coefs_df$Estimate))
# 
# coefs_df <- coefs_df[order(coefs_df$AbsEstimate, decreasing = TRUE), ]
# 
# coefs_df$Color <- ifelse(coefs_df$Estimate > 0, "positive", "negative")
# 
# coefs_df$Significance <- ifelse(coefs_df$`Pr(>|z|)` < 0.001, "***",
#                                 ifelse(coefs_df$`Pr(>|z|)` < 0.01, "**",
#                                        ifelse(coefs_df$`Pr(>|z|)` < 0.05, "*", "")))
# 
# # Plot with significance indicators
# importance <- ggplot(coefs_df, aes(x = reorder(Variable, AbsEstimate), y = AbsEstimate, fill = Color)) +
#   geom_bar(stat = "identity") +
#   geom_text(aes(label = Significance, hjust = ifelse(AbsEstimate < 0, 1.2, -0.2))) +
#   coord_flip() +
#   scale_fill_manual(values = c("positive" = "steelblue", "negative" = "firebrick")) +
#   theme_minimal() +
#   labs(#title = "Feature Importance in Logistic Regression Model",
#     #subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
#     x = "Features",
#     y = "log(Absolute Coefficient Value)",
#     fill = "Coefficient Direction") +
#   theme(legend.position = "bottom")
# 
# ggsave(paste0("feat_importance_", site, ".png"), importance, bg = "white", dpi = 300, height = 5, width = 8)
# 
# 
# 
# 
# # Predict probabilities on the test (holdout) set
# preds_prob <- predict(fit_IBD, newdata = df_test_IBD, type = "prob")[, "R"]
# 
# # Define the range of decision_thresholds
# decision_thresholds <- seq(0, 1, by = 0.05)
# 
# # Initialize an empty results dataframe
# results <- data.frame(eCOI_pairs = character(),
#                       decision_threshold = numeric(),
#                       sensitivity = numeric(),
#                       specificity = numeric(),
#                       R_pairs = numeric(),
#                       NI_pairs = numeric(),
#                       stringsAsFactors = FALSE)
# 
# #decision_thresholds <- 0.5 # if wanting to use only 0.5 decision threshold for all pair types...
# 
# # Loop through each decision_threshold
# for (thresh in decision_thresholds) {
#   
#   # Convert probabilities to binary predictions at the current decision_threshold
#   preds <- ifelse(preds_prob >= thresh, "R", "NI")
#   
#   for (strain_comb in unique(TEST_META$eCOI_pairs)) {
#     
#     subset_indices <- TEST_META$eCOI_pairs == strain_comb
#     subset_TEST_labels <- TEST_labels[subset_indices]
#     
#     r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
#     ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
#     
#     # Get the corresponding predictions for the current subset
#     subset_preds <- preds[subset_indices]
#     
#     # Evaluate the confusion matrix for the current subset
#     cm <- confusionMatrix(as.factor(subset_preds), as.factor(subset_TEST_labels), positive = "R")
#     
#     sens <- cm$byClass["Sensitivity"]
#     spec <- cm$byClass["Specificity"]
#     
#     results <- rbind(results, data.frame(eCOI_pairs = strain_comb,
#                                          decision_threshold = thresh,
#                                          sensitivity = sens,
#                                          specificity = spec,
#                                          R_pairs = r,
#                                          NI_pairs = ni,
#                                          stringsAsFactors = FALSE))
#   }
# }
# 
# 
# # Select the best decision_threshold per eCOI_pair based on balance between sensitivity and specificity
# best_decision_thresholds <- results %>%
#   mutate(youden_j = sensitivity + specificity - 1) %>%  # Compute Youden’s J
#   group_by(eCOI_pairs) %>%
#   slice_max(youden_j) %>%  # Select rows with the largest Youden's J
#   #slice_min(decision_threshold, with_ties = FALSE) %>%  # If ties, pick the lowest decision_threshold
#   slice_min(abs(decision_threshold - 0.5), with_ties = FALSE) %>%  # Pick decision_threshold closest to 0.5
#   ungroup() 
# 
# 
# # Reshape data to long format for easy plotting
# best_decision_thresholds_long <- best_decision_thresholds %>%
#   select(eCOI_pairs, sensitivity, specificity) %>%
#   pivot_longer(cols = c(sensitivity, specificity), 
#                names_to = "Metric", 
#                values_to = "Value")
# 
# # Create bar plot
# metrics <- ggplot(best_decision_thresholds_long, aes(x = eCOI_pairs, y = Value, fill = Metric)) +
#   geom_bar(stat = "identity", position = "dodge") +  # Dodge separates bars for clarity
#   labs(title = "",
#        x = "Pair Type",
#        y = "Value") +
#   theme_minimal() +
#   scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
#   geom_hline(yintercept = 0.9, linetype = "solid", color = "black")+
#   geom_hline(yintercept = 0.80, linetype = "dashed", color = "black")
# 
# 
# # save model and model results
# saveRDS(fit_IBD, paste0("LogReg_model_", site, ".RDS"))
# write.csv(best_decision_thresholds, paste0("training_Results_LogReg_", site, ".csv"), row.names = F)
# ggsave(paste0(site, "_model_results_LogReg.png"), metrics, bg = "white", dpi = 300, height = 5, width = 7)


# --------------------------
# 1) SETUP
# --------------------------
ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

df_train_IBD$label <- relevel(df_train_IBD$label, ref = "R")

# --------------------------
# 2) TRAINING WITH GLMNET
# --------------------------
fit_IBD <- train(
  label ~ .,
  data = df_train_IBD,
  method = "glmnet",
  family = "binomial",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 10
)

print(fit_IBD)

# --------------------------
# 3) FEATURE IMPORTANCE (GLMNET)
# --------------------------
best_lambda <- fit_IBD$bestTune$lambda
coefs <- coef(fit_IBD$finalModel, s = best_lambda)
coefs_df <- as.data.frame(as.matrix(coefs))
colnames(coefs_df) <- "Estimate"
coefs_df$Variable <- rownames(coefs_df)

coefs_df <- coefs_df[coefs_df$Variable != "(Intercept)", ]
coefs_df$AbsEstimate <- log(abs(coefs_df$Estimate) + 1e-8)  # avoid log(0)
coefs_df$Color <- ifelse(coefs_df$Estimate > 0, "positive", "negative")

importance <- ggplot(coefs_df, aes(x = reorder(Variable, Estimate), y = Estimate, fill = Color)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("positive" = "steelblue", "negative" = "firebrick")) +
  theme_minimal() +
  labs(x = "Features", y = "Coefficient Value", fill = "Coefficient Direction") +
  theme(legend.position = "bottom")

ggsave(paste0("feat_importance_", site, ".png"), importance, bg = "white", dpi = 300, height = 5, width = 8)

# --------------------------
# 4) PREDICT ON TEST SET
# --------------------------
preds_prob <- predict(fit_IBD, newdata = df_test_IBD, type = "prob")[, "R"]

decision_thresholds <- seq(0, 1, by = 0.05)
results <- data.frame(eCOI_pairs = character(),
                      decision_threshold = numeric(),
                      sensitivity = numeric(),
                      specificity = numeric(),
                      R_pairs = numeric(),
                      NI_pairs = numeric(),
                      stringsAsFactors = FALSE)

for (thresh in decision_thresholds) {
  preds <- ifelse(preds_prob >= thresh, "R", "NI")
  
  for (strain_comb in unique(TEST_META$eCOI_pairs)) {
    subset_indices <- TEST_META$eCOI_pairs == strain_comb
    subset_TEST_labels <- TEST_labels[subset_indices]
    
    r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
    ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
    
    subset_preds <- preds[subset_indices]
    cm <- confusionMatrix(as.factor(subset_preds), as.factor(subset_TEST_labels), positive = "R")
    
    sens <- cm$byClass["Sensitivity"]
    spec <- cm$byClass["Specificity"]
    
    results <- rbind(results, data.frame(eCOI_pairs = strain_comb,
                                         decision_threshold = thresh,
                                         sensitivity = sens,
                                         specificity = spec,
                                         R_pairs = r,
                                         NI_pairs = ni,
                                         stringsAsFactors = FALSE))
  }
}

# --------------------------
# 5) BEST THRESHOLDS (YOUDEN'S J)
# --------------------------
best_decision_thresholds <- results %>%
  mutate(youden_j = sensitivity + specificity - 1) %>%
  group_by(eCOI_pairs) %>%
  slice_max(youden_j) %>%
  slice_min(abs(decision_threshold - 0.5), with_ties = FALSE) %>%
  ungroup()

best_decision_thresholds_long <- best_decision_thresholds %>%
  select(eCOI_pairs, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity),
               names_to = "Metric",
               values_to = "Value")

metrics <- ggplot(best_decision_thresholds_long, aes(x = eCOI_pairs, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Pair Type", y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0.9, linetype = "solid", color = "black") +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "black")

# --------------------------
# 6) SAVE RESULTS
# --------------------------
saveRDS(fit_IBD, paste0("GLMNET_model_", site, ".RDS"))
write.csv(best_decision_thresholds, paste0("training_Results_GLMNET_", site, ".csv"), row.names = F)
ggsave(paste0(site, "_model_results_GLMNET.png"), metrics, bg = "white", dpi = 300, height = 5, width = 7)



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
  geom_vline(data = best_decision_thresholds, aes(xintercept = decision_threshold), color = "red", linetype = "solid") +
  facet_wrap(~eCOI_pairs) +
  labs(title = "",
       x = "Decision Threshold",
       y = "Value") +
  theme_minimal()


##### 6) TEST MODEL ON REAL DATA USING IBD ONLY (BEST FEATURE) --------

#add threhold data
best_decision_thresholds <- best_decision_thresholds %>% rename(pair_type = eCOI_pairs)
REAL_DATA <- left_join(REAL_DATA, best_decision_thresholds[c("pair_type", "decision_threshold")], by = c("pair_type"))

REAL_DATA$prediction_prob <- NA
REAL_DATA$predictions <- NA

# Loop over each row in REAL_DATA to apply the corresponding decision_threshold
for (i in 1:nrow(REAL_DATA)) {
  
  decision_threshold <- REAL_DATA$decision_threshold[i]
  
  feats <- REAL_DATA %>% select(all_of(features_to_use))
  
  newdata <- feats[i, , drop = FALSE]
  
  prediction_prob <- predict(fit_IBD, newdata = newdata, type = "prob")[, "R"]
  
  # Classify using the current (best) decision_threshold (instead of 0.5)
  prediction_class <- ifelse(prediction_prob >= decision_threshold, "R", "NI")
  
  REAL_DATA$prediction_prob[i] <- prediction_prob
  REAL_DATA$predictions[i] <- prediction_class
}

REAL_DATA <- REAL_DATA %>% select(PairsID, NIDA1, NIDA2, pair_type, c(features_to_use), decision_threshold, prediction_prob, predictions) %>% arrange(PairsID)


## OUTPUT RESULTS 
write.csv(REAL_DATA, paste0(site, "_REAL_DATA_PREDICTIONS.csv"), row.names = F)
ggsave(paste0(site, "_REAL_DATA_THRESHOLDS_PLOT.png"), sens_spec_plot, bg = "white", dpi = 300, height = 9, width = 12)

paste0("n_amps used: ", length(unique(data_all$locus)))
