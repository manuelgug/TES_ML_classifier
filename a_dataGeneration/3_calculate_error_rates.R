
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)

# inputs: 
# 1) TRUTH_CLEAN.csv
# 2) LAB_CONTROLS_CLEAN.csv 
# 3) control_metadata.csv (strain content of each sampleID)

TRUTH <- read.csv("TRUTH_CLEAN.csv")
TRUTH$allele <- gsub("__", "___", TRUTH$allele) # needed an extra underscore
LAB_CONTROLS <- read.csv("LAB_CONTROLS_CLEAN.csv")

# #format sampleID
# LAB_CONTROLS <- LAB_CONTROLS %>%
#   mutate(sampleID = if_else(!grepl("^N", sampleID), paste0("N", sampleID), sampleID))

# make the strain content a vector
LAB_CONTROLS_METADATA <- read.csv("control_metadata.csv")
LAB_CONTROLS_METADATA <- LAB_CONTROLS_METADATA %>%
  mutate(
    strain_content = pmap(
      select(., matches("^strain[1-9]$")),
      ~ na.omit(c(...)) %>% trimws() %>% .[. != ""]
    )
  ) %>% 
  select(NIDA, strain_content, PROPORTIONS) %>%
  rename(sampleID=NIDA)


#check:
# Get unique IDs
meta_ids <- unique(LAB_CONTROLS_METADATA$sampleID)
control_ids <- unique(LAB_CONTROLS$sampleID)

# Intersection
shared_ids <- intersect(meta_ids, control_ids)

# Exclusive to metadata: THESE ARE THE ONES FROM MANHICA AND ISG CONTAMINATED RUNS I DIDN'T USE (MANUALLY REMOVED FROM THE FILTERED_DATA_FROM_CLUSTER FOLDER); ALSO TRUTH STRAINS
only_in_meta <- setdiff(meta_ids, control_ids)

# Exclusive to controls: SHOULD BE ZERO
only_in_controls <- setdiff(control_ids, meta_ids)
only_in_controls


# MERGE
LAB_CONTROLS <- merge(LAB_CONTROLS, LAB_CONTROLS_METADATA, by = "sampleID")
length(unique(LAB_CONTROLS$sampleID))



############## calculate FPR and FNR ###############

SAMPLE_IDS <- unique(LAB_CONTROLS$sampleID)

# Initialize empty list to store results
results_list <- list()

for (i in seq_along(SAMPLE_IDS)) {
  SID <- SAMPLE_IDS[i]
  
  # Subset LAB_CONTROLS for the current sample ID
  subset_data <- LAB_CONTROLS[LAB_CONTROLS$sampleID == SID, ]
  
  # Extract strain vector
  strains <- unlist(unique(subset_data$strain_content))
  
  # Subset TRUTH for the relevant strains
  subset_truth <- TRUTH[TRUTH$Strain %in% strains, ]
  
  # extract shared loci for a fair comparison
  shared_loci <- intersect(subset_data$locus, subset_truth$locus)
  subset_data <- subset_data[subset_data$locus %in% shared_loci,]
  subset_truth <- subset_truth[subset_truth$locus %in% shared_loci,]
  
  # Get allele sets
  observed_alleles <- subset_data$allele
  expected_alleles <- unique(subset_truth$allele)
  
  # Calculate metrics: proportion of alleles relative to the total alleles shared (expected)
  FP <- round(length(setdiff(observed_alleles, expected_alleles)) / length(expected_alleles), 2)
  FN <- round(length(setdiff(expected_alleles, observed_alleles)) / length(expected_alleles), 2)
  COI <- length(strains)
  
  # Save as a named list
  results_list[[i]] <- data.frame(
    sampleID = SID,
    strains = paste(strains, collapse = ","),
    COI = COI,
    FP = FP,
    FN = FN,
    stringsAsFactors = FALSE
  )
}

# Combine all into one data frame
allele_evaluation <- do.call(rbind, results_list)

# add strain proportions
allele_evaluation <- left_join(allele_evaluation, LAB_CONTROLS[c("sampleID", "PROPORTIONS")], by = "sampleID")

allele_evaluation <- unique(allele_evaluation)

# add minoruty clone column
allele_evaluation <- allele_evaluation %>%
  mutate(
    minority_clone = if_else(
      str_detect(PROPORTIONS, "_"),
      str_extract(PROPORTIONS, "[^_]+$"),  # Extract content after last "_"
      PROPORTIONS                          # If no "_", use entire string
    )
  )

allele_evaluation$minority_clone <- as.numeric(allele_evaluation$minority_clone)


plot_fp_fn_histograms <- function(data, group_var) {
  group_var_sym <- rlang::sym(group_var)
  
  # Calculate FP stats
  fp_stats <- data %>%
    group_by(!!group_var_sym) %>%
    summarise(
      mean_FP = mean(FP, na.rm = TRUE),
      median_FP = median(FP, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0("Median: ", median_FP)
    )
  
  # Calculate FN stats
  fn_stats <- data %>%
    group_by(!!group_var_sym) %>%
    summarise(
      mean_FN = mean(FN, na.rm = TRUE),
      median_FN = median(FN, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0("Median: ", median_FN)
    )
  
  # FP Plot
  fp_plot <- ggplot(data, aes(x = FP)) +
    geom_histogram(fill = "#4E79A7", color = "black", bins = 20) +
    facet_wrap(vars(!!group_var_sym), scales = "free_y") +
    geom_text(
      data = fp_stats,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.1, vjust = 1.5,
      inherit.aes = FALSE,
      size = 3.5
    ) +
    labs(
      title = paste("False Positives by", group_var),
      x = "False Positives (FP)",
      y = "Count"
    ) +
    theme_minimal()
  
  # FN Plot
  fn_plot <- ggplot(data, aes(x = FN)) +
    geom_histogram(fill = "#4E79A7", color = "black", bins = 20) +
    facet_wrap(vars(!!group_var_sym), scales = "free_y") +
    geom_text(
      data = fn_stats,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.1, vjust = 1.5,
      inherit.aes = FALSE,
      size = 3.5
    ) +
    labs(
      title = paste("False Negatives by", group_var),
      x = "False Negatives (FN)",
      y = "Count"
    ) +
    theme_minimal()
  
  # Save plots
  ggsave(paste0("FP_by_", group_var, "_histogram.png"), fp_plot, dpi = 300, width = 8, height = 6, bg = "white")
  ggsave(paste0("FN_by_", group_var, "_histogram.png"), fn_plot, dpi = 300, width = 8, height = 6, bg = "white")
}

# Example usage:
plot_fp_fn_histograms(allele_evaluation, "PROPORTIONS")
plot_fp_fn_histograms(allele_evaluation, "minority_clone")





## for each lab_control mix
## extract the clones present in the lab_control mix
## generate a mix of the truth samples using the clones from the lab_control mix
## calculate False positive rate: alleles found in the lab_controls that are not in the truth mix
## calculate false negative rate: alleles found in the truth mix that are not in the lab_controls mix
## append results to a df with col0= name of the lab_controls mix sample, col1= vector of strains present in the lab_controls mix, col2= n_strains present in the lab_controls mix, col3= FNR, col4= FPR

## ((samples with HIGH FPR are probably mislabeled, remove and rerun from script 1_clean_truth_and_lab_controls_data.R))

## plot the distribution of FNR and calculate summary stats; STRATIFIED BY m_Strains