# ------------------------------------------------------------------------------
# 0. Load Libraries
# ------------------------------------------------------------------------------
library(parallel)
library(progress)
library(dplyr)
library(data.table)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)

# ------------------------------------------------------------------------------
# 1. Parameters & Inputs
# ------------------------------------------------------------------------------
site <- "Tete"
DATA_TYPE <- "TRAINING_DATA"  # for naming consistency

# Columns and labels
selected_columns <- c("PairsID", "labels")
labels <- fread(paste0(site, "_TRAINING_DATA.csv"), select = selected_columns)

# Genomic and metadata
allele_data <- readRDS(paste0("PAIRS_GENOMIC_", site, ".RDS"))
METADATA    <- readRDS(paste0("PAIRS_METADATA_", site, ".RDS"))

# Keep only needed columns
allele_data <- allele_data %>%
  select(PairsID, time_point, locus, allele) %>%
  as.data.table()

gc()

# Input features testing data for subsetting
test_data <- readRDS(paste0(site, "_test_data.RDS"))
TEST_META <- test_data %>%
  mutate(PairsID = as.character(PairsID),
         labels = factor(labels, levels = c("NI", "R"))) %>%
  select(PairsID, labels, pair_type = eCOI_pairs) %>%  # Renamed to pair_type for consistency
  as.data.table()

# Subset test data for fair comparison 
METADATA <- METADATA[METADATA$PairsID %in% TEST_META$PairsID,]
allele_data <- allele_data[allele_data$PairsID %in% TEST_META$PairsID,]

# ------------------------------------------------------------------------------
# 2. Pivot Allele Data for Matching Comparison (in chunks)
# ------------------------------------------------------------------------------
num_cores   <- detectCores()
chunk_size  <- 20000
output_file <- paste0(site, "_allele_data_pivoted_for_allele_matching_comparison.csv")

pb <- progress_bar$new(
  format = "Processing [:bar] :percent :elapsedfull",
  total  = ceiling(length(unique(allele_data$PairsID)) / chunk_size),
  width  = 100
)

process_pair <- function(pair_id) {
  sub <- allele_data[PairsID == pair_id]
  sub <- sub[, .(alleles = list(unique(allele))), by = .(PairsID, time_point, locus)]
  wide <- dcast(sub, PairsID + locus ~ time_point, value.var = "alleles", fill = NA)
  setnames(wide, c("D0", "Dx"), c("D0_alleles", "Dx_alleles"))
  wide
}

write_chunk <- function(df, path, append = FALSE) {
  fwrite(df, path, append = append)
}

unique_pairs <- unique(allele_data$PairsID)

for (i in seq(1, length(unique_pairs), by = chunk_size)) {
  chunk_ids <- unique_pairs[i:min(i + chunk_size - 1, length(unique_pairs))]
  cat("Processing PairsID", first(chunk_ids), "to", last(chunk_ids), "...\n")
  
  results <- mclapply(chunk_ids, process_pair, mc.cores = num_cores)
  results <- Filter(Negate(is.null), results)
  
  if (length(results) > 0) {
    df_chunk <- rbindlist(results, fill = TRUE)
    write_chunk(df_chunk, output_file, append = (i > 1))
  }
  
  pb$tick()
}
gc()

# Reload pivoted data
allele_data <- fread(output_file) %>% as.data.table()

# ------------------------------------------------------------------------------
# 3. Compute Allele Matching & WHO Criteria
# ------------------------------------------------------------------------------
allele_data[, alleles_match := as.integer(map2_lgl(D0_alleles, Dx_alleles, ~ any(.x %in% .y)))]

# Remove allele lists to save memory
allele_data[, c("D0_alleles", "Dx_alleles") := NULL]
gc()

# Number of loci for WHO criteria
n_loci      <- length(unique(allele_data$locus))
two_thirds  <- round(n_loci * 2/3)

# WHO 2/3 and 3/3 classification
who_classification <- allele_data %>%
  group_by(PairsID) %>%
  summarise(
    classif_2_3 = ifelse(sum(alleles_match) >= two_thirds, "R", "NI"),
    classif_3_3 = ifelse(sum(alleles_match) == n_loci,      "R", "NI")
  ) %>%
  ungroup()

# xN/N classification for each x = 1…n_loci
xN_classification <- allele_data %>%
  group_by(PairsID) %>%
  summarise(alleles_match = sum(alleles_match)) %>%
  ungroup()

for (x in seq_len(n_loci)) {
  col <- paste0("xN_", x, "_", n_loci)
  xN_classification[[col]] <- ifelse(xN_classification$alleles_match >= x, "R", "NI")
}

xN_classification$alleles_match <- NULL
gc()

# ------------------------------------------------------------------------------
# 4. Merge Classifications & Add Metadata
# ------------------------------------------------------------------------------
all_classifs <- merge(labels, xN_classification, by = "PairsID")
all_classifs <- left_join(all_classifs, METADATA[, c("PairsID", "pair_type")], by = "PairsID")
all_classifs <- as.data.frame(all_classifs)

write.csv(all_classifs,
          paste0(site, "_allele_matching_classifications_TRAINING_DATA.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Compute Sensitivity, Specificity & Youden’s J by pair_type
# ------------------------------------------------------------------------------
# Helper function
calc_metrics <- function(pred, actual) {
  TP <- sum(pred == "R"  & actual == "R")
  TN <- sum(pred == "NI" & actual == "NI")
  FP <- sum(pred == "R"  & actual == "NI")
  FN <- sum(pred == "NI" & actual == "R")
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  youd <- sens + spec - 1
  c(Sensitivity = sens, Specificity = spec, YoudensJ = youd)
}

# Reshape to long for group-wise metrics
long_df <- all_classifs %>%
  pivot_longer(cols = -c(PairsID, labels, pair_type),
               names_to = "Classifier", values_to = "Prediction")

results_grouped <- long_df %>%
  group_by(pair_type, Classifier) %>%
  summarise(
    list_Prediction = list(Prediction),
    list_labels     = list(labels),
    .groups = "drop"
  ) %>%
  mutate(metrics = map2(list_Prediction, list_labels, ~ calc_metrics(.x, .y))) %>%
  unnest_wider(metrics)

# Reshape for plotting
results_long <- results_grouped %>%
  pivot_longer(cols = c(Sensitivity, Specificity, YoudensJ),
               names_to = "Metric", values_to = "Value") %>%
  mutate(
    Order = as.numeric(str_extract(Classifier, "(?<=xN_)[0-9]+"))
  ) %>%
  arrange(pair_type, Order) %>%
  group_by(pair_type) %>%
  mutate(Classifier = factor(Classifier, levels = unique(Classifier))) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 6. Line Plot: Sensitivity & Specificity Faceted by pair_type
# ------------------------------------------------------------------------------
plot_data <- results_long %>% filter(Metric %in% c("Sensitivity", "Specificity"))

best_classifiers <- results_long %>%
  filter(Metric == "YoudensJ") %>%
  group_by(pair_type) %>%
  filter(Value == max(Value)) %>%
  select(pair_type, Classifier) %>%
  distinct() %>%
  mutate(Order = as.numeric(str_extract(Classifier, "(?<=xN_)[0-9]+")))

thresholds <- ggplot(plot_data, aes(x = Order, y = Value, linetype = Metric)) +
  geom_line() +
  geom_vline(data = best_classifiers, aes(xintercept = Order),
             color = "red", linetype = "solid") +
  facet_wrap(~pair_type, scales = "free_x") +
  labs(x = "Decision Threshold Index (1/n loci to classsify as R)", y = "Metric Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

thresholds

ggsave(filename = paste0(site, "_allele_matching_THRESHOLDS_PLOT.png"),
       plot = thresholds,
       bg = "white", dpi = 300, height = 9, width = 12)

# ------------------------------------------------------------------------------
# 7. Bar Plot: Best Sensitivity & Specificity per pair_type
# ------------------------------------------------------------------------------
best_metrics <- results_long %>%
  filter(Metric == "YoudensJ") %>%
  group_by(pair_type) %>%
  slice_max(Value) %>%
  left_join(results_long, by = c("pair_type", "Classifier")) %>%
  filter(Metric.y %in% c("Sensitivity", "Specificity")) %>%
  select(pair_type, Metric = Metric.y, Value = Value.y)

metrics_plot <- ggplot(best_metrics, aes(x = pair_type, y = Value, fill = Metric)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Sensitivity" = "#008080",
                               "Specificity" = "orange")) +
  geom_hline(yintercept = 0.9, linetype = "solid") +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  labs(x = "Pair Type", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

metrics_plot

ggsave(filename = paste0(site, "_allele_matching_results.png"),
       plot = metrics_plot,
       bg = "white", dpi = 300, height = 5, width = 7)
