
library(dplyr)
library(tidyr)

# de los datos problema, necesito saber:

# 1) COI
# 2) PROPORTIONS

labcontrols_metadata<- readRDS("Lab_Controls_metadata.RDS")
labcontrols_genomic <- readRDS("Lab_Controls_genomic.RDS")

# 1) COI
# coi (in this case it's strain count)
combined_vec <- unique(paste0(labcontrols_metadata$D0_sample, "_____" ,labcontrols_metadata$D0_nstrains))
coi_df <- as.data.frame(do.call(rbind, strsplit(combined_vec, "_____")))
colnames(coi_df) <- c("NIDA", "COI")

# 2) PROPORTIONS
# for proportions, coi should be the numver of strains, so i'm picking all loci that have n_alleles == coi for each mix and then averaging out the norm.reads.locus

allele_count_per_locus <- labcontrols_genomic %>% 
  group_by(NIDA, run, locus) %>%
  summarise(n_alleles = length(unique(allele)))

# Merge the two data frames by NIDA
merged_df <- left_join(allele_count_per_locus, coi_df, by = "NIDA")
merged_df <- merged_df[merged_df$COI > 1,] # no COI = 1

# Filter rows where COI matches n_alleles: here are the loci that will be used for the proportions calculation for each mix (NIDA)
matched_df <- merged_df %>%
  filter(COI == n_alleles)

# Subset labcontrols_genomic by rowwise matching of NIDA and locus from matched_df
subset_labcontrols_genomic <- semi_join(labcontrols_genomic, matched_df, by = c("NIDA", "locus"))

# remove columns that generate redundancies
subset_labcontrols_genomic <- subset_labcontrols_genomic %>% select(-PairsID, -time_point)

# keep uniques
subset_labcontrols_genomic <- unique(subset_labcontrols_genomic)


# Step 1: Rank norm.reads.locus descendingly within each NIDA and locus
ranked_df <- subset_labcontrols_genomic %>%
  group_by(NIDA, locus) %>%
  arrange(desc(norm.reads.locus), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Step 2: Compute the average of norm.reads.locus by rank per NIDA
average_by_rank <- ranked_df %>%
  group_by(NIDA, rank) %>%
  summarise(avg_norm_reads = mean(norm.reads.locus, na.rm = TRUE), .groups = "drop") %>%
  arrange(NIDA, rank)

average_by_rank_wide <- average_by_rank %>%
  pivot_wider(
    names_from = rank,
    values_from = avg_norm_reads,
    names_prefix = "rank_"
  )


## subset wors with any clone in < proportion of < 0.06
subset_below_threshold <- average_by_rank_wide %>%
  filter(if_any(starts_with("rank_"), ~ . <= 0.1))


# readd coi
subset_below_threshold <- left_join(subset_below_threshold, coi_df[c("NIDA", "COI")], by = "NIDA")


# # Define FNR thresholds: 3 TO 4 WAS PULLED OUT OF MY ASS
# FNR_minclone_2_to_3 <- 0.14
# FNR_minclone_3_to_4 <- 0.08
# FNR_minclone_4_to_5 <- 0.03
# 
# subset_below_threshold <- subset_below_threshold %>%
#   rowwise() %>%
#   mutate(
#     last_rank_value = last(na.omit(c_across(starts_with("rank_")))),
#     FNR = case_when(
#       last_rank_value >= 0 & last_rank_value < 0.03 ~ FNR_minclone_2_to_3,
#       last_rank_value >= 0.03 & last_rank_value < 0.04 ~ FNR_minclone_3_to_4,
#       last_rank_value >= 0.04 & last_rank_value < 0.05 ~ FNR_minclone_4_to_5,
#       TRUE ~ NA_real_
#     )
#   ) %>%
#   ungroup()

## FINALLY, FALSE NEGATIVE RATES TO APPLY:
write.csv(subset_below_threshold, "FNRs.csv", row.names = F )

# 
# 
# ## FINALLY, FALSE NEGATIVE RATES TO APPLY:
# FNRs_vec <- unique(paste0(subset_below_threshold$COI, "__", subset_below_threshold$FNR))
# FNRs <- as.data.frame(do.call(rbind, strsplit(FNRs_vec, "__")))
# colnames(FNRs) <- c("COI", "FNR")
# 
# FNRs

