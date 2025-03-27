

library(dplyr)
library(tidyr)
library(ggplot2)


site <- "Cabo_Delgado"
top_n_amps <- 50 


clones_genomic <- read.csv(paste0("clones_genomic_data_",site,"_top_",top_n_amps,"_amps.csv"), stringsAsFactors = FALSE, colClasses = c(sampleID = "character"))
metadata_updated <- read.csv(paste0("metadata_updated_", site, "_top_", top_n_amps,"_amps.csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))




##### 1) DETERMINE PAIRS OF MIXES BASED ON THE DATA'S COIs ------

metadata_updated$post_effective_coi_med <- round(metadata_updated$post_effective_coi_med)

metadata_updated_wide <- metadata_updated %>%
  pivot_wider(
    id_cols = PairsID, 
    names_from = time_point, 
    values_from = c(NIDA, post_effective_coi_med), 
    names_glue = "{.value}_{time_point}"
  )

unique_combos <- metadata_updated_wide %>%
  distinct(post_effective_coi_med_D0, post_effective_coi_med_Dx)

unique_combos$post_effective_coi_med_D0 <- paste0("mix", unique_combos$post_effective_coi_med_D0)
unique_combos$post_effective_coi_med_Dx <- paste0("mix", unique_combos$post_effective_coi_med_Dx)

unique_combos <- unique_combos %>% arrange(post_effective_coi_med_D0, post_effective_coi_med_Dx)


## 3) CREATE PAIRS------------------------------

MIXES_METADATA <- readRDS(paste0("MIXES_METADATA_",site,"_top_", top_n_amps,"amps_", N_CLONES, "_clones.RDS"))
MIXES_GENOMIC <- readRDS(paste0("MIXES_GENOMIC_" ,site,"_top_", top_n_amps, "amps_", N_CLONES, "_clones.RDS"))

nidas_all <- MIXES_METADATA$NIDA

# Create all possible pairs of sampleIDs
pairs_df <- expand.grid(NIDA1 = nidas_all, NIDA2 = nidas_all, stringsAsFactors = FALSE)

# # remove self pairs
# pairs_df <- pairs_df[pairs_df$NIDA1 != pairs_df$NIDA2,]

# subsample pairs based on unique_combos
pairs_df <- pairs_df[rowSums(sapply(1:nrow(unique_combos), function(i) {
  grepl(unique_combos$post_effective_coi_med_D0[i], pairs_df$NIDA1) & 
    grepl(unique_combos$post_effective_coi_med_Dx[i], pairs_df$NIDA2)
})) > 0, ]


# Create a unique PairsID for each pair
pairs_df <- pairs_df %>%
  mutate(PairsID = row_number()) %>%
  pivot_longer(cols = c(NIDA1, NIDA2), names_to = "time_point", values_to = "NIDA")

# Assign time_point based on the column names
pairs_df <- pairs_df %>%
  mutate(time_point = ifelse(time_point == "NIDA1", "D0", "Dx"))

# Rearrange columns to match the desired format
pairs_df <- pairs_df %>%
  select(PairsID, NIDA, time_point)

#merge control_pairs and CONTROLS_ALL by NIDA
colnames(MIXES_GENOMIC)[colnames(MIXES_GENOMIC) == "mixID"] <- "NIDA"
merged_dfs <-  left_join(pairs_df, MIXES_GENOMIC, by = "NIDA")

length(unique(merged_dfs$PairsID))

gc()

saveRDS(merged_dfs, paste0("PAIRS_GENOMIC_",site,"_top_",top_n_amps,"_amps.RDS"))


############################################
### COMPARE ALLELE CONTENT OF PAIRS

# compare alleles shared between mixes:::   HASTA AQUÍ BIEN!!!!!!!!!! 13/MARZO/2025
alleles <- merged_dfs %>%
  group_by(PairsID,NIDA, time_point) %>%
  summarize(alleles = list(allele))

alleles <- alleles %>% arrange(PairsID, time_point)

gc()

# Calculate the percentage of shared alleles for each pair
alleles_shared_prop <- alleles%>%
  group_by(PairsID) %>%
  summarize(
    
    NIDA1 = NIDA[1],
    NIDA2 = NIDA[2],
    
    # Calculate intersection and union for shared allele percentage
    shared_count = length(intersect(alleles[[1]], alleles[[2]])),
    union_count = length(union(alleles[[1]], alleles[[2]])),
    shared_prop = (shared_count / union_count)
  ) %>%
  ungroup()


alleles_shared_prop <- alleles_shared_prop %>%
  mutate(
    NIDA1_trimmed = sub("_.*", "", NIDA1),
    NIDA2_trimmed = sub("_.*", "", NIDA2),
    # Paste the trimmed values together
    pair_type = paste0(NIDA1_trimmed, "_", NIDA2_trimmed)
  ) %>%
  select(-NIDA1_trimmed, -NIDA2_trimmed)


shared_alleles_distribution <- ggplot(alleles_shared_prop, aes(x = pair_type, y = shared_prop)) +
  geom_violin(fill = "lightblue", color = "black", alpha = 0.5) +  # Adds the distribution shape
  geom_boxplot(width = 0.2, color = "black", outliers = F, outlier.color = "black", outlier.shape = 16) + 
  labs(
    title = "",
    x = "Pair Type",
    y = "Proportion of Alleles Shared between Mixes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

shared_alleles_distribution

ggsave(paste0("Prop_Alleles_Shared_by_Mix_Type_",site,"_top_",top_n_amps,"_amps.png"), shared_alleles_distribution, height = 7, width = 10, bg = "white", dpi = 300)



#############################################
## LABEL DATA 

# Replace all "" with NA in the dataframe
MIXES_METADATA <- MIXES_METADATA %>%
  mutate(across(where(is.character), ~ na_if(., "")))

# create strain vectors
MIXES_METADATA <- MIXES_METADATA %>%
  rowwise() %>%
  mutate(
    STRAINS = list(na.omit(
      c_across(matches("^strain\\d?$"))  # Matches "strain" or "strain" + single digit, excludes "strain_"
    ))
  ) %>%
  ungroup()

#add strain content of each control
PAIRS <- merge(pairs_df, MIXES_METADATA[c("NIDA", "STRAINS")], by = "NIDA")
PAIRS <- PAIRS%>%
  arrange(PairsID, time_point)


#label data 
labels <- PAIRS %>%
  group_by(PairsID) %>%
  summarise(
    
    # NEW_CONDITION: If none of the strains from Dx are found in D0, labels will be "NI", else R (same as If any strain from Dx is found in D0, labels will be "R", else NI.)
    labels = ifelse(any(unlist(STRAINS[time_point == "Dx"]) %in% unlist(STRAINS[time_point == "D0"])), "R", "NI"),  
    
  )


#add metadata:
PAIRS <- PAIRS %>%
  rowwise() %>%
  mutate(nstrains = length(STRAINS), #nstrains
         sample = paste(STRAINS, collapse= "_")) #strain content

PAIRS_wide <- PAIRS %>%
  pivot_wider(
    names_from = time_point,  # Reshape based on the time_point
    values_from = c(sample, nstrains),  # Use these columns for D0 and Dx
    names_glue = "{time_point}_{.value}"  # Create names like D0_sample, Dx_sample, etc.
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

#add labels
PAIRS_metadata <- inner_join(PAIRS_metadata, labels, by = "PairsID")
#add shared alleles for stratification
PAIRS_metadata <- inner_join(PAIRS_metadata, alleles_shared_prop, by = "PairsID")


saveRDS(PAIRS_metadata, paste0("PAIRS_METADATA_",site,"_top_",top_n_amps,"_amps.RDS"))



####################################
### PARIS METADATA mini EDA

PAIRS_metadata

PAIRS_summary <- PAIRS_metadata %>%
  group_by(pair_type) %>%
  summarise(
    
    n_pairs = length(PairsID),
    
    # Proportion of each label
    NI_prop = mean(labels == "NI"),  # proportion of "NI" labels in pair_type
    R_prop = mean(labels == "R"),  # repeat for each specific label as needed
    
    # Median shared_pct
    median_shared_prop = median(shared_prop, na.rm = TRUE),
    
    #NUMBERS
    NI_size = NI_prop * n_pairs,
    R_size = R_prop * n_pairs
    
  )

# add n of unique mixes to summary
unique_mixes_summary <- PAIRS_metadata %>%
  group_by(pair_type) %>%
  summarise(unique_mixes = n_distinct(Dx_sample))

PAIRS_summary <- PAIRS_summary %>%
  left_join(unique_mixes_summary, by = "pair_type")

PAIRS_summary

write.csv(PAIRS_summary, paste0("PAIRS_SUMMARY_",site,"_top_",top_n_amps,"_amps.csv"), row.names = F)

