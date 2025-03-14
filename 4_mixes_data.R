library(dplyr)
library(tidyr)


clones_genomic <- read.csv("clones_genomic_data.csv")
metadata_updated <- read.csv("metadata_updated.csv", stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))



###### 1) IDENTIFY THE PAIRS DATA REQUIREMENTS ------

min_coi <- round(min(metadata_updated$post_effective_coi_med))
max_coi <- round(max(metadata_updated$post_effective_coi_med))



###### 2) CREATE ALL POSSIBLE MIXES FROM 1 TO 5 STRAINS -----

# Load unique strains into vector
nidas <- unique(clones_genomic$sampleID)
nidas_df <- data.frame(strain_1 = nidas)

# Function to create combinations of a specific size and convert to a data frame
create_combinations_df <- function(vector, combination_size) {
  comb_matrix <- combn(vector, combination_size)
  comb_df <- as.data.frame(t(comb_matrix))
  
  # Rename columns dynamically based on the combination size
  colnames(comb_df) <- paste0("strain_", 1:combination_size)
  
  return(comb_df)
}

# Create a list to store each combination set
strain_mixes <- lapply(seq(min_coi + 1, max_coi), function(n) {
  create_combinations_df(nidas, n)
})

# Optionally, name the list elements dynamically
names(strain_mixes) <- paste0("mix_", seq(min_coi + 1, max_coi))

# Add the nidas data frame as the first element in the list
strain_mixes <- c(list("mix_1" = nidas_df), strain_mixes)
strain_mixes <- lapply(strain_mixes, as.data.frame)

# Order the list alphabetically by the names
strain_mixes <- strain_mixes[sort(names(strain_mixes))]



####### 3) SUBSAMPLE THE MIXES -------

# Define the initial sample size from the first data frame
initial_sample_size <- nrow(strain_mixes$mix_1)
sample_size <- initial_sample_size


# Loop over each data frame in strain_mixes and apply the subsampling
strain_mixes_subsampled <- list()

n_strains<- length(strain_mixes)

for (i in 1:n_strains) {

  mix_name <- paste0("mix_", i)
  

  if (mix_name %in% names(strain_mixes)) {

    sample_size <- initial_sample_size  # initial_sample_size would give a more balanced training data. however, could add " * (n_strains - i + n_strains)" to get more mixes of lower strain content or maybe other approach if needed

    strain_mixes_subsampled[[mix_name]] <- strain_mixes[[mix_name]] %>% 
      sample_n(min(nrow(strain_mixes[[mix_name]]), sample_size))
  }
}

# check
lapply(strain_mixes_subsampled, nrow)



# ###### 4) OUTPUT MIXES------
# 
# saveRDS(strain_mixes_subsampled, "mixes_list.RDS")


##### 4) CREATE PAIRS ------

# Create empty df dynamically for metadata
strain_cols <- paste0("strain", seq(min_coi, max_coi))
strain_prop_cols <- paste0("strain_prop", seq(min_coi, max_coi))

MIXES_METADATA <- data.frame(
  NIDA = character(0),
  as.data.frame(matrix(ncol = length(strain_cols), nrow = 0, dimnames = list(NULL, strain_cols))),
  as.data.frame(matrix(ncol = length(strain_prop_cols), nrow = 0, dimnames = list(NULL, strain_prop_cols)))
)

MIXES_METADATA[strain_cols] <- lapply(MIXES_METADATA[strain_cols], as.character)
MIXES_METADATA[strain_prop_cols] <- lapply(MIXES_METADATA[strain_prop_cols], as.numeric)

#create genomic df
MIXES_GENOMIC <- data.frame(mixID=character(0), locus=character(0), allele=character(0), read_counts=numeric(0), norm.reads.locus=numeric(0), n.alleles=numeric(0))


# Loop through each mix type
metadata_list <- list()
genomic_list <- list()

for (mix_num in min_coi:max_coi) {  
  
  # Dynamically refer to each mix column in `strain_mixes_subsampled`
  mix_column <- paste0("mix_", mix_num)
  
  print(paste0("Processing mixes of ", mix_column))
  
  for (numbar in 1:nrow(strain_mixes_subsampled[[mix_column]])){
    
    # Subset genomic data only once
    selected_genomic <- clones_genomic[clones_genomic$sampleID %in% strain_mixes_subsampled[[mix_column]][numbar,],]
    total_reads <- sum(selected_genomic$reads)
    
    # Calculate strain proportions
    prop_strains <- selected_genomic %>%
      group_by(sampleID) %>%
      summarise(prop_reads = sum(reads) / total_reads, .groups = 'drop')
    
    # Get unique strains and mix ID
    nidas <- unique(selected_genomic$sampleID)
    mixID <- paste0("mix", length(nidas), "_ID", numbar)
    
    # Summarize genomic reads data for the combination
    rows_genomic <- selected_genomic %>%
      group_by(locus, allele) %>%
      summarise(read_counts = sum(reads), .groups = 'drop') %>%
      mutate(mixID = mixID) %>%
      group_by(mixID, locus) %>%
      mutate(norm.reads.locus = read_counts / sum(read_counts),
             n.alleles = n_distinct(allele)) %>%
      ungroup()
    
    ROW_METADATA <- c(
      mixID,
      nidas,
      rep("", max_coi - length(nidas)),
      prop_strains$prop_reads,
      rep("", max_coi - length(nidas))
    )
    
    ROW_METADATA <- as.data.frame(t(ROW_METADATA), stringsAsFactors = FALSE)
    colnames(ROW_METADATA) <- colnames(MIXES_METADATA)
    
    metadata_list[[length(metadata_list) + 1]] <- ROW_METADATA
    genomic_list[[length(genomic_list) + 1]] <- rows_genomic
  }
}

MIXES_METADATA <- do.call(rbind, metadata_list)
MIXES_GENOMIC <- do.call(rbind, genomic_list)


## CHECK: all mixes present in both metadata and genomic data?
nrow(MIXES_METADATA) == length(unique(MIXES_GENOMIC$mixID))
all(MIXES_METADATA$NIDA %in% unique(MIXES_GENOMIC$mixID))


saveRDS(MIXES_METADATA, "MIXES_METADATA.RDS")
saveRDS(MIXES_GENOMIC, "MIXES_GENOMIC.RDS")



###### 3) EDA OF MIXES --------

#allelic richness per mix
alleles_per_mix <- MIXES_GENOMIC %>%
  group_by(mixID) %>%
  summarize(alleles_per_mix = length(unique(allele)), .groups = "drop")

alleles_per_mix <- alleles_per_mix %>%
  separate(mixID, into = c("mix_type", "ID"), sep = "_", remove = FALSE)

allele_counts_per_mix <- ggplot(alleles_per_mix, aes(x = alleles_per_mix, fill = mix_type)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.5) +
  labs(
    title = "",
    x = "Number of Unique  Alleles",
    y = "Count",
    fill = "Mix Type"
  ) +
  theme_minimal()+
  xlim(0,NA)

allele_counts_per_mix

ggsave("mixes_EDA.png", allele_counts_per_mix, width = 9, height = 6, dpi = 300, bg = "white")

