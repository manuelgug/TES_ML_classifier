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

    sample_size <- initial_sample_size * (n_strains - i + 5) # Calculate the sample size based on the iteration

    strain_mixes_subsampled[[mix_name]] <- strain_mixes[[mix_name]] %>% 
      sample_n(min(nrow(strain_mixes[[mix_name]]), sample_size))
  }
}

# check
lapply(strain_mixes_subsampled, nrow)



###### 4) OUTPUT MIXES------

saveRDS(strain_mixes_subsampled, "mixes_list.RDS")

