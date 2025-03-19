
library(parallel)
library(progress)
library(dplyr)
library(data.table)
library(tidyr)
library(purrr)

# objective: test with WHO's 2/3 and 3/3 criteria:

#   3/3: matching alleles present in all of the loci
#   2/3: matching alleles present in at least 2/3 of the loci


### inputs

selected_columns <- c("PairsID", "labels")
labels <- fread("TRAINING_DATA_CLEAN.csv", select = selected_columns)

allele_data <- readRDS("PAIRS_GENOMIC.RDS")
allele_data <- allele_data[allele_data$PairsID %in% unique(labels$PairsID),]

setDT(allele_data)

selected_columns2 <- c("PairsID", "time_point", "locus", "allele")
allele_data <- allele_data %>% select(selected_columns2)

gc()

# Setup for parallel processing and output
num_cores <- detectCores() - 0
chunk_size <- 20000  # Adjust based on available memory
output_file <- "allele_data_pivoted.csv"

# Initialize the progress bar
pb <- progress_bar$new(
  format = "Processing [:bar] :percent :elapsedfull",
  total = ceiling(length(unique(allele_data$PairsID)) / chunk_size),
  width = 100
)

# Function to process each pair in a chunk
process_pair <- function(pair_id) {
  # Filter data for the given pair_id
  allele_subset <- allele_data[PairsID == pair_id]
  
  # Aggregate alleles by unique values
  allele_subset <- allele_subset[, .(alleles = list(unique(allele))), by = .(PairsID, time_point, locus)]
  
  # Pivot the data
  allele_subset <- dcast(allele_subset, PairsID + locus ~ time_point, value.var = "alleles", fill = NA)
  
  # Rename columns for clarity
  setnames(allele_subset, old = c("D0", "Dx"), new = c("D0_alleles", "Dx_alleles"))
  
  return(allele_subset)
}

# Write results to disk in chunks
write_chunk_to_csv <- function(data, file_path, append = FALSE) {
  fwrite(data, file_path, append = append)
}

# Process data in chunks by PairsID
unique_pairs <- unique(allele_data$PairsID)
for (i in seq(1, length(unique_pairs), by = chunk_size)) {
  # Define the current chunk of pair IDs
  chunk_pair_ids <- unique_pairs[i:min(i + chunk_size - 1, length(unique_pairs))]
  
  print(paste0("Processing PairsID ", first(chunk_pair_ids), " to ", last(chunk_pair_ids), "..."))
  
  # Parallel processing of the chunk
  results <- mclapply(chunk_pair_ids, process_pair, mc.cores = num_cores)
  
  # Filter out any NULL results (if any)
  metrics_list_chunk <- Filter(Negate(is.null), results)
  
  # Combine results and write to disk incrementally
  if (length(metrics_list_chunk) > 0) {
    allele_data_df_chunk <- rbindlist(metrics_list_chunk, fill = TRUE)
    
    # Append the chunk results to the output file
    write_chunk_to_csv(allele_data_df_chunk, output_file, append = (i > 1))
  }
  
  pb$tick()
}

gc()


allele_data <- read.csv("allele_data_pivoted.csv")
setDT(allele_data)


# Add a new column indicating if D0 and Dx alleles match (1 if any matches, 0 if none), per PairsID per locus
allele_data <- allele_data %>%
  mutate(alleles_match = ifelse(map2_lgl(D0_alleles, Dx_alleles, ~any(.x %in% .y)), 1, 0))

allele_data <- allele_data %>%
  select(-D0_alleles, -Dx_alleles)

gc()


### WHO 2/3 and 3/3 criteria
n_loci <- length(unique(allele_data$locus))
two_thirds <- round(n_loci*2/3,0)

who_classification <- allele_data %>%
  group_by(PairsID) %>%
  summarise(classif_2_3 = ifelse(sum(alleles_match) >= two_thirds, "R", "NI"),
            classif_3_3 = ifelse(sum(alleles_match) == n_loci, "R", "NI")) %>%
  ungroup()

gc()


#TEST THIS NEW CLASSIF!!
# Loop over each 'x' and create new columns for the x/N classification
xN_classification <- allele_data %>%
  group_by(PairsID) %>%
  summarise(across(starts_with("alleles_match"), sum))  # Optional: Summarize alleles_match per group

for (x in 1:n_loci) {
  column_name <- paste0("xN_", x, "_", n_loci)
  
  xN_classification <- xN_classification %>%
    mutate(!!sym(column_name) := ifelse(alleles_match >= x, "R", "NI"))
}

xN_classification <- xN_classification %>% ungroup()


gc()


all_classifs <- merge(labels, who_classification, by = "PairsID")
all_classifs <- merge(all_classifs, xN_classification, by = "PairsID")
all_classifs <- as.data.frame(all_classifs)

all_classifs <- all_classifs %>%
  select(-alleles_match)


### output
write.csv(all_classifs, "WHO_xN_classifications_TRAINING_DATA.csv", row.names = F)

