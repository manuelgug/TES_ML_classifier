
## I'M USING TRUTH STRAINS (NO FALSE NEGATIVES NOR POSITIVES) PROVIDED BY KATHRYN TO CREATE THE TRAINING DATA: MIXES, PAIRS, FEATURES
# THEN, CREATE THE MODEL
# THEN, I'M USING LAB CONTROLS (MIXES) TO CREATE THE TESTING DATA: PAIRS, FEATURES
# THEN, I'M USING THE MODEL TRAINED ON TRUTH DATA TO PREDICT THE LAB CONTROLS (WHICH SHOULD HAVE ERRORS)

library(progress)
library(purrr)
library(dplyr)
library(fs)
library(stringr)
library(tidyr)
library(moire)
library(data.table)


# set min allele freq threshold
MAF <- 0.02


#####################################################################################3
#### 1. CREATE CONTROL DATA ------

# ISG / CISM

directory_path <- "../FILTERED_DATA_FROM_CLUSTER/"
data_all <- data.frame()

pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = length(dir_ls(path = directory_path, regexp = "_FILTERED", ignore.case = TRUE))
)

for (folder_path in dir_ls(path = directory_path, regexp = "_FILTERED", ignore.case = TRUE)) {
  pb$tick() 
  
  folder_name <- path_file(folder_path)
  file_path <- file.path(folder_path, "allele_data_global_max_0_filtered.csv")
  
  cat("\n")
  print(folder_name)
  
  if (file.exists(file_path)) {
    data <- read.csv(file_path)
    data$run <- sub("^.*/([^/]+)/[^/]+$", "\\1", file_path)
    data_all <- rbind(data_all, data)
  }
}

data_all$run <- sub("_RESULTS_v0.1.8_FILTERED$", "", data_all$run)

#explore control names
#print(sort(unique(data_all$sampleID), decreasing = F), max = 10000)

# subset controls
control_patterns <- c("dd2", "3d7", "hb3", "d10", "ds2", "ds4", "qs1", "qs2", "ts2")

data_controls <- data_all %>%
  filter(str_detect(tolower(sampleID), paste(control_patterns, collapse = "|")))

data_controls <- data_controls%>%
  select(-Category)

colnames(data_controls)

#how many isg/cism controls?
unique(paste0(data_controls$sampleID, "__", data_controls$run))


##########################################################################################################################################


# UCSF

controls_ucsf <- read.csv("combined_mixture_control_allele_data_2023_11_29.csv")

controls_ucsf$run <- "UCSF" 

controls_ucsf <- controls_ucsf %>%
  select(-X, -Barcode, -qpcr, -qpcr_flag)

colnames(controls_ucsf)[1] <- "sampleID"

#rearrange columns
controls_ucsf <- controls_ucsf %>%
  select(sampleID, locus, asv, reads, allele, pseudo_cigar, norm.reads.locus, n.alleles, run)

colnames(controls_ucsf)

unique(paste0(controls_ucsf$sampleID, "__", controls_ucsf$run))


##########################################################################################################################################


# UCSF KATRHYN

library(readxl)

# Get the sheet names
sheets <- excel_sheets("FORMAT_final_pool_1A_table_separate_mixes.xlsx")

# Read each sheet into a separate data frame inside a list
list_of_dfs <- lapply(sheets, function(sheet) {
  read_xlsx("FORMAT_final_pool_1A_table_separate_mixes.xlsx", sheet = sheet)
})

names(list_of_dfs) <- sheets


# Function to clean and transform each data frame
clean_and_transform_df <- function(df) {
  
  # Step 1: Remove rows where the "class" column contains "fn"
  if ("class" %in% colnames(df)) {
    df <- df[df$class != "fn", ]
  }
  
  # 1.5 remove pseudopool since it was not used before
  if ("pool" %in% colnames(df)) {
    df <- df[df$pool == "pseudo" & df$omega == 120, ]
  }
  
  # Step 2: Remove rows where any column with "Percentage" in the name contains 99
  percentage_cols <- grep("Percentage", colnames(df), value = TRUE)  # Find "Percentage" columns
  if (length(percentage_cols) > 0) {
    for (col in percentage_cols) {
      df <- df[df[[col]] != 99, ]
    }
  }
  
  # Step 3: Paste contents of sampleID, pool, omega, Strain and Percentage columns into sampleID
  strain_cols <- grep("Strain", colnames(df), value = TRUE)  # Find "Strain" columns
  combined_cols <- c("sampleID", "pool", "omega", strain_cols, percentage_cols)
  combined_cols <- combined_cols[combined_cols %in% colnames(df)]  # Ensure columns exist
  
  if (length(combined_cols) > 1) {  # If more than just sampleID exists
    df$sampleID <- apply(df[, combined_cols], 1, function(row) paste(na.omit(row), collapse = "_"))
  }
  
  # Step 4: Remove pool, omega, Strain, Percentage columns, and reads.locus
  cols_to_remove <- c("pool", "omega", strain_cols, percentage_cols, "reads.locus")
  df <- df[ , !(colnames(df) %in% cols_to_remove)]
  
  # Step 5: Create a column named "asv" and fill it with NA
  df$asv <- NA
  
  #Step 6: create allele column
  df$allele <-paste0(df$locus, "__", df$pseudo_cigar) 
  
  #Step 7: create n.alleles column
  df <- df %>%
    group_by(sampleID, locus) %>%
    mutate(n.alleles = n_distinct(allele)) %>%
    ungroup()
  
  # Step 8: Order the columns according to the required order
  ordered_cols <- c("sampleID", "locus", "asv", "reads", "allele", "pseudo_cigar", "norm.reads.locus", "n.alleles", "run")
  df <- df[ , ordered_cols[ordered_cols %in% colnames(df)]]
  
  
  return(df)
}

# Apply the function to each data frame in the list
cleaned_transformed_list_of_dfs <- lapply(list_of_dfs, clean_and_transform_df)

controls_ucsf_kathryn <- bind_rows(cleaned_transformed_list_of_dfs)


unique(paste0(controls_ucsf_kathryn$sampleID, "__", controls_ucsf_kathryn$run))


##########################################################################################################################################


# MERGE CONTROLS

CONTROLS_ALL <- rbind(data_controls, controls_ucsf)
CONTROLS_ALL <- rbind(CONTROLS_ALL, controls_ucsf_kathryn)

# remove masking and aggregate
CONTROLS_ALL$pseudo_cigar <- gsub("\\d+\\+[^N]*N", "", CONTROLS_ALL$pseudo_cigar) # Remove masking
CONTROLS_ALL$pseudo_cigar <- ifelse(CONTROLS_ALL$pseudo_cigar == "" | is.na(CONTROLS_ALL$pseudo_cigar), ".", CONTROLS_ALL$pseudo_cigar) # If empty, add "." since it was reference
CONTROLS_ALL <- CONTROLS_ALL %>% group_by(sampleID, locus, allele, pseudo_cigar, run) %>% summarise(reads = (sum(reads)),
                                                                      norm.reads.locus = sum(norm.reads.locus)) %>% ungroup()

#create allele column (locus +  pseudocigar)
CONTROLS_ALL$allele <-paste0(CONTROLS_ALL$locus, "__", CONTROLS_ALL$pseudo_cigar) 

#add run name to sampleID to differentiate across runs
CONTROLS_ALL$sampleID <- paste(CONTROLS_ALL$sampleID, CONTROLS_ALL$run, sep ="___")

#check for identical controls and keep one of each to avoid repeated pairs later on
# Group by sampleID and create a concatenated string of sorted alleles for each sampleID
allele_groups <- CONTROLS_ALL %>%
  group_by(sampleID) %>%
  summarise(allele_content = paste(sort(unique(allele)), collapse = ",")) %>%
  ungroup()

# Group by the concatenated allele content and find all sampleIDs that have the same allele content
duplicate_samples <- allele_groups %>%
  group_by(allele_content) %>%
  summarise(samples = list(sampleID)) %>%
  filter(lengths(samples) > 1)

# Extract all samples that are not the first in each group
non_first_samples <- duplicate_samples %>%
  mutate(non_first = map(samples, ~ .x[-1])) %>%  # Remove the first element from each group
  pull(non_first) %>%
  unlist()

# Filter CONTROLS_ALL to remove these non-first samples
CONTROLS_ALL <- CONTROLS_ALL %>%
  filter(!sampleID %in% non_first_samples)


#remove 99% - 1% mixes since MAF filter is 0.01  anyways
remove_mixes_pattern <- c("_99_S", "_1_S", "_1_99_", "_99_1_", "DC11-", "DC1-", "DC11_", "DC1_")

# Filter out rows where sampleID matches any of the patterns
CONTROLS_ALL <- CONTROLS_ALL %>%
  filter(!str_detect(tolower(sampleID), paste(tolower(remove_mixes_pattern), collapse = "|")))

length(unique(paste(CONTROLS_ALL$sampleID, CONTROLS_ALL$run)))

unique(paste(CONTROLS_ALL$sampleID, CONTROLS_ALL$run))

#keep only diversity amplicons "-1A"
CONTROLS_ALL <- CONTROLS_ALL[grepl("-1A$", CONTROLS_ALL$locus),]

#remove spaces if any
CONTROLS_ALL$sampleID <- gsub(" ", "_", CONTROLS_ALL$sampleID)

#### allele name formatting ####
CONTROLS_ALL$allele <- paste0(CONTROLS_ALL$locus, "___", CONTROLS_ALL$pseudo_cigar)


CONTROLS_ALL <- CONTROLS_ALL %>% select(sampleID, locus, allele, run, reads, norm.reads.locus)



#####################################################################################3
#### 4. SELECT SAMPLES BASED ON REGULAR CRITERIA (min loci, min read counts) ------

# What's the  minimum amount of loci sequenced for a sample to be considered good quality ?
n_loci_per_sample <- CONTROLS_ALL %>%
  group_by(sampleID) %>%
  summarise(n_loci = length(unique(locus)))

#hist(n_loci_per_sample$n_loci)

# Filter and identify unique NIDA with at least 50 loci having >= 100 reads (NANNA'S AND SIMONE'S FILTER)
unique_nida <- CONTROLS_ALL %>%
  # Filter rows where reads >= 100
  filter(reads >= 100) %>%
  # Group by NIDA
  group_by(sampleID) %>%
  # Count distinct loci for each NIDA
  summarise(n_loci = n_distinct(locus)) %>%
  # Filter NIDA where n_loci >= 50
  filter(n_loci >= 50) %>%
  # Select the unique NIDA values
  pull(sampleID)

# keep HIGH QUALITY SAMPLES
CONTROLS_ALL <- CONTROLS_ALL %>%
  filter(sampleID %in% unique_nida)

gc()


#####################################################################################
####### 7. ALLELE FILTERING #######--------------------

# Convert merged_dfs to a data.table for faster operations
CONTROLS_ALL <- as.data.table(CONTROLS_ALL)
setDT(CONTROLS_ALL)

# Apply MAF filtering and removes rows with alleles containing 'I=' or 'D='
CONTROLS_ALL <- CONTROLS_ALL[
  norm.reads.locus >= MAF &          # Filter based on MAF threshold
    reads >=10 &
    !grepl("I=|D=", allele)            # Remove alleles with 'I=' (insertion) or 'D=' (deletion) in one step
]





#### KEEP PFPHAST LOCI ONLY ###

madhito_amps <- read.csv("../../madhito_20amps.csv")

CONTROLS_ALL <- CONTROLS_ALL[CONTROLS_ALL$locus %in% madhito_amps$locus,]





#####################################################################################3
### 2. CREATE CONTROL METADATA (PAIRS) ----

sampleIDs <- unique(CONTROLS_ALL$sampleID)

# Create all possible pairs of sampleIDs
pairs_df <- expand.grid(NIDA1 = sampleIDs, NIDA2 = sampleIDs, stringsAsFactors = FALSE)

# #just remove pairs of the same sample
# pairs_df <- pairs_df %>%
#   filter(NIDA1 != NIDA2)


# Create a unique PairsID for each pair
pairs_df <- pairs_df %>%
  mutate(PairsID = row_number()) %>%
  pivot_longer(cols = c(NIDA1, NIDA2), names_to = "time_point", values_to = "NIDA")

# Assign time_point based on the column names
pairs_df <- pairs_df %>%
  mutate(time_point = ifelse(time_point == "NIDA1", "D0", "Dx"))

# Rearrange columns to match the desired format
control_pairs <- pairs_df %>%
  select(PairsID, NIDA, time_point)

# add metadata to the same df

# Step 1: Extract 'run' from NIDA column
control_pairs <- control_pairs %>%
  mutate(run = str_extract(NIDA, "(?<=___).*"))



#####################################################################################3
####### 3. FORMATTING #######------------------

#paired_data import and formatting
control_pairs$NIDA <- gsub("\\.", "_", control_pairs$NIDA) #change "." for "_" on NIDA column
control_pairs$NIDA <- ifelse(startsWith(control_pairs$NIDA, "N"), 
                             control_pairs$NIDA, 
                             paste0("N", control_pairs$NIDA)) #add initial N to all NIDAs


# Remove "_S*" from sampleID only if "_S" is not present in control_pairs$NIDA
contains_S <- any(grepl("_S", control_pairs$NIDA))
if (!contains_S) {
  CONTROLS_ALL$sampleID <- sub("_S.*$", "", CONTROLS_ALL$sampleID)
}
colnames(CONTROLS_ALL)[1] <- "NIDA"
CONTROLS_ALL$NIDA <- gsub("\\.", "_", CONTROLS_ALL$NIDA) #change "." for "_" on NIDA column
CONTROLS_ALL$NIDA <- ifelse(startsWith(CONTROLS_ALL$NIDA, "N"), 
                            CONTROLS_ALL$NIDA, 
                            paste0("N", CONTROLS_ALL$NIDA)) #add initial N to all NIDAs

#are all the NIDAs in paired_data present in CONTROLS_ALL?
control_pairs_NIDA <- unique(control_pairs$NIDA)
CONTROLS_ALL_NIDA <- unique(CONTROLS_ALL$NIDA)
missing_NIDAs <- control_pairs_NIDA[!control_pairs_NIDA %in% CONTROLS_ALL_NIDA]

if (length(missing_NIDAs) > 0) {
  print("These NIDAs from control_pairs are not present in CONTROLS_ALL:")
  print(missing_NIDAs)
  
  #remove paired samples with missing NIDAs
  missing_NIDAs_PairsID<-subset(control_pairs, NIDA %in% missing_NIDAs)[2]
  control_pairs <- control_pairs[!(control_pairs$PairsID %in% missing_NIDAs_PairsID$PairsID), ]
  
  
} else {
  print("All NIDAs from control_pairs are present in CONTROLS_ALL.")
}

# ### remove asv, pseudo_cigar and run columns because ugh
# if ("run" %in% colnames(CONTROLS_ALL)) {
#   # If present, select all columns except 'asv', 'pseudo_cigar', and 'run'
#   CONTROLS_ALL <- CONTROLS_ALL %>%
#     select(-pseudo_cigar, -run, -n.alleles)
# } else {
#   # If not present, keep only the 'run' column (if it exists)
#   CONTROLS_ALL <- CONTROLS_ALL %>%
#     select(-pseudo_cigar, -n.alleles)
# }
# 
# if ("run" %in% colnames(control_pairs)) {
#   # If present, remove 'run' column
#   control_pairs <- control_pairs %>%
#     select(-run)
# }

CONTROLS_ALL <- CONTROLS_ALL %>% select(-run)

#merge control_pairs and CONTROLS_ALL by NIDA
merged_dfs <-  inner_join(control_pairs, CONTROLS_ALL, by = "NIDA")

### remove asv column because ugh
merged_dfs <- merged_dfs %>%
  arrange(PairsID, time_point)


#check
all(unique(control_pairs$NIDA) %in% unique(merged_dfs$NIDA))
all(unique(control_pairs$PairsID) %in% unique(merged_dfs$PairsID))

all(unique(merged_dfs$NIDA) %in% unique(control_pairs$NIDA))
all(unique(merged_dfs$PairsID) %in% unique(control_pairs$PairsID))


# #### RAM DUMP
# rm(controls_ucsf, controls_ucsf_kathryn, data_all, data)
# gc()

merged_dfs <- merged_dfs[merged_dfs$PairsID %in% control_pairs$PairsID,]

#check
all(unique(control_pairs$NIDA) %in% unique(merged_dfs$NIDA))
all(unique(control_pairs$PairsID) %in% unique(merged_dfs$PairsID))

all(unique(merged_dfs$NIDA) %in% unique(control_pairs$NIDA))
all(unique(merged_dfs$PairsID) %in% unique(control_pairs$PairsID))





## format metadata. add labels. add nstrains for each time

# 1) import control metadata
meta <- read.csv("control_metadata.csv")

# 2) count strains for each mix
library(stringr)

meta <- meta %>%
  rowwise() %>%
  mutate(
    n_strains = sum(c_across(contains("strain")) != "")
  ) %>%
  ungroup()

# 3) merge with control_pairs to get the n_strains
control_pairs <- left_join(control_pairs, meta[c("n_strains", "NIDA")], by = "NIDA")

# 4) reformat
library(tidyr)

# Pivot the data wider by time_point (D0 vs Dx)
formatted_pairs <- control_pairs %>%
  mutate(row = row_number()) %>%
  pivot_wider(
    id_cols = PairsID,
    names_from = time_point,
    values_from = c(NIDA, n_strains),
    names_sep = "_"
  ) %>%
  # Rename columns as requested
  rename(
    D0_nstrains = n_strains_D0,
    Dx_nstrains = n_strains_Dx,
    D0_sample = NIDA_D0,
    Dx_sample = NIDA_Dx) %>%
  select(PairsID,  D0_sample, Dx_sample, D0_nstrains, Dx_nstrains)

formatted_pairs$NIDA1 <- formatted_pairs$D0_sample
formatted_pairs$NIDA2 <- formatted_pairs$Dx_sample


## 4) assign labels
library(purrr)

# Step 1: Get all unique NIDA values
nida_list <- unique(meta$NIDA)

# Step 2: Create all pairwise combinations (including both orders if needed)
nida_pairs <- expand.grid(NIDA1 = nida_list, NIDA2 = nida_list) %>%
  filter(NIDA1 != NIDA2)

# # Optional: Remove duplicate reversed pairs
# nida_pairs <- nida_pairs %>%
#   rowwise() %>%
#   mutate(pair_id = paste(sort(c(NIDA1, NIDA2)), collapse = "_")) %>%
#   distinct(pair_id, .keep_all = TRUE) %>%
#   select(-pair_id) %>%
#   ungroup()

# Step 3: Create a helper function to get strains from a given NIDA
get_strains <- function(nida_name) {
  meta %>%
    filter(NIDA == nida_name) %>%
    select(starts_with("strain")) %>%
    unlist(use.names = FALSE) %>%
    na.omit() %>%
    trimws() %>%
    .[. != ""]
}

# Step 4: Compare strain sets for each pair
compare_strains <- nida_pairs %>%
  rowwise() %>%
  mutate(
    shared = length(intersect(get_strains(NIDA1), get_strains(NIDA2))) > 0,
    labels = ifelse(shared, "R", "NI")
  ) %>%
  ungroup() %>%
  select(NIDA1, NIDA2, labels)

# View result
print(compare_strains)

### 5) merge all to create final metadata
final_metadata <- formatted_pairs %>%
  left_join(compare_strains, by = c("NIDA1", "NIDA2"))

final_metadata$labels <- ifelse(final_metadata$NIDA1 == final_metadata$NIDA2, "R", final_metadata$labels)

# 6) final touches
# reoving some pairs that don't appear in the metadata. this is becauase metadata was originally done manually and i added more runs later. all is good. no need to add new runs anyways
final_metadata <- final_metadata[!is.na(final_metadata$labels),]

# subset genomic
merged_dfs <- merged_dfs[merged_dfs$PairsID %in% final_metadata$PairsID,]

# 1. Check PairsID matches both ways
all(merged_dfs$PairsID %in% final_metadata$PairsID)     # A in B
all(final_metadata$PairsID %in% merged_dfs$PairsID)     # B in A

# 2. Check NIDA matches both ways (against NIDA1 and NIDA2)
all(merged_dfs$NIDA %in% c(final_metadata$NIDA1, final_metadata$NIDA2))  # A in B
all(c(final_metadata$NIDA1, final_metadata$NIDA2) %in% merged_dfs$NIDA)  # B in A


################################## export

saveRDS(merged_dfs,  "Lab_Controls_genomic.RDS")
saveRDS(final_metadata,  "Lab_Controls_metadata.RDS")


#################################################################################### AQUÍ VOY!!! 11/JUN/2025


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


  # for the training data:
  PAIRS_METADATA <- readRDS(paste0("Lab_Controls_metadata.RDS"))
  PAIRS_GENOMIC <- readRDS(paste0("Lab_Controls_genomic.RDS"))
  
  PAIRS_GENOMIC <- as.data.table(PAIRS_GENOMIC)
  PAIRS_METADATA <- as.data.table(PAIRS_METADATA)
  
  PAIRS_GENOMIC <- left_join(PAIRS_GENOMIC, PAIRS_METADATA, by = "PairsID") 



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
unique_pairs <- unique(PAIRS_METADATA$PairsID)

# Create a named vector for mapping PairsID to D0 and Dx for faster access
pairs_to_d0 <- PAIRS_METADATA$NIDA1[match(unique_pairs, PAIRS_METADATA$PairsID)]
pairs_to_dx <- PAIRS_METADATA$NIDA2[match(unique_pairs, PAIRS_METADATA$PairsID)]

# Preprocess merged_dfs_filtered into a list for fast access
merged_dfs_filtered <- split(PAIRS_GENOMIC, PAIRS_GENOMIC$NIDA)
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
output_file <- paste0("delta_features_Lab_Controls_TRAINING_DATA.csv")

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

delta_metrics_df_final <- read.csv(paste0("delta_features_Lab_Controls_TRAINING_DATA.csv"))

delta_metrics_df_final <- delta_metrics_df_final %>%
  select(PairsID, everything())

write.csv(delta_metrics_df_final, paste0("delta_features_Lab_Controls_TRAINING_DATA.csv"), row.names = F)



##### MERGE FEATURES AND OUTPUT ------

# FEATURES <- left_join(dres0_long_final_summarized, delta_metrics_df_final, by = "PairsID")

FEATURES <- left_join(delta_metrics_df_final, PAIRS_METADATA, by = "PairsID")

# ECOI CHANGE
FEATURES$coi_change <- FEATURES$Dx_nstrains - FEATURES$D0_nstrains

# some NAs for whatever reason...
FEATURES <- FEATURES[!is.na(FEATURES$coi_change),]


write.csv(FEATURES, paste0("Lab_Controls_TRAINING_DATA.csv"), row.names = F)



####################### EDA @@@@@@@@@@@@@@@@@@

library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(uwot)
library(ggpubr)
library(rstatix)
library(cluster)

#----------------------------------------------------------
# Data Import & Preparation
#----------------------------------------------------------
library(corrplot)

TRAINING_DATA <- read.csv("Lab_Controls_TRAINING_DATA.csv", row.names = 1)
TRAINING_DATA$eCOI_pairs <- paste0(TRAINING_DATA$D0_nstrains, "__", TRAINING_DATA$Dx_nstrains)
LABELS <- data.frame(labels = TRAINING_DATA$labels)
LABELS$labels <- as.factor(LABELS$labels)

#REAL_DATA <- read.csv(paste0(site, "_REAL_DATA.csv"), stringsAsFactors = FALSE, colClasses = c(NIDA1 = "character", NIDA2= "character")) 

features_to_use <- colnames(TRAINING_DATA)[!colnames(TRAINING_DATA) %in% c("PairsID", "NIDA1", "NIDA2", "pair_type", "IBD_estimate",
                                                                           "offset_naive_coi_D0", "offset_naive_coi_Dx", 
                                                                           "replacement_pattern_score", "locus_discordance_rate", 
                                                                           "D0_sample", "Dx_sample","D0_nstrains", "Dx_nstrains", "eCOI_pairs", "labels", "shared_count","union_count", "shared_prop")]

#----------------------------------------------------------
# 1) FEATURE CORRELATIONS
#----------------------------------------------------------
corr_data <- TRAINING_DATA %>% select(all_of(features_to_use))
corrplot(cor(corr_data, use = "complete.obs"), method = "pie", tl.cex = 0.7)

#----------------------------------------------------------
# 2) FEATURE DISTRIBUTIONS: GLOBAL & STRATIFIED BY eCOI_pairs
#----------------------------------------------------------
# Convert to long format
training_long <- TRAINING_DATA %>%
  select(all_of(features_to_use), eCOI_pairs, labels) %>%
  pivot_longer(cols = all_of(features_to_use), names_to = "feature", values_to = "value")

# Global distributions (density plots)
ggplot(training_long, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  facet_wrap(~ feature, scales = "free") +
  theme_minimal() +
  labs(title = "", x = "Value", y = "Density")

# Stratified by eCOI_pairs
distros_strat <- ggplot(training_long, aes(x = value, fill = labels, color = labels)) +
  geom_density(alpha = 0.5) +
  facet_grid(feature~ eCOI_pairs, scales = "free") +
  theme_minimal() +
  labs(title = "", x = "Value", y = "Density")+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8) )

ggsave(paste0("distributions_stratified_Lab_Controls.png"), distros_strat, height = 12, width = 14, dpi = 300, bg= "white")

#----------------------------------------------------------
# 3) BOXPLOTS OF FEATURES: GLOBAL & STRATIFIED BY eCOI_pairs WITH WILCOX TEST P-VALUES
#----------------------------------------------------------
# Boxplots by label (global)
p1 <- ggplot(training_long, aes(x = labels, y = value, fill = labels)) +
  geom_boxplot(alpha = 1, color = "black", outlier.alpha = 0.05) +  # hide outliers; outlier.shape = NA
  facet_grid(feature ~ eCOI_pairs, scales = "free") +
  # stat_compare_means(method = "wilcox.test", label = "p.signif", label.y.npc = "top") +
  #scale_fill_brewer(palette = "Set3") +
  theme_minimal(base_size = 14) +
  labs(title = "", 
       x = "Label", 
       y = "Value") +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )

ggsave(paste0("boxplots_stratified_Lab_Controls.png"), p1, height = 12, width = 17, dpi = 300, bg= "white")


# p2 <- ggplot(training_long, aes(x = eCOI_pairs, y = value, fill = eCOI_pairs)) +
#   geom_boxplot(alpha = 0.8, color = "black", outlier.alpha = 0.05) +  # hide outliers; outlier.shape = NA
#   facet_wrap(~ feature, scales = "free") +
#   #stat_compare_means(method = "wilcox.test", label = "p.signif", label.y.npc = "top") +
#   scale_fill_brewer(palette = "Set3") +
#   theme_minimal(base_size = 14) +
#   labs(title = "", 
#        x = "Pair Type", 
#        y = "Value") +
#   theme(legend.position = "none",
#         strip.background = element_rect(fill = "lightgrey", color = "black"),
#         strip.text = element_text(face = "bold"),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.border = element_rect(color = "black", fill = NA, size = 0.8) )
# 
# ggsave(paste0("boxplots_global_", site, ".png"), p2, height = 8, width = 13, dpi = 300, bg= "white")



#----------------------------------------------------------
# 4) UMAP OF FEATURES: COLORED BY LABEL, SHAPED BY eCOI_pairs
#----------------------------------------------------------
# Run UMAP using uwot
umap_res <- uwot::umap(TRAINING_DATA %>% select(all_of(features_to_use)),
                       n_neighbors = 15,
                       seed = 420, 
                       n_threads = 20, 
                       n_components = 3,
                       verbose = T)

# Convert the result into a data frame and rename columns
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2", "UMAP3")

# Add label and eCOI_pairs information
umap_df$label <- LABELS$labels
umap_df$eCOI_pairs <- TRAINING_DATA$eCOI_pairs

# Plot UMAP embedding colored by label and shaped by eCOI_pairs
umap_all <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, size = UMAP3, color = label)) +
  geom_point(size = 3, alpha = 0.3) +
  theme_minimal() +
  labs(title = "", x = "UMAP1", y = "UMAP2", 
       color = "Label", shape = "eCOI_pairs")

#umap_all

ggsave(paste0("umap_global_Lab_Controls.png"), umap_all, height = 7, width = 8, dpi = 300, bg= "white")




### PREDICT WITH THE MODEL CREATED FROM TRUTH STRAINS!

model <- readRDS("LogReg_model_TRUTH_CLEAN.RDS")

thresholds <- read.csv("training_Results_LogReg_TRUTH_CLEAN.csv")
thresholds <- thresholds %>% rename(pair_type = eCOI_pairs)


# metatada
PAIRS_METADATA <- readRDS("Lab_Controls_metadata.RDS")

# Input features testing data for subsetting
test_data <- read.csv("Lab_Controls_TRAINING_DATA.csv")
test_data$pair_type <- paste0(test_data$D0_nstrains, "__", test_data$Dx_nstrains)
TEST_META <- test_data %>%
  mutate(PairsID = as.character(PairsID),
         labels = factor(labels, levels = c("NI", "R"))) %>%
  select(PairsID, labels, pair_type) %>%  # Renamed to pair_type for consistency
  as.data.table()

# subset metadata to only account for test data
PAIRS_METADATA <- PAIRS_METADATA[PAIRS_METADATA$PairsID %in% TEST_META$PairsID,]


TEST_META_threshs <- left_join(TEST_META, thresholds[c("pair_type", "decision_threshold")], by = "pair_type")
TEST_META_threshs$preds_prob <- predict(model, newdata = test_data, type = "prob")[, "R"]
TEST_META_threshs$preds <- ifelse(TEST_META_threshs$preds_prob >= TEST_META_threshs$decision_threshold, "R", "NI")


# Function to calculate performance metrics
calculate_metrics <- function(TEST_META_threshs) {
  # Initialize an empty results dataframe
  results <- data.frame(pair_type = character(),
                        sensitivity = numeric(),
                        specificity = numeric(),
                        R_pairs = numeric(),
                        NI_pairs = numeric(),
                        stringsAsFactors = FALSE)
  
  # Loop through each unique pair_type combination
  for (strain_comb in unique(TEST_META_threshs$pair_type)) {
    
    # Subset the TEST and TEST_labels based on the current combination
    subset_indices <- TEST_META_threshs$pair_type == strain_comb
    subset_TEST_labels <- TEST_META_threshs$labels[subset_indices]
    
    # Count R and NI pairs
    r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
    ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
    
    # Get the corresponding predictions for the current subset
    subset_preds <- TEST_META_threshs$preds[subset_indices]
    
    # Evaluate the confusion matrix for the current subset
    cm <- caret::confusionMatrix(as.factor(subset_preds), as.factor(subset_TEST_labels), positive = "R")
    
    sens <- cm$byClass["Sensitivity"]
    spec <- cm$byClass["Specificity"]
    
    # Append results
    results <- rbind(results, data.frame(pair_type = strain_comb,
                                         sensitivity = sens,
                                         specificity = spec,
                                         R_pairs = r,
                                         NI_pairs = ni,
                                         stringsAsFactors = FALSE))
  }
  
  # Join with thresholds
  results <- left_join(results, thresholds[c("pair_type", "decision_threshold")], by = "pair_type")
  return(results)
}

# Calculate metrics
results <- calculate_metrics(TEST_META_threshs)


# Calculate metrics
results <- calculate_metrics(TEST_META_threshs)


# Reshape data to long format for easy plotting
best_decision_thresholds_long <- results %>%
  select(pair_type, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity), 
               names_to = "Metric", 
               values_to = "Value")

# Create bar plot
metrics <- ggplot(best_decision_thresholds_long, aes(x = pair_type, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +  # Dodge separates bars for clarity
  labs(title = "",
       x = "Pair Type",
       y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept = 0.9, linetype = "solid", color = "black")+
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "black")

metrics


# TO DO:
## TRIPLE CHECK BOTH THE LAB CONTROLS AND THE TRUTH DATA. MAYBE THERE ARE SOME SAMPLES/MIXES IN THE LAB CONTROLS THAT MAY BE BETTER OFF NOT INCLUDING DUE TO UNCERTAIN COMPOSITION
## IS THE AMOUNT OF PAIRS TOO LOW FOR THE TRUTH DATA?
## MAYBE DON'T SPLIT THE TRUTH DATA THIS TIME?? IS NOT NEEDED GIVEN THAT THE TSTING WILL BE DONE WOTH THE LAB CONTROLS...
# REMOVE MIXES WITH DD2k? AND W2? APPARENTLY BOTH ARE THE SAME AS DD2. THEY ARE TAGGED AS DD2, BUT MAYBE THERE COULD BE SOME TAGGING ISSUES IN SOME LABS...
# MAYBE RESTRICT TO >0.98 MIXES??? AT THE MOMENT IS THOUGHT AS IF IT WERE >0.99

