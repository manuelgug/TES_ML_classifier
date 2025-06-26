
## I'M USING TRUTH STRAINS (NO FALSE NEGATIVES NOR POSITIVES) PROVIDED BY KATHRYN TO CREATE THE TRAINING DATA: MIXES, PAIRS, FEATURES
# THEN, CREATE THE MODEL
# THEN, I'M USING LAB CONTROLS (MIXES) TO CREATE THE TESTING DATA: PAIRS, FEATURES
# THEN, I'M USING THE MODEL TRAINED ON TRUTH DATA TO PREDICT THE LAB CONTROLS (WHICH SHOULD HAVE ERRORS)
# HOW ROBUST IS A MODEL TRAINED ON TRUTH AND TESTED ON ERROR-PRONE DATA?


##############################################################################
################### TRUTH ##########################
##############################################################################

library(dplyr)
library(ggdendro)
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(viridis)


TRUTH <- read.csv("pool1A_truth.csv")


## subset pfPHAST amplicons
madhito_amps <- read.csv("../../madhito_20amps.csv")
TRUTH <- TRUTH[TRUTH$locus %in% madhito_amps$locus,]


# Ignore masking, turn it into ref (.)
TRUTH$pseudo_cigar <- gsub("\\d+\\+[^N]*N", "", TRUTH$pseudo_cigar) # Remove masking
TRUTH$pseudo_cigar <- ifelse(TRUTH$pseudo_cigar == "" | is.na(TRUTH$pseudo_cigar), ".", TRUTH$pseudo_cigar) # If empty, add "." since it was reference


# Aggregate unmasked sequences
TRUTH <- unique(TRUTH)


# Create allele column
TRUTH$allele <- paste0(TRUTH$locus, "__", TRUTH$pseudo_cigar)


# Remove indels
TRUTH <- TRUTH[!grepl("I=", TRUTH$allele), ] # Remove alleles with I (insertion)
TRUTH <- TRUTH[!grepl("D=", TRUTH$allele), ] # Remove alleles with D (deletion)



# loci check
TRUTH %>% group_by(Strain) %>% summarise(length(unique(locus)))

# allele check
TRUTH %>% group_by(Strain) %>% summarise(length(allele))


##################################################################

# SHARED ALLELES COMPARISON

# Step 1: Create the presence/absence matrix
allele_matrix <- TRUTH %>%
  distinct(Strain, allele) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Strain, values_from = present, values_fill = 0)

# Step 2: Extract matrix (remove allele column)
allele_data <- as.data.frame(allele_matrix)
allele_only <- as.matrix(allele_data[ , -1])  # remove 'allele' column

# Step 3: Get strain names from column names
strain_names <- colnames(allele_only)

# Step 4: Compute pairwise shared alleles matrix
shared_matrix <- t(allele_only) %*% allele_only  # dot product

# Step 5: Assign row/column names (based on matrix dimensions)
rownames(shared_matrix) <- strain_names
colnames(shared_matrix) <- strain_names

# Sort strain names by total shared alleles (row sums)
strain_order <- shared_matrix %>%
  rowSums() %>%
  sort(decreasing = TRUE) %>%
  names()

# Convert to long format
shared_df <- as.data.frame(as.table(shared_matrix))
colnames(shared_df) <- c("Strain1", "Strain2", "SharedAlleles")

# Convert strain factors to sorted order
shared_df <- shared_df %>%
  mutate(Strain1 = factor(Strain1, levels = strain_order),
         Strain2 = factor(Strain2, levels = strain_order))

# Normalize shared alleles to similarity [0,1]
similarity_matrix <- shared_matrix / max(shared_matrix)

# Convert similarity to distance
distance_matrix <- 1 - similarity_matrix

# Convert to 'dist' object for clustering
dist_obj <- as.dist(distance_matrix)


# Define a viridis color ramp like ggplot2
viridis_colors <- colorRamp2(
  seq(0, max(shared_matrix), length.out = 100),
  viridis(100)
)

heatmap <- Heatmap(shared_matrix,
        name = "Shared Alleles",
        col = viridis_colors,  # equivalent to scale_fill_viridis_c()
        clustering_distance_rows = dist_obj,
        clustering_distance_columns = dist_obj,
        clustering_method_rows = "average",
        clustering_method_columns = "average",
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        row_names_side = "left",
        column_names_rot = 45)


png("TRUTH_heatmap.png", width = 8, height = 7, units = "in", res = 300, bg = "white")
draw(heatmap)
dev.off()


### REMOVE GENETICALLY IDENTICAL CLONES
remove_clones <- c("W2", "DD2k")
TRUTH_CLEAN <- TRUTH[!TRUTH$Strain %in% remove_clones,]


# checks:
unique(TRUTH_CLEAN$Strain)


##  EXPORT CLEAN DATASET
write.csv(TRUTH_CLEAN, "TRUTH_CLEAN.csv", row.names = F)



##############################################################################
################### lab_controls ##########################
##############################################################################

library(progress)
library(purrr)
library(dplyr)
library(fs)
library(stringr)
library(tidyr)
library(moire)
library(data.table)


# set min allele freq threshold
MAF <- 0.01


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


#remove 99% - 1% mixes since MAF filter is 0.02  anyways
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


### keep only what's in control_metadata.csv (no newer runs! avoid extra work here). this was done by hand (strain compositions)
control_metadata<- read.csv("control_metadata.csv")

CONTROLS_ALL <- CONTROLS_ALL %>%
  mutate(sampleID = if_else(!grepl("^N", sampleID), paste0("N", sampleID), sampleID)) #format sampleID

keep_samples <- intersect((CONTROLS_ALL$sampleID), unique(control_metadata$NIDA))

length(keep_samples)

CONTROLS_ALL <- CONTROLS_ALL[CONTROLS_ALL$sampleID %in% keep_samples,]


## export
write.csv(CONTROLS_ALL, "LAB_CONTROLS_CLEAN.csv", row.names = F)




