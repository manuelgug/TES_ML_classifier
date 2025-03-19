
library(dplyr)
library(stringr)



####### 0) PARAMETERS------------------

main_dir <- "." # directory where all the filtered runs data are located
metadata_file <- "Final_list_TES_Envio_ISG_INPUT_TABLE.csv" # should have the nidas, the timepoints and the pairs to which samples belong to
maf_filter <- 0.02
min_allele_read_count <- 10 
top_n_amps <- 50 # top most variable amps to use


####### 1) IMPORT DATA-------------------

#import metadata
metadata <- read.csv(metadata_file)
metadata$NIDA <- as.character(metadata$NIDA)


# import genomic data
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



####### 2) FORMAT DATA------------

# format nidas
merged_dfs <- merged_dfs %>%
  mutate(sampleID = gsub("N|_S.*$", "", sampleID),  # Remove "N" and everything after "_S"
         sampleID = if_else(str_detect(sampleID, "_"), sampleID, paste0(sampleID, ".0")),
         sampleID = gsub("_", ".", sampleID))  # change "_" for "."

metadata <- metadata %>%
  mutate(NIDA = gsub("N|_S.*$", "", NIDA),  # Remove "N" and everything after "_S"
         NIDA = if_else(!str_detect(NIDA, "\\."), paste0(NIDA, ".0"), NIDA))  # If no ".", add ".0"


# subset genomic data
merged_dfs <- merged_dfs[merged_dfs$sampleID %in% metadata$NIDA,]

# add run name to NIDA in case there is a duplicate sample from many runs
merged_dfs$sampleID <- paste0(merged_dfs$sampleID, "__", merged_dfs$run)



####### 3) CLEAN DATA --------

# identify sampleID with >=50 loci having >= 100 reads (NANNA'S AND SIMONE'S FILTER)
good_sampleID <- merged_dfs %>%
  
  group_by(sampleID, locus) %>%
  summarize(reads = sum(reads)) %>%
  
  # Filter rows where reads >= 100
  filter(reads >= 100) %>%   
  
  # Filter NIDA where n_loci >= 50
  summarise(n_loci = n_distinct(locus)) %>%
  filter(n_loci >= 50) %>%
  
  pull(sampleID)

# Keep good quality samples
merged_dfs <- merged_dfs[merged_dfs$sampleID %in% good_sampleID, ]

# keep samples with > 10000 reads: ANDRÉS
read_counts_above_10K_reads <- merged_dfs %>% group_by(sampleID) %>% summarise(total_reads = sum(reads)) 
read_counts_above_10K_reads <- read_counts_above_10K_reads[read_counts_above_10K_reads$total_reads > 10000,]
merged_dfs <- merged_dfs[merged_dfs$sampleID %in% read_counts_above_10K_reads$sampleID, ]


### keep pool 1A
merged_dfs <- merged_dfs[grepl("-1A$", merged_dfs$locus),]


# ignore masking, turn it into ref (.)
merged_dfs$pseudo_cigar <-  gsub("\\d+\\+[^N]*N", "", merged_dfs$pseudo_cigar) #remove masking
merged_dfs$pseudo_cigar <- ifelse(merged_dfs$pseudo_cigar == "" | is.na(merged_dfs$pseudo_cigar), ".", merged_dfs$pseudo_cigar) # if empty, add "." since it was reference

#aggregate unmasked sequences
merged_dfs_agg <- merged_dfs %>% group_by(sampleID, locus, pseudo_cigar) %>% summarise(reads = sum(reads), norm.reads.locus = sum(norm.reads.locus))


#create allele column
merged_dfs_agg$allele <- paste0(merged_dfs_agg$locus, "__", merged_dfs_agg$pseudo_cigar)

# remove indels
merged_dfs_agg <- merged_dfs_agg[!grepl("I=", merged_dfs_agg$allele),] #remove alleles with I (insertion)
merged_dfs_agg <- merged_dfs_agg[!grepl("D=", merged_dfs_agg$allele),] #remove alleles with D (deletion)


# apply MAF filter
merged_dfs_agg <- merged_dfs_agg[merged_dfs_agg$norm.reads.locus > maf_filter,] #from contaminants study


# remove alleles with low read counts
merged_dfs_agg <- merged_dfs_agg[merged_dfs_agg$reads > min_allele_read_count,]


# remove pseudocigar column
merged_dfs_agg <- merged_dfs_agg %>% select(-pseudo_cigar)



###### 4) KEEP SAMPLES WITH PAIRS

# keep samples that have a pair
ids <- data.frame(NIDA = unique(merged_dfs_agg$sampleID))
ids$NIDA <- gsub("__.*$", "", ids$NIDA)
ids_merged <- merge(ids, metadata, by =  "NIDA")

#select pairs to keep
pairs_counts <- ids_merged %>% group_by(PairsID) %>% summarise(n_samples = length(NIDA))
pairs_counts <- pairs_counts[pairs_counts$n_samples == 2,]

#update metadata
metadata_updated <- metadata[metadata$PairsID %in% pairs_counts$PairsID,]

#resubset data 
merged_dfs_agg <- merged_dfs_agg[gsub("__.*$", "", merged_dfs_agg$sampleID) %in% metadata_updated$NIDA,]



####### 4) SHARED AMPLICONS --------

#amount of samples on which each amp was sequenced
amps_data <- merged_dfs_agg %>% group_by(locus) %>% summarise(in_n_samples = length(unique(sampleID))) %>% arrange(-in_n_samples)

hist(amps_data$in_n_samples)

#selected based on the max amount of amplicons shared. COULD BE ADJUSTED!
selected_amps <- amps_data[amps_data$in_n_samples == max(amps_data$in_n_samples),]$locus

# subset data based on amplicons
merged_dfs_agg <- merged_dfs_agg[merged_dfs_agg$locus %in% selected_amps,]

# calculate heterozygosity of loci and remove those with He = 0
heterozygosity_per_sample <- merged_dfs_agg %>% #Compute Hₑ per sample per locus
  group_by(sampleID, locus) %>%
  mutate(freq = norm.reads.locus / sum(norm.reads.locus, na.rm = TRUE)) %>%
  summarise(
    heterozygosity = 1 - sum(freq^2, na.rm = TRUE), 
    .groups = "drop"
  )

variability_per_locus <- heterozygosity_per_sample %>% # Compute variability (SD, CV) across samples for each locus
  group_by(locus) %>%
  summarise(
    mean_He = mean(heterozygosity, na.rm = TRUE),
    sd_He = sd(heterozygosity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(-mean_He)  # Rank by SD (or use desc(cv_He) for relative variability)

hist(variability_per_locus$mean_He)

# keep amplicons with mean He above 10%
# top_var_amps <- variability_per_locus[variability_per_locus$mean_He > 0.1,]$locus 

# keep top 50 most variable amps
top_var_amps <- variability_per_locus[1:top_n_amps,]$locus

# subset data based on amplicons
merged_dfs_agg <- merged_dfs_agg[merged_dfs_agg$locus %in% top_var_amps,]


####### 5) CHECKS --------

length(unique(merged_dfs_agg$sampleID)) # nidas
length(unique(merged_dfs_agg$sampleID))/2 # pairs
length(unique(merged_dfs_agg$locus)) # amplicons shared by all nidas


###### 6) OUTPUTS ---------

write.csv(merged_dfs_agg, paste0("genomic_updated_top_", top_n_amps, "_amps.csv"), row.names = F)
write.csv(metadata_updated, paste0("metadata_updated_", top_n_amps, "_amps.csv"), row.names = F)

