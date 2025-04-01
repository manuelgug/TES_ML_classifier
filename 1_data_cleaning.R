
library(dplyr)
library(stringr)



####### 0) PARAMETERS------------------

site <- "Zambezia"

cum_curve_threshold <- 0.99

main_dir <- "." # directory where all the filtered runs data are located
metadata_file <- paste0("metadata_tes_", site, ".csv") # should have the nidas, the timepoints and the pairs to which samples belong to
maf_filter <- 0.02
min_allele_read_count <- 10 



####### 1) IMPORT DATA-------------------

#import metadata
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))


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

merged_dfs <- merged_dfs %>% select(sampleID, locus, pseudo_cigar, reads, norm.reads.locus, run)



####### 3) CLEAN TES DATA --------

clean_data <- function(merged_dfs, maf_filter, min_allele_read_count) {
  library(dplyr)
  
  # Identify sampleID with >=50 loci having >= 100 reads (NANNA'S AND SIMONE'S FILTER)
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
  
  # Keep samples with > 10000 reads: ANDRÉS
  read_counts_above_10K_reads <- merged_dfs %>%
    group_by(sampleID) %>%
    summarise(total_reads = sum(reads)) 
  
  read_counts_above_10K_reads <- read_counts_above_10K_reads[read_counts_above_10K_reads$total_reads > 10000,]
  merged_dfs <- merged_dfs[merged_dfs$sampleID %in% read_counts_above_10K_reads$sampleID, ]
  
  # Keep pool 1A
  merged_dfs <- merged_dfs[grepl("-1A$", merged_dfs$locus), ]
  
  # Ignore masking, turn it into ref (.)
  merged_dfs$pseudo_cigar <- gsub("\\d+\\+[^N]*N", "", merged_dfs$pseudo_cigar) # Remove masking
  merged_dfs$pseudo_cigar <- ifelse(merged_dfs$pseudo_cigar == "" | is.na(merged_dfs$pseudo_cigar), ".", merged_dfs$pseudo_cigar) # If empty, add "." since it was reference
  
  # Aggregate unmasked sequences
  merged_dfs_agg <- merged_dfs %>%
    group_by(sampleID, locus, pseudo_cigar) %>%
    summarise(reads = sum(reads), norm.reads.locus = sum(norm.reads.locus))
  
  # Create allele column
  merged_dfs_agg$allele <- paste0(merged_dfs_agg$locus, "__", merged_dfs_agg$pseudo_cigar)
  
  # Remove indels
  merged_dfs_agg <- merged_dfs_agg[!grepl("I=", merged_dfs_agg$allele), ] # Remove alleles with I (insertion)
  merged_dfs_agg <- merged_dfs_agg[!grepl("D=", merged_dfs_agg$allele), ] # Remove alleles with D (deletion)
  
  # Apply MAF filter
  merged_dfs_agg <- merged_dfs_agg[merged_dfs_agg$norm.reads.locus > maf_filter, ] # From contaminants study
  
  # Remove alleles with low read counts
  merged_dfs_agg <- merged_dfs_agg[merged_dfs_agg$reads > min_allele_read_count, ]
  
  # Remove pseudo_cigar column
  merged_dfs_agg <- merged_dfs_agg %>% select(-pseudo_cigar)
  
  return(merged_dfs_agg)
}

# clean data
merged_dfs_agg <- clean_data(merged_dfs, maf_filter, min_allele_read_count) 




###### 4) KEEP SAMPLES WITH PAIRS -----

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




###### 5) CLEAN SITE DATA -----

genomic_site <- read.csv(paste0("genomic_site_",site,".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA2 = "character"))
genomic_site <- genomic_site %>% rename(sampleID = NIDA2)

#clean data
genomic_site_agg <- clean_data(genomic_site, maf_filter, min_allele_read_count)

#remove tes data from site data 
genomic_site_agg <- genomic_site_agg[!sapply(genomic_site_agg$sampleID, function(x) any(grepl(x, merged_dfs_agg$sampleID, fixed = TRUE))), ]



##### 5) ADD DATA_TYPE COLUMN AND MERGE
merged_dfs_agg$data_type <- "tes"
genomic_site_agg$data_type <- "site"

data_all <- rbind(merged_dfs_agg, genomic_site_agg)




####### 6) SHARED AMPLICONS --------

#amount of samples on which each amp was sequenced
amps_data <- data_all %>% group_by(locus) %>% summarise(in_n_samples = length(unique(sampleID))) %>% arrange(-in_n_samples)

hist(amps_data$in_n_samples)

#selected based on the max amount of amplicons shared. COULD BE ADJUSTED!
selected_amps <- amps_data[amps_data$in_n_samples == max(amps_data$in_n_samples),]$locus

# subset data based on amplicons
data_all <- data_all[data_all$locus %in% selected_amps,]

# calculate heterozygosity of loci and remove those with He = 0
heterozygosity_per_sample <- data_all %>% #Compute Hₑ per sample per locus
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


variability_per_locus$locus <- as.factor(variability_per_locus$locus)

amps_mean_He <- ggplot(variability_per_locus, aes(x = reorder(locus, mean_He), y = mean_He)) +
  geom_bar(stat = "identity", color = "black") +
  # geom_errorbar(aes(ymin = mean_He - sd_He, ymax = mean_He + sd_He), 
  #               width = 0.4, color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "",
       x = "Locus",
       y = "Mean He")+
  coord_flip()

amps_mean_He

#ggsave(paste0("amps_mean_He_", site, ".png"), amps_mean_He, dpi = 300, height = 18, width = 10, bg = "white")


# variability_per_locus must be sorted descending by mean_He
variability_per_locus <- variability_per_locus %>% 
  arrange(desc(mean_He)) %>%
  mutate(prob_identity = 1 - mean_He)

# Initialize a vector to hold cumulative multilocus heterozygosity
n_loci <- nrow(variability_per_locus)
cumulative_He <- numeric(n_loci)

# Calculate cumulative multilocus heterozygosity using the product of (1 - He)
for (i in 1:n_loci) {
  prod_identity <- prod(variability_per_locus$prob_identity[1:i])
  cumulative_He[i] <- 1 - prod_identity
}

# Create a data frame for plotting
cum_curve <- data.frame(
  loci_included = 1:n_loci,
  multilocus_He = cumulative_He
)


amps_var_data <- cbind(variability_per_locus, cum_curve)

write.csv(amps_var_data, paste0("amps_variation_", site, ".csv"), row.names = F)


# Plot the cumulative diversity curve
cum_curve_plot <- ggplot(cum_curve, aes(x = loci_included, y = multilocus_He)) +
  geom_line(color = "black") +
  geom_point(color = "red") +
  labs(
    x = "Number of Loci (ordered by descending mean He)",
    y = "Cumulative Multilocus Heterozygosity",
    title = paste0("Cumulative Diversity Curve for ", site)
  ) +
  theme_minimal()

cum_curve_plot

ggsave(paste0("cumulative_diversity_curve_", site, ".png"), cum_curve_plot, dpi = 300, height = 7, width = 7, bg = "white")



##################### SELECT TOP N AMPS USING THE CUMULATIVE CURVE #############################3

# select top n amps
top_n_amps <- max(cum_curve[cum_curve$multilocus_He < 0.99,]$loci_included)

# keep top 50 most variable amps
top_var_amps <- variability_per_locus[1:top_n_amps,]$locus

# subset data based on amplicons
data_all <- data_all[data_all$locus %in% top_var_amps,]


####### 5) CHECKS --------

length(unique(data_all$sampleID)) # nidas
length(unique(data_all$locus)) # amplicons shared by all nidas

length(unique(data_all[data_all$data_type == "tes",]$sampleID))/2 # pairs



###### 6) OUTPUTS ---------

write.csv(data_all, paste0("genomic_updated_", site, ".csv"), row.names = F)
write.csv(metadata_updated, paste0("metadata_updated_", site, ".csv"), row.names = F)

