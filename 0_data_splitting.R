

library(dplyr)
library(stringr)


## 1) split metadata by geographic location

metadata_file <- "Final_list_TES_Envio_ISG_INPUT_TABLE.csv" # should have the nidas, the timepoints and the pairs to which samples belong to

#import metadata
metadata <- read.csv(metadata_file)
metadata$NIDA <- as.character(metadata$NIDA)

metadata <- metadata %>%
  mutate(NIDA = gsub("N|_S.*$", "", NIDA),  # Remove "N" and everything after "_S"
         NIDA = if_else(!str_detect(NIDA, "\\."), paste0(NIDA, ".0"), NIDA))  # If no ".", add ".0"


# Get unique province names
unique_provinces <- unique(metadata$Province)

# Loop through each province and save as a CSV file
for (province in unique_provinces) {
  subset_data <- metadata %>% filter(Province == province)
  write.csv(subset_data, paste0("metadata_tes_", province, ".csv"), row.names = FALSE)
  }



## 2) collect and import 2022 data from simone's paper

main_dir = "SITE_DATA/results_v0.1.8_RESMARKERS_FIX/"

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

genomic_site <- bind_rows(list_of_dfs)

# format nidas
genomic_site <- genomic_site %>%
  mutate(sampleID = gsub("N|_S.*$", "", sampleID),  # Remove "N" and everything after "_S"
         sampleID = gsub("-", "_", sampleID), # Change - for _
         sampleID = if_else(str_detect(sampleID, "_"), sampleID, paste0(sampleID, ".0")),
         sampleID = gsub("_", ".", sampleID))  # change "_" for "."

# site metadata
site_db <- haven::read_dta(paste0(dirname(main_dir), "/",  "DrugRes_2021_2022_DB_ALLDATA_22Jun2024.dta"))
site_db <- as_tibble(site_db)
site_db <- site_db[site_db$year == 2022,] 
site_db$province <- gsub(" ", "_", site_db$province)
site_db <- site_db %>% rename(Province = province)

site_db <- site_db %>%
  mutate(NIDA2 = gsub("N|_S.*$", "", NIDA2),  # Remove "N" and everything after "_S"
         NIDA2 = if_else(!str_detect(NIDA2, "\\."), paste0(NIDA2, ".0"), NIDA2))  # If no ".", add ".0"


# merge metadata with genomic
genomic_metadata <- left_join(site_db, genomic_site, by = c("NIDA2" = "sampleID")) %>% select(NIDA2, Province, locus, pseudo_cigar, reads, norm.reads.locus, run)

#specify data type = site
genomic_metadata$data_type <- "site"


# Loop through each province and save as a CSV file
for (province in unique_provinces) {
  subset_nidas <- site_db %>% filter(Province == province) %>% select(NIDA2)
  subset_genomic <- genomic_metadata[genomic_metadata$NIDA2 %in% subset_nidas$NIDA2 & genomic_metadata$Province == province,]
  write.csv(subset_genomic, paste0("genomic_site_", province, ".csv"), row.names = FALSE)
}
