
library(dplyr)
library(reshape2)
library(ggplot2)


site <- "Tete"


metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
data <- read.csv(paste0("genomic_updated_",site, ".csv"), stringsAsFactors = FALSE, colClasses = c(sampleID = "character"))
data <- data[data$data_type == "tes",]
coi_stats <- read.csv(paste0("coi_stats_", site, ".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))

# edit nidas
data$sampleID <- gsub("__.*","", data$sampleID)
coi_stats$NIDA <- gsub("__.*","", coi_stats$NIDA)

# keep D0 samples only
data <- left_join(data, metadata_updated[c("NIDA", "time_point")], by =c("sampleID" = "NIDA"))
#data <- data[data$time_point == "D0",]

coi_stats <- left_join(coi_stats, metadata_updated[c("NIDA", "time_point")], by =c("NIDA"))
coi_stats <- coi_stats[coi_stats$time_point == "D0",]

#allele frequencies
AF <- moire::summarize_allele_freqs(readRDS(paste0("coi_mcmc_", site, ".RDS")))
AF <- AF %>% select(locus, allele, post_allele_freqs_mean)



###### 1) SEPARATE MONOCLONAL INFECTIONS from D0 -----

clones <- coi_stats[coi_stats$naive_coi < 1.1,]$NIDA #ecoi < 1.1 #  !grepl("__", coi_stats$NIDA) avoids clones from the tes study
#clones <- coi_stats[coi_stats$naive_coi == 1,]$NIDA #naive coi = 1

clones_genomic <- data[data$sampleID %in% clones,]



###### 1) CREATE ALL POSSIBLE CLONES FROM POLYCLONAL INFECTIONS from D0 -----

N_CLONES = 1000 - length(clones) # aim for 100 clones, minus the ones already isolated above

data_polyclonal_D0 <- data[!data$sampleID %in% clones,] #subset polyclonal data

polyclonal_samples <- unique(data_polyclonal_D0$sampleID)
iterations <- ceiling(N_CLONES / length(polyclonal_samples)) # how many clones to draw from each sample


# unction to perform random allele draws for each locus using AF as weights
random_draw <- function(locus_data) {
  slice_sample(locus_data, n = 1, weight_by = post_allele_freqs_mean) %>% pull(allele)
}

# Loop over each polyclonal sample
sampled_monoclonals <- data.frame()

for (sample_id in polyclonal_samples) {

  samp <- data_polyclonal_D0[data_polyclonal_D0$sampleID == sample_id,]
  
  # Step 1: Join allele frequencies to loci in the sample dataframe (samp)
  allele_combinations <- samp %>%
    select(locus, allele) %>%
    distinct() %>%
    left_join(AF, by = c("locus", "allele"))
  
  # Step 3: Loop to generate random combinations (1 per clone) for the current sample
  set.seed(42069)  # Set the seed for reproducibility
  
  for (i in 1:iterations) {

    clone_id <- paste(sample_id, "__clone_", i, sep = "")
    
    # Step 4: Randomly draw one allele for each locus for this clone
    clone_result <- allele_combinations %>%
      group_by(locus) %>%
      summarise(allele = random_draw(pick(everything())), .groups = "drop") %>%  # Use `summarise` instead of `mutate`
      mutate(sampleID = clone_id) %>%
      select(sampleID, locus, allele)
    
    clone_result <- left_join(clone_result, samp[c("allele", "reads", "norm.reads.locus", "data_type", "time_point")], by = c("allele"))
    
    # Step 5: Bind the result to the main dataframe
    sampled_monoclonals <- bind_rows(sampled_monoclonals, clone_result)
  }
}

# merge sampled monoclonals with clones
all_clones <- rbind(sampled_monoclonals, clones_genomic)

length(unique(all_clones$sampleID))


# # check if clones are indeed monoallelic
# check <- all_clones %>% group_by(sampleID, locus) %>% summarise(length(unique(allele)))
# all(check$`length(unique(allele))`== 1)


###### 2) CALCULATE PAIRWISE PROPORTION OF SHARED ALLELES -----

# compare alleles shared between coi = 1 samples
alleles <- all_clones %>%
  group_by(sampleID) %>%
  summarize(alleles = list(allele))

# Create pairwise comparison matrix
n <- nrow(alleles)
comparison_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(alleles$sampleID, alleles$sampleID))

# Calculate pairwise proportions of shared alleles 
for (i in 1:n) {
  for (j in i:n) {
    alleles_i <- alleles$alleles[[i]]
    alleles_j <- alleles$alleles[[j]]
    shared_alleles <- length(intersect(alleles_i, alleles_j))
    total_unique_alleles <- length(unique(c(alleles_i, alleles_j)))
    shared_percentage <- (shared_alleles / total_unique_alleles)
    comparison_matrix[i, j] <- shared_percentage
    comparison_matrix[j, i] <- shared_percentage # matrix is symmetric
  }
}

# Convert the matrix to a dataframe before melting
comparison_df <- as.data.frame(comparison_matrix)

# Ensure row names are preserved as a column
comparison_df$RowName <- rownames(comparison_matrix)

# Melt while keeping row and column names
comparison_long <- reshape2::melt(comparison_df, id.vars = "RowName", variable.name = "ColumnName", value.name = "value")
colnames(comparison_long) <- c("Var1", "Var2", "value")

comparison_long <- comparison_long %>%
  arrange(value)

# remove the same comparison in the opposite order and self comparisons
comparison_long <- comparison_long[as.character(comparison_long$Var1) < as.character(comparison_long$Var2),]

#plot a histogram
hist<- ggplot(comparison_long, aes(x = value)) +
  geom_histogram(fill = "skyblue", color = "black", bins = 30) +
  labs(title = "", x = "Proportion of Shared Alleles", y = "Monoclonal Pairwise Comparisons") +
  theme_minimal()+
  xlim(0, 1)

hist

ggsave(paste0("hist_shared_alleles_",site,".png"), hist, height = 5, width = 8, bg = "white", dpi = 300)



###### 3) CLEAN MONOCLONAL DATA ------

# if 2 samples have all of their alleles shared, might as well just keep one of those to avoid repeated data 
same_clones <- comparison_long[comparison_long$value > 0.7,] # separate comparisons that are the same clone

# group samples that are the same clone into a single group
groups <- list()

for (i in seq_len(nrow(same_clones))) {
  v1 <- same_clones$Var1[i]
  v2 <- same_clones$Var2[i]
  
  found <- FALSE
  for (j in seq_along(groups)) {
    if (v1 %in% groups[[j]] || v2 %in% groups[[j]]) {
      groups[[j]] <- unique(c(groups[[j]], v1, v2))
      found <- TRUE
      break
    }
  }
  
  if (!found) {
    groups <- append(groups, list(c(v1, v2)))
  }
}

# Remove the last element from each group, this is not gonna be removed from the data (it's the representative sample from each group of identical clones)
groups <- lapply(groups, function(x) x[1])

# these are the redundant clones to remove from the clone data
clones_to_remove <- unlist(groups)

#remove redundant clones
all_clones <- all_clones[!all_clones$sampleID %in% clones_to_remove,]

length(unique(all_clones$sampleID))


###### 4) EXPORT CLONE DATA

write.csv(all_clones, paste0("clones_genomic_data_",site,".csv"), row.names = F)

