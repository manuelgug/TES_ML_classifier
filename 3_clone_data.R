
library(dplyr)
library(reshape2)
library(ggplot2)



coi_stats <- read.csv("coi_stats.csv")
data <- read.csv("genomic_updated.csv")



###### 1) SEPARATE MONOCLONAL INFECTIONS -----

clones <- coi_stats[coi_stats$post_effective_coi_med < 1.1,]$NIDA #ecoi < 1.1
#clones <- coi_stats[coi_stats$naive_coi == 1,]$NIDA #naive coi = 1

clones_genomic <- data[data$sampleID %in% clones,]


###### 2) CALCULATE PAIRWISE PROPORTION OF SHARED ALLELES -----

# compare alleles shared between coi = 1 samples
alleles <- clones_genomic %>%
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

comparison_long <- melt(comparison_matrix)
comparison_long <- comparison_long %>%
  arrange(value)

# remove the same comparison in the opposite order and self comparisons
comparison_long <- comparison_long[as.character(comparison_long$Var1) < as.character(comparison_long$Var2),]

#plot a histogram
hist<- ggplot(comparison_long, aes(x = value)) +
  geom_histogram(fill = "skyblue", color = "black", bins = 100) +
  labs(title = "", x = "Proportion of Shared Alleles", y = "Monoclonal Pairwise Comparisons") +
  theme_minimal()+
  xlim(0, 1)

hist

ggsave("hist_shared_alleles.png", hist, height = 5, width = 8, bg = "white", dpi = 300)



###### 3) CLEAN MONOCLONAL DATA ------

# if 2 samples have all of their alleles shared, might as well just keep one of those to avoid repeated data 
same_clones <- comparison_long[comparison_long$value == 1,] # separate comparisons that are the same clone

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

# Remove the first element from each group, this is not gonna be removed from the data (it's the representative sample from each group of identical clones)
groups <- lapply(groups, function(x) x[-1])

# these are the redundant clones to remove from the clone data
clones_to_remove <- unlist(groups)

#remove redundant clones
clones_genomic <- clones_genomic[!clones_genomic$sampleID %in% clones_to_remove,]



###### 4) EXPORT CLONE DATA

write.csv(clones_genomic, "clones_genomic_data.csv", row.names = F)
