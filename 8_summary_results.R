


library(dplyr)
library(tidyr)
library(png)
library(grid)
library(gridExtra)



sites = c("Cabo_Delgado", "Inhambane", "Tete", "Zambezia")



##### 1) SUMMARY RESULTS -----

summary_Results <- data.frame()

for (site in sites){
  
  data <- read.csv(paste0("genomic_updated_",site, ".csv"), stringsAsFactors = FALSE, colClasses = c(sampleID = "character"))
  clones_genomic <- read.csv(paste0("clones_genomic_data_",site,".csv"), stringsAsFactors = FALSE, colClasses = c(sampleID = "character"))
  metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
  PAIRS_SUMMARY <- read.csv(paste0("PAIRS_SUMMARY_",site,".csv"))
  
  
  # amount of regional samples used
  N_REGIONAL <- length(unique(data[data$data_type == "site",]$sampleID))
  
  # amount of amplicons that had multilocus He up to 0.99
  N_LOCI <- length(unique(data$locus))
  
  # amount of clones across tes and site(regional) data used in the study
  N_CLONES <- length(unique(clones_genomic$sampleID)) 
  
  # max (rounded to int) ecoi across tes samples
  MAX_COI <- round(max(metadata_updated$post_effective_coi_med))
  
  # amount of mixes used in the study
  N_MIXES <- sum(PAIRS_SUMMARY$unique_mixes)
  
  # amount of simulated pairs
  N_PAIRS <- sum(PAIRS_SUMMARY$n_pairs)
  
  # amount of unique pair types
  N_PAIR_TYPES <- nrow(PAIRS_SUMMARY)
  
  # amount of simulated pairs classified as NI
  NI_PROP <- sum(PAIRS_SUMMARY$NI_size) / N_PAIRS
  
  # amount of simulated pairs classified as NI
  R_PROP <- sum(PAIRS_SUMMARY$R_size) / N_PAIRS
  
  # row 
  data_row <- data.frame(site, N_REGIONAL, N_LOCI, N_CLONES, MAX_COI, N_MIXES, N_PAIRS, N_PAIR_TYPES, NI_PROP, R_PROP)
  
  summary_Results<- rbind(summary_Results, data_row)
  
}


write.csv(summary_Results, "SUMMARY_RESULTS_ACROSS_SITES.csv", row.names = F)



##### 2) MODEL PERFORMANCE -----

labels <- c("A", "B", "C", "D")

image_list <- lapply(1:length(sites), function(i) {
  img_path <- paste0(sites[i], "_model_results_LogReg.png")
  img <- rasterGrob(readPNG(img_path), interpolate = TRUE)
  
  gTree(children = gList(
    img,
    textGrob(labels[i], x = unit(0.05, "npc"), y = unit(0.95, "npc"),
             gp = gpar(fontsize = 20, fontface = "bold"))
  ))
})

# Save the 2x2 panel as a PNG
png("MODEL_PERFORMANCE_RESULTS_ACROSS_SITES.png", width = 8, height = 8, units = "in", res = 300)

grid.arrange(grobs = image_list, ncol = 2, nrow = 2)

dev.off()



##### 3) DUMMY MODEL PERFORMANCE -----

labels <- c("A", "B", "C", "D")

image_list <- lapply(1:length(sites), function(i) {
  img_path <- paste0(sites[i], "_sensitivity_dummy_model_comparison.png")
  img <- rasterGrob(readPNG(img_path), interpolate = TRUE)
  
  gTree(children = gList(
    img,
    textGrob(labels[i], x = unit(0.05, "npc"), y = unit(0.95, "npc"),
             gp = gpar(fontsize = 20, fontface = "bold"))
  ))
})

# Save the 2x2 panel as a PNG
png("DUMMY_MODEL_PERFORMANCE_RESULTS_ACROSS_SITES.png", width = 8, height = 8, units = "in", res = 300)

grid.arrange(grobs = image_list, ncol = 2, nrow = 2)

dev.off()




##### 3) PREDICTION RESULTS -----

predictions_Results <- data.frame()

for (site in sites){
  
  partial_df <- read.csv(paste0(site, "_REAL_DATA_PREDICTIONS.csv"), stringsAsFactors = FALSE, colClasses = c(NIDA1 = "character", NIDA2="character"))

  partial_df$site <- site

  metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))
  
  partial_df <- left_join(partial_df, metadata_updated[c("PairsID", "Study")], by = c("PairsID"))
  
  predictions_Results <- rbind(predictions_Results, partial_df)

}


predictions_Results <- predictions_Results %>% rename(pair_type = eCOI_pairs) %>% select(PairsID, Study, site, everything(), -NIDA1, -NIDA2) %>% arrange(PairsID) %>% distinct()
predictions_Results$IBD_estimate <- as.character(predictions_Results$IBD_estimate)

write.csv(predictions_Results, "PREDICTION_RESULTS_ACROSS_SITES.csv", row.names = F)


