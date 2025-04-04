
library(moire)
library(dplyr)


site <- "Tete"

# ECOI = FALSE

data <- read.csv(paste0("genomic_updated_", site, ".csv"))
data <- data %>% rename(sample_id = sampleID)

data <- data[data$data_type == "tes",] # ONLY TES DATA!!!

metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))



# if(ECOI){
  
  # set MOIRE parameters
  dat_filter <- moire::load_long_form_data(data)
  burnin <- 1e4
  num_samples <- 1e4
  pt_chains <- seq(1, .5, length.out = 20)
  
  # run moire
  mcmc_results <- moire::run_mcmc(
    dat_filter, is_missing = dat_filter$is_missing,
    verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
    pt_chains = pt_chains, pt_num_threads = length(pt_chains),
    thin = 10)
  
  saveRDS(mcmc_results, paste0("coi_mcmc_", site,".RDS")) # save checkpoint ffs
  
  
  # extract ecoi and naive coi
  mcmc_results <- readRDS(paste0("coi_mcmc_", site,".RDS"))
  
  coi_stats <- merge(summarize_effective_coi(mcmc_results), summarize_coi(mcmc_results), by = "sample_id")
  coi_stats <- coi_stats %>% rename(NIDA = sample_id)
  
  write.csv(coi_stats, paste0("coi_stats_", site, ".csv"), row.names = F)
  
  
  # update metadata
  coi_stats$NIDA <- gsub("__.*", "", coi_stats$NIDA)
  metadata_updated2 <- merge(metadata_updated, coi_stats, by = c("NIDA"))
  metadata_updated2 <- metadata_updated2 %>% arrange(SampleID)
  
# }else{
#   
#   # CALCULATE ONLY NAIVE COI
#   
#   naive_coi <- data %>% group_by(sample_id, locus) %>% summarise(n_alleles=length(unique(allele))) %>% group_by(sample_id) %>% summarise(naive_coi = max(n_alleles)) %>% rename(NIDA= sample_id)
#   naive_coi$NIDA <- gsub("__.*", "", naive_coi$NIDA)
#   metadata_updated2 <- merge(metadata_updated, naive_coi, by = c("NIDA"))
#   metadata_updated2 <- metadata_updated2 %>% arrange(SampleID)
# 
# }
  
write.csv(metadata_updated2, paste0("metadata_updated_", site, ".csv"), row.names = F)
