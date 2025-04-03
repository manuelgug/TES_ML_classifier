
library(moire)
library(dplyr)


site <- "Zambezia"


data <- read.csv(paste0("genomic_updated_", site, ".csv"))
data <- data %>% rename(sample_id = sampleID)

metadata_updated <- read.csv(paste0("metadata_updated_", site, ".csv"), stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))


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
metadata_updated <- merge(metadata_updated, coi_stats, by = c("NIDA"))
metadata_updated <- metadata_updated %>% arrange(SampleID)

write.csv(metadata_updated, paste0("metadata_updated_", site, ".csv"), row.names = F)

