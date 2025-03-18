
library(moire)
library(dplyr)


data <- read.csv("genomic_updated.csv")
data <- data %>% rename(sample_id = sampleID)

metadata_updated <- read.csv("metadata_updated.csv", stringsAsFactors = FALSE, colClasses = c(NIDA = "character"))


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

saveRDS(mcmc_results, "coi_mcmc.RDS") # save checkpoint ffs


# extract ecoi and naive coi
mcmc_results <- readRDS("coi_mcmc.RDS")

coi_stats <- merge(summarize_effective_coi(mcmc_results), summarize_coi(mcmc_results), by = "sample_id")
coi_stats <- coi_stats %>% rename(NIDA = sample_id)

write.csv(coi_stats, "coi_stats.csv", row.names = F)


# update metadata
coi_stats$NIDA <- gsub("__.*", "", coi_stats$NIDA)
metadata_updated <- merge(metadata_updated, coi_stats, by = c("NIDA"))
metadata_updated <- metadata_updated %>% arrange(SampleID)

write.csv(metadata_updated, "metadata_updated.csv", row.names = F)

