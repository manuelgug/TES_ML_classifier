library(moire)
library(dplyr)


data <- read.csv("genomic_updated.csv")

data <- data %>% rename(sample_id = sampleID)


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


# extract ecoi

mcmc_results <- readRDS("coi_mcmc.RDS")

coi_stats <- merge(summarize_effective_coi(mcmc_results), summarize_coi(mcmc_results), by = "sample_id")

write.csv(coi_stats, "coi_stats.csv", row.names = F)
