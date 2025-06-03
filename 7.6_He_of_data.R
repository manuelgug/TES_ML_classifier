library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(ggridges)

sites <- c("Cabo_Delgado", "Inhambane", "Tete", "Zambezia")

# Function to process one site
calc_he_for_site <- function(site) {
  # Load data
  real_data <- read_csv(paste0("genomic_updated_", site, ".csv"))
  metadata_updated <- read_csv(
    paste0("metadata_updated_", site, ".csv"),
    col_types = cols(NIDA = col_character())
  )
  
  # Fix missing time points
  metadata_updated$time_point <- ifelse(
    is.na(metadata_updated$time_point), "D0", metadata_updated$time_point
  )
  
  # Subset D0 samples
  D0_nidas <- metadata_updated %>%
    filter(time_point == "D0") %>%
    pull(NIDA)
  
  real_data <- real_data %>%
    filter(sampleID %in% D0_nidas) %>%
    select(sampleID, locus, allele)
  
  # Calculate He per sampleID/locus
  he <- real_data %>%
    group_by(sampleID, locus, allele) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(freq = count / sum(count)) %>%
    summarise(He = 1 - sum(freq^2), .groups = "drop") %>%
    mutate(site = site)  # tag site for plotting
  
  return(he)
}

# Run for all sites and combine
he_all_sites <- map_dfr(sites, calc_he_for_site)

he_all_sites <- he_all_sites %>% group_by(site, locus) %>% summarise(mean_He = mean(He))

# Compute mean He per site for ordering
site_order <- he_all_sites %>%
  group_by(site) %>%
  summarise(mesian_He = median(mean_He, na.rm = TRUE)) %>%
  arrange(desc(mesian_He)) %>%
  pull(site)

# Convert site to factor with desired order
he_all_sites$site <- factor(he_all_sites$site, levels = site_order)

# Ridge plot
library(viridis)

ridgeplot <- ggplot(he_all_sites, aes(x = mean_He, y = site, fill = site)) +
  ggridges::geom_density_ridges(scale = 2, rel_min_height = 0.01, alpha = 0.7) +
  scale_fill_viridis_d(option = "D") +  # Clean categorical color palette
  labs(
    title = "",
    x = "Genome-wide He",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(paste0("ridgeplot_all_sites", site, ".png"), ridgeplot, height = 6, width = 8, dpi = 300, bg= "white")


# 
# # individual mixes that create the pairs and train/test data
# mixes <- readRDS(paste0("MIXES_GENOMIC_", site, ".RDS"))
# mixes <- mixes %>% select(mixID, locus, allele)
# 
# he_by_locus_mix <- mixes %>%
#   group_by(mixID, locus, allele) %>%
#   summarise(count = n(), .groups = "drop_last") %>%
#   mutate(freq = count / sum(count)) %>%
#   summarise(
#     He = 1 - sum(freq^2),
#     .groups = "drop"
#   )
# 
# he_by_locus_mix %>% group_by(locus) %>% summarise(mean(He))
