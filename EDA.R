library(dplyr)
library(tidyr)
library(ggplot2)
library(corrplot)
library(uwot)
library(ggpubr)
library(rstatix)
library(cluster)

#----------------------------------------------------------
# Data Import & Preparation
#----------------------------------------------------------
site <- "Inhambane"

TRAINING_DATA <- read.csv(paste0(site, "_TRAINING_DATA.csv"), row.names = 1)
TRAINING_DATA$eCOI_pairs <- paste0(TRAINING_DATA$D0_nstrains, "__", TRAINING_DATA$Dx_nstrains)
LABELS <- data.frame(labels = TRAINING_DATA$labels)
LABELS$labels <- as.factor(LABELS$labels)

REAL_DATA <- read.csv(paste0(site, "_REAL_DATA.csv"),
                      stringsAsFactors = FALSE,
                      colClasses = c(NIDA1 = "character", NIDA2 = "character"))

features_to_use <- colnames(REAL_DATA)[!colnames(REAL_DATA) %in% 
                                         c("PairsID", "NIDA1", "NIDA2", "eCOI_pairs", 
                                           "locus_concordance_rate", "naive_coi_D0", "naive_coi_Dx", "pair_type" )]

#----------------------------------------------------------
# 1) FEATURE CORRELATIONS
#----------------------------------------------------------
corr_data <- TRAINING_DATA %>% select(all_of(features_to_use))
corrplot(cor(corr_data, use = "complete.obs"), method = "pie", tl.cex = 0.7)

#----------------------------------------------------------
# 2) FEATURE DISTRIBUTIONS: GLOBAL & STRATIFIED BY eCOI_pairs
#----------------------------------------------------------
# Convert to long format
training_long <- TRAINING_DATA %>%
  select(all_of(features_to_use), eCOI_pairs, labels) %>%
  pivot_longer(cols = all_of(features_to_use), names_to = "feature", values_to = "value")

# Global distributions (density plots)
ggplot(training_long, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  facet_wrap(~ feature, scales = "free") +
  theme_minimal() +
  labs(title = "", x = "Value", y = "Density")

# Stratified by eCOI_pairs
distros_strat <- ggplot(training_long, aes(x = value, fill = labels, color = labels)) +
  geom_density(alpha = 0.5) +
  facet_grid(feature~ eCOI_pairs, scales = "free") +
  theme_minimal() +
  labs(title = "", x = "Value", y = "Density")+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8) )

ggsave(paste0("distributions_stratified_", site, ".png"), distros_strat, height = 12, width = 14, dpi = 300, bg= "white")

#----------------------------------------------------------
# 3) BOXPLOTS OF FEATURES: GLOBAL & STRATIFIED BY eCOI_pairs WITH WILCOX TEST P-VALUES
#----------------------------------------------------------
# Boxplots by label (global)
p1 <- ggplot(training_long, aes(x = labels, y = value, fill = labels)) +
  geom_boxplot(alpha = 1, color = "black", outlier.alpha = 0.05) +  # hide outliers; outlier.shape = NA
  facet_grid(feature ~ eCOI_pairs, scales = "free") +
  # stat_compare_means(method = "wilcox.test", label = "p.signif", label.y.npc = "top") +
  #scale_fill_brewer(palette = "Set3") +
  theme_minimal(base_size = 14) +
  labs(title = "", 
       x = "Label", 
       y = "Value") +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "lightgrey", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8)
  )

ggsave(paste0("boxplots_stratified_", site, ".png"), p1, height = 12, width = 17, dpi = 300, bg= "white")


# p2 <- ggplot(training_long, aes(x = eCOI_pairs, y = value, fill = eCOI_pairs)) +
#   geom_boxplot(alpha = 0.8, color = "black", outlier.alpha = 0.05) +  # hide outliers; outlier.shape = NA
#   facet_wrap(~ feature, scales = "free") +
#   #stat_compare_means(method = "wilcox.test", label = "p.signif", label.y.npc = "top") +
#   scale_fill_brewer(palette = "Set3") +
#   theme_minimal(base_size = 14) +
#   labs(title = "", 
#        x = "Pair Type", 
#        y = "Value") +
#   theme(legend.position = "none",
#         strip.background = element_rect(fill = "lightgrey", color = "black"),
#         strip.text = element_text(face = "bold"),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.border = element_rect(color = "black", fill = NA, size = 0.8) )
# 
# ggsave(paste0("boxplots_global_", site, ".png"), p2, height = 8, width = 13, dpi = 300, bg= "white")



#----------------------------------------------------------
# 4) UMAP OF FEATURES: COLORED BY LABEL, SHAPED BY eCOI_pairs
#----------------------------------------------------------
# Run UMAP using uwot
umap_res <- uwot::umap(TRAINING_DATA %>% select(all_of(features_to_use)),
                       n_neighbors = 15,
                       seed = 420, 
                       n_threads = 20, 
                       n_components = 3,
                       verbose = T)

# Convert the result into a data frame and rename columns
umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2", "UMAP3")

# Add label and eCOI_pairs information
umap_df$label <- LABELS$labels
umap_df$eCOI_pairs <- TRAINING_DATA$eCOI_pairs

# Plot UMAP embedding colored by label and shaped by eCOI_pairs
umap_all <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, size = UMAP3, color = label)) +
  geom_point(size = 3, alpha = 0.3) +
  theme_minimal() +
  labs(title = "", x = "UMAP1", y = "UMAP2", 
       color = "Label", shape = "eCOI_pairs")

#umap_all

ggsave(paste0("umap_global_", site, ".png"), umap_all, height = 7, width = 8, dpi = 300, bg= "white")

#----------------------------------------------------------
# 5) UMAP OF FEATURES: COLORED BY EACH FEATURE VALUE (One Panel per Feature)
#----------------------------------------------------------
# Add each feature's value to the UMAP dataframe
for(feat in features_to_use){
  umap_df[[feat]] <- TRAINING_DATA[[feat]]
}

umap_long <- umap_df %>%
  pivot_longer(cols = all_of(features_to_use), names_to = "feature", values_to = "value")

umap_strat <- ggplot(umap_long, aes(x = UMAP1, y = UMAP2, size = UMAP3, color = value)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_grid(feature~ eCOI_pairs, scales = "free") +
  scale_color_gradient(low = "black", high = "red") +
  theme_minimal() +
  labs(title = "", x = "UMAP1", y = "UMAP2", color = "Value")+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8) )+
  ylim(min(umap_long$UMAP2),max(umap_long$UMAP2))+ xlim(min(umap_long$UMAP1),max(umap_long$UMAP1))

#umap_strat

ggsave(paste0("umap_stratified_", site, ".png"), umap_strat, height = 13, width = 15, dpi = 300, bg= "white")


# #----------------------------------------------------------
# # 6) K-MEANS CLUSTERING OF FEATURES, PLOTTED WITH eCOI_pairs AS SHAPE
# #----------------------------------------------------------
# set.seed(123)
# k_clusters <- 2   # choose number of clusters (e.g., via elbow method)
# kmeans_res <- kmeans(TRAINING_DATA %>% select(all_of(features_to_use)), 
#                      centers = k_clusters, nstart = 25)
# TRAINING_DATA$cluster <- as.factor(kmeans_res$cluster)
# 
# # Convert UMAP results to a data frame and add grouping variables
# umap_df <- as.data.frame(umap_res)
# colnames(umap_df) <- c("UMAP1", "UMAP2")
# umap_df$eCOI_pairs <- TRAINING_DATA$eCOI_pairs
# umap_df$cluster <- TRAINING_DATA$cluster
# 
# # Plot UMAP embedding with clusters and eCOI_pairs as shape
# ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
#   geom_point(size = 2, alpha = 0.8) +
#   theme_minimal() +
#   labs(title = "", 
#        x = "UMAP1", y = "UMAP2", 
#        color = "Cluster", 
#        shape = "Pair Type")

