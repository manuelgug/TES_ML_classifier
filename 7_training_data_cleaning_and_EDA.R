

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(parallel)
library(corrplot)
library(patchwork)
library(ggsignif)



###### 1) INPUT AND FORMAT ----------

TRAINING_DATA <- read.csv("TRAINING_DATA_top_50_amps_30_clones.csv", row.names = 1) 
LABELS <- data.frame(labels = TRAINING_DATA$labels)
#LABELS$labels <- ifelse(FEATS_LABELS_FINAL$labels == 'NI', 0, ifelse(FEATS_LABELS_FINAL$labels == 'R', 1, FEATS_LABELS_FINAL$labels))

TRAINING_DATA_EDA <-TRAINING_DATA %>% select(1:which(names(TRAINING_DATA) == "IBD_estimate"))



######  2) CORRPLOT AND REMOVAL OF CORELATED VARIABLES ---------

# check correlations between variables
cor_var <-cor(TRAINING_DATA_EDA, method = "spearman")
corrplot_Feats <- corrplot(cor_var, method = "pie")


# Subset by correlation??
threshold <- 0.95 # 1 = all features pass

# Find pairs of variables with high correlation (both positive and negative)
high_correlations <- which(abs(cor_var) > threshold, arr.ind = TRUE)
high_correlations <- high_correlations[high_correlations[,1] < high_correlations[,2], ]

# Create a data frame of high correlations
high_correlations_df <- data.frame(
  var1 = colnames(cor_var)[high_correlations[,1]],
  var2 = colnames(cor_var)[high_correlations[,2]],
  correlation = cor_var[high_correlations]
)

# Select variables based on the high correlation data frame
TRAINING_DATA_EDA <- TRAINING_DATA_EDA[, !colnames(TRAINING_DATA_EDA) %in% high_correlations_df$var2] # COULD BE EITHER VAR1 OR VAR2, IT'S THE SAME, JUST KEEP ONE



###### 3) DISTRIBUTIONS -------------------

# Reshape the data for faceted histograms
TRAINING_DATA_EDA_long <- TRAINING_DATA_EDA %>%
  pivot_longer(cols = everything(), names_to = "Feature", values_to = "Value")

# Adjust the histogram creation to better handle distributions
distributions <- ggplot(TRAINING_DATA_EDA_long, aes(x = Value)) +
  geom_histogram(fill = "steelblue", color = "black", bins = 30) + # Adjust bin count for better resolution
  facet_wrap(~ Feature, scales = "free", ncol = 3) + # "free" scales for better visualization
  labs(title = "", x = "Value", y = "Count") +
  theme_minimal() +
  theme(strip.text = element_text(size = 8)) # Smaller facet labels for better readability

distributions

ggsave("training_data_EDA_distributions_top_50_amps_30_clones.png", distributions, bg = "white", height = 8, width = 8, dpi = 300)



###### 4) CONVERGENCE ANALYSIS -------------------

num_cores <- 17

# Function to perform bootstrapping with varying sample sizes in parallel
bootstrap_convergence <- TRAINING_DATA_EDA %>%
  pivot_longer(cols = everything(), names_to = "Feature", values_to = "Value") %>%
  group_by(Feature) %>%
  summarize(
    bootstrap_medians_mean_by_size = list(
      mclapply(
        seq(0, nrow(TRAINING_DATA_EDA), 2500),
        FUN = function(sample_size) {
          tibble(
            sample_size = sample_size,
            median_value = median(replicate(1000, median(sample(Value, size = sample_size, replace = TRUE)))),
            mean_value = mean(replicate(1000, mean(sample(Value, size = sample_size, replace = TRUE))))
          )
        },
        mc.cores = num_cores,  # Number of cores to use
        mc.preschedule = FALSE  # Optionally disable prescheduling
      )
    ),
    .groups = "drop"
  ) %>%
  unnest(bootstrap_medians_mean_by_size)

bootstrap_convergence_df <- bootstrap_convergence %>%
  unnest(bootstrap_medians_mean_by_size)

saveRDS(bootstrap_convergence, "training_data_EDA_median_bootstrap_convergence_top_50_amps_30_clones.RDS")


# Create the plot
p <- ggplot(bootstrap_convergence_df, aes(x = sample_size, y = mean_value)) +
  geom_line() +
  facet_wrap(~Feature, scales = "free") +
  labs(
    title = "",
    x = "Sample Size",
    y = "Bootstrapped Mean"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )

p

# Save the plot
ggsave("training_data_EDA_median_bootstrap_convergence_top_50_amps_30_clones.png", plot = p, width = 12, height = 10, dpi = 300, bg = "white")



### 5) DIFFERENCES IN FEATURES BETWEEN CLASSES ------

# Function to perform Wilcoxon test and calculate medians
perform_w_test_and_medians <- function(variable) {
  w_test_result <- wilcox.test(variable ~ LABELS$labels)
  median1 <- median(variable[LABELS$labels == "R"], na.rm = TRUE)
  median2 <- median(variable[LABELS$labels == "NI"], na.rm = TRUE)
  return(c(P_Value = w_test_result$p.value, Median1 = median1, Median2 = median2))
}

# Apply the function to each variable
p_values_and_medians <- TRAINING_DATA_EDA %>% 
  map_dfr(perform_w_test_and_medians)

wilcox_results <- data.frame(Feature = colnames(TRAINING_DATA_EDA), p_values_and_medians)

plot_list <- list()

for (feature in wilcox_results$Feature) {
  p_value <- wilcox_results$P_Value[wilcox_results$Feature == feature]
  
  # Determine significance level based on p-value
  if (p_value < 0.001) {
    signif_label <- "***"
  } else if (p_value < 0.01) {
    signif_label <- "**"
  } else if (p_value < 0.05) {
    signif_label <- "*"
  } else {
    signif_label <- "ns"
  }
  
  p <- ggplot(cbind(TRAINING_DATA_EDA, LABELS), aes_string(x = "labels", y = feature, fill = "labels")) +
    geom_boxplot() +
    labs(title = "",
         x = "",
         y = feature) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("NI" = "green", "R" = "red")) +
    guides(legend = "none", fill = "none") +
    
    # Add Wilcoxon significance bracket with asterisks
    geom_signif(comparisons = list(c("NI", "R")), 
                annotations = signif_label,
                y_position = max(TRAINING_DATA_EDA[[feature]], na.rm = TRUE) * 1.1,
                tip_length = 0.03) +
    
    # Increase ylim to ensure space for annotations
    ylim(0, max(TRAINING_DATA_EDA[[feature]], na.rm = TRUE) * 1.3)  # Increase y-axis limit for space
  
  plot_list[[feature]] <- p
}

var_boxplots <- wrap_plots(plot_list, ncol = round(sqrt(length(wilcox_results$Feature)), 0)) 
var_boxplots

ggsave("training_data_EDA_boxplots_top_50_amps_30_clones.png", var_boxplots, width = 12, height = 12, bg = "white")

# remove non significant features
non_significant_features <- wilcox_results[wilcox_results$P_Value > 0.05,]$Feature

TRAINING_DATA_EDA <- TRAINING_DATA_EDA[,!colnames(TRAINING_DATA_EDA) %in% non_significant_features]



###### 6) UMAP --------

# scale data
data_scaled <- scale(TRAINING_DATA_EDA)

#umap
set.seed(42069)
umap_result <- uwot::umap(data_scaled, n_components = 3, verbose = T, n_neighbors = 15, n_threads = 20)

umap_df <- data.frame(UMAP1 = umap_result[,1],
                      UMAP2 = umap_result[,2],
                      UMAP3 = umap_result[,3])

#add classifications
umap_df <- cbind(umap_df, labels = as.factor(LABELS$labels))


# Create a new dataframe with different combinations of UMAP axes
umap_combined <- bind_rows(
  umap_df %>% mutate(x = UMAP1, y = UMAP2, z = UMAP3, panel = "UMAP1 vs UMAP2"),
  umap_df %>% mutate(x = UMAP1, y = UMAP3, z = UMAP2, panel = "UMAP1 vs UMAP3"),
  umap_df %>% mutate(x = UMAP2, y = UMAP3, z = UMAP1, panel = "UMAP2 vs UMAP3")
)

p <- ggplot(umap_combined, aes(x = x, y = y, color = labels, shape = labels)) +
  geom_jitter(aes(size = z), alpha = 0.1, width=0, height =0) +
  scale_size_continuous(range = c(1, 5)) +
  labs(x = "Axis 1",
       y = "Axis 2") +
  theme_minimal() +
  scale_color_manual(values = c("NI" = "green", "R" = "red")) +
  scale_shape_manual(values = c("NI" = 19, "R" = 19))+
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.box = "horizontal") +
  guides(size = "none") +
  facet_wrap(~ panel, scales = "free", ncol=2)
p

ggsave("training_data_EDA_umap_top_50_amps_30_clones.png", p, width = 10, height = 8, bg = "white")



#per variable coloring
variables <- colnames(TRAINING_DATA_EDA)

plot_list <- list()

for (var in variables) {
  pl <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes_string(color = var, size = "UMAP3"), alpha = 0.25, data = cbind(umap_df, TRAINING_DATA_EDA[var])) +
    scale_color_gradient(low = "green", high = "red",) +
    scale_size_continuous(range = c(1, 3)) +
    labs(title = var,
         x = "UMAP1",
         y = "UMAP2",
         color = var) +
    theme_minimal() +
    theme(legend.position = "right",
          legend.title = element_blank(),  # This removes the legend title
          legend.text = element_text(size = 10),
          legend.box = "horizontal") +
    guides(size = "none")
  
  plot_list[[var]] <- pl
}

final_plot <- wrap_plots(plot_list, ncol = round(sqrt(length(variables)), 0))

ggsave("training_data_EDA_umap_feats_top_50_amps_30_clones.png", final_plot, width = 18, height = 12, bg = "white")




###### 7) K-MEANS ----------

K_RANGE <- c(2, 3, 4)

for (k in K_RANGE){
  
  print(paste0("Processing k = ", k))
  
  # Apply k-means clustering
  set.seed(69)
  kmeans_result <- kmeans(data_scaled, centers = k, nstart = 10, iter.max = 500)
  
  umap_df$Cluster <- as.factor(kmeans_result$cluster)
  # umap_df$labels <- LABELS$labels
  
  # # describe the clusters
  clusters_df <- cbind(Cluster = umap_df$Cluster, TRAINING_DATA_EDA)
  
  # Convert data to long format for easier plotting
  long_df <- clusters_df %>%
    pivot_longer(cols = -Cluster, names_to = "Variable", values_to = "Value")
  
  cluster_description <- ggplot(long_df, aes(x = as.factor(Cluster), y = Value, fill = as.factor(Cluster) )) +
    #geom_jitter(width = 0.1, aes(color = as.factor(Cluster)))+
    geom_boxplot(alpha = 0.5, outliers = F) +
    facet_wrap(~ Variable, scales = "free_y") +  # Create a facet for each variable
    labs(title = "Cluster summary",
         x = "Cluster",
         y = "Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),  
      strip.text = element_text(size = 7, face = "bold", hjust = 0),  
      strip.placement = "top",  
      panel.background = element_rect(fill = "white", color = "black"), 
      panel.grid.major = element_line(color = "grey90"), 
      panel.grid.minor = element_line(color = "grey90"),
      panel.border = element_rect(color = "black", fill = NA),
      legend.position = "none"
    )
  
  # Visualize the clusters on the UMAP projection
  umap_cluster <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster, size = UMAP3, shape = as.factor(labels))) + 
    geom_point(alpha = 0.1) +
    scale_size_continuous(range = c(2, 5)) +
    labs(title = paste("k =", k),
         x = "UMAP1",
         y = "UMAP2",
         color = "Cluster") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 10), 
          legend.text = element_text(size = 10),
          legend.box = "horizontal") +   
    guides(size = "none", 
           color = guide_legend(title = "Cluster", 
                                title.position = "top",  
                                override.aes = list(size = 7, shape = 19)))
  
  # output plots
  panels <- gridExtra::arrangeGrob(umap_cluster, cluster_description, ncol = 2)
  ggsave(paste0("training_data_EDA_KMEANS_CLEAN_k_", k,"_top_50_amps_30_clones.png"), panels, width = 15, height = 7.5, bg = "white")  
  
}
  


###### 8) CLEAN EXPORTS ----------

before_cols <- colnames(TRAINING_DATA[,1:which(names(TRAINING_DATA) == "IBD_estimate")]) # assumes ibd_Estimate is kept!
after_cols <- colnames(TRAINING_DATA_EDA[,1:which(names(TRAINING_DATA_EDA) == "IBD_estimate")]) # assumes ibd_Estimate is kept!

cols_to_remove <- setdiff(before_cols, after_cols)

#remove shit features
TRAINING_DATA <- TRAINING_DATA %>% select(-cols_to_remove)

#order columns
TRAINING_DATA$PairsID <- c(1:nrow(TRAINING_DATA))

TRAINING_DATA <- TRAINING_DATA %>% select(PairsID, everything())

write.csv(TRAINING_DATA, "TRAINING_DATA_CLEAN_top_50_amps_30_clones.csv", row.names = F)

