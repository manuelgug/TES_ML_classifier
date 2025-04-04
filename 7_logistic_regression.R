
library(caret)    
library(dplyr)    
library(tidyr)    
library(ggplot2)  
library(broom)


site <- "Tete"



### 1) IMPORT TRAINING AND REAL DATA ----------

TRAINING_DATA <- read.csv(paste0(site, "_TRAINING_DATA.csv"), row.names = 1) 
TRAINING_DATA$eCOI_pairs <- paste0(TRAINING_DATA$D0_nstrains, "__", TRAINING_DATA$Dx_nstrains) # Create a new variable by combining 'D0nstrains' and 'Dxnstrains'
LABELS <- data.frame(labels = TRAINING_DATA$labels)
LABELS$labels <- as.factor(LABELS$labels)

REAL_DATA <- read.csv(paste0(site, "_REAL_DATA.csv"), stringsAsFactors = FALSE, colClasses = c(NIDA1 = "character", NIDA2= "character")) 



### 2) SPLIT DATA ------------

# Step 1: Stratified Sampling by `pair_type`
train_indices <- createDataPartition(TRAINING_DATA$eCOI_pairs, p = 0.7, list = FALSE)
train_data <- TRAINING_DATA[train_indices, ]
test_data <- TRAINING_DATA[-train_indices, ]

# Step 2: Extract Metadata and Labels
TRAIN_META <- train_data %>% select(-IBD_estimate, -Jaccard)
TRAIN <- train_data %>% select(IBD_estimate, Jaccard)
TRAIN_labels <- LABELS[train_indices, ]

TEST_META <- test_data %>% select(-IBD_estimate, -Jaccard)
TEST <- test_data %>% select(IBD_estimate, Jaccard)
TEST_labels <- LABELS[-train_indices, ]

# Step 3: Verify Stratification
prop.table(table(TRAIN_META$eCOI_pairs))
prop.table(table(TRAIN_labels))
table(TRAIN_labels)

prop.table(table(TEST_META$eCOI_pairs))
prop.table(table(TEST_labels))
table(TEST_labels)



##### 4) TEST MODEL USING IBD ONLY (BEST FEATURE)--------

# Create training and test data frames for the current feature
df_train_IBD <- data.frame(TRAIN[c("IBD_estimate", "Jaccard")], label = as.factor(TRAIN_labels))
df_test_IBD <- data.frame(TEST[c("IBD_estimate", "Jaccard")], label = as.factor(TEST_labels))

# Fit a simple logistic regression on the training set
fit_IBD <- glm(label ~ ., data = df_train_IBD, family = binomial)
summary(fit_IBD)

# Predict probabilities on the test (holdout) set
preds_prob <- predict(fit_IBD, newdata = df_test_IBD, type = "response")

# Define the range of decision_thresholds
decision_thresholds <- seq(0, 1, by = 0.05)

# Initialize an empty results dataframe
results <- data.frame(eCOI_pairs = character(),
                      decision_threshold = numeric(),
                      sensitivity = numeric(),
                      specificity = numeric(),
                      R_pairs = numeric(),
                      NI_pairs = numeric(),
                      stringsAsFactors = FALSE)

#decision_thresholds <- 0.5 # if wanting to use only 0.5 decision threshold for all pair types...

# Loop through each decision_threshold
for (thresh in decision_thresholds) {
  
  # Convert probabilities to binary predictions at the current decision_threshold
  preds <- ifelse(preds_prob >= thresh, "R", "NI")
  
  # Loop through each unique eCOI_pairs combination
  for (strain_comb in unique(TEST_META$eCOI_pairs)) {
    
    # Subset the TEST and TEST_labels based on the current combination
    subset_indices <- TEST_META$eCOI_pairs == strain_comb
    subset_TEST_labels <- TEST_labels[subset_indices]
    
    # Count R and NI pairs
    r <- sum(subset_TEST_labels == "R", na.rm = TRUE)
    ni <- sum(subset_TEST_labels == "NI", na.rm = TRUE)
    
    # Get the corresponding predictions for the current subset
    subset_preds <- preds[subset_indices]
    
    # Evaluate the confusion matrix for the current subset
    cm <- confusionMatrix(as.factor(subset_preds), as.factor(subset_TEST_labels), positive = "R")
    
    sens <- cm$byClass["Sensitivity"]
    spec <- cm$byClass["Specificity"]
    
    # Append results
    results <- rbind(results, data.frame(eCOI_pairs = strain_comb,
                                         decision_threshold = thresh,
                                         sensitivity = sens,
                                         specificity = spec,
                                         R_pairs = r,
                                         NI_pairs = ni,
                                         stringsAsFactors = FALSE))
  }
}


# Select the best decision_threshold per eCOI_pair based on balance between sensitivity and specificity
best_decision_thresholds <- results %>%
  mutate(youden_j = sensitivity + specificity - 1) %>%  # Compute Youden’s J
  group_by(eCOI_pairs) %>%
  slice_max(youden_j) %>%  # Select rows with the largest Youden's J
  #slice_min(decision_threshold, with_ties = FALSE) %>%  # If ties, pick the lowest decision_threshold
  slice_min(abs(decision_threshold - 0.5), with_ties = FALSE) %>%  # Pick decision_threshold closest to 0.5
ungroup() 


print(best_decision_thresholds)


# Reshape data to long format for easy plotting
best_decision_thresholds_long <- best_decision_thresholds %>%
  select(eCOI_pairs, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity), 
               names_to = "Metric", 
               values_to = "Value")

# Create bar plot
metrics <- ggplot(best_decision_thresholds_long, aes(x = eCOI_pairs, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +  # Dodge separates bars for clarity
  labs(title = "",
       x = "eCOI Pairs",
       y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black")

metrics


# save model and model results
saveRDS(fit_IBD, paste0("LogReg_model_", site, ".RDS"))
write.csv(best_decision_thresholds, paste0("training_Results_LogReg_", site, ".csv"), row.names = F)
ggsave(paste0(site, "_model_results_LogReg.png"), metrics, bg = "white", dpi = 300, height = 5, width = 7)



# Reshape data into long format for plotting
results_long <- results %>%
  pivot_longer(cols = c(sensitivity, specificity), 
               names_to = "Metric", 
               values_to = "Value")

# Merge best_decision_thresholds to get the specific decision_threshold for each eCOI_pair
best_decision_thresholds_long <- best_decision_thresholds %>%
  select(eCOI_pairs, decision_threshold, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity), 
               names_to = "Metric", 
               values_to = "Value")

# Plot both Sensitivity and Specificity in the same graph
sens_spec_plot <- ggplot(results_long, aes(x = decision_threshold, y = log(Value), linetype = Metric)) +
  geom_line() +
  #geom_point(data = best_decision_thresholds_long, aes(x = decision_threshold, y = log(Value)), shape = 19, size = 3, stroke = 1.5) + 
  geom_vline(data = best_decision_thresholds, aes(xintercept = decision_threshold), color = "red", linetype = "solid") +
  facet_wrap(~eCOI_pairs) +
  labs(title = "",
       x = "Decision Threshold",
       y = "Log(Value)") +
  theme_minimal()

sens_spec_plot



###### 5) BENCHMAKR AGAINST A DUMMY RANDOM CLASSIFIER ----------------

set.seed(420) 

# Create a dummy classifier that randomly assigns labels based on class distribution
class_probs <- prop.table(table(TRAIN_labels))
TEST$dummy_random <- sample(names(class_probs), size = nrow(TEST), replace = TRUE, prob = class_probs)

# Evaluate performance of both dummy classifiers per eCOI_pair
dummy_results <- data.frame(eCOI_pairs = character(),
                            model = character(),
                            sensitivity = numeric(),
                            specificity = numeric())

for (strain_comb in unique(TEST_META$eCOI_pairs)) {
  
  subset_indices <- TEST_META$eCOI_pairs == strain_comb
  subset_TEST_labels <- TEST_labels[subset_indices]
  
  preds <- TEST$dummy_random[subset_indices]
  
  cm <- confusionMatrix(as.factor(preds), as.factor(subset_TEST_labels), positive = "R")
  
  dummy_results <- rbind(dummy_results, data.frame(
    eCOI_pairs = strain_comb,
    model = "dummy",
    sensitivity = cm$byClass["Sensitivity"],
    specificity = cm$byClass["Specificity"]
  ))
}

# # Combine dummy and LR model results
# best_decision_thresholds$model <- "LogReg"
# comparison_results <- bind_rows(best_decision_thresholds[c("eCOI_pairs", "sensitivity", "specificity", "model")], dummy_results)
# comparison_results <- comparison_results[comparison_results$eCOI_pairs %in% unique(TRAINING_DATA$eCOI_pairs),]
# 
# # Plot comparison
# dummy_comparison <- ggplot(comparison_results, aes(x = eCOI_pairs, y = sensitivity, fill = model)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(title = "",
#        x = "eCOI Pairs",
#        y = "Sensitivity") +
#   theme_minimal() +
#   scale_fill_manual(values = c("LogReg" = "#008080", "dummy_random" = "#002080")) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_hline(yintercept = 0.9, linetype = "dashed", color = "black")
# 
# dummy_comparison
# 
# ggsave(paste0(site, "_sensitivity_dummy_model_comparison.png"), dummy_comparison, bg = "white", dpi = 300, height = 5, width = 7)

# Reshape data to long format for easy plotting
dummy_results_long <- dummy_results %>%
  select(eCOI_pairs, sensitivity, specificity) %>%
  pivot_longer(cols = c(sensitivity, specificity), 
               names_to = "Metric", 
               values_to = "Value")

dummy_comparison <- ggplot(dummy_results_long, aes(x = eCOI_pairs, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +  # Dodge separates bars for clarity
  labs(title = "",
       x = "eCOI Pairs",
       y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = c("sensitivity" = "#008080", "specificity" = "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black")+
  ylim(0,1)

dummy_comparison

ggsave(paste0(site, "_sensitivity_dummy_model_comparison.png"), dummy_comparison, bg = "white", dpi = 300, height = 5, width = 7)



##### 6) TEST MODEL ON REAL DATA USING IBD ONLY (BEST FEATURE) --------

# cover both directions 1__2 and 2__1
best_decision_thresholds_rev <- best_decision_thresholds %>%
  mutate(eCOI_pairs = sapply(strsplit(eCOI_pairs, "__"), function(x) paste0(rev(x), collapse = "__")))

best_decision_thresholds <- unique(rbind(best_decision_thresholds, best_decision_thresholds_rev))

#best_decision_thresholds$decision_threshold <- 0.5 # if wanting to use 0.5 for all real samples

#add threhold data
REAL_DATA <- left_join(REAL_DATA, best_decision_thresholds[c("eCOI_pairs", "decision_threshold")], by = c("eCOI_pairs"))

# Initialize a vector to store predictions
REAL_DATA$prediction_prob <- NA
REAL_DATA$predictions <- NA

# Loop over each row in REAL_DATA to apply the corresponding decision_threshold
for (i in 1:nrow(REAL_DATA)) {
  # Subset the current row's IBD_estimate and decision_threshold
  ibd_value <- REAL_DATA$IBD_estimate[i]
  jaccard <- REAL_DATA$Jaccard[i]
  decision_threshold <- REAL_DATA$decision_threshold[i]
  
  # Create a new data frame for prediction (based on the current decision_threshold)
  newdata <- data.frame(IBD_estimate = ibd_value, Jaccard = jaccard)
  
  # Predict probabilities using the fitted logistic regression model
  prediction_prob <- predict(fit_IBD, newdata = newdata, type = "response")
  
  # Classify using the current decision_threshold (instead of 0.5)
  prediction_class <- ifelse(prediction_prob >= decision_threshold, "R", "NI")
  
  # Store the prediction in the REAL_DATA data frame
  REAL_DATA$prediction_prob[i] <- prediction_prob
  REAL_DATA$predictions[i] <- prediction_class

}

REAL_DATA <- REAL_DATA %>% select(PairsID, NIDA1, NIDA2, eCOI_pairs, everything()) %>% arrange(PairsID)

print(REAL_DATA)


## OUTPUT RESULTS 

write.csv(REAL_DATA, paste0(site, "_REAL_DATA_PREDICTIONS.csv"), row.names = F)

ggsave(paste0(site, "_REAL_DATA_THRESHOLDS_PLOT.png"), sens_spec_plot, bg = "white", dpi = 300, height = 9, width = 12)

