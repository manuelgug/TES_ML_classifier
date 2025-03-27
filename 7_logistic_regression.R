


library(caret)    
library(dplyr)    
library(tidyr)    
library(ggplot2)  
library(broom)



site <- "Cabo_Delgado"
top_n_amps <- 50 



### 1) IMPORT TRAINING AND REAL DATA ----------

TRAINING_DATA <- read.csv(paste0(site, "_TRAINING_DATA_top_",top_n_amps,"_amps.csv"), row.names = 1) 
TRAINING_DATA$eCOI_pairs <- paste0(TRAINING_DATA$D0_nstrains, "__", TRAINING_DATA$Dx_nstrains) # Create a new variable by combining 'D0nstrains' and 'Dxnstrains'
LABELS <- data.frame(labels = TRAINING_DATA$labels)
LABELS$labels <- as.factor(LABELS$labels)

REAL_DATA <- read.csv(paste0(site, "_REAL_DATA_top_",top_n_amps,"_amps.csv"), stringsAsFactors = FALSE, colClasses = c(NIDA1 = "character", NIDA2= "character")) 



### 2) SPLIT DATA ------------

# Assuming TRAINING_DATA and LABELS have the same number of rows and aligned order
set.seed(420)

n <- nrow(TRAINING_DATA)
train_indices <- sample(1:n, size = floor(0.7 * n))
test_indices <- setdiff(1:n, train_indices)

# Split the data
TRAIN <- TRAINING_DATA[train_indices, ]
TRAIN_META <- TRAIN %>% select(-(1:which(names(TRAIN) == "IBD_estimate")))
TRAIN <-TRAIN %>% select(1:which(names(TRAIN) == "IBD_estimate"))
TRAIN_labels <- LABELS[train_indices, ]
table(TRAIN_labels)

TEST <- TRAINING_DATA[test_indices, ]
TEST_META <- TEST %>% select(-(1:which(names(TEST) == "IBD_estimate")))
TEST <-TEST %>% select(1:which(names(TEST) == "IBD_estimate"))
TEST_labels <- LABELS[test_indices, ]
table(TEST_labels)



##### 4) TEST MODEL USING IBD ONLY (BEST FEATURE)--------

# Create training and test data frames for the current feature
df_train_IBD <- data.frame(x = TRAIN[["IBD_estimate"]], label = as.factor(TRAIN_labels))
df_test_IBD <- data.frame(x = TEST[["IBD_estimate"]], label = as.factor(TEST_labels))

# Fit a simple logistic regression on the training set
fit_IBD <- glm(label ~ x, data = df_train_IBD, family = binomial)
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


##### 5) TEST MODEL ON REAL DATA USING IBD ONLY (BEST FEATURE) [[NEEDS COI FROM PAIRS TO APPLY DECIDION THRESHOLS!!]]--------

#add threhold data
REAL_DATA <- merge(REAL_DATA, best_decision_thresholds[c("eCOI_pairs", "decision_threshold")], by = c("eCOI_pairs"))

# Initialize a vector to store predictions
REAL_DATA$prediction_prob <- NA
REAL_DATA$predictions <- NA

# Loop over each row in REAL_DATA to apply the corresponding decision_threshold
for (i in 1:nrow(REAL_DATA)) {
  # Subset the current row's IBD_estimate and decision_threshold
  ibd_value <- REAL_DATA$IBD_estimate[i]
  decision_threshold <- REAL_DATA$decision_threshold[i]
  
  # Create a new data frame for prediction (based on the current decision_threshold)
  newdata <- data.frame(x = ibd_value)
  
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

write.csv(REAL_DATA, paste0(site, "_REAL_DATA_top_",top_n_amps,"_amps_PREDICTIONS.csv"), row.names = F)

ggsave(paste0(site, "_REAL_DATA_top_",top_n_amps,"_amps_THRESHOLDS_PLOT.png"), sens_spec_plot, bg = "white", dpi = 300, height = 9, width = 12)
