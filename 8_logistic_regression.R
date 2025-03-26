


library(caret)    
library(dplyr)    
library(tidyr)    
library(ggplot2)  
library(broom)



site <- "Cabo_Delgado"
top_n_amps <- 50 



### 1) IMPORT TRAINING AND REAL DATA ----------

TRAINING_DATA <- read.csv(paste0(site, "_TRAINING_DATA_top_",top_n_amps,"_amps.csv"), row.names = 1) 
LABELS <- data.frame(labels = TRAINING_DATA$labels)
LABELS$labels <- as.factor(LABELS$labels)

TRAINING_DATA <-TRAINING_DATA %>% select(1:which(names(TRAINING_DATA) == "IBD_estimate"))

REAL_DATA <- read.csv(paste0(site, "_REAL_DATA_top_",top_n_amps,"_amps.csv"), row.names = 1) 



### 2) SPLIT DATA ------------

# Assuming TRAINING_DATA and LABELS have the same number of rows and aligned order
set.seed(420)

n <- nrow(TRAINING_DATA)
train_indices <- sample(1:n, size = floor(0.7 * n))
test_indices <- setdiff(1:n, train_indices)

# Split the data
TRAIN <- TRAINING_DATA[train_indices, ]
TRAIN_labels <- LABELS[train_indices, ]
table(TRAIN_labels)

TEST <- TRAINING_DATA[test_indices, ]
TEST_labels <- LABELS[test_indices, ]
table(TEST_labels)



### 3) FULL LOGISTIC REGRESSION MODEL (ALL FEATURES) --------

# Train logistic regression using all features
full_model <- glm(TRAIN_labels ~ ., data = TRAIN, family = binomial)
summary(full_model)

# Predict probabilities on the test set
test_probs <- predict(full_model, newdata = TEST, type = "response")

# Convert probabilities to binary class predictions (threshold = 0.5)
test_preds <- ifelse(test_probs >= 0.5, "R", "NI")

# Evaluate the model using a confusion matrix
cm <- confusionMatrix(as.factor(test_preds), TEST_labels, positive = "R")

print(cm$byClass[c("Sensitivity", "Specificity")])


# Extract model summary into a tidy format
coef_df <- tidy(full_model)

# Plot coefficients with confidence intervals
ggplot(coef_df, aes(x = estimate, y = term)) +
  geom_point(color = "blue", size = 3) +  # Plot coefficients
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error, 
                     xmax = estimate + 1.96 * std.error), 
                 height = 0.2, color = "red") +  # Add confidence intervals
  theme_minimal() +
  labs(title = "Logistic Regression Coefficients",
       x = "Estimate (Effect Size)",
       y = "Feature")


coef_df$odds_ratio <- exp(coef_df$estimate)  # Convert to odds ratios

ggplot(coef_df, aes(x = odds_ratio, y = term)) +
  geom_point(color = "blue", size = 3) + 
  geom_errorbarh(aes(xmin = exp(estimate - 1.96 * std.error), 
                     xmax = exp(estimate + 1.96 * std.error)), 
                 height = 0.2, color = "red") + 
  scale_x_log10() +  # Log scale for better visualization
  theme_minimal() +
  labs(title = "Odds Ratios from Logistic Regression",
       x = "Odds Ratio",
       y = "Feature")




### 3.1) TEST FULL MODEL ON REAL DATA ------

predictions <- predict(full_model, newdata = REAL_DATA, type = "response")
preds <- ifelse(predictions >= 0.5, "R", "NI")
print(as.data.frame(cbind(PairsID = rownames(REAL_DATA), prediction = preds, R_prob = predictions)))




### 4) INDIVIDUAL LOGISTIC REGRESSION MODELS (EACH FEATURE) --------

results <- data.frame(feature = character(),
                      sensitivity = numeric(),
                      specificity = numeric(),
                      stringsAsFactors = FALSE)

# Loop over each feature in the training data
for(feature in colnames(TRAIN)) {
  
  # Create training and test data frames for the current feature
  df_train <- data.frame(x = TRAIN[[feature]], label = as.factor(TRAIN_labels))
  df_test <- data.frame(x = TEST[[feature]], label = as.factor(TEST_labels))
  
  # Fit a simple logistic regression on the training set
  fit <- glm(label ~ x, data = df_train, family = binomial)
  
  # Predict probabilities on the test (holdout) set
  preds_prob <- predict(fit, newdata = df_test, type = "response")
  
  # Convert probabilities to binary predictions using threshold 0.5
  preds <- ifelse(preds_prob >= 0.5, "R", "NI")
  
  # Calculate the confusion matrix (assuming positive class is "1")
  cm <- confusionMatrix(as.factor(preds), df_test$label, positive = "R")
  
  # Extract sensitivity and specificity
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  # Append results
  results <- rbind(results, data.frame(feature = feature,
                                       sensitivity = sens,
                                       specificity = spec,
                                       stringsAsFactors = FALSE))
}

results

# Reshape the results for plotting
results_long <- results %>%
  pivot_longer(cols = c("sensitivity", "specificity"),
               names_to = "metric",
               values_to = "value")

# Create an informative bar plot
ggplot(results_long, aes(x = feature, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "",
       x = "Feature",
       y = "Metric Value") +
  theme_minimal()



##### 4) TEST MODEL ON REAL DATA USING IBD ONLY (BEST FEATURE)--------

# Create training and test data frames for the current feature
df_train <- data.frame(x = TRAIN[["IBD_estimate"]], label = as.factor(TRAIN_labels))
df_test <- data.frame(x = TEST[["IBD_estimate"]], label = as.factor(TEST_labels))

# Fit a simple logistic regression on the training set
fit_IBD <- glm(label ~ x, data = df_train, family = binomial)
summary(fit_IBD)

#subset ibd from real data
newdata <- data.frame(x = REAL_DATA$IBD_estimate)

#predict
predictions <- predict(fit_IBD, newdata = newdata, type = "response")
preds <- ifelse(predictions >= 0.5, "R", "NI")

print(as.data.frame(cbind(PairsID = rownames(REAL_DATA), prediction = preds, R_prob = predictions)))
