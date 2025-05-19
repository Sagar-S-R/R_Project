# Brain_Tumor_Analysis_Project/scripts/brain_tumor_survival_modeling.R

# Load required libraries
library(tidyverse)
library(caret)
library(rpart)
library(rpart.plot)
library(pROC)
library(glmnet)  # Added for regularized regression

# Load cleaned data
brain_tumor_data <- read.csv("data/processed/brain_tumor_data_clean.csv", stringsAsFactors = TRUE)

# Inspect the structure of Survival_Binary
print("Levels in Survival_Binary:")
print(levels(brain_tumor_data$Survival_Binary))

# Create treatment combination column
brain_tumor_data$treatment_combo <- paste(
  brain_tumor_data$Radiation_Treatment,
  brain_tumor_data$Surgery_Performed,
  brain_tumor_data$Chemotherapy,
  sep = "_"
)

# Check if we need to create a binary outcome (since current one only has one level)
if(length(levels(brain_tumor_data$Survival_Binary)) < 2) {
  print("WARNING: Survival_Binary has only one level. Creating a new binary outcome based on Survival_Rate")
  
  # Create a new binary outcome based on Survival_Rate
  # Use median as threshold to split the data roughly equally
  survival_threshold <- median(brain_tumor_data$Survival_Rate)
  brain_tumor_data$Survival_Binary_New <- factor(
    ifelse(brain_tumor_data$Survival_Rate >= survival_threshold, "High", "Low"),
    levels = c("Low", "High")
  )
  
  print(paste("Created new binary outcome using threshold:", survival_threshold))
  print(table(brain_tumor_data$Survival_Binary_New))
} else {
  # Use existing binary outcome
  brain_tumor_data$Survival_Binary_New <- brain_tumor_data$Survival_Binary
}

# Split data into training and test sets
set.seed(123)
train_indices <- createDataPartition(brain_tumor_data$Survival_Binary_New, p = 0.8, list = FALSE)
train_data <- brain_tumor_data[train_indices, ]
test_data <- brain_tumor_data[-train_indices, ]

# ---- ANALYSIS 5: Survival Rate Analysis ----
survival_dist_plot <- ggplot(brain_tumor_data, aes(x = Survival_Rate)) +
  geom_histogram(binwidth = 5, fill = "#4287f5", color = "black", alpha = 0.8) +
  labs(title = "Distribution of Survival Rates for Brain Tumor Patients", x = "Survival Rate (%)", y = "Count") +
  theme_classic() + 
  theme(panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white"))
ggsave("plots/survival_modeling/survival_rate_distribution.png", survival_dist_plot, width = 8, height = 6)

survival_by_type <- ggplot(brain_tumor_data, aes(x = Tumor_Type, y = Survival_Rate, fill = Tumor_Type)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Survival Rate by Brain Tumor Type", x = "Tumor Type", y = "Survival Rate (%)") +
  theme_classic() + 
  theme(panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")
ggsave("plots/survival_modeling/survival_by_tumor_type.png", survival_by_type, width = 10, height = 6)

survival_by_treatment <- ggplot(brain_tumor_data, aes(x = treatment_combo, y = Survival_Rate, fill = treatment_combo)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(title = "Survival Rate by Treatment Combination", x = "Treatment Combination", y = "Survival Rate (%)") +
  theme_classic() + 
  theme(panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")
ggsave("plots/survival_modeling/survival_by_treatment.png", survival_by_treatment, width = 10, height = 6)

# ---- ANALYSIS 6: Linear Regression - Survival Rate Prediction ----
survival_formula <- as.formula(paste("Survival_Rate ~", 
                                     paste(c("Age", "Tumor_Size"),
                                           collapse = " + ")))

lm_model <- lm(survival_formula, data = train_data)
summary_lm <- summary(lm_model)

sink("output/survival_modeling/linear_model_summary.txt")
print(summary_lm)
sink()

test_data$predicted_survival <- predict(lm_model, newdata = test_data)

lm_rmse <- sqrt(mean((test_data$Survival_Rate - test_data$predicted_survival)^2))
lm_mae <- mean(abs(test_data$Survival_Rate - test_data$predicted_survival))
lm_r2 <- cor(test_data$Survival_Rate, test_data$predicted_survival)^2

lm_metrics <- data.frame(
  Metric = c("RMSE", "MAE", "R-squared"),
  Value = c(lm_rmse, lm_mae, lm_r2)
)
write.csv(lm_metrics, "output/survival_modeling/linear_model_metrics.csv", row.names = FALSE)

lm_plot <- ggplot(test_data, aes(x = Survival_Rate, y = predicted_survival)) +
  geom_point(alpha = 0.6, color = "#E94B3C", size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "#2D7D90", linetype = "dashed", linewidth = 1.2) +
  labs(title = "Linear Model: Actual vs Predicted Survival Rate",
       x = "Actual Survival Rate",
       y = "Predicted Survival Rate") +
  theme_classic() + 
  theme(panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white"))
ggsave("plots/survival_modeling/linear_model_predictions.png", lm_plot, width = 8, height = 6)

# ---- ANALYSIS 7: Logistic Regression - Binary Survival Prediction ----
# Now using new binary outcome
print("Checking Survival_Binary_New levels before modeling:")
print(table(train_data$Survival_Binary_New))

# Make sure we have both levels in training and test data
if(length(levels(train_data$Survival_Binary_New)) < 2 || 
   length(unique(train_data$Survival_Binary_New)) < 2) {
  print("ERROR: Training data still only has one level for binary outcome. Cannot fit logistic regression.")
} else {
  # Simplified logistic regression model with fewer predictors
  logistic_formula <- as.formula(paste("Survival_Binary_New ~", 
                                       paste(c("Age", "Tumor_Size"),
                                             collapse = " + ")))
  
  # Fit logistic regression
  logistic_model <- glm(logistic_formula, 
                       data = train_data, 
                       family = "binomial")
  
  summary_logistic <- summary(logistic_model)
  
  sink("output/survival_modeling/logistic_model_summary.txt")
  print(summary_logistic)
  sink()
  
  test_data$predicted_prob <- predict(logistic_model, newdata = test_data, type = "response")
  test_data$predicted_class <- factor(ifelse(test_data$predicted_prob >= 0.5, 
                                            "High", "Low"),
                                     levels = c("Low", "High"))
  
  # Now create confusion matrix
  conf_matrix <- confusionMatrix(test_data$predicted_class, test_data$Survival_Binary_New)
  
  sink("output/survival_modeling/logistic_model_confusion_matrix.txt")
  print(conf_matrix)
  sink()
  
  # ROC curve
  roc_obj <- roc(as.numeric(test_data$Survival_Binary_New) - 1, test_data$predicted_prob)
  auc_value <- auc(roc_obj)
  
  # Add a colorful ROC curve
  png("plots/survival_modeling/logistic_roc_curve.png", width = 800, height = 800)
  par(bg = "white")
  plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc_value, 3), ")"),
       col = "#4287f5", lwd = 3)
  abline(a = 0, b = 1, lty = 2, col = "darkgray")
  dev.off()
  
  # Additional plot: Probability Distribution by Actual Class
  prob_dist_plot <- ggplot(test_data, aes(x = predicted_prob, fill = Survival_Binary_New)) +
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = c("#E74C3C", "#3498DB"), name = "Survival Class") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +
    labs(title = "Distribution of Predicted Probabilities by Actual Class",
         x = "Predicted Probability", y = "Density") +
    theme_classic() +
    theme(panel.background = element_rect(fill = "white", color = "black"),
          plot.background = element_rect(fill = "white"))
  ggsave("plots/survival_modeling/probability_distribution.png", prob_dist_plot, width = 8, height = 6)
}

# ---- ANALYSIS 8: Decision Tree for Survival Classification ----
if(length(unique(train_data$Survival_Binary_New)) < 2) {
  print("ERROR: Cannot fit decision tree with only one class level")
} else {
  # Simplified decision tree with fewer predictors to avoid overfitting
  tree_formula <- as.formula(paste("Survival_Binary_New ~", 
                                   paste(c("Age", "Tumor_Size", "Tumor_Growth_Rate", 
                                           "Gender", "Tumor_Type"),
                                         collapse = " + ")))
  
  dt_model <- rpart(tree_formula, data = train_data, method = "class",
                    cp = 0.01, minsplit = 20)
  
  # Create a more colorful decision tree plot
  png("plots/survival_modeling/decision_tree.png", width = 1200, height = 800, bg = "white")
  rpart.plot(dt_model, extra = 104, fallen.leaves = TRUE, 
             main = "Decision Tree for Brain Tumor Survival Classification",
             box.palette = "Blues", shadow.col = "gray", 
             branch.lty = 3, branch.lwd = 2)
  dev.off()
  
  dt_predictions <- predict(dt_model, newdata = test_data, type = "class")
  
  # Create confusion matrix
  dt_conf_matrix <- confusionMatrix(dt_predictions, test_data$Survival_Binary_New)
  
  sink("output/survival_modeling/decision_tree_evaluation.txt")
  print(dt_conf_matrix)
  sink()
  
  # Add feature importance plot
  importance <- dt_model$variable.importance
  if(length(importance) > 0) {
    imp_df <- data.frame(
      Feature = names(importance),
      Importance = as.numeric(importance)
    )
    imp_df <- imp_df[order(-imp_df$Importance),]
    
    var_imp_plot <- ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance, fill = Importance)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "#9BD770", high = "#3A5F0B") +
      coord_flip() +
      labs(title = "Variable Importance in Decision Tree Model",
           x = "Feature",
           y = "Importance Score") +
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            plot.background = element_rect(fill = "white"))
    
    ggsave("plots/survival_modeling/variable_importance.png", var_imp_plot, width = 8, height = 6)
  }
}

# ---- ANALYSIS 9: Additional Plot - Survival Rate by Age Group ----
# Create age groups
brain_tumor_data$Age_Group <- cut(brain_tumor_data$Age, 
                                 breaks = c(0, 20, 40, 60, 80, 100),
                                 labels = c("0-20", "21-40", "41-60", "61-80", "80+"))

age_survival_plot <- ggplot(brain_tumor_data, aes(x = Age_Group, y = Survival_Rate, fill = Age_Group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_brewer(palette = "Spectral") +
  labs(title = "Survival Rate by Age Group", x = "Age Group", y = "Survival Rate (%)") +
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white")) +
  guides(fill = "none")
ggsave("plots/survival_modeling/survival_by_age.png", age_survival_plot, width = 8, height = 6)

# ---- ANALYSIS 10: Correlation Plot for Numerical Variables ----
numerical_vars <- c("Age", "Tumor_Size", "Tumor_Growth_Rate", "Survival_Rate")
cor_data <- brain_tumor_data[, numerical_vars]
cor_matrix <- cor(cor_data, use = "complete.obs")

# Create correlation plot with viridis color scale
library(corrplot)
png("plots/survival_modeling/correlation_plot.png", width = 800, height = 800, bg = "white")
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust",
         col = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
         addCoef.col = "black", tl.col = "black", tl.srt = 45,
         diag = FALSE, title = "Correlation Plot of Numerical Variables",
         mar = c(0, 0, 2, 0))
dev.off()