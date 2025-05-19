# Brain_Tumor_Analysis_Project/scripts/brain_tumor_survival_modeling.R

# Load required libraries
library(tidyverse)
library(caret)
library(rpart)
library(rpart.plot)
library(pROC)

# Load cleaned data
brain_tumor_data <- read.csv("data/processed/brain_tumor_data_clean.csv", stringsAsFactors = TRUE)

# Create treatment combination column
brain_tumor_data$treatment_combo <- paste(
  brain_tumor_data$Radiation_Treatment,
  brain_tumor_data$Surgery_Performed,
  brain_tumor_data$Chemotherapy,
  sep = "_"
)

# Split data into training and test sets
set.seed(123)
train_indices <- createDataPartition(brain_tumor_data$Survival_Binary, p = 0.8, list = FALSE)
train_data <- brain_tumor_data[train_indices, ]
test_data <- brain_tumor_data[-train_indices, ]

# ---- ANALYSIS 5: Survival Rate Analysis ----
survival_dist_plot <- ggplot(brain_tumor_data, aes(x = Survival_Rate)) +
  geom_histogram(binwidth = 5, fill = "purple", color = "black") +
  labs(title = "Distribution of Survival Rates for Brain Tumor Patients", x = "Survival Rate (%)", y = "Count") +
  theme_minimal()
ggsave("plots/survival_modeling/survival_rate_distribution.png", survival_dist_plot, width = 8, height = 6)

survival_by_type <- ggplot(brain_tumor_data, aes(x = Tumor_Type, y = Survival_Rate, fill = Tumor_Type)) +
  geom_boxplot() +
  labs(title = "Survival Rate by Brain Tumor Type", x = "Tumor Type", y = "Survival Rate (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")
ggsave("plots/survival_modeling/survival_by_tumor_type.png", survival_by_type, width = 10, height = 6)

survival_by_treatment <- ggplot(brain_tumor_data, aes(x = treatment_combo, y = Survival_Rate, fill = treatment_combo)) +
  geom_boxplot() +
  labs(title = "Survival Rate by Treatment Combination", x = "Treatment Combination", y = "Survival Rate (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")
ggsave("plots/survival_modeling/survival_by_treatment.png", survival_by_treatment, width = 10, height = 6)

# ---- ANALYSIS 6: Linear Regression - Survival Rate Prediction ----
survival_formula <- as.formula(paste("Survival_Rate ~", 
                                     paste(c("Age", "Tumor_Size", "Tumor_Growth_Rate", 
                                             "Gender", "Tumor_Type", "Stage",
                                             "Radiation_Treatment", "Surgery_Performed", 
                                             "Chemotherapy", "Family_History"),
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
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Linear Model: Actual vs Predicted Survival Rate",
       x = "Actual Survival Rate",
       y = "Predicted Survival Rate") +
  theme_minimal()
ggsave("plots/survival_modeling/linear_model_predictions.png", lm_plot, width = 8, height = 6)

# ---- ANALYSIS 7: Logistic Regression - Binary Survival Prediction ----
logistic_formula <- as.formula(paste("Survival_Binary ~", 
                                     paste(c("Age", "Tumor_Size", "Tumor_Growth_Rate", 
                                             "Gender", "Tumor_Type", "Stage",
                                             "Radiation_Treatment", "Surgery_Performed", 
                                             "Chemotherapy", "Family_History"),
                                           collapse = " + ")))

logistic_model <- glm(logistic_formula, data = train_data, family = "binomial")
summary_logistic <- summary(logistic_model)

sink("output/survival_modeling/logistic_model_summary.txt")
print(summary_logistic)
sink()

test_data$predicted_prob <- predict(logistic_model, newdata = test_data, type = "response")
test_data$predicted_class <- ifelse(test_data$predicted_prob >= 0.5, "Good", "Poor")

conf_matrix <- confusionMatrix(factor(test_data$predicted_class, levels = c("Poor", "Good")), 
                               test_data$Survival_Binary)

sink("output/survival_modeling/logistic_model_confusion_matrix.txt")
print(conf_matrix)
sink()

roc_obj <- roc(test_data$Survival_Binary, test_data$predicted_prob)
auc_value <- auc(roc_obj)

png("plots/survival_modeling/logistic_roc_curve.png", width = 800, height = 800)
plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc_value, 3), ")"))
dev.off()

# ---- ANALYSIS 8: Decision Tree for Survival Classification ----
tree_formula <- as.formula(paste("Survival_Binary ~", 
                                 paste(c("Age", "Tumor_Size", "Tumor_Growth_Rate", 
                                         "Gender", "Tumor_Type", "Stage",
                                         "Radiation_Treatment", "Surgery_Performed", 
                                         "Chemotherapy", "Family_History"),
                                       collapse = " + ")))

dt_model <- rpart(tree_formula, data = train_data, method = "class",
                  cp = 0.01, minsplit = 20)

png("plots/survival_modeling/decision_tree.png", width = 1200, height = 800)
rpart.plot(dt_model, extra = 104, fallen.leaves = TRUE, 
           main = "Decision Tree for Brain Tumor Survival Classification")
dev.off()

dt_predictions <- predict(dt_model, newdata = test_data, type = "class")

dt_conf_matrix <- confusionMatrix(dt_predictions, test_data$Survival_Binary)

sink("output/survival_modeling/decision_tree_evaluation.txt")
print(dt_conf_matrix)
sink()
