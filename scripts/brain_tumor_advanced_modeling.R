# Load required libraries
library(tidyverse)
library(caret)
library(randomForest)
library(ROSE)
library(pROC)
library(gridExtra)
library(corrplot)
library(car)  # For VIF in linear regression

# Create output directories
dir.create("output/advanced_modeling", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/advanced_modeling", recursive = TRUE, showWarnings = FALSE)


# Redirect output for Random Forest
sink("output/advanced_modeling/advanced_modeling_results.txt")

# Load and preprocess data
data <- read.csv("data/processed/brain_tumor_data_clean.csv") %>%
  dplyr::select(-Patient_ID) %>%
  filter(complete.cases(.), Age >= 0, Age <= 100, Survival_Rate >= 0, Survival_Rate <= 100, Tumor_Size > 0)

# Convert categorical variables to factors
data[, c("Gender", "Tumor_Type", "Location", "Histology", "Stage", "Radiation_Treatment",
         "Surgery_Performed", "Chemotherapy", "Family_History", "MRI_Result",
         "Follow_Up_Required", "Survival_Binary")] <- lapply(data[, c("Gender", "Tumor_Type", "Location",
                                                                     "Histology", "Stage", "Radiation_Treatment",
                                                                     "Surgery_Performed", "Chemotherapy",
                                                                     "Family_History", "MRI_Result",
                                                                     "Follow_Up_Required", "Survival_Binary")], as.factor)

# Handle Survival_Binary with single class
if (length(unique(data$Survival_Binary)) < 2) {
  cat("WARNING: Survival_Binary has one level. Creating Survival_Binary_New.\n")
  threshold <- median(data$Survival_Rate, na.rm = TRUE)
  data$Survival_Binary_New <- factor(ifelse(data$Survival_Rate >= threshold, "1", "0"), levels = c("0", "1"))
  cat("Threshold:", threshold, "\n")
  print(table(data$Survival_Binary_New))
} else {
  data$Survival_Binary_New <- data$Survival_Binary
}

# Create copies of the data for each analysis
data_rf <- data  # For Random Forest (normalized)
data_lr <- data  # For Linear Regression (will be unnormalized)

# Normalize numeric columns for Random Forest
numeric_cols <- c("Age", "Tumor_Size", "Tumor_Growth_Rate")
data_rf[numeric_cols] <- lapply(data_rf[numeric_cols], function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))

# Unnormalize numeric columns for Linear Regression
data_lr$Age <- data_lr$Age * (79 - 20) + 20  # Original range: 20 to 79
data_lr$Tumor_Size <- data_lr$Tumor_Size * (9.3591 - 0.3394) + 0.3394  # Original range: 0.3394 to 9.3591
data_lr$Tumor_Growth_Rate <- data_lr$Tumor_Growth_Rate * (5.0936 - (-0.1596)) + (-0.1596)  # Original range: -0.1596 to 5.0936

# --- Random Forest Analysis ---
# Machine Learning for Random Forest
set.seed(123)
features <- c("Age", "Gender", "Tumor_Type", "Tumor_Size", "Location", "Histology", "Stage",
              "Tumor_Growth_Rate", "Radiation_Treatment", "Surgery_Performed", "Chemotherapy", "Family_History")
data_model_rf <- data_rf[, c(features, "Survival_Binary_New")]

# Remove outliers for Random Forest
for (col in numeric_cols) {
  q <- quantile(data_model_rf[[col]], c(0.05, 0.95), na.rm = TRUE)
  data_model_rf <- data_model_rf[data_model_rf[[col]] >= q[1] & data_model_rf[[col]] <= q[2], ]
}
data_model_rf <- na.omit(data_model_rf)

# Split data into training and testing sets
train_idx <- createDataPartition(data_model_rf$Survival_Binary_New, p = 0.7, list = FALSE)
train_data <- data_model_rf[train_idx, ]
test_data <- data_model_rf[-train_idx, ]

# Check training class distribution and balance if necessary
cat("Training class distribution:\n")
print(table(train_data$Survival_Binary_New))
if (length(unique(train_data$Survival_Binary_New)) >= 2) {
  train_data <- ovun.sample(Survival_Binary_New ~ ., data = train_data, method = "both", p = 0.5, seed = 123)$data
} else {
  cat("ERROR: Cannot balance classes; only one class in training data.\n")
  sink()
  stop("Single class in training data")
}

# Train Random Forest model
rf <- randomForest(Survival_Binary_New ~ ., data = train_data, ntree = 50)

# Feature Importance
imp <- importance(rf)
imp_df <- data.frame(Feature = rownames(imp), Importance = imp[, "MeanDecreaseGini"])[1:8, ]
imp_df <- imp_df[order(-imp_df$Importance), ]

# Output top 8 features to file
sink()
sink("output/advanced_modeling/feature_importance.txt", append = TRUE)
cat("Top 8 Features:\n")
print(imp_df$Feature)
sink()
sink("output/advanced_modeling/advanced_modeling_results.txt", append = TRUE)

# Plot feature importance
png("plots/advanced_modeling/feature_importance.png", width = 800, height = 600)
ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance", x = "Feature") +
  theme_minimal()
dev.off()

# Close all graphics devices
graphics.off()

# Close sink for Random Forest
sink()

# --- Linear Regression Analysis ---
# Redirect output for Linear Regression
sink("output/advanced_modeling/survival_prediction_results.txt")

# Remove outliers for Linear Regression (in original scale)
# Load required library for VIF
library(car)

# Remove outliers for Linear Regression
data_model_lr <- data_lr
numeric_cols <- c("Age", "Tumor_Size", "Tumor_Growth_Rate")
for (col in numeric_cols) {
  q <- quantile(data_model_lr[[col]], c(0.05, 0.95), na.rm = TRUE)
  data_model_lr <- data_model_lr[data_model_lr[[col]] >= q[1] & data_model_lr[[col]] <= q[2], ]
}
data_model_lr <- na.omit(data_model_lr)

# Fit a simplified linear regression model (already dropped Tumor_Growth_Rate, now dropping Age)
model <- lm(Survival_Rate ~ Tumor_Size + Tumor_Type, data = data_model_lr)

# Check for multicollinearity
cat("\nVariance Inflation Factors (VIF):\n")
print(vif(model))

# Create prediction grid without Age
size_range <- seq(min(data_model_lr$Tumor_Size), max(data_model_lr$Tumor_Size), length.out = 20)  # 1.2 to 8.2 (after outlier removal)
prediction_grid <- expand.grid(
  Tumor_Size = size_range,
  Tumor_Type = factor("Malignant", levels = levels(data_model_lr$Tumor_Type))
)

# Predict survival rates
prediction_grid$Predicted_Survival <- predict(model, newdata = prediction_grid)

# Ensure predictions are within the valid range [67.10, 100]
prediction_grid$Predicted_Survival <- pmin(pmax(prediction_grid$Predicted_Survival, 67.10), 100)

# Create a line plot (since Age is dropped, we only vary Tumor_Size)
p7 <- ggplot(prediction_grid, aes(x = Tumor_Size, y = Predicted_Survival)) +
  geom_line(color = "steelblue", size = 1.5) +
  geom_point(color = "steelblue", size = 3) +
  labs(
    title = "Predicted Survival Rate by Tumor Size",
    subtitle = "Based on linear regression model (Tumor Type: Malignant)",
    x = "Tumor Size (cm)",
    y = "Predicted Survival Rate (%)"
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#F8F8F8", color = NA),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  )

# Save the plot
ggsave("plots/advanced_modeling/survival_prediction_lineplot.png", p7, width = 10, height = 8, bg = "white")

# Print model summary
cat("\nSurvival Prediction Model Summary:\n")
print(summary(model))

# Close sink for Linear Regression
sink()