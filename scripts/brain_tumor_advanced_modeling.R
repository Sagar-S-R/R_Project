# Load required libraries
library(tidyverse)
library(caret)
library(randomForest)
library(ROSE)
library(glmnet)
library(pROC)
library(moments)
library(gridExtra)
library(corrplot)
library(psych)

# Create output directories
dir.create("output/advanced_modeling", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/advanced_modeling", recursive = TRUE, showWarnings = FALSE)

# Redirect output
sink("output/advanced_modeling/advanced_modeling_results.txt")

# Load and preprocess data
data <- read.csv("data/processed/brain_tumor_data_clean.csv") %>%
  dplyr::select(-Patient_ID) %>%
  filter(complete.cases(.), Age >= 0, Age <= 100, Survival_Rate >= 0, Survival_Rate <= 100, Tumor_Size > 0)
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

# Normalize numeric columns
numeric_cols <- c("Age", "Tumor_Size", "Tumor_Growth_Rate")
data[numeric_cols] <- lapply(data[numeric_cols], function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))

# Univariate Analysis
describe_numeric <- function(x) {
  data.frame(N = length(na.omit(x)), Mean = mean(x, na.rm = TRUE), Median = median(x, na.rm = TRUE),
             SD = sd(x, na.rm = TRUE), Min = min(x, na.rm = TRUE), Max = max(x, na.rm = TRUE),
             Range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE),
             Skewness = skewness(x, na.rm = TRUE), Kurtosis = kurtosis(x, na.rm = TRUE))
}
describe_categorical <- function(x) {
  freq <- table(x)
  data.frame(Category = names(freq), Frequency = as.numeric(freq),
             Percentage = round(as.numeric(prop.table(freq)) * 100, 2))
}
numeric_stats <- do.call(rbind, lapply(data[numeric_cols], describe_numeric))
rownames(numeric_stats) <- numeric_cols
sink(); sink("output/advanced_modeling/numeric_descriptive_stats.txt", append = TRUE)
cat("Numeric Descriptive Statistics:\n")
print(numeric_stats)
sink(); sink("output/advanced_modeling/advanced_modeling_results.txt", append = TRUE)
categorical_cols <- c("Gender", "Tumor_Type", "Location", "Histology", "Stage", "Radiation_Treatment",
                      "Surgery_Performed", "Chemotherapy", "Family_History", "MRI_Result", "Follow_Up_Required", "Survival_Binary")
cat_stats <- lapply(data[categorical_cols], describe_categorical)
sink(); sink("output/advanced_modeling/categorical_descriptive_stats.txt", append = TRUE)
cat("Categorical Descriptive Statistics:\n")
for (i in seq_along(cat_stats)) {
  cat(names(cat_stats)[i], ":\n")
  print(cat_stats[[i]])
}
sink(); sink("output/advanced_modeling/advanced_modeling_results.txt", append = TRUE)

# Numeric visualizations
for (col in numeric_cols) {
  p1 <- ggplot(data, aes(x = !!sym(col))) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
    geom_density(alpha = 0.5, fill = "orange") + ggtitle(paste("Histogram of", col)) + theme_minimal()
  p2 <- ggplot(data, aes(y = !!sym(col))) +
    geom_boxplot(fill = "lightgreen", color = "black") + ggtitle(paste("Boxplot of", col)) + theme_minimal()
  p3 <- ggplot(data, aes(sample = !!sym(col))) +
    stat_qq() + stat_qq_line() + ggtitle(paste("Q-Q Plot of", col)) + theme_minimal()
  combined_plot <- grid.arrange(p1, p2, p3, ncol = 3)
  ggsave(paste0("plots/advanced_modeling/", col, "_distribution.png"), combined_plot, width = 15, height = 5)
}

# Multivariate Analysis
cor_matrix <- cor(data[, numeric_cols], use = "complete.obs")
sink(); sink("output/advanced_modeling/correlation_matrix.txt", append = TRUE)
cat("Correlation Matrix:\n")
print(cor_matrix)
sink(); sink("output/advanced_modeling/advanced_modeling_results.txt", append = TRUE)
png("plots/advanced_modeling/correlation_heatmap.png", width = 800, height = 600)
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8, addCoef.col = "black")
dev.off()
png("plots/advanced_modeling/scatterplot_matrix.png", width = 1000, height = 1000)
pairs(data[, numeric_cols], main = "Scatterplot Matrix")
dev.off()
sink(); sink("output/advanced_modeling/stats_by_gender.txt", append = TRUE)
cat("Stats by Gender:\n")
print(describeBy(data[numeric_cols], group = data$Gender))
sink(); sink("output/advanced_modeling/advanced_modeling_results.txt", append = TRUE)
sink(); sink("output/advanced_modeling/stats_by_survival.txt", append = TRUE)
cat("Stats by Survival:\n")
print(describeBy(data[numeric_cols], group = data$Survival_Binary))
sink(); sink("output/advanced_modeling/advanced_modeling_results.txt", append = TRUE)

# Machine Learning
set.seed(123)
features <- c("Age", "Gender", "Tumor_Type", "Tumor_Size", "Location", "Histology", "Stage",
              "Tumor_Growth_Rate", "Radiation_Treatment", "Surgery_Performed", "Chemotherapy", "Family_History")
data_model <- data[, c(features, "Survival_Binary_New")]
for (col in numeric_cols) {
  q <- quantile(data_model[[col]], c(0.05, 0.95), na.rm = TRUE)
  data_model <- data_model[data_model[[col]] >= q[1] & data_model[[col]] <= q[2], ]
}
data_model <- na.omit(data_model)
train_idx <- createDataPartition(data_model$Survival_Binary_New, p = 0.7, list = FALSE)
train_data <- data_model[train_idx, ]
test_data <- data_model[-train_idx, ]
cat("Training class distribution:\n")
print(table(train_data$Survival_Binary_New))
if (length(unique(train_data$Survival_Binary_New)) >= 2) {
  train_data <- ovun.sample(Survival_Binary_New ~ ., data = train_data, method = "both", p = 0.5, seed = 123)$data
} else {
  cat("ERROR: Cannot balance classes; only one class in training data.\n")
  sink()
  stop("Single class in training data")
}
rf <- randomForest(Survival_Binary_New ~ ., data = train_data, ntree = 50)
imp <- importance(rf)
imp_df <- data.frame(Feature = rownames(imp), Importance = imp[, "MeanDecreaseGini"])[1:8, ]
imp_df <- imp_df[order(-imp_df$Importance), ]
sink(); sink("output/advanced_modeling/feature_importance.txt", append = TRUE)
cat("Top 8 Features:\n")
print(imp_df$Feature)
sink(); sink("output/advanced_modeling/advanced_modeling_results.txt", append = TRUE)
png("plots/advanced_modeling/feature_importance.png", width = 800, height = 600)
ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") + coord_flip() +
  labs(title = "Feature Importance", x = "Feature") + theme_minimal()
dev.off()
train_subset <- train_data[, c(imp_df$Feature, "Survival_Binary_New")]
test_subset <- test_data[, c(imp_df$Feature, "Survival_Binary_New")]
x_train <- model.matrix(Survival_Binary_New ~ ., train_subset)[, -1]
x_test <- model.matrix(Survival_Binary_New ~ ., test_subset)[, -1]
y_train <- train_subset$Survival_Binary_New
log_model <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, nfolds = 5)
log_pred_prob <- as.numeric(predict(log_model, x_test, type = "response", s = "lambda.min"))
roc_obj <- roc(test_subset$Survival_Binary_New, log_pred_prob)
threshold <- coords(roc_obj, "best", ret = "threshold")$threshold
log_pred <- factor(ifelse(log_pred_prob > threshold, "1", "0"), levels = c("0", "1"))
cm <- confusionMatrix(log_pred, test_subset$Survival_Binary_New)
metrics <- data.frame(Metric = c("Accuracy", "F1-Score", "AUC-ROC"),
                      Value = c(mean(log_pred == test_subset$Survival_Binary_New), cm$byClass["F1"], auc(roc_obj)))
sink(); sink("output/advanced_modeling/logistic_regression_metrics.txt", append = TRUE)
print(metrics)
sink(); sink("output/advanced_modeling/advanced_modeling_results.txt", append = TRUE)
cm_table <- cm$table
cm_percent <- round((cm_table / nrow(test_subset)) * 100, 2)
cm_combined <- matrix(paste0(cm_table, " (", cm_percent, "%)"), nrow = 2)
dimnames(cm_combined) <- dimnames(cm_table)
sink(); sink("output/advanced_modeling/logistic_regression_confusion_matrix.txt", append = TRUE)
print(cm_combined)
sink(); sink("output/advanced_modeling/advanced_modeling_results.txt", append = TRUE)
png("plots/advanced_modeling/logistic_regression_confusion_matrix.png", width = 800, height = 600)
ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() + geom_text(aes(label = Freq), color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Logistic Regression Confusion Matrix", x = "Actual", y = "Predicted") + theme_minimal()
dev.off()

# Close all graphics devices
graphics.off()

# Close sink
sink()