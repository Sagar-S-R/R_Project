# Brain Tumor Analysis Script
# This script performs comprehensive analysis on brain tumor data

# Create necessary directories if they don't exist
dir.create("output/advanced_modeling", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/advanced_modeling", recursive = TRUE, showWarnings = FALSE)

# Install and Load the packages
# --------- Set CRAN Mirror ---------
options(repos = c(CRAN = "https://cloud.r-project.org"))

# --------- Install Packages (only if not already installed) ---------
packages <- c(
  "tidyverse", "psych", "moments", "gridExtra", "corrplot", "ggpubr", "skimr",
  "caret", "randomForest", "class", "ROSE", "glmnet", "pROC",
  "rcompanion", "httr", "lubridate", "plotly", "rio", "MASS",
  "shiny", "ggcorrplot"
)

installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

# --------- Load Libraries ---------
library(tidyverse)
library(psych)
library(moments)
library(gridExtra)
library(corrplot)
library(ggpubr)
library(skimr)
library(caret)
library(randomForest)
library(class)
library(ROSE)
library(glmnet)
library(pROC)
library(rcompanion)
library(httr)
library(lubridate)
library(plotly)
library(rio)
library(MASS)
library(shiny)
library(ggcorrplot)


# Load the dataset
data <- read.csv("data/processed/brain_tumor_data_clean.csv")

#*******************************************************************************
## Data Preprocessing

# Step 1: Read the data (using brain tumor dataset)
data <- read.csv("data/processed/brain_tumor_data_clean.csv")

# Step 2: Basic preprocessing
# Remove rows with NA values
data <- data %>% filter(complete.cases(.))

# Step 3: Convert categorical variables to factors
data$Gender <- as.factor(data$Gender)
data$Tumor_Type <- as.factor(data$Tumor_Type)
data$Location <- as.factor(data$Location)
data$Histology <- as.factor(data$Histology)
data$Stage <- as.factor(data$Stage)
data$Symptom_1 <- as.factor(data$Symptom_1)
data$Symptom_2 <- as.factor(data$Symptom_2)
data$Symptom_3 <- as.factor(data$Symptom_3)
data$Radiation_Treatment <- as.factor(data$Radiation_Treatment)
data$Surgery_Performed <- as.factor(data$Surgery_Performed)
data$Chemotherapy <- as.factor(data$Chemotherapy)
data$Family_History <- as.factor(data$Family_History)
data$MRI_Result <- as.factor(data$MRI_Result)
data$Follow_Up_Required <- as.factor(data$Follow_Up_Required)
data$Survival_Binary <- as.factor(data$Survival_Binary)

# Step 4: Remove unnecessary or redundant columns
# Remove Patient_ID as it's not needed for analysis
data <- data %>% dplyr::select(-Patient_ID)

# Step 5: Remove extreme outliers for numeric columns
# Define function to remove outliers based on IQR
remove_outliers <- function(x) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- IQR(x, na.rm = TRUE)
  lower_bound <- qnt[1] - 1.5 * iqr
  upper_bound <- qnt[2] + 1.5 * iqr
  x >= lower_bound & x <= upper_bound
}

# Apply outlier removal to numeric columns
numeric_cols <- c("Age", "Tumor_Size", "Survival_Rate", "Tumor_Growth_Rate")
for (col in numeric_cols) {
  data <- data %>% filter(remove_outliers(!!sym(col)))
}

# Step 6: Additional plausibility checks
# Ensure age is within a reasonable range (e.g., 0 to 100)
data <- data %>% filter(Age >= 0 & Age <= 100)

# Ensure Survival_Rate is between 0 and 100
data <- data %>% filter(Survival_Rate >= 0 & Survival_Rate <= 100)

# Ensure Tumor_Size is within a reasonable range 
data <- data %>% filter(Tumor_Size > 0 & Tumor_Size <= 200)

# Step 7: Summary of preprocessed data
preprocessing_summary <- summary(data)
# Save summary to file
sink("output/advanced_modeling/preprocessing_summary.txt")
print("Preprocessing Summary:")
print(preprocessing_summary)
sink()

# Step 8: Normalizing data
# Define the min-max normalization function
min_max_normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Create a copy of the dataset
normalized_data <- data

# Apply normalization to numeric columns
for (col in numeric_cols) {
  normalized_data[[col]] <- min_max_normalize(normalized_data[[col]])
}

# Verify the normalization
normalized_summary <- summary(normalized_data[numeric_cols])
# Save normalized summary to file
sink("output/advanced_modeling/normalized_data_summary.txt")
print("Normalized Data Summary:")
print(normalized_summary)
sink()

#*******************************************************************************
# Data visualization using R

## Scatter Plots

# Scatter Plot of Age vs Survival Rate grouped by Survival_Binary
survival_scatter <- ggplot(data, aes(x = Age, y = Survival_Rate, color = Survival_Binary)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("red", "blue"), labels = c("Deceased", "Survived")) +
  labs(title = "Survival Rate Across Ages by Survival Outcome",
       x = "Age (years)",
       y = "Survival Rate (%)",
       color = "Survival Status") +
  theme_minimal()
print(survival_scatter)
ggsave("plots/advanced_modeling/age_vs_survival_rate.png", survival_scatter)

# Correlation between age and Survival_Rate
age_survival_cor <- cor(data$Age, data$Survival_Rate, use = "complete.obs")
sink("output/advanced_modeling/age_survival_correlation.txt")
print(paste("Correlation between Age and Survival Rate:", age_survival_cor))
sink()

# Logistic regression to assess the impact of age and tumor size on survival
age_tumorsize_model <- glm(Survival_Binary ~ Age + Tumor_Size, data = data, family = "binomial")
age_tumorsize_summary <- summary(age_tumorsize_model)
sink("output/advanced_modeling/age_tumorsize_logistic_model.txt")
print("Logistic Regression: Impact of Age and Tumor Size on Survival")
print(age_tumorsize_summary)
sink()

#*******************************************************************************

# Scatter Plot of Age vs Tumor Size grouped by Survival_Binary
tumor_scatter <- ggplot(data, aes(x = Age, y = Tumor_Size, color = Survival_Binary)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("red", "blue"), labels = c("Deceased", "Survived")) +
  labs(title = "Tumor Size vs Age by Survival Status",
       x = "Age (years)",
       y = "Tumor Size (mm)",
       color = "Survival Status") +
  theme_minimal()
print(tumor_scatter)
ggsave("plots/advanced_modeling/age_vs_tumor_size.png", tumor_scatter)

# Correlation between Age and Tumor Size
age_tumor_cor <- cor(data$Age, data$Tumor_Size, use = "complete.obs")
sink("output/advanced_modeling/age_tumor_correlation.txt")
print(paste("Correlation between Age and Tumor Size:", age_tumor_cor))
sink()

#*******************************************************************************

# Scatter Plot of Tumor Size vs Tumor Growth Rate grouped by Survival_Binary
growth_scatter <- ggplot(data, aes(x = Tumor_Size, y = Tumor_Growth_Rate, color = Survival_Binary)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("red", "blue"), labels = c("Deceased", "Survived")) +
  labs(title = "Tumor Growth Rate vs Tumor Size by Survival Status",
       x = "Tumor Size (mm)",
       y = "Tumor Growth Rate (mm/month)",
       color = "Survival Status") +
  theme_minimal()
print(growth_scatter)
ggsave("plots/advanced_modeling/tumor_size_vs_growth_rate.png", growth_scatter)

# Correlation between Tumor Size and Growth Rate
tumor_growth_cor <- cor(data$Tumor_Size, data$Tumor_Growth_Rate, use = "complete.obs")
sink("output/advanced_modeling/tumor_growth_correlation.txt")
print(paste("Correlation between Tumor Size and Growth Rate:", tumor_growth_cor))
sink()

#*******************************************************************************

# Scatter Plot of Tumor Type vs Survival Rate grouped by Survival_Binary
type_survival_plot <- ggplot(data, aes(x = Tumor_Type, y = Survival_Rate, color = Survival_Binary)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  scale_color_manual(values = c("red", "blue"), labels = c("Deceased", "Survived")) +
  labs(title = "Survival Rate by Tumor Type",
       x = "Tumor Type",
       y = "Survival Rate (%)",
       color = "Survival Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(type_survival_plot)
ggsave("plots/advanced_modeling/tumor_type_vs_survival.png", type_survival_plot)

#*******************************************************************************

# Scatter Plot of Stage vs Survival Rate grouped by Survival_Binary
stage_survival_plot <- ggplot(data, aes(x = Stage, y = Survival_Rate, color = Survival_Binary)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
  scale_color_manual(values = c("red", "blue"), labels = c("Deceased", "Survived")) +
  labs(title = "Survival Rate by Cancer Stage",
       x = "Cancer Stage",
       y = "Survival Rate (%)",
       color = "Survival Status") +
  theme_minimal()
print(stage_survival_plot)
ggsave("plots/advanced_modeling/stage_vs_survival.png", stage_survival_plot)

# Logistic regression to assess the impact of stage on survival
stage_model <- glm(Survival_Binary ~ Stage, data = data, family = "binomial")
stage_summary <- summary(stage_model)
sink("output/advanced_modeling/stage_logistic_model.txt")
print("Logistic Regression: Impact of Cancer Stage on Survival")
print(stage_summary)
sink()

#*******************************************************************************

## Box Plots

# Box Plot of Tumor Size grouped by Survival_Binary
size_box <- ggplot(data, aes(x = Survival_Binary, y = Tumor_Size, fill = Survival_Binary)) +
  geom_boxplot() +
  scale_fill_manual(values = c("salmon", "darkcyan"), labels = c("Deceased", "Survived")) +
  labs(title = "Tumor Size Distribution by Survival Status",
       x = "Survival Status",
       y = "Tumor Size (mm)",
       fill = "Survival Status") +
  theme_minimal()
print(size_box)
ggsave("plots/advanced_modeling/tumor_size_boxplot.png", size_box)

# Summary statistics for tumor size by survival status
tumor_size_summary <- data %>%
  group_by(Survival_Binary) %>%
  summarise(
    Mean = mean(Tumor_Size, na.rm = TRUE),
    Median = median(Tumor_Size, na.rm = TRUE),
    SD = sd(Tumor_Size, na.rm = TRUE),
    Min = min(Tumor_Size, na.rm = TRUE),
    Max = max(Tumor_Size, na.rm = TRUE)
  )
sink("output/advanced_modeling/tumor_size_by_survival.txt")
print("Tumor Size Statistics by Survival Status:")
print(tumor_size_summary)
sink()

#*******************************************************************************

# Box Plot of Age grouped by Survival_Binary
age_box <- ggplot(data, aes(x = Survival_Binary, y = Age, fill = Survival_Binary)) +
  geom_boxplot() +
  scale_fill_manual(values = c("salmon", "darkcyan"), labels = c("Deceased", "Survived")) +
  labs(title = "Age Distribution by Survival Status",
       x = "Survival Status",
       y = "Age (years)",
       fill = "Survival Status") +
  theme_minimal()
print(age_box)
ggsave("plots/advanced_modeling/age_boxplot.png", age_box)

# Summary statistics for age by survival status
age_summary <- data %>%
  group_by(Survival_Binary) %>%
  summarise(
    Mean = mean(Age, na.rm = TRUE),
    Median = median(Age, na.rm = TRUE),
    SD = sd(Age, na.rm = TRUE),
    Min = min(Age, na.rm = TRUE),
    Max = max(Age, na.rm = TRUE)
  )
sink("output/advanced_modeling/age_by_survival.txt")
print("Age Statistics by Survival Status:")
print(age_summary)
sink()

#*******************************************************************************

# Box Plot of Tumor Growth Rate grouped by Treatment Status
# First, create a treatment indicator (received at least one treatment)
data$Treatment_Status <- as.factor(ifelse(
  data$Radiation_Treatment == "Yes" | 
  data$Surgery_Performed == "Yes" | 
  data$Chemotherapy == "Yes", 
  "Treated", "Untreated"))

growth_treatment_box <- ggplot(data, aes(x = Treatment_Status, y = Tumor_Growth_Rate, fill = Treatment_Status)) +
  geom_boxplot() +
  labs(title = "Tumor Growth Rate by Treatment Status",
       x = "Treatment Status",
       y = "Tumor Growth Rate (mm/month)",
       fill = "Treatment Status") +
  theme_minimal()
print(growth_treatment_box)
ggsave("plots/advanced_modeling/growth_rate_by_treatment.png", growth_treatment_box)

#*******************************************************************************

# Box Plot of Survival Rate grouped by Stage
survival_stage_box <- ggplot(data, aes(x = Stage, y = Survival_Rate, fill = Stage)) +
  geom_boxplot() +
  labs(title = "Survival Rate Distribution by Cancer Stage",
       x = "Cancer Stage",
       y = "Survival Rate (%)",
       fill = "Stage") +
  theme_minimal()
print(survival_stage_box)
ggsave("plots/advanced_modeling/survival_by_stage_boxplot.png", survival_stage_box)

#*******************************************************************************

# Descriptive statistics

cat("\nDetailed skim summary:\n")
skim_result <- skim(data)
sink("output/advanced_modeling/skim_summary.txt")
print(skim_result)
sink()

# Select columns
numeric_cols <- c("Age", "Tumor_Size", "Survival_Rate", "Tumor_Growth_Rate")
categorical_cols <- c("Gender", "Tumor_Type", "Location", "Histology", "Stage", 
                     "Radiation_Treatment", "Surgery_Performed", "Chemotherapy", 
                     "Family_History", "MRI_Result", "Follow_Up_Required", "Survival_Binary")

#*******************************************************************************
# --- UNIVARIATE ANALYSIS ---

# Function: Numeric stats
describe_numeric <- function(x) {
  stats <- data.frame(
    N = length(na.omit(x)),
    Mean = mean(x, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    SD = sd(x, na.rm = TRUE),
    Min = min(x, na.rm = TRUE),
    Max = max(x, na.rm = TRUE),
    Range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE),
    IQR = IQR(x, na.rm = TRUE),
    Skewness = skewness(x, na.rm = TRUE),
    Kurtosis = kurtosis(x, na.rm = TRUE)
  )
  return(stats)
}

# Function: Categorical stats
describe_categorical <- function(x) {
  freq <- table(x)
  prop <- prop.table(freq)
  stats <- data.frame(
    Category = names(freq),
    Frequency = as.numeric(freq),
    Percentage = round(as.numeric(prop) * 100, 2)
  )
  return(stats)
}

# Descriptive stats for numeric variables
numeric_stats <- lapply(data[numeric_cols], describe_numeric)
numeric_stats_df <- do.call(rbind, numeric_stats)
rownames(numeric_stats_df) <- numeric_cols
sink("output/advanced_modeling/numeric_descriptive_stats.txt")
print("Descriptive Statistics for Numeric Variables:")
print(numeric_stats_df)
sink()

# Descriptive stats for categorical variables
cat_stats <- lapply(data[categorical_cols], describe_categorical)
names(cat_stats) <- categorical_cols
sink("output/advanced_modeling/categorical_descriptive_stats.txt")
print("Descriptive Statistics for Categorical Variables:")
for (i in 1:length(cat_stats)) {
  cat(paste("\n", names(cat_stats)[i], ":\n", sep=""))
  print(cat_stats[[i]])
}
sink()

# Calculate the mode
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Visualizations: numeric variables
for (col in numeric_cols) {
  # Histogram with density
  p1 <- ggplot(data, aes(x = !!sym(col))) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
    geom_density(alpha = 0.5, fill = "orange") +
    ggtitle(paste("Histogram of", col)) +
    theme_minimal()
  
  # Boxplot
  p2 <- ggplot(data, aes(y = !!sym(col))) +
    geom_boxplot(fill = "lightgreen", color = "black") +
    ggtitle(paste("Boxplot of", col)) +
    theme_minimal()
  
  # Q-Q Plot
  p3 <- ggplot(data, aes(sample = !!sym(col))) +
    stat_qq() +
    stat_qq_line() +
    ggtitle(paste("Q-Q Plot of", col)) +
    theme_minimal()
  
  combined_plot <- grid.arrange(p1, p2, p3, ncol = 3)
  ggsave(paste0("plots/advanced_modeling/", col, "_distribution.png"), combined_plot, width = 15, height = 5)
}

#*******************************************************************************
# --- MULTIVARIATE ANALYSIS ---

# Correlation matrix
# Simplified selection of numeric columns without all_of()
cor_data <- data[, numeric_cols, drop = FALSE]
cor_matrix <- cor(cor_data, use = "complete.obs")

sink("output/advanced_modeling/correlation_matrix.txt")
print("Correlation Matrix:")
print(cor_matrix)
sink()

# Correlation heatmap
corr_plot <- corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8, addCoef.col = "black")
png("plots/advanced_modeling/correlation_heatmap.png", width = 800, height = 600)
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8, addCoef.col = "black")
dev.off()

# Scatterplot matrix
png("plots/advanced_modeling/scatterplot_matrix.png", width = 1000, height = 1000)
pairs(cor_data, main = "Scatterplot Matrix of Numeric Variables")
dev.off()

# --- DESCRIPTIVE STATS BY GROUP ---

sink("output/advanced_modeling/stats_by_gender.txt")
cat("Descriptive Stats by Gender:\n")
gender_desc <- describeBy(data[numeric_cols], group = data$Gender)
print(gender_desc)
sink()

sink("output/advanced_modeling/stats_by_survival.txt")
cat("Descriptive Stats by Survival Status:\n")
survival_desc <- describeBy(data[numeric_cols], group = data$Survival_Binary)
print(survival_desc)
sink()

#*******************************************************************************

# Hypothesis Testing

# 1. T-test: Compare Tumor Size between patients who survived and didn't survive
# Subsetting data for Survival_Binary = 0 and Survival_Binary = 1
size_survival_0 <- data$Tumor_Size[data$Survival_Binary == 0]
size_survival_1 <- data$Tumor_Size[data$Survival_Binary == 1]

# Performing the T-test
t_test_result <- t.test(size_survival_0, size_survival_1, var.equal = FALSE)
sink("output/advanced_modeling/t_test_tumor_size.txt")
print("T-test Result (Tumor Size vs Survival):")
print(t_test_result)
sink()

# Plotting boxplot for T-test
data$survival_factor <- factor(data$Survival_Binary, labels = c("Deceased", "Survived"))
survival_boxplot <- ggplot(data, aes(x = survival_factor, y = Tumor_Size, fill = survival_factor)) +
  geom_boxplot() +
  labs(title = "Tumor Size Distribution by Survival Status",
       x = "Survival Status",
       y = "Tumor Size (mm)") +
  theme_minimal()
print(survival_boxplot)
ggsave("plots/advanced_modeling/tumor_size_survival_boxplot.png", survival_boxplot)

#*******************************************************************************

# 2. Correlation: Between Tumor Size and Survival Rate
cor_test_result <- cor.test(data$Tumor_Size, data$Survival_Rate, method = "pearson")
sink("output/advanced_modeling/cor_test_size_survival.txt")
print("Correlation Test Result (Tumor Size vs Survival Rate):")
print(cor_test_result)
sink()

# Plotting scatterplot for correlation
size_survival_scatter <- ggplot(data, aes(x = Tumor_Size, y = Survival_Rate, color = survival_factor)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Scatterplot of Tumor Size vs Survival Rate",
       x = "Tumor Size (mm)",
       y = "Survival Rate (%)") +
  theme_minimal()
print(size_survival_scatter)
ggsave("plots/advanced_modeling/size_vs_survival_scatterplot.png", size_survival_scatter)

#*******************************************************************************

# 3. Chi-Square Test: Association between Radiation Treatment and Survival
# Creating a contingency table
contingency_table <- table(data$Radiation_Treatment, data$Survival_Binary)
rownames(contingency_table) <- c("No Radiation", "Radiation")
colnames(contingency_table) <- c("Deceased", "Survived")

# Performing Chi-Square test
chi_square_result <- chisq.test(contingency_table)
sink("output/advanced_modeling/chi_square_radiation_survival.txt")
print("Chi-Square Test Result (Radiation Treatment vs Survival):")
print(chi_square_result)
print("Contingency Table:")
print(contingency_table)
sink()

# Plotting bar plot for Chi-Square test
radiation_survival_df <- as.data.frame(contingency_table)
colnames(radiation_survival_df) <- c("Radiation", "Survival", "Count")
radiation_survival_df$Radiation <- factor(radiation_survival_df$Radiation, 
                                       labels = c("No Radiation", "Radiation"))
radiation_survival_df$Survival <- factor(radiation_survival_df$Survival, 
                                      labels = c("Deceased", "Survived"))

radiation_barplot <- ggplot(radiation_survival_df, aes(x = Radiation, y = Count, fill = Survival)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Radiation Treatment vs Survival Status",
       x = "Radiation Treatment",
       y = "Count") +
  theme_minimal()
print(radiation_barplot)
ggsave("plots/advanced_modeling/radiation_vs_survival_barplot.png", radiation_barplot)

#*******************************************************************************

# 4. ANOVA: Compare Tumor Size across cancer stages
# Performing ANOVA
anova_result <- aov(Tumor_Size ~ Stage, data = data)
sink("output/advanced_modeling/anova_tumor_size_stage.txt")
print("ANOVA Result (Tumor Size vs Cancer Stage):")
print(summary(anova_result))
sink()

# Plotting boxplot for ANOVA
stage_boxplot <- ggplot(data, aes(x = Stage, y = Tumor_Size, fill = Stage)) +
  geom_boxplot() +
  labs(title = "Tumor Size Distribution by Cancer Stage",
       x = "Cancer Stage",
       y = "Tumor Size (mm)") +
  theme_minimal()
print(stage_boxplot)
ggsave("plots/advanced_modeling/tumor_size_by_stage_boxplot.png", stage_boxplot)

#*******************************************************************************

# Machine Learning

# Set seed for reproducibility
set.seed(123)

# Selecting features and target
features <- c("Age", "Gender", "Tumor_Type", "Tumor_Size", "Location", 
              "Histology", "Stage", "Tumor_Growth_Rate", "Radiation_Treatment", 
              "Surgery_Performed", "Chemotherapy", "Family_History")
data_model <- data[, c(features, "Survival_Binary")]

# Data Preprocessing: Handle outliers using IQR for numeric columns
numeric_model_cols <- c("Age", "Tumor_Size", "Tumor_Growth_Rate")
for (col in numeric_model_cols) {
  outlier_filter <- function(x) {
    iqr <- IQR(x)
    q <- quantile(x, probs = c(0.25, 0.75))
    lower <- q[1] - 1.5 * iqr
    upper <- q[2] + 1.5 * iqr
    x[x < lower | x > upper] <- NA
    return(x)
  }
  data_model[[col]] <- outlier_filter(data_model[[col]])
}
data_model <- na.omit(data_model)

# Converting categorical variables to factors (if not already done)
for (col in setdiff(features, numeric_model_cols)) {
  if (!is.factor(data_model[[col]])) {
    data_model[[col]] <- as.factor(data_model[[col]])
  }
}

# Train-test split (70% train, 30% test)
trainIndex <- createDataPartition(data_model$Survival_Binary, p = 0.7, list = FALSE)
train_data <- data_model[trainIndex, ]
test_data <- data_model[-trainIndex, ]

# Address class imbalance using ROSE (SMOTE-like oversampling/undersampling)
train_data <- ovun.sample(Survival_Binary ~ ., data = train_data, method = "both", 
                          p = 0.5, seed = 123)$data

# Feature Selection using Random Forest (initial model to get importance)
initial_rf <- randomForest(Survival_Binary ~ ., data = train_data, ntree = 50)
importance_scores <- importance(initial_rf)
important_features <- rownames(importance_scores[order(importance_scores[, "MeanDecreaseGini"], decreasing = TRUE), ])[1:8]
sink("output/advanced_modeling/feature_importance.txt")
print("Top 8 Important Features:")
print(important_features)
sink()

# Plot feature importance
importance_df <- data.frame(
  Feature = rownames(importance_scores),
  Importance = importance_scores[, "MeanDecreaseGini"]
)
importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]
feature_plot <- ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance", x = "Feature", y = "Importance") +
  theme_minimal()
print(feature_plot)
ggsave("plots/advanced_modeling/feature_importance.png", feature_plot)

# Prepare data with selected features
train_data_subset <- train_data[, c(important_features, "Survival_Binary")]
test_data_subset <- test_data[, c(important_features, "Survival_Binary")]

# Create interaction terms
# Assuming tumor size and stage are important features and present in the top 8
if ("Tumor_Size" %in% important_features && "Stage" %in% important_features) {
  # Convert Stage to numeric for interaction
  stage_levels <- levels(train_data_subset$Stage)
  train_data_subset$Stage_num <- as.numeric(train_data_subset$Stage)
  test_data_subset$Stage_num <- as.numeric(test_data_subset$Stage)
  
  train_data_subset$Size_Stage <- train_data_subset$Tumor_Size * train_data_subset$Stage_num
  test_data_subset$Size_Stage <- test_data_subset$Tumor_Size * test_data_subset$Stage_num
  
  # Add interaction term to important features
  important_features <- c(important_features, "Size_Stage")
}

# Total number of test predictions (for percentage calculation)
total_test <- nrow(test_data)

# Define cross-validation control
train_control <- trainControl(method = "cv", number = 5, summaryFunction = twoClassSummary, classProbs = TRUE)

#*******************************************************************************
# 1. Logistic Regression with glmnet (regularization) and threshold optimization

# Create model matrix for glmnet
x_train <- model.matrix(Survival_Binary ~ ., data = train_data_subset)[, -1]
x_test <- model.matrix(Survival_Binary ~ ., data = test_data_subset)[, -1]
y_train <- train_data_subset$Survival_Binary

# Train logistic regression model with regularization
log_model <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, nfolds = 5)
log_pred_prob <- predict(log_model, newx = x_test, type = "response", s = "lambda.min")

# Fix for roc(): Convert matrix to numeric vector
log_pred_prob_vec <- as.numeric(log_pred_prob)
roc_obj <- roc(test_data_subset$Survival_Binary, log_pred_prob_vec)
best_threshold <- coords(roc_obj, "best", ret = "threshold")$threshold
log_pred <- ifelse(log_pred_prob_vec > best_threshold, "1", "0")
log_pred <- factor(log_pred, levels = c("0", "1"))

# Calculate metrics
log_cm <- confusionMatrix(log_pred, test_data_subset$Survival_Binary)
log_accuracy <- mean(log_pred == test_data_subset$Survival_Binary)
log_f1 <- log_cm$byClass["F1"]
log_auc <- auc(roc_obj)

# Print metrics
sink("output/advanced_modeling/logistic_regression_metrics.txt")
print("Logistic Regression Metrics:")
print(paste("Accuracy:", round(log_accuracy, 3)))
print(paste("F1-Score:", round(log_f1, 3)))
print(paste("AUC-ROC:", round(log_auc, 3)))
sink()

# Confusion Matrix for Logistic Regression with percentages
sink("output/advanced_modeling/logistic_regression_confusion_matrix.txt")
print("Logistic Regression Confusion Matrix (Counts and Percentages):")
log_cm_table <- log_cm$table
log_cm_percent <- round((log_cm_table / total_test) * 100, 2)
log_cm_combined <- matrix(paste0(log_cm_table, " (", log_cm_percent, "%)"), nrow = 2)
dimnames(log_cm_combined) <- dimnames(log_cm_table)
print(log_cm_combined)
sink()

# Plotting Confusion Matrix for Logistic Regression
log_cm_df <- as.data.frame(as.table(log_cm$table))
colnames(log_cm_df) <- c("Predicted", "Actual", "Freq")
log_cm_plot <- ggplot(log_cm_df, aes(x = Actual, y = Predicted, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "white") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Confusion Matrix: Logistic Regression",
       x = "Actual Survival Status",
       y = "Predicted Survival Status") +
  theme_minimal()
print(log_cm_plot)
ggsave("plots/advanced_modeling/logistic_regression_confusion_matrix.png", log_cm_plot)

#*******************************************************************************