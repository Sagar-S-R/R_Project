# Create directory if it doesn't exist
if (!dir.exists("output/advanced_analyses")) {
  dir.create("output/advanced_analyses", recursive = TRUE)
}

# Start logging to a file
sink("output/advanced_analyses/statistical_analysis_report.txt", split = TRUE)

# --- INSTALLATION (optional run once) ---
pkgs <- c("tidyverse", "psych", "moments", "gridExtra", "corrplot", "ggpubr", "skimr")
installed <- pkgs %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(pkgs[!installed], repos = "https://cran.rstudio.com/")
}

# --- LIBRARY LOADING ---
library(tidyverse)
library(psych)
library(moments)
library(gridExtra)
library(corrplot)
library(ggpubr)
library(skimr)

# --- DATA LOADING ---
data <- read.csv("data/processed/brain_tumor_data_clean.csv")

# --- BASIC STRUCTURE ---
cat("Structure of the dataset:\n")
str(data)

cat("\nSummary of the dataset:\n")
summary(data)

cat("\nDetailed skim summary:\n")
print(skim(data))

# --- CONVERT CATEGORICAL VARIABLES TO FACTORS ---
data <- data %>%
  mutate(
    Gender = factor(Gender, levels = c("Male", "Female")),
    Tumor_Type = factor(Tumor_Type, levels = c("Benign", "Malignant")),
    Location = factor(Location),
    Histology = factor(Histology),
    Stage = factor(Stage, levels = c("I", "II", "III", "IV")),
    Symptom_1 = factor(Symptom_1),
    Symptom_2 = factor(Symptom_2),
    Symptom_3 = factor(Symptom_3),
    Radiation_Treatment = factor(Radiation_Treatment, levels = c(0, 1), labels = c("No", "Yes")),
    Surgery_Performed = factor(Surgery_Performed, levels = c(0, 1), labels = c("No", "Yes")),
    Chemotherapy = factor(Chemotherapy, levels = c(0, 1), labels = c("No", "Yes")),
    Family_History = factor(Family_History, levels = c(0, 1), labels = c("No", "Yes")),
    Follow_Up_Required = factor(Follow_Up_Required, levels = c(0, 1), labels = c("No", "Yes")),
    Survival_Binary = factor(Survival_Binary, levels = c("Poor", "Good"))
    # Note: MRI_Result is numeric, not converted to factor
    # If Survival_Binary needs thresholding (e.g., Survival_Rate < 80 = "Poor"), uncomment:
    # Survival_Binary = factor(ifelse(Survival_Rate >= 80, "Good", "Poor"), levels = c("Poor", "Good"))
  )

# Define numeric and categorical columns
numeric_cols <- c("Age", "Tumor_Size", "Survival_Rate", "Tumor_Growth_Rate", "MRI_Result")
categorical_cols <- c("Gender", "Tumor_Type", "Location", "Histology", "Stage", 
                     "Symptom_1", "Symptom_2", "Symptom_3", "Radiation_Treatment", 
                     "Surgery_Performed", "Chemotherapy", "Family_History", 
                     "Follow_Up_Required", "Survival_Binary")

# --- DESCRIPTIVE STATISTICS FUNCTIONS ---
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

# --- NUMERIC VARIABLES STATS ---
cat("\nDescriptive Statistics for Numeric Variables:\n")
numeric_stats <- lapply(data[numeric_cols], describe_numeric)
numeric_stats_df <- do.call(rbind, numeric_stats)
rownames(numeric_stats_df) <- numeric_cols
print(numeric_stats_df)

# --- CATEGORICAL VARIABLES STATS ---
cat("\nDescriptive Statistics for Categorical Variables:\n")
cat_stats <- lapply(data[categorical_cols], describe_categorical)
names(cat_stats) <- categorical_cols
print(cat_stats)

# --- CORRELATION MATRIX ---
cor_data <- data %>% select(all_of(numeric_cols))
cor_matrix <- cor(cor_data, use = "complete.obs")

cat("\nCorrelation Matrix:\n")
print(cor_matrix)

# --- DESCRIPTIVE STATS BY GROUP ---
cat("\nDescriptive Stats by Tumor Type (Malignant vs Benign):\n")
tumor_type_desc <- describeBy(data[numeric_cols], group = data$Tumor_Type)
print(tumor_type_desc)

cat("\nDescriptive Stats by Survival Binary (Good vs Poor):\n")
survival_desc <- describeBy(data[numeric_cols], group = data$Survival_Binary)
print(survival_desc)

# Stop writing to file
sink()



# Create output directory for plots if it doesn't exist
if (!dir.exists("plots/advanced_analyses")) {
  dir.create("plots/advanced_analyses", recursive = TRUE)
}

# Define a colorful palette
color_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")  # Blue, Orange, Green, Red, Purple

# Save histogram and boxplot for each numeric variable
for (i in seq_along(numeric_cols)) {
  col <- numeric_cols[i]
  col_color <- color_palette[i %% length(color_palette) + 1]  # Cycle through colors
  
  # Histogram + Density
  p1 <- ggplot(data, aes_string(x = col)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = col_color, color = "black", alpha = 0.7) +
    geom_density(alpha = 0.3, fill = col_color, color = col_color) +
    ggtitle(paste("Histogram of", col)) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          axis.text = element_text(color = "black"),
)
  
  ggsave(filename = paste0("plots/advanced_analyses/", col, "_hist_density.png"),
         plot = p1, width = 6, height = 4, dpi = 300)

  # Boxplot
  p2 <- ggplot(data, aes_string(y = col)) +
    geom_boxplot(fill = col_color, color = "black", alpha = 0.7) +
    ggtitle(paste("Boxplot of", col)) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          axis.text = element_text(color = "black"))
  
  ggsave(filename = paste0("plots/advanced_analyses/", col, "_boxplot.png"),
         plot = p2, width = 4, height = 4, dpi = 300)
}

# Save correlation matrix plot
png("plots/advanced_analyses/correlation_matrix.png", width = 800, height = 600, bg = "white")
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8, addCoef.col = "black",
         col = colorRampPalette(c("#d62728", "#ffffff", "#1f77b4"))(200))
dev.off()

# Rest of the script (sink() closing, etc.) remains the same
sink()