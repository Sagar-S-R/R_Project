if (!require("tidyverse")) install.packages("tidyverse")
if (!require("corrplot")) install.packages("corrplot")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

library(tidyverse)
library(corrplot)
library(ggplot2)
library(dplyr)

# Create output directories
output_dir <- "output/eda"
plots_dir <- "plots/eda"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Load cleaned data
brain_tumor_data <- read.csv("data/processed/brain_tumor_data_clean.csv", stringsAsFactors = TRUE)

# Ensure factor levels are consistent
brain_tumor_data$Tumor_Type <- factor(brain_tumor_data$Tumor_Type, levels = c("Benign", "Malignant"))
brain_tumor_data$Stage <- factor(brain_tumor_data$Stage, levels = c("I", "II", "III", "IV"))
brain_tumor_data$Survival_Binary <- factor(brain_tumor_data$Survival_Binary, levels = c("Poor", "Good"))

# Check for missing data
if (anyNA(brain_tumor_data)) {
  cat("Warning: Missing data detected in brain_tumor_data.\n")
  cat("Number of missing values per column:\n")
  print(colSums(is.na(brain_tumor_data)))
  cat("\n")
}
# Write text output to a single file
output_file <- file.path(output_dir, "eda_results.txt")
sink(output_file)

cat("=== Exploratory Data Analysis Results ===\n\n")

# ---- ANALYSIS 1: Patient Demographics ----
cat("1. Patient Demographics\n")
cat("-----------------------\n")
demographics <- brain_tumor_data %>%
  summarise(
    count = n(),
    mean_age = mean(Age, na.rm = TRUE),
    median_age = median(Age, na.rm = TRUE),
    min_age = min(Age, na.rm = TRUE),
    max_age = max(Age, na.rm = TRUE),
    gender_ratio = sum(Gender == "Male") / sum(Gender == "Female")
  )
cat("Demographics Summary:\n")
print(demographics)
cat("\n")
write.csv(demographics, file.path(output_dir, "demographics_summary.csv"), row.names = FALSE)

# Age distribution
age_plot <- ggplot2::ggplot(brain_tumor_data, ggplot2::aes(x = Age)) +
  ggplot2::geom_histogram(binwidth = 5, fill = "#1f77b4", color = "black") +
  ggplot2::labs(title = "Age Distribution of Brain Tumor Patients", x = "Age", y = "Count") +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA),
                 panel.background = ggplot2::element_rect(fill = "white", color = NA))
ggplot2::ggsave(file.path(plots_dir, "age_distribution.png"), age_plot, width = 8, height = 6)
cat("Age distribution plot saved as: plots/eda/age_distribution.png\n\n")

# Gender distribution
gender_plot <- ggplot2::ggplot(brain_tumor_data, ggplot2::aes(x = Gender, fill = Gender)) +
  ggplot2::geom_bar() +
  ggplot2::labs(title = "Gender Distribution", x = "Gender", y = "Count") +
  ggplot2::scale_fill_manual(values = c("Male" = "#ff7f0e", "Female" = "#d62728")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA),
                 panel.background = ggplot2::element_rect(fill = "white", color = NA))
ggplot2::ggsave(file.path(plots_dir, "gender_distribution.png"), gender_plot, width = 8, height = 6)
cat("Gender distribution plot saved as: plots/eda/gender_distribution.png\n\n")

# ---- ANALYSIS 2: Tumor Characteristics ----
cat("2. Tumor Characteristics\n")
cat("-----------------------\n")
# Tumor type distribution
tumor_type_plot <- ggplot2::ggplot(brain_tumor_data, ggplot2::aes(x = Tumor_Type, fill = Tumor_Type)) +
  ggplot2::geom_bar() +
  ggplot2::labs(title = "Distribution of Brain Tumor Types", x = "Tumor Type", y = "Count") +
  ggplot2::scale_fill_manual(values = c("Benign" = "#1f77b4", "Malignant" = "#ff7f0e")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 plot.background = ggplot2::element_rect(fill = "white", color = NA),
                 panel.background = ggplot2::element_rect(fill = "white", color = NA))
ggplot2::ggsave(file.path(plots_dir, "tumor_type_distribution.png"), tumor_type_plot, width = 10, height = 6)
cat("Tumor type distribution plot saved as: plots/eda/tumor_type_distribution.png\n\n")

# Tumor size distribution
tumor_size_plot <- ggplot2::ggplot(brain_tumor_data, ggplot2::aes(x = Tumor_Size)) +
  ggplot2::geom_histogram(binwidth = 0.5, fill = "#2ca02c", color = "black") +
  ggplot2::labs(title = "Brain Tumor Size Distribution", x = "Tumor Size (cm)", y = "Count") +
  ggplot2::theme_minimal() +
  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA),
                 panel.background = ggplot2::element_rect(fill = "white", color = NA))
ggplot2::ggsave(file.path(plots_dir, "tumor_size_distribution.png"), tumor_size_plot, width = 8, height = 6)
cat("Tumor size distribution plot saved as: plots/eda/tumor_size_distribution.png\n\n")

# Tumor location distribution
location_plot <- ggplot2::ggplot(brain_tumor_data, ggplot2::aes(x = Location, fill = Location)) +
  ggplot2::geom_bar() +
  ggplot2::labs(title = "Brain Tumor Location Distribution", x = "Location", y = "Count") +
  ggplot2::scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b")) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 plot.background = ggplot2::element_rect(fill = "white", color = NA),
                 panel.background = ggplot2::element_rect(fill = "white", color = NA))
ggplot2::ggsave(file.path(plots_dir, "location_distribution.png"), location_plot, width = 10, height = 6)
cat("Tumor location distribution plot saved as: plots/eda/location_distribution.png\n\n")

# ---- ANALYSIS 3: Statistical Tests ----
cat("3. Statistical Tests\n")
cat("-----------------------\n")

# T-test: Tumor Size vs Tumor Type (Malignant/Benign)
cat("3.1 T-test: Tumor Size vs Tumor Type (Malignant/Benign)\n")
tumor_size_malignant <- brain_tumor_data$Tumor_Size[brain_tumor_data$Tumor_Type == "Malignant"]
tumor_size_benign <- brain_tumor_data$Tumor_Size[brain_tumor_data$Tumor_Type == "Benign"]
t_test_result <- stats::t.test(tumor_size_malignant, tumor_size_benign, var.equal = FALSE)
cat("T-test Result (Tumor_Size vs Tumor_Type):\n")
print(t_test_result)
cat("\n")
tumor_size_plot <- ggplot2::ggplot(brain_tumor_data, ggplot2::aes(x = Tumor_Type, y = Tumor_Size, fill = Tumor_Type)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = c("Benign" = "#1f77b4", "Malignant" = "#ff7f0e")) +
  ggplot2::labs(title = "Tumor Size Distribution by Tumor Type",
                x = "Tumor Type",
                y = "Tumor Size (cm)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", color = "white"),
    plot.background = ggplot2::element_rect(fill = "white"),
    legend.title = ggplot2::element_blank()
  )
ggplot2::ggsave(file.path(plots_dir, "tumor_size_vs_tumor_type_boxplot.png"), tumor_size_plot, width = 8, height = 6)
cat("Boxplot saved as: plots/eda/tumor_size_vs_tumor_type_boxplot.png\n\n")

# Chi-Square Test: Radiation Treatment vs Tumor Type
cat("3.2 Chi-Square Test: Radiation Treatment vs Tumor Type\n")
contingency_table <- table(brain_tumor_data$Radiation_Treatment, brain_tumor_data$Tumor_Type)
rownames(contingency_table) <- c("No Radiation", "Radiation")
colnames(contingency_table) <- c("Benign", "Malignant")
chi_square_result <- stats::chisq.test(contingency_table)
cat("Chi-Square Test Result (Radiation_Treatment vs Tumor_Type):\n")
print(chi_square_result)
cat("\n")
radiation_tumor_df <- as.data.frame(contingency_table)
colnames(radiation_tumor_df) <- c("Radiation", "Tumor_Type", "Count")
radiation_tumor_df$Radiation <- factor(radiation_tumor_df$Radiation, labels = c("No Radiation", "Radiation"))
radiation_tumor_df$Tumor_Type <- factor(radiation_tumor_df$Tumor_Type, labels = c("Benign", "Malignant"))
radiation_tumor_plot <- ggplot2::ggplot(radiation_tumor_df, ggplot2::aes(x = Radiation, y = Count, fill = Tumor_Type)) +
  ggplot2::geom_bar(stat = "identity", position = "dodge") +
  ggplot2::scale_fill_manual(values = c("Benign" = "#1f77b4", "Malignant" = "#ff7f0e")) +
  ggplot2::labs(title = "Radiation Treatment vs Tumor Type",
                x = "Radiation Treatment",
                y = "Count") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", color = "white"),
    plot.background = ggplot2::element_rect(fill = "white")
  )
ggplot2::ggsave(file.path(plots_dir, "radiation_vs_tumor_type_barplot.png"), radiation_tumor_plot, width = 8, height = 6)
cat("Bar plot saved as: plots/eda/radiation_vs_tumor_type_barplot.png\n\n")

# ANOVA: Survival Rate vs Tumor Stage
cat("3.3 ANOVA: Survival Rate vs Tumor Stage\n")
anova_result <- stats::aov(Survival_Rate ~ Stage, data = brain_tumor_data)
cat("ANOVA Result (Survival_Rate vs Stage):\n")
print(summary(anova_result))
cat("\n")
survival_stage_plot <- ggplot2::ggplot(brain_tumor_data, ggplot2::aes(x = Stage, y = Survival_Rate, fill = Stage)) +
  ggplot2::geom_boxplot() +
  ggplot2::scale_fill_manual(values = c("I" = "#1f77b4", "II" = "#ff7f0e", "III" = "#2ca02c", "IV" = "#d62728")) +
  ggplot2::labs(title = "Survival Rate Distribution by Tumor Stage",
                x = "Tumor Stage",
                y = "Survival Rate (%)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", color = "white"),
    plot.background = ggplot2::element_rect(fill = "white"),
    legend.title = ggplot2::element_blank()
  )
ggplot2::ggsave(file.path(plots_dir, "survival_rate_vs_stage_boxplot.png"), survival_stage_plot, width = 8, height = 6)
cat("Boxplot saved as: plots/eda/survival_rate_vs_stage_boxplot.png\n\n")

# ---- ANALYSIS 4: Correlation Analysis ----
cat("4. Correlation Analysis\n")
cat("-----------------------\n")
numeric_vars <- brain_tumor_data %>%
  dplyr::select(Age, Tumor_Size, Survival_Rate, Tumor_Growth_Rate)
correlation_matrix <- cor(numeric_vars, use = "complete.obs")
cat("Correlation Matrix:\n")
print(correlation_matrix)
cat("\n")
write.csv(correlation_matrix, file.path(output_dir, "correlation_matrix.csv"), row.names = TRUE)
png(file.path(plots_dir, "correlation_plot.png"), width = 800, height = 800)
corrplot::corrplot(correlation_matrix, method = "circle", type = "upper", order = "hclust",
                   tl.col = "black", tl.srt = 45, col = colorRampPalette(c("#d62728", "#1f77b4"))(200))
dev.off()
cat("Correlation plot saved as: plots/eda/correlation_plot.png\n\n")

# Correlation: Tumor Size vs Survival Rate
cor_test_result <- stats::cor.test(brain_tumor_data$Tumor_Size, brain_tumor_data$Survival_Rate, method = "pearson")
cat("Correlation Test Result (Tumor_Size vs Survival_Rate):\n")
print(cor_test_result)
cat("\n")
tumor_survival_plot <- ggplot2::ggplot(brain_tumor_data, ggplot2::aes(x = Tumor_Size, y = Survival_Rate)) +
  ggplot2::geom_point(alpha = 0.7, color = "#2ca02c", size = 3) +
  ggplot2::geom_smooth(method = "lm", color = "#d62728", fill = alpha("#d62728", 0.2)) +
  ggplot2::labs(title = "Scatterplot of Tumor Size vs Survival Rate",
                x = "Tumor Size (cm)",
                y = "Survival Rate (%)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", color = "white"),
    plot.background = ggplot2::element_rect(fill = "white")
  )
ggplot2::ggsave(file.path(plots_dir, "tumor_size_vs_survival_rate_scatterplot.png"), tumor_survival_plot, width = 8, height = 6)
cat("Scatterplot saved as: plots/eda/tumor_size_vs_survival_rate_scatterplot.png\n\n")

cat("=== End of Exploratory Data Analysis Results ===\n")
sink()