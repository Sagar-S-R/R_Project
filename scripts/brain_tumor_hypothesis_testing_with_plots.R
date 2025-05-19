# Installing required packages if not already installed
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

# Loading the libraries
library(ggplot2)
library(dplyr)

# Importing the dataset
data <- read.csv("data/processed/brain_tumor_data_clean.csv")

# Create the directory structure for saving plots if it doesn't exist
output_dir <- "output/relationShip_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Output and Results Section
cat("=== Output and Results ===\n\n")

# 1. T-test: Compare Tumor_Size between Malignant and Benign tumors
cat("1. T-test: Tumor Size vs Tumor Type (Malignant/Benign)\n")
cat("-------------------------------------------------\n")

# Subsetting data for Tumor_Type = Malignant and Benign
tumor_size_malignant <- data$Tumor_Size[data$Tumor_Type == "Malignant"]
tumor_size_benign <- data$Tumor_Size[data$Tumor_Type == "Benign"]

# Performing the T-test
t_test_result <- t.test(tumor_size_malignant, tumor_size_benign, var.equal = FALSE)
cat("T-test Result (Tumor_Size vs Tumor_Type):\n")
print(t_test_result)
cat("\n")

# Plotting boxplot for T-test
data$Tumor_Type_factor <- factor(data$Tumor_Type, levels = c("Benign", "Malignant"))
tumor_size_plot <- ggplot(data, aes(x = Tumor_Type_factor, y = Tumor_Size, fill = Tumor_Type_factor)) +
  geom_boxplot() +
  labs(title = "Tumor Size Distribution by Tumor Type",
       x = "Tumor Type",
       y = "Tumor Size (cm)") +
  theme_minimal()
print(tumor_size_plot)
ggsave(file.path(output_dir, "tumor_size_vs_tumor_type_boxplot.png"))
cat("Boxplot saved as: output/relationShip_analysis/tumor_size_vs_tumor_type_boxplot.png\n\n")

# 2. Paired T-test: Comparing Tumor_Size and Tumor_Growth_Rate
cat("2. Paired T-test: Tumor Size vs Tumor Growth Rate\n")
cat("-------------------------------------------------\n")

# Note: This is for demonstration, as Tumor_Size and Tumor_Growth_Rate aren't true paired measures
paired_t_test_result <- t.test(data$Tumor_Size, data$Tumor_Growth_Rate, paired = TRUE)
cat("Paired T-test Result (Tumor_Size vs Tumor_Growth_Rate):\n")
print(paired_t_test_result)
cat("\n")
# No plot for paired T-test, as the original code didn't include one

# 3. Correlation: Between Tumor_Size and Survival_Rate
cat("3. Correlation: Tumor Size vs Survival Rate\n")
cat("-------------------------------------------------\n")

# Performing correlation test
cor_test_result <- cor.test(data$Tumor_Size, data$Survival_Rate, method = "pearson")
cat("Correlation Test Result (Tumor_Size vs Survival_Rate):\n")
print(cor_test_result)
cat("\n")

# Plotting scatterplot for correlation
tumor_survival_plot <- ggplot(data, aes(x = Tumor_Size, y = Survival_Rate)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "blue") +
  labs(title = "Scatterplot of Tumor Size vs Survival Rate",
       x = "Tumor Size (cm)",
       y = "Survival Rate (%)") +
  theme_minimal()
print(tumor_survival_plot)
ggsave(file.path(output_dir, "tumor_size_vs_survival_rate_scatterplot.png"))
cat("Scatterplot saved as: output/relationShip_analysis/tumor_size_vs_survival_rate_scatterplot.png\n\n")

# 4. Chi-Square Test: Association between Radiation_Treatment and Tumor_Type
cat("4. Chi-Square Test: Radiation Treatment vs Tumor Type\n")
cat("-------------------------------------------------\n")

# Creating a contingency table
contingency_table <- table(data$Radiation_Treatment, data$Tumor_Type)
rownames(contingency_table) <- c("No Radiation", "Radiation")
colnames(contingency_table) <- c("Benign", "Malignant")

# Performing Chi-Square test
chi_square_result <- chisq.test(contingency_table)
cat("Chi-Square Test Result (Radiation_Treatment vs Tumor_Type):\n")
print(chi_square_result)
cat("\n")

# Plotting bar plot for Chi-Square test
radiation_tumor_df <- as.data.frame(contingency_table)
colnames(radiation_tumor_df) <- c("Radiation", "Tumor_Type", "Count")
radiation_tumor_df$Radiation <- factor(radiation_tumor_df$Radiation, labels = c("No Radiation", "Radiation"))
radiation_tumor_df$Tumor_Type <- factor(radiation_tumor_df$Tumor_Type, labels = c("Benign", "Malignant"))

radiation_tumor_plot <- ggplot(radiation_tumor_df, aes(x = Radiation, y = Count, fill = Tumor_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Radiation Treatment vs Tumor Type",
       x = "Radiation Treatment",
       y = "Count") +
  theme_minimal()
print(radiation_tumor_plot)
ggsave(file.path(output_dir, "radiation_vs_tumor_type_barplot.png"))
cat("Bar plot saved as: output/relationShip_analysis/radiation_vs_tumor_type_barplot.png\n\n")

# 5. ANOVA: Compare Survival_Rate across tumor Stage levels
cat("5. ANOVA: Survival Rate vs Tumor Stage\n")
cat("-------------------------------------------------\n")

# Converting Stage to factor
data$Stage_factor <- factor(data$Stage, levels = c("I", "II", "III", "IV"))

# Performing ANOVA
anova_result <- aov(Survival_Rate ~ Stage_factor, data = data)
cat("ANOVA Result (Survival_Rate vs Stage):\n")
print(summary(anova_result))
cat("\n")

# Plotting boxplot for ANOVA
survival_stage_plot <- ggplot(data, aes(x = Stage_factor, y = Survival_Rate, fill = Stage_factor)) +
  geom_boxplot() +
  labs(title = "Survival Rate Distribution by Tumor Stage",
       x = "Tumor Stage",
       y = "Survival Rate (%)") +
  theme_minimal()
print(survival_stage_plot)
ggsave(file.path(output_dir, "survival_rate_vs_stage_boxplot.png"))
cat("Boxplot saved as: output/relationShip_analysis/survival_rate_vs_stage_boxplot.png\n\n")

cat("=== End of Output and Results ===\n")