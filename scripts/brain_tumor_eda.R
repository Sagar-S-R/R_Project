# Brain_Tumor_Analysis_Project/scripts/brain_tumor_eda.R

# Load required libraries
library(tidyverse)
library(corrplot)

# Load cleaned data
brain_tumor_data <- read.csv("data/processed/brain_tumor_data_clean.csv", stringsAsFactors = TRUE)

# ---- ANALYSIS 1: Basic Patient Demographics ----
demographics <- brain_tumor_data %>%
  summarise(
    count = n(),
    mean_age = mean(Age),
    median_age = median(Age),
    min_age = min(Age),
    max_age = max(Age),
    gender_ratio = sum(Gender == "Male") / sum(Gender == "Female")
  )

write.csv(demographics, "output/eda/demographics_summary.csv", row.names = FALSE)

# Age distribution
age_plot <- ggplot(brain_tumor_data, aes(x = Age)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black") +
  labs(title = "Age Distribution of Brain Tumor Patients", x = "Age", y = "Count") +
  theme_minimal()
ggsave("plots/eda/age_distribution.png", age_plot, width = 8, height = 6)

# Gender distribution
gender_plot <- ggplot(brain_tumor_data, aes(x = Gender, fill = Gender)) +
  geom_bar() +
  labs(title = "Gender Distribution", x = "Gender", y = "Count") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")
ggsave("plots/eda/gender_distribution.png", gender_plot, width = 8, height = 6)

# ---- ANALYSIS 2: Tumor Characteristics Analysis ----
tumor_type_plot <- ggplot(brain_tumor_data, aes(x = Tumor_Type, fill = Tumor_Type)) +
  geom_bar() +
  labs(title = "Distribution of Brain Tumor Types", x = "Tumor Type", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")
ggsave("plots/eda/tumor_type_distribution.png", tumor_type_plot, width = 10, height = 6)

tumor_size_plot <- ggplot(brain_tumor_data, aes(x = Tumor_Size)) +
  geom_histogram(binwidth = 0.5, fill = "darkgreen", color = "black") +
  labs(title = "Brain Tumor Size Distribution", x = "Tumor Size (cm)", y = "Count") +
  theme_minimal()
ggsave("plots/eda/tumor_size_distribution.png", tumor_size_plot, width = 8, height = 6)

location_plot <- ggplot(brain_tumor_data, aes(x = Location, fill = Location)) +
  geom_bar() +
  labs(title = "Brain Tumor Location Distribution", x = "Location", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")
ggsave("plots/eda/location_distribution.png", location_plot, width = 10, height = 6)

# ---- ANALYSIS 3: Treatment Analysis ----
brain_tumor_data$treatment_combo <- paste(
  ifelse(brain_tumor_data$Radiation_Treatment == "Yes", "Rad", "NoRad"),
  ifelse(brain_tumor_data$Surgery_Performed == "Yes", "Surg", "NoSurg"),
  ifelse(brain_tumor_data$Chemotherapy == "Yes", "Chemo", "NoChemo"),
  sep = "_"
)

treatment_combo_plot <- ggplot(brain_tumor_data, aes(x = treatment_combo, fill = treatment_combo)) +
  geom_bar() +
  labs(title = "Treatment Combinations for Brain Tumors", x = "Treatment Combination", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = FALSE)
ggsave("plots/eda/treatment_combinations.png", treatment_combo_plot, width = 10, height = 6)

treatment_by_type <- ggplot(brain_tumor_data, aes(x = Tumor_Type, fill = treatment_combo)) +
  geom_bar(position = "fill") +
  labs(title = "Treatment Approaches by Brain Tumor Type", x = "Tumor Type", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/eda/treatment_by_tumor_type.png", treatment_by_type, width = 10, height = 6)

# ---- ANALYSIS 4: Correlation Analysis ----
numeric_vars <- brain_tumor_data %>%
  select(Age, Tumor_Size, Survival_Rate, Tumor_Growth_Rate)
correlation_matrix <- cor(numeric_vars, use = "complete.obs")
write.csv(correlation_matrix, "output/eda/correlation_matrix.csv")

png("plots/eda/correlation_plot.png", width = 800, height = 800)
corrplot(correlation_matrix, method = "circle", type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45)
dev.off()