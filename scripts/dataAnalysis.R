# Load required libraries
pkgs <- c("tidyverse", "ggplot2", "viridis", "RColorBrewer", "gridExtra", "GGally")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cran.rstudio.com/")
} else {
  message("✅ All required packages are already installed.")
}
invisible(lapply(pkgs, library, character.only = TRUE))

# Create directories
if (!dir.exists("plots/dataAnalysis")) {
  dir.create("plots/dataAnalysis", recursive = TRUE)
}
if (!dir.exists("output/dataAnalysis")) {
  dir.create("output/dataAnalysis", recursive = TRUE)
}

# Redirect output to file
sink("output/dataAnalysis/analysis_output.txt")

# Set a vibrant color palette
vibrant_colors <- c("#FF5733", "#33FF57", "#3357FF", "#FF33F5", "#F5FF33", "#33FFF5")
pastel_colors <- c("#FFB6C1", "#98FB98", "#ADD8E6", "#FFA07A", "#FFFACD", "#E6E6FA")
bold_colors <- c("#FF0000", "#00FF00", "#0000FF", "#FF00FF", "#FFFF00", "#00FFFF")

# Set a seed for reproducibility
set.seed(42)

# Load the dataset
data <- read.csv("data/processed/brain_tumor_data_clean.csv")

# Print basic summary
cat("Original Data Summary:\n")
cat("Dimensions:", dim(data), "\n")
cat("First few rows:\n")
print(head(data))

# --- Data Preprocessing ---
cat("\nCount of NA values per column:\n")
print(colSums(is.na(data)))

# Handle missing values
if (any(is.na(data))) {
  numeric_cols <- c("Age", "Tumor_Size", "Survival_Rate", "Tumor_Growth_Rate")
  for (col in numeric_cols) {
    if (any(is.na(data[[col]]))) {
      data[[col]][is.na(data[[col]])] <- median(data[[col]], na.rm = TRUE)
      cat("Filled NA values in", col, "with median\n")
    }
  }
  categorical_cols <- setdiff(names(data), c("Patient_ID", numeric_cols))
  for (col in categorical_cols) {
    if (any(is.na(data[[col]]))) {
      mode_val <- names(sort(table(data[[col]]), decreasing = TRUE))[1]
      data[[col]][is.na(data[[col]])] <- mode_val
      cat("Filled NA values in", col, "with mode:", mode_val, "\n")
    }
  }
}

# Check Survival_Binary values
cat("\nUnique values in Survival_Binary:", unique(data$Survival_Binary), "\n")

# Convert categorical variables to factors
data <- data %>%
  mutate(
    Survival_Binary = case_when(
      Survival_Binary == "0" ~ "Poor",
      Survival_Binary == "1" ~ "Good",
      TRUE ~ as.character(Survival_Binary)
    ),
    Survival_Binary = factor(Survival_Binary),
    Gender = factor(Gender),
    Tumor_Type = factor(Tumor_Type),
    Location = factor(Location),
    Histology = factor(Histology),
    Radiation_Treatment = ifelse(Radiation_Treatment == "1", "Yes", "No"),
    Surgery_Performed = ifelse(Surgery_Performed == "1", "Yes", "No"),
    Chemotherapy = ifelse(Chemotherapy == "1", "Yes", "No"),
    Family_History = ifelse(Family_History == "1", "Yes", "No"),
    Follow_Up_Required = ifelse(Follow_Up_Required == "1", "Yes", "No"),
    MRI_Result = factor(MRI_Result),
    Symptom_1 = factor(Symptom_1),
    Symptom_2 = factor(Symptom_2),
    Symptom_3 = factor(Symptom_3)
  )

cat("\nSurvival_Binary levels after conversion:", levels(data$Survival_Binary), "\n")

# Remove Patient_ID
data <- data %>% select(-Patient_ID)

# Print summary after preprocessing
cat("\nData After Preprocessing:\n")
cat("Dimensions:", dim(data), "\n")
print(summary(data))

# Sample data for visualization
max_points <- 1000
if (nrow(data) > max_points) {
  cat("\nSampling", max_points, "rows for clearer visualization\n")
  data_sample <- data %>% sample_n(max_points)
} else {
  data_sample <- data
}

# --- Visualization Theme ---
create_vibrant_theme <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "#F8F8F8", color = NA),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      plot.title = element_text(face = "bold", size = 16, color = "#333333"),
      plot.subtitle = element_text(size = 12, color = "#666666"),
      axis.title = element_text(face = "bold", size = 12, color = "#333333"),
      axis.text = element_text(size = 10, color = "#333333"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#E0E0E0"),
      strip.text = element_text(face = "bold", color = "#333333")
    )
}

# --- Visualizations ---

# 1. Age vs Tumor Size by Survival Outcome
p1 <- ggplot(data_sample, aes(x = Age, y = Tumor_Size, color = Survival_Binary)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Poor" = "#FF5733", "Good" = "#3357FF")) +
  labs(
    title = "Age vs Tumor Size by Survival Outcome",
    subtitle = "How tumor size relates to patient age and survival",
    x = "Age (years)",
    y = "Tumor Size (cm)",
    color = "Survival Outcome"
  ) +
  create_vibrant_theme()
ggsave("plots/dataAnalysis/age_vs_tumor_size_survival.png", p1, width = 10, height = 7, bg = "white")

# 2. Survival Rate Distribution by Tumor Type
p2 <- ggplot(data_sample, aes(x = Tumor_Type, y = Survival_Rate, fill = Tumor_Type)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, alpha = 0.9, fill = "white") +
  scale_fill_manual(values = c("Benign" = "#33FF57", "Malignant" = "#FF33F5")) +
  labs(
    title = "Survival Rate Distribution by Tumor Type",
    subtitle = "Comparing survival rates between benign and malignant tumors",
    x = "Tumor Type",
    y = "Survival Rate (%)"
  ) +
  create_vibrant_theme()
ggsave("plots/dataAnalysis/survival_rate_by_tumor_type.png", p2, width = 10, height = 7, bg = "white")

# 3. Mean Survival Rate by Treatment Combination
treatment_data <- data %>%
  group_by(Surgery_Performed, Radiation_Treatment, Chemotherapy) %>%
  summarize(
    Mean_Survival = mean(Survival_Rate, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  ) %>%
  filter(Count >= 5)
treatment_data$Treatment_Combo <- paste(
  treatment_data$Surgery_Performed,
  treatment_data$Radiation_Treatment,
  treatment_data$Chemotherapy,
  sep = " + "
)
p3 <- ggplot(treatment_data, aes(x = reorder(Treatment_Combo, Mean_Survival), y = Mean_Survival, fill = Mean_Survival)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    title = "Mean Survival Rate by Treatment Combination",
    subtitle = "Comparing effectiveness of different treatment approaches",
    x = "Treatment Combination",
    y = "Mean Survival Rate (%)",
    fill = "Mean Survival"
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#F8F8F8", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  )
ggsave("plots/dataAnalysis/treatment_effectiveness.png", p3, width = 12, height = 7, bg = "white")

# 4. Tumor Growth Rate vs Survival Rate by Tumor Type
p4 <- ggplot(data_sample, aes(x = Tumor_Growth_Rate, y = Survival_Rate, color = Tumor_Type)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Benign" = "#33FF57", "Malignant" = "#FF33F5")) +
  labs(
    title = "Tumor Growth Rate vs Survival Rate",
    subtitle = "Impact of tumor growth rate on patient survival",
    x = "Tumor Growth Rate (cm/year)",
    y = "Survival Rate (%)",
    color = "Tumor Type"
  ) +
  create_vibrant_theme()
ggsave("plots/dataAnalysis/growth_vs_survival_by_type.png", p4, width = 10, height = 7, bg = "white")

# 5. Tumor Size Distribution by Location
p5 <- ggplot(data_sample, aes(x = reorder(Location, Tumor_Size, FUN = median), y = Tumor_Size, fill = Location)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Tumor Size Distribution by Location",
    subtitle = "Comparing tumor sizes across different brain regions",
    x = "Tumor Location",
    y = "Tumor Size (cm)"
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#F8F8F8", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  )
ggsave("plots/dataAnalysis/tumor_size_by_location.png", p5, width = 12, height = 7, bg = "white")

# 6. Age Distribution by Gender and Survival Outcome
p6 <- ggplot(data_sample, aes(x = Age, fill = Survival_Binary)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~Gender) +
  scale_fill_manual(values = c("Poor" = "#FF5733", "Good" = "#3357FF")) +
  labs(
    title = "Age Distribution by Gender and Survival Outcome",
    subtitle = "How age and gender relate to survival",
    x = "Age (years)",
    y = "Density",
    fill = "Survival Outcome"
  ) +
  create_vibrant_theme()
ggsave("plots/dataAnalysis/age_distribution_by_gender_survival.png", p6, width = 10, height = 7, bg = "white")

# 7. Predicted Survival Rate by Age and Tumor Size
tryCatch({
  model <- lm(Survival_Rate ~ Age + Tumor_Size + Tumor_Growth_Rate + Tumor_Type, data = data)
  age_range <- seq(min(data$Age), max(data$Age), length.out = 20)
  size_range <- seq(min(data$Tumor_Size), max(data$Tumor_Size), length.out = 20)
  prediction_grid <- expand.grid(
    Age = age_range,
    Tumor_Size = size_range,
    Tumor_Growth_Rate = median(data$Tumor_Growth_Rate),
    Tumor_Type = factor("Malignant", levels = levels(data$Tumor_Type))
  )
  prediction_grid$Predicted_Survival <- predict(model, newdata = prediction_grid)
  p7 <- ggplot(prediction_grid, aes(x = Age, y = Tumor_Size, fill = Predicted_Survival)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", name = "Predicted\nSurvival Rate") +
    labs(
      title = "Predicted Survival Rate by Age and Tumor Size",
      subtitle = "Based on linear regression model",
      x = "Age (years)",
      y = "Tumor Size (cm)"
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "#F8F8F8", color = NA),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12)
    )
  ggsave("plots/dataAnalysis/survival_prediction_heatmap.png", p7, width = 10, height = 8, bg = "white")
  cat("\nSurvival Prediction Model Summary:\n")
  print(summary(model))
}, error = function(e) {
  cat("Could not create survival prediction model:", e$message, "\n")
})

# 8. Mean Survival Rate by Symptom and Tumor Type
symptom_counts <- c(table(data$Symptom_1), table(data$Symptom_2), table(data$Symptom_3))
symptom_counts <- sort(symptom_counts, decreasing = TRUE)
top_symptoms <- names(symptom_counts)[1:min(5, length(symptom_counts))]
symptom_data <- data.frame(
  Symptom = c(as.character(data$Symptom_1), as.character(data$Symptom_2), as.character(data$Symptom_3)),
  Tumor_Type = rep(as.character(data$Tumor_Type), 3),
  Survival_Rate = rep(data$Survival_Rate, 3)
)
symptom_data <- symptom_data %>%
  filter(Symptom %in% top_symptoms & !is.na(Symptom))
symptom_summary <- symptom_data %>%
  group_by(Symptom, Tumor_Type) %>%
  summarize(
    Mean_Survival = mean(Survival_Rate, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  ) %>%
  filter(Count >= 5)
p8 <- ggplot(symptom_summary, aes(x = reorder(Symptom, Mean_Survival), y = Mean_Survival, fill = Tumor_Type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("Benign" = "#33FF57", "Malignant" = "#FF33F5")) +
  labs(
    title = "Mean Survival Rate by Symptom and Tumor Type",
    subtitle = "Impact of common symptoms on survival",
    x = "Symptom",
    y = "Mean Survival Rate (%)",
    fill = "Tumor Type"
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#F8F8F8", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12)
  )
ggsave("plots/dataAnalysis/symptom_impact_analysis.png", p8, width = 12, height = 7, bg = "white")

# 9. Numeric Variables Relationships
tryCatch({
  pair_data <- data_sample %>%
    select(Age, Tumor_Size, Survival_Rate, Tumor_Growth_Rate, Survival_Binary)
  p9 <- ggpairs(pair_data,
                columns = 1:4,
                aes(color = Survival_Binary),
                upper = list(continuous = "density", combo = "box_no_facet"),
                lower = list(continuous = "points", combo = "dot_no_facet"),
                diag = list(continuous = "densityDiag"),
                axisLabels = "show") +
    scale_color_manual(values = c("Poor" = "#FF5733", "Good" = "#3357FF")) +
    scale_fill_manual(values = c("Poor" = "#FF5733", "Good" = "#3357FF")) +
    labs(title = "Relationships Between Key Numeric Variables",
         subtitle = "Colored by survival outcome") +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "#F8F8F8", color = NA),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12)
    )
  ggsave("plots/dataAnalysis/numeric_variables_relationships.png", p9, width = 12, height = 10, bg = "white")
}, error = function(e) {
  cat("Could not create pairplot:", e$message, "\n")
})

# --- Key Correlations ---
cat("\n--- KEY CORRELATIONS ---\n")
safe_cor <- function(x, y, group = NULL) {
  tryCatch({
    if (is.null(group)) {
      result <- cor(x, y, use = "pairwise.complete.obs")
    } else {
      result <- cor(x[group], y[group], use = "pairwise.complete.obs")
    }
    return(result)
  }, error = function(e) {
    return(NA)
  })
}
cat("Age vs Tumor Size: ", safe_cor(data$Age, data$Tumor_Size), "\n")
cat("Age vs Survival Rate: ", safe_cor(data$Age, data$Survival_Rate), "\n")
cat("Tumor Size vs Survival Rate: ", safe_cor(data$Tumor_Size, data$Survival_Rate), "\n")
cat("Tumor Growth Rate vs Survival Rate: ", safe_cor(data$Tumor_Growth_Rate, data$Survival_Rate), "\n")
if (sum(data$Survival_Binary == "Poor", na.rm = TRUE) > 5) {
  poor_indices <- which(data$Survival_Binary == "Poor")
  cat("\nCorrelations for Poor Survival group:\n")
  cat("Age vs Tumor Size: ", safe_cor(data$Age, data$Tumor_Size, poor_indices), "\n")
  cat("Tumor Size vs Survival Rate: ", safe_cor(data$Tumor_Size, data$Survival_Rate, poor_indices), "\n")
}
if (sum(data$Survival_Binary == "Good", na.rm = TRUE) > 5) {
  good_indices <- which(data$Survival_Binary == "Good")
  cat("\nCorrelations for Good Survival group:\n")
  cat("Age vs Tumor Size: ", safe_cor(data$Age, data$Tumor_Size, good_indices), "\n")
  cat("Tumor Size vs Survival Rate: ", safe_cor(data$Tumor_Size, data$Survival_Rate, good_indices), "\n")
}

cat("\nAll visualizations saved to plots/dataAnalysis directory.\n")

# Close the sink
sink()