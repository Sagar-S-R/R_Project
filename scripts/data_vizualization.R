# Load required libraries
# List of required packages
pkgs <- c("tidyverse", "ggplot2", "viridis", "RColorBrewer", "gridExtra", "GGally")

# Install only the packages that are not already installed
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cran.rstudio.com/")
} else {
  message("✅ All required packages are already installed.")
}

# Load all required packages
invisible(lapply(pkgs, library, character.only = TRUE))

# Create the plots directory if it doesn't exist
if (!dir.exists("plots/dataAnalysis")) {
  dir.create("plots/dataAnalysis", recursive = TRUE)
}

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

# Check for NA values before attempting to filter
cat("\nCount of NA values per column:\n")
print(colSums(is.na(data)))

# Handle missing values properly - instead of removing rows, let's fill NA values
if(any(is.na(data))) {
  # For numeric columns, fill with median
  numeric_cols <- c("Age", "Tumor_Size", "Survival_Rate", "Tumor_Growth_Rate")
  for(col in numeric_cols) {
    if(any(is.na(data[[col]]))) {
      data[[col]][is.na(data[[col]])] <- median(data[[col]], na.rm = TRUE)
      cat("Filled NA values in", col, "with median\n")
    }
  }
  
  # For categorical columns, fill with mode
  categorical_cols <- setdiff(names(data), c("Patient_ID", numeric_cols))
  for(col in categorical_cols) {
    if(any(is.na(data[[col]]))) {
      mode_val <- names(sort(table(data[[col]]), decreasing = TRUE))[1]
      data[[col]][is.na(data[[col]])] <- mode_val
      cat("Filled NA values in", col, "with mode:", mode_val, "\n")
    }
  }
}

# Check Survival_Binary values before conversion
cat("\nUnique values in Survival_Binary:", unique(data$Survival_Binary), "\n")

# Convert categorical variables to factors with appropriate levels
data <- data %>% 
  mutate(
    # Try to handle the case where Survival_Binary might have unexpected values
    Survival_Binary = case_when(
      Survival_Binary == "0" ~ "Poor",
      Survival_Binary == "1" ~ "Good",
      TRUE ~ as.character(Survival_Binary)  # Keep other values as is
    ),
    Survival_Binary = factor(Survival_Binary),
    
    Gender = factor(Gender),
    Tumor_Type = factor(Tumor_Type),
    Location = factor(Location),
    Histology = factor(Histology),
    
    # Convert other binary variables to Yes/No format
    Radiation_Treatment = ifelse(Radiation_Treatment == "1", "Yes", "No"),
    Surgery_Performed = ifelse(Surgery_Performed == "1", "Yes", "No"),
    Chemotherapy = ifelse(Chemotherapy == "1", "Yes", "No"),
    Family_History = ifelse(Family_History == "1", "Yes", "No"),
    Follow_Up_Required = ifelse(Follow_Up_Required == "1", "Yes", "No"),
    
    # Convert MRI_Result to factor
    MRI_Result = factor(MRI_Result)
  )

# Print updated factor levels
cat("\nSurvival_Binary levels after conversion:", levels(data$Survival_Binary), "\n")

# Remove Patient_ID
data <- data %>% select(-Patient_ID)

# Print summary after preprocessing
cat("\nData After Preprocessing:\n")
cat("Dimensions:", dim(data), "\n")
print(summary(data))

# --- Ensure a manageable dataset size for visualization ---

# If the dataset is too large, sample a subset
max_points <- 1000  # Maximum number of points for clear visualization
if(nrow(data) > max_points) {
  cat("\nSampling", max_points, "rows for clearer visualization\n")
  data_sample <- data %>% sample_n(max_points)
} else {
  data_sample <- data
}

# --- Create Enhanced Visualizations ---

# Function to create a consistent theme for all plots
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

# 1. Age vs Tumor Size by Survival Outcome - Improved Scatter Plot
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

# 2. Survival Rate Distribution by Tumor Type - Violin Plot
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

# 3. Treatment Effectiveness - Bar Chart
treatment_data <- data %>%
  group_by(Surgery_Performed, Radiation_Treatment, Chemotherapy) %>%
  summarize(
    Mean_Survival = mean(Survival_Rate, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  ) %>%
  filter(Count >= 5)  # Filter out small groups

treatment_data$Treatment_Combo <- paste(
  treatment_data$Surgery_Performed,
  treatment_data$Radiation_Treatment,
  treatment_data$Chemotherapy,
  sep = " + "
)

p3 <- ggplot(treatment_data, aes(x = reorder(Treatment_Combo, Mean_Survival), 
                                y = Mean_Survival, 
                                fill = Mean_Survival)) +
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

# 4. Tumor Growth Rate vs Survival Rate by Tumor Type - Scatter Plot with Density
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

# 5. Tumor Size Distribution by Location - Box Plot
p5 <- ggplot(data_sample, aes(x = reorder(Location, Tumor_Size, FUN = median), 
                             y = Tumor_Size, 
                             fill = Location)) +
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

# 6. Age Distribution by Gender with Survival Outcome - Density Plot
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

# 7. Treatment Effect on Survival by Tumor Type - Bar Chart
# Create a combined treatment variable
data$Treatment <- paste(
  ifelse(data$Surgery_Performed == "Yes", "Surgery", "No Surgery"),
  ifelse(data$Radiation_Treatment == "Yes", "Radiation", "No Radiation"),
  ifelse(data$Chemotherapy == "Yes", "Chemo", "No Chemo"),
  sep = " + "
)

treatment_effect <- data %>%
  group_by(Treatment, Tumor_Type) %>%
  summarize(
    Mean_Survival = mean(Survival_Rate, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  ) %>%
  filter(Count >= 5)  # Filter out rare combinations

p7 <- ggplot(treatment_effect, aes(x = reorder(Treatment, Mean_Survival), 
                                  y = Mean_Survival, 
                                  fill = Tumor_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Benign" = "#33FF57", "Malignant" = "#FF33F5")) +
  labs(
    title = "Mean Survival Rate by Treatment and Tumor Type",
    subtitle = "Comparing treatment effectiveness across tumor types",
    x = "Treatment Combination",
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

ggsave("plots/dataAnalysis/treatment_effect_by_tumor_type.png", p7, width = 14, height = 7, bg = "white")

# 8. Survival Prediction Heatmap
# Create a basic linear model for survival prediction
tryCatch({
  model <- lm(Survival_Rate ~ Age + Tumor_Size + Tumor_Growth_Rate + Tumor_Type, data = data)
  
  # Create prediction grid
  age_range <- seq(min(data$Age), max(data$Age), length.out = 20)
  size_range <- seq(min(data$Tumor_Size), max(data$Tumor_Size), length.out = 20)
  
  prediction_grid <- expand.grid(
    Age = age_range,
    Tumor_Size = size_range,
    Tumor_Growth_Rate = median(data$Tumor_Growth_Rate),
    Tumor_Type = factor("Malignant", levels = levels(data$Tumor_Type))
  )
  
  # Add predictions
  prediction_grid$Predicted_Survival <- predict(model, newdata = prediction_grid)
  
  # Create heatmap
  p8 <- ggplot(prediction_grid, aes(x = Age, y = Tumor_Size, fill = Predicted_Survival)) +
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
  
  ggsave("plots/dataAnalysis/survival_prediction_heatmap.png", p8, width = 10, height = 8, bg = "white")
  
  # Print model summary
  cat("\nSurvival Prediction Model Summary:\n")
  print(summary(model))
  
}, error = function(e) {
  cat("Could not create survival prediction model:", e$message, "\n")
})

# 9. Symptom Analysis
# Extract the top 5 most common symptoms
symptom_counts <- c(
  table(data$Symptom_1),
  table(data$Symptom_2),
  table(data$Symptom_3)
)
symptom_counts <- sort(symptom_counts, decreasing = TRUE)
top_symptoms <- names(symptom_counts)[1:min(5, length(symptom_counts))]

# Create a dataset for symptom analysis
symptom_data <- data.frame(
  Symptom = c(as.character(data$Symptom_1), 
              as.character(data$Symptom_2), 
              as.character(data$Symptom_3)),
  Tumor_Type = rep(as.character(data$Tumor_Type), 3),
  Survival_Rate = rep(data$Survival_Rate, 3)
)

# Filter for top symptoms
symptom_data <- symptom_data %>%
  filter(Symptom %in% top_symptoms & !is.na(Symptom))

# Calculate average survival rate by symptom and tumor type
symptom_summary <- symptom_data %>%
  group_by(Symptom, Tumor_Type) %>%
  summarize(
    Mean_Survival = mean(Survival_Rate, na.rm = TRUE),
    Count = n(),
    .groups = "drop"
  ) %>%
  filter(Count >= 5)  # Filter out rare combinations

# Create plot
p9 <- ggplot(symptom_summary, aes(x = reorder(Symptom, Mean_Survival), 
                                y = Mean_Survival, 
                                fill = Tumor_Type)) +
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

ggsave("plots/dataAnalysis/symptom_impact_analysis.png", p9, width = 12, height = 7, bg = "white")

# 10. Advanced Pair Plot for Numeric Variables
# Try to create a pairplot for better variable relationships
tryCatch({
  # Select numeric variables and add Survival_Binary for color
  pair_data <- data_sample %>%
    select(Age, Tumor_Size, Survival_Rate, Tumor_Growth_Rate, Survival_Binary)
  
  # Create pairplot
  p10 <- ggpairs(pair_data,
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
  
  ggsave("plots/dataAnalysis/numeric_variables_relationships.png", p10, width = 12, height = 10, bg = "white")
}, error = function(e) {
  cat("Could not create pairplot:", e$message, "\n")
})

# Optional: Create a multi-panel dashboard of selected plots
# Try to combine 4 of the most informative plots
tryCatch({
  # Create a 2x2 grid of selected plots
  combined_plot <- grid.arrange(
    p1, p2, p4, p5,
    ncol = 2,
    top = grid::textGrob("Brain Tumor Analysis Dashboard", 
                         gp = grid::gpar(fontsize = 20, fontface = "bold"))
  )
  
  ggsave("plots/dataAnalysis/brain_tumor_dashboard.png", combined_plot, width = 16, height = 12, bg = "white")
}, error = function(e) {
  cat("Could not create dashboard:", e$message, "\n")
})

# Print key correlations - with better error handling
cat("\n--- KEY CORRELATIONS ---\n")

# Safely calculate correlations
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

# Overall correlations
cat("Age vs Tumor Size: ", safe_cor(data$Age, data$Tumor_Size), "\n")
cat("Age vs Survival Rate: ", safe_cor(data$Age, data$Survival_Rate), "\n")
cat("Tumor Size vs Survival Rate: ", safe_cor(data$Tumor_Size, data$Survival_Rate), "\n")
cat("Tumor Growth Rate vs Survival Rate: ", safe_cor(data$Tumor_Growth_Rate, data$Survival_Rate), "\n")

# Correlations by groups - only if we have data for both groups
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