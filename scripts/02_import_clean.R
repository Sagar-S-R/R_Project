# Brain_Tumor_Analysis_Project/scripts/02_import_clean.R

library(mice)
# Import the dataset
brain_tumor_data <- read.csv("data/raw/brain_tumor.csv")

# Check for missing values
missing_values <- colSums(is.na(brain_tumor_data))
print(missing_values)

# Handle missing values
# Function to get mode of a vector
get_mode <- function(x) {
  unique_x <- unique(x[!is.na(x)])
  if (length(unique_x) == 0) return(NA)
  return(unique_x[which.max(tabulate(match(x, unique_x)))])
}

# Identify numeric and categorical columns
numeric_cols <- sapply(brain_tumor_data, is.numeric)
categorical_cols <- !numeric_cols & colnames(brain_tumor_data) != "Patient_ID"

# Impute missing values in numeric columns using MICE
if (sum(sapply(brain_tumor_data[,numeric_cols], function(x) sum(is.na(x)))) > 0) {
  imp_model <- mice(brain_tumor_data[,numeric_cols], m = 5, method = "pmm")
  brain_tumor_data[,numeric_cols] <- complete(imp_model)
}

# Impute categorical columns with mode
for (col in colnames(brain_tumor_data)[categorical_cols]) {
  if (sum(is.na(brain_tumor_data[[col]])) > 0) {
    brain_tumor_data[[col]][is.na(brain_tumor_data[[col]])] <- get_mode(brain_tumor_data[[col]])
  }
}

# Convert categorical variables to factors
factor_columns <- c("Gender", "Tumor_Type", "Location", "Histology", "Stage",
                   "Symptom_1", "Symptom_2", "Symptom_3", 
                   "Radiation_Treatment", "Surgery_Performed", "Chemotherapy",
                   "Family_History", "MRI_Result", "Follow_Up_Required")

for (col in factor_columns) {
  if (col %in% colnames(brain_tumor_data)) {
    brain_tumor_data[[col]] <- as.factor(brain_tumor_data[[col]])
  }
}

# Create binary target variable for modeling
brain_tumor_data$Survival_Binary <- ifelse(brain_tumor_data$Survival_Rate >= 50, 1, 0)
brain_tumor_data$Survival_Binary <- factor(brain_tumor_data$Survival_Binary, 
                                         levels = c(0, 1), 
                                         labels = c("Poor", "Good"))

# Save cleaned data
write.csv(brain_tumor_data, "data/processed/brain_tumor_data_clean.csv", row.names = FALSE)