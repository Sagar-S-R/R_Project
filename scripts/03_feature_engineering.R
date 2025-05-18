# Brain_Tumor_Analysis_Project/scripts/03_feature_engineering.R

# Load cleaned data
brain_tumor_data <- read.csv("data/processed/brain_tumor_data_clean.csv", stringsAsFactors = TRUE)

# Create age groups
brain_tumor_data$Age_Group <- cut(brain_tumor_data$Age, 
                                 breaks = c(0, 20, 40, 60, 80, 100),
                                 labels = c("0-20", "21-40", "41-60", "61-80", "81+"),
                                 include.lowest = TRUE)

# Create tumor size categories
brain_tumor_data$Size_Category <- cut(brain_tumor_data$Tumor_Size,
                                    breaks = quantile(brain_tumor_data$Tumor_Size, probs = seq(0, 1, 0.25)),
                                    labels = c("Small", "Medium", "Large", "Very Large"),
                                    include.lowest = TRUE)

# Create treatment intensity score (0-3)
brain_tumor_data$Treatment_Intensity <- (brain_tumor_data$Radiation_Treatment == "Yes") +
                                       (brain_tumor_data$Surgery_Performed == "Yes") +
                                       (brain_tumor_data$Chemotherapy == "Yes")

# Create symptom count
brain_tumor_data$Symptom_Count <- (brain_tumor_data$Symptom_1 != "None") +
                                 (brain_tumor_data$Symptom_2 != "None") +
                                 (brain_tumor_data$Symptom_3 != "None")

# Create risk score based on multiple factors
brain_tumor_data$Risk_Score <- scale(brain_tumor_data$Tumor_Size) + 
                              scale(brain_tumor_data$Tumor_Growth_Rate) +
                              (brain_tumor_data$Family_History == "Yes") - 
                              scale(brain_tumor_data$Age) / 2

# Save the feature-engineered dataset
write.csv(brain_tumor_data, "data/processed/brain_tumor_data_engineered.csv", row.names = FALSE)