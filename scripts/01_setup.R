required_packages <- c("tidyverse", "caret", "corrplot", "rpart", "rpart.plot", 
                      "randomForest", "e1071", "glmnet", "gbm", "survival", 
                      "survminer", "cluster", "factoextra", "skimr", "MASS",
                      "pROC", "mice", "pROC")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
  library(pkg, character.only = TRUE)
}

# Create project directories
dirs <- c("data/raw", "data/processed", "scripts", 
         "plots/part1", "plots/part2", "plots/part3", "plots/part4",
         "output/part1", "output/part2", "output/part3", "output/part4")
for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}