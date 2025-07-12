# Brain Tumor Analysis Project

This project provides a comprehensive analysis of a brain tumor dataset, including data cleaning, feature engineering, exploratory data analysis (EDA), statistical testing, advanced modeling, and survival analysis. The workflow is implemented in R and produces both visualizations and detailed statistical reports.

## Project Structure

```
.
├── .gitignore
├── README.md
├── data/
│   ├── raw/ 
│   │   ├── brain_tumor.csv
│   │   ├── brain_tumor_old.csv # Raw input data 
│   └── processed/              # Cleaned and feature-engineered datasets
│       ├── brain_tumor_data_clean.csv
│       └── brain_tumor_data_engineered.csv
├── output/
│   ├── advanced_analyses/    # Statistical analysis reports
│   ├── advanced_modeling/    # Advanced modeling results
│   ├── dataAnalysis/         # EDA and basic analysis outputs
│   ├── eda/                  # EDA text outputs
│   └── survival_modeling/    # Survival modeling results
├── plots/
│   ├── advanced_analyses/    # Plots from advanced analyses
│   ├── advanced_modeling/    # Plots from advanced modeling
│   ├── dataAnalysis/         # EDA plots
│   ├── eda/                  # EDA plots
│   └── survival_modeling/    # Survival modeling plots
└── scripts/
    ├── 01_setup.R                    # Install packages and setup directories
    ├── 02_import_clean.R             # Data import and cleaning
    ├── 03_feature_engineering.R      # Feature engineering
    ├── brain_tumor_eda.R             # Exploratory data analysis
    ├── dataAnalysis.R                # Data analysis and visualization
    ├── brain_tumor_advanced_analyses.R   # Advanced statistical analyses
    ├── brain_tumor_advanced_modeling.R   # Advanced machine learning models
    └── brain_tumor_survival_modeling.R   # Survival analysis and modeling
```

## Workflow Overview

1. **Setup**  
   Run [`scripts/01_setup.R`](scripts/01_setup.R) to install required R packages and create necessary directories.

2. **Data Import & Cleaning**  
   Use [`scripts/02_import_clean.R`](scripts/02_import_clean.R) to import raw data, handle missing values, and save a cleaned dataset.

3. **Feature Engineering**  
   Run [`scripts/03_feature_engineering.R`](scripts/03_feature_engineering.R) to create new features and save the engineered dataset.

4. **Exploratory Data Analysis (EDA)**  
   - [`scripts/brain_tumor_eda.R`](scripts/brain_tumor_eda.R): Generates summary statistics, visualizations, and statistical tests.
   - Outputs: `output/eda/eda_results.txt`, plots in `plots/eda/`.

5. **Data Analysis & Visualization**  
   - [`scripts/dataAnalysis.R`](scripts/dataAnalysis.R): Produces advanced visualizations and correlation analysis.
   - Outputs: `output/dataAnalysis/analysis_output.txt`, plots in `plots/dataAnalysis/`.

6. **Advanced Statistical Analyses**  
   - [`scripts/brain_tumor_advanced_analyses.R`](scripts/brain_tumor_advanced_analyses.R): Performs detailed descriptive and group statistics.
   - Outputs: `output/advanced_analyses/statistical_analysis_report.txt`, plots in `plots/advanced_analyses/`.

7. **Advanced Modeling**  
   - [`scripts/brain_tumor_advanced_modeling.R`](scripts/brain_tumor_advanced_modeling.R): Random Forest, linear regression, and feature importance.
   - Outputs: `output/advanced_modeling/`, plots in `plots/advanced_modeling/`.

8. **Survival Modeling**  
   - [`scripts/brain_tumor_survival_modeling.R`](scripts/brain_tumor_survival_modeling.R): Linear/logistic regression, decision trees, ROC analysis.
   - Outputs: `output/survival_modeling/`, plots in `plots/survival_modeling/`.

## How to Run

1. Open the project in RStudio or VS Code.
2. Run the scripts in order as described above.
3. Check the `output/` and `plots/` directories for results and visualizations.

## Requirements

- R (>= 4.0)
- R packages: `tidyverse`, `caret`, `randomForest`, `corrplot`, `rpart`, `rpart.plot`, `pROC`, `mice`, `skimr`, `GGally`, `ROSE`, `car`, `moments`, `ggpubr`, `gridExtra`, `viridis`, `RColorBrewer`, `glmnet`, `e1071`, `gbm`, `survival`, `survminer`, `cluster`, `factoextra`, `MASS`

All required packages are installed by [`scripts/01_setup.R`](scripts/01_setup.R).

## Outputs

- **Text Reports:** Statistical summaries and model results in `output/`
- **Plots:** All visualizations in `plots/` subfolders

## License

This project is for academic and research purposes.

---

*For questions or contributions, please open an issue or submit a pull request.*