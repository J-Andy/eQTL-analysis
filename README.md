# Gene Knockout and Mutation Analysis

## Overview

This project analyzes the association between specific cancer mutations and gene knockout experiments using Welch's t-test in R. The analysis is encapsulated in a Singularity container, which includes all necessary dependencies and R scripts.

## Comment: Extending the Analysis with Advanced Models
While Welch’s t-test is suitable for this toy example, more advanced methods would be necessary for larger datasets or when controlling for additional covariates. If we had more data or multiple covariates, I would be using linear regression or ANOVA. Linear regression would help model the relationship between mutation status and knockout outcomes while accounting for factors such as tissue type or growth conditions.
For efficiently handling large-scale data with many gene-mutation pairs, I would use the MatrixEQTL package in R (https://github.com/andreyshabalin/MatrixEQTL). MatrixEQTL allows the use of linear models and is optimized for large-scale genomic data, automatically performs False Discovery Rate (FDR) correction and has some model customization options.


## Singularity Container

The container runs the R script `main.R` which:
- Loads mutation data (`Mutations.tsv`) and gene knockout data (`Gene_KOs.tsv`). 
- Performs Welch’s t-tests for each mutation-gene pair.
- Outputs results with False Discovery Rate (FDR) correction.
- Generates Q-Q plots, histograms, and volcano plots.
- Saves significant results and all plots in the output directory.

## Installation Instructions

### 1. Clone the repository

Clone the repository containing the `Singularity` definition file, R scripts, and data files:

```bash
git clone https://github.com/J-Andy/eQTL-analysis.git
cd eQTL-analysis
```

### 2. Build the Singularity image
You can build the Singularity image locally using the following command:
```bash
singularity build r-analysis.sif Singularity
```

If you don't have root access on the system, like me, you can use the remote build service:
```bash
singularity build --remote r-analysis.sif Singularity
```
### 3. Running the Analysis
After building the container, run the analysis using this script: 

```./run_analysis.sh```

**Make sure the `run_analysis.sh` script is executable**:
```bash
  chmod +x run_analysis.sh
```


The shell script will:
 - Create an output directory in your current directory.
 - Run the Singularity container and bind the local output directory to the container’s /usr/src/app/output.

### 4. Running Unit Tests
To run the unit tests, use ```run_tests.sh ``` script.

***Unit Test Details***

1. Data Alignment Test:
Ensures that the mutation data (Mutations.tsv) and gene knockout data (Gene_KOs.tsv) are properly aligned by common models (cell lines).
The test verifies that the aligned datasets share the same number of common cell lines.

2. T-test Validation:
Checks that Welch’s t-tests return valid p-values when comparing mutation presence with gene knockout outcomes.
The test simulates mutation and fold change data and ensures the p-value is numeric and less than 1.
Test Output:

If all tests pass, the script will print a success message: "All tests passed successfully!"

### 5. Outputs
The results and plots generated by the analysis are saved in the output directory. These include:

- significant_mutation_gene_associations.tsv: A table of significant mutation-gene associations after FDR correction.
- qqplot_pvalues.png: A Q-Q plot for visualizing the distribution of p-values.
- histogram_pvalues.png: A histogram displaying the p-value distribution.
- volcano_plot.png: A volcano plot of mutation-gene associations.
- Boxplots: Individual PNG files for each significant mutation-gene association.


### Download the Prebuilt Singularity Image
If you prefer not to build the container yourself, you can download the prebuilt Singularity image (`r-analysis.sif`) from the following link: 

https://drive.google.com/file/d/17zTK5rTlSquPqw5A8xGvv1kwUP1szQDB/view?usp=sharing 

---

###

---


# Machine Learning Models for Predicting CRISPR Knockout Effects

### Classical ML Approach

For each gene:
1. Define the features as the binary mutations in the cell lines.
2. Define the target: the target is a continuous variable "knock-out effect".
3. Train a separate Random Forest regression model for each gene: each model is trained on the same mutation data (features), but the target variable changes based on the gene whose knock-out effect we are trying to predict.
4. I would then test several other models that people previously used in the literature, like Elastic Net Regression and Gradient Boosting Machines. 

### Multi-Output Models
In the current approach, I train separate models for each gene. However, there is an opportunity to take advantage of the correlation between gene knockouts by predicting multiple gene knockout effects simultaneously using a Multi-Output Random Forest or indeed multi-task learning with Convolutional Neural Networks (becasue neural network can be set up to have multiple output neurons, one for each target variable).

### Multi-Task Linear Regression
An even more advanced approach would be to apply a Bayesian multi-task linear regression model or multi-task learning with CNNs to solve this problem. For example, see here: https://www.sciencedirect.com/science/article/abs/pii/S1046202323002128?via%3Dihub

### Implementation
I implemented a simple RF model, for a single gene, in the script ```random_forest_model_for_gene_knockout_prediction_using_mutation_data.py```
You can also run it on Colab: https://colab.research.google.com/drive/1aIjK0wWs3UwiVI_bn-OcuUrHaQCpnZZo?usp=sharing 


# Generalization to Unseen Data

## Initial Model Build
The generalization of the model will primarily depend the diversity and size of training data and regularization techniques (i.e. feature selection) we used. We can assess the initial generalizability with k-fold cross-validation and bootstrapping.

***Pan-Cancer vs. Tissue-Specific Models***
In a broader context, models could be trained in two ways:

•	Pan-cancer models pool data from multiple cancer types, treating the knockouts as a universal phenomenon across all cell lines. This assumes mutation-response relationships are similar across different tissue types.

•	Tissue-specific models account for lineage-specific effects by training separate models for individual cancer types. 

## Model in Production
Once the models are deployed to production two issues are likely to occur **data drift** (changing mutation patterns) and **concept drift** (shifting biological relationships/interpretation).

**1. Data Drift:**
This will occur when the distribution of mutations in the new cell lines is significantly different from the training data. For instance, if the model was trained on cell lines that mostly contained mutations in specific genes, but new cell lines have mutations in entirely different genes or exhibit new mutation patterns, the model's performance will degrade.

***Suggested Solution:*** To detect data drift, we could monitor the distribution of mutations (input features) in new cell lines and compare them with the original training data. We can use statistical tests like Kolmogorov-Smirnov test or domain classifiers (which distinguish between old and new data) to detect these changes. A significant drift will require model retraining.


**2. Concept Drift:**
This I believe will be less of a concern. Concept drift refers to a change in the underlying relationship between the input features and the target variable over time. In our context this could happen when we introduce new experimental conditions, such as different cell environments or evolve resistant cell lines.

***Suggested Solution:*** I think detecting this type of drift in our context would be difficult and would have to be informed by wet-lab scientists (are they doing anything different?). We could monitor the model’s predictions and their accuracy on new cell lines. If the KO effect predictions are consistently inaccurate for new data even after retraining that to me would indicate concept drift. 

