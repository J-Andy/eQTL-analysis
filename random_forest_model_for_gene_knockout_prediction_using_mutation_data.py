# -*- coding: utf-8 -*-
"""Random Forest Model for Gene Knockout Prediction Using Mutation Data.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1aIjK0wWs3UwiVI_bn-OcuUrHaQCpnZZo
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import seaborn as sns

mutations_url = "https://raw.githubusercontent.com/J-Andy/eQTL-analysis/main/Mutations.tsv"
gene_kos_url = "https://raw.githubusercontent.com/J-Andy/eQTL-analysis/main/Gene_KOs.tsv"

# Load the mutation and gene KO datasets
mutations = pd.read_csv(mutations_url, sep='\t', index_col=0)
gene_kos = pd.read_csv(gene_kos_url, sep='\t', index_col=0)
print(mutations.head())
print(gene_kos.head())

# Transpose the mutations dataset to align with gene_kos
mutations_t = mutations.T
# Merge the datasets on index
combined_data = mutations_t.merge(gene_kos, left_index=True, right_index=True)

# Separate features (mutations) and target (GeneA_KO)
features = combined_data.iloc[:, :mutations_t.shape[1]]
target = combined_data['GeneA_KO']

# Split the data into training and test sets
train_features, test_features, train_target, test_target = train_test_split(
    features, target, test_size=0.2, random_state=321)

# Train the Random Forest model
rf_model = RandomForestRegressor(n_estimators=500, max_features=5, random_state=321)
rf_model.fit(train_features, train_target)

# Make predictions on the test set
predictions = rf_model.predict(test_features)

# Calculate Mean Squared Error (MSE)
mse = mean_squared_error(test_target, predictions)

# Calculate R-squared
rsq = r2_score(test_target, predictions)

# Print the performance metrics
print(f"Mean Squared Error (MSE): {mse}")
print(f"R-squared (R²): {rsq}")

# Feature importance
importances = rf_model.feature_importances_
feature_names = features.columns

# Create a DataFrame for feature importance
feat_importances = pd.DataFrame({'Feature': feature_names, 'Importance': importances})
feat_importances = feat_importances.sort_values(by='Importance', ascending=False)

# Plot feature importance
plt.figure(figsize=(10, 8))
sns.barplot(x='Importance', y='Feature', data=feat_importances)
plt.title('Feature Importance')
plt.show()

# Set up cross-validation with 10 folds
cv_scores = cross_val_score(rf_model, train_features, train_target, cv=10, scoring='neg_mean_squared_error')

# Print cross-validated results
print(f"Cross-validated MSE: {np.mean(-cv_scores)}")

# Hyperparameter tuning using GridSearchCV (to tune the mtry equivalent in Random Forest)
param_grid = {
    'max_features': [3, 5, 7, 9],  # Equivalent to R's mtry
    'n_estimators': [200, 500]     # Number of trees
}

# Grid Search with cross-validation
grid_search = GridSearchCV(estimator=rf_model, param_grid=param_grid,
                           cv=10, scoring='neg_mean_squared_error', verbose=1)

grid_search.fit(train_features, train_target)

# Best mtry (max_features in Python)
best_mtry = grid_search.best_params_['max_features']

# Print the best hyperparameters
print(f"Best max_features (mtry): {best_mtry}")