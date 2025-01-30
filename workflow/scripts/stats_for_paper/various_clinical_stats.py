import pandas as pd
import numpy as np
import os 
import re
from scipy.stats import ttest_ind
from lifelines import KaplanMeierFitter
from datetime import datetime
from lifelines.statistics import logrank_test
from scipy.stats import fisher_exact
from lifelines import KaplanMeierFitter
from scipy.stats import mannwhitneyu
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr
from scipy.stats import chi2_contingency

# LOAD DATASETS
DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"
PATH_bladder_clinical_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"
PATH_kidney_clinical_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/RCC clinical data - mRCC clinical Data.csv"

source_functions = "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/visualization/UTILITIES_make_chip_plots.py"
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# CH mutations
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[(all_vars_chip["Dependent"] == False) & (all_vars_chip["Timepoint"] == "Baseline")]
bladder_chip = all_vars_chip[all_vars_chip["Diagnosis"] == "Bladder"]
kidney_chip = all_vars_chip[all_vars_chip["Diagnosis"] == "Kidney"]

bladder_clin_main = pd.read_csv(PATH_bladder_clinical_data)
bladder_clin_main = bladder_clin_main[bladder_clin_main["First sample?"]]

# 1. Do BCG exposed patients have higher CH burden?
bcg = bladder_clin_main[["Patient_id", "Disease stage at initial diagnosis", "Age at blood draw", "Type of intravesical therapy administered for NMIBC"]]
bcg = bcg[bcg["Disease stage at initial diagnosis"] != "Metastatic"]

# Are mibc patients older than localized nmibc?
mibc_age = bcg[bcg["Disease stage at initial diagnosis"] == "Localized MIBC"]["Age at blood draw"]
nmibc_age = bcg[bcg["Disease stage at initial diagnosis"] == "Localized NMIBC"]["Age at blood draw"]
mannwhitneyu(mibc_age, nmibc_age)
np.median(mibc_age)
np.median(nmibc_age)

# Get CH status
ch_status = annotate_mutation_status(bladder_chip[bladder_chip["VAF_n"] > 2], "Bladder", PATH_sample_information, annotate_what = "CHIP")
ch_status = ch_status[ch_status["Timepoint"] == "Baseline"].reset_index(drop = True)[["Patient_id", "CHIP status"]].merge(bcg[["Patient_id", "Disease stage at initial diagnosis"]])

contingency_table = pd.crosstab(ch_status['Disease stage at initial diagnosis'], ch_status['CHIP status'])
print(contingency_table)
odds_ratio, p_value = fisher_exact(contingency_table)
print(f"Odds Ratio: {odds_ratio}")
print(f"P-value: {p_value}")

# Patients with CH and ctDNA mutations in the same gene
PATH_kidney_clinical = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/CHIP Supplementary tables - Clinical_data_mRCC.csv"
PATH_clinical_bladder = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"

kidney_clin = pd.read_csv(PATH_kidney_clinical)
bladder_clin = pd.read_csv(PATH_clinical_bladder)
bladder_clin = bladder_clin[bladder_clin["First sample?"]]
bladder_clin = bladder_clin.rename(columns = {"Age at blood draw": "Age at GUBB draw"})

age_df = pd.concat([kidney_clin[["Patient_id", "Age at GUBB draw"]], bladder_clin[["Patient_id", "Age at GUBB draw"]]]).reset_index(drop = True)
bins = [20, 40, 50, 60, 70, 80, 90]
labels = [f'{bins[i]}-{bins[i+1]}' for i in range(len(bins) - 1)]
age_df['Age_Bin'] = pd.cut(age_df['Age at GUBB draw'], bins=bins, labels=labels, right=False)
age_df['Age_Bin'] = pd.Categorical(age_df['Age_Bin'], categories=labels, ordered=True)

mut_status = annotate_mutation_status(all_vars_chip, "Both", PATH_sample_information, annotate_what = "CHIP")
mut_status = mut_status[mut_status["Timepoint"] == "Baseline"]
age_df = mut_status.merge(age_df, how = "inner")

# Step 1: Calculate fraction of positives in each age bin
age_bins = []
fractions = []

for i, group in age_df.groupby("Age_Bin"):
    n_pos = group[group["CHIP status"] == "Positive"].shape[0]
    n_neg = group[group["CHIP status"] == "Negative"].shape[0]
    fr_pos = n_pos / (n_pos + n_neg)  # Fraction of Positive
    age_bins.append(i)
    fractions.append(fr_pos)

# Convert age bins and fractions to numpy arrays
# Create a mapping of age bins to numeric labels
age_bin_mapping = {'20-40': 0, '40-50': 1, '50-60': 2, '60-70': 3, '70-80': 4, '80-90': 5}

# Map the age bins in your DataFrame to numeric labels
age_bins_numeric = np.array([age_bin_mapping[bin_label] for bin_label in age_bins])

# Now use these numeric labels for fitting the polynomial
coefficients = np.polyfit(age_bins_numeric, fractions, 2)
polynomial = np.poly1d(coefficients)

# Plotting will remain the same
x_vals = np.linspace(age_bins_numeric.min(), age_bins_numeric.max(), 100)
y_vals = polynomial(x_vals)

# Step 4: Plot the data and the fitted polynomial
# Step 3: Create figure and axis from scratch
fig, ax = plt.subplots(figsize=(8, 6))

# Step 4: Plot data points
ax.scatter(age_bins_numeric, fractions, color='blue', label='Data Points')

# Step 5: Plot the polynomial curve
ax.plot(x_vals, y_vals, color='red', label='2nd Degree Polynomial Fit')

plt.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/test.png")
plt.show()




import statsmodels.api as sm

# Convert 'CHIP status' to binary (0 = Negative, 1 = Positive)
age_df['CHIP_binary'] = age_df['CHIP status'].apply(lambda x: 1 if x == 'Positive' else 0)

# If 'Age_Bin' is already categorical, convert it to numeric codes (ordered)
age_df['Age_Bin'] = pd.Categorical(age_df['Age_Bin'], ordered=True)
age_df['Age_Bin_numeric'] = age_df['Age_Bin'].cat.codes  # Convert Age Bin to numeric

# Step 1: Prepare X (predictor) and y (outcome)
X = age_df['Age_Bin_numeric']  # Age bins as predictor
y = age_df['CHIP_binary']      # Binary CHIP status as outcome

# Step 2: Add constant to X for the intercept
X = sm.add_constant(X)

# Step 3: Fit logistic regression model
logit_model = sm.Logit(y, X)
result = logit_model.fit()

# Step 4: Print summary of the model
print(result.summary())

# Step 5: Calculate odds ratios
odds_ratios = pd.DataFrame({
    "Odds Ratio": np.exp(result.params),
    "p-value": result.pvalues
})
print(odds_ratios)



















# Step 1: Generate predicted probabilities from the logistic regression model
age_bins_numeric_range = np.linspace(age_df['Age_Bin_numeric'].min(), age_df['Age_Bin_numeric'].max(), 100)
X_plot = sm.add_constant(age_bins_numeric_range)
y_pred_probs = result.predict(X_plot)

# Step 2: Create the figure and axis
fig, ax = plt.subplots(figsize=(8, 6))

# Step 3: Plot the actual fraction of CH-positive individuals by age bin
ax.scatter(sorted(age_df['Age_Bin_numeric'].unique()), fractions, color='blue', label='Actual Fraction Positive')

# Step 4: Plot the logistic regression predicted probabilities
ax.plot(age_bins_numeric_range, y_pred_probs, color='red', label='Logistic Regression Fit')

# Step 5: Customize the plot
ax.set_xlabel('Age Bin')
ax.set_ylabel('Probability of being CH-positive')
ax.set_title('Probability of CH-positive by Age Bin')
ax.set_xticks(sorted(age_df['Age_Bin_numeric'].unique()))
ax.set_xticklabels(age_df['Age_Bin'].cat.categories, rotation=45)  # Use the original age bin labels
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/test.png")

# Step 6: Show the plot
plt.tight_layout()

