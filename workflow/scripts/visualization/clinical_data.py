import numpy as np 
import pandas as pd
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import statsmodels.api as sm
import matplotlib.patches as mpatches
from scipy.stats import kruskal
import scipy.stats as stats
from scipy.stats import linregress



PATH_clinical = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/kidney.tsv"
PATH_muts = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined.csv"

clinical = pd.read_csv(PATH_clinical, "\t")
muts = pd.read_csv(PATH_muts)

muts = muts.merge(clinical, "left")
muts = muts.dropna(subset=['Age_at_draw']) # Filter rows where 'Age_at_draw' is not NaN
muts = muts[muts["n_callers"] != 1].reset_index(drop = True)




# muts["Gene_colors"] = muts["Gene"].map(gene_colors)
# muts["Gene_colors"] = muts["Gene_colors"].fillna("Black")

#==================================================================
# Number of mutations and age 
filtered_muts = muts.drop_duplicates(["Patient_id", "Protein_annotation"]).reset_index(drop = True)
filtered_muts = filtered_muts.drop_duplicates(["Patient_id", "Protein_annotation"]).reset_index(drop = True)
filtered_muts= filtered_muts[filtered_muts["n_callers"] != 1].reset_index(drop = True)

# 2. Group by 'Patient_id' and count the number of mutations for each patient
mutation_counts = filtered_muts.groupby('Patient_id').size().reset_index(name='Mutation_Count')
filtered_muts = filtered_muts.merge(mutation_counts, "left")

bins = list(range(9))  # Define the mutation count bins as needed
filtered_muts['Mutation_Bins'] = pd.cut(filtered_muts['Mutation_Count'], bins)
age_by_mutation_bins = filtered_muts.groupby('Mutation_Bins')['Age_at_draw'].mean()  # or median()
filtered_muts["Bins_enumerated"] = filtered_muts['Mutation_Bins'].map({sample: i / (len(filtered_muts['Mutation_Bins'].unique()) - 1) for i, sample in enumerate(filtered_muts['Mutation_Bins'].unique())})
filtered_muts= filtered_muts.drop_duplicates(["Patient_id", "Age_at_draw", "Bins_enumerated"])

filtered_muts.to_csv("/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/figures/test.csv")

# Calculate significance and add asterisks
significance_threshold = 0.05  # Set the significance threshold

# Perform Kruskal-Wallis test for all unique groups
group_values = filtered_muts.groupby('Bins_enumerated')['Age_at_draw'].apply(list)
_, p_value = kruskal(*group_values.values)

fig, ax = plt.subplots(figsize=(5, 5))

# Plot box plots
boxplot_data = [group['Age_at_draw'] for _, group in filtered_muts.groupby('Bins_enumerated')]
boxprops = dict(facecolor='lightgray', color='black', linewidth=1.5)
medianprops = dict(color='black', linewidth=2)  
ax.boxplot(boxplot_data, positions=np.sort(filtered_muts["Bins_enumerated"].unique()), patch_artist=True, medianprops=medianprops, boxprops=boxprops, showfliers=False, widths=0.07)

# Plot scatter points
for bin_value, group in filtered_muts.groupby('Bins_enumerated'):
    jitter_x = np.random.uniform(-0.02, 0.02, len(group))  # Add jitter within the range [-0.02, 0.02]
    ax.scatter(group["Bins_enumerated"] + jitter_x, group['Age_at_draw'], s = 14, alpha=1, color="black", zorder = 10)

font_size = 12
ax.set_xticks(np.sort(filtered_muts["Bins_enumerated"].unique()))
ax.set_xticklabels(['1', '2', '3', '4', '6', '7', '8'], fontsize=font_size)
ax.set_yticks([20, 30, 40, 50, 60, 70, 80, 90])
ax.set_yticklabels(["20", "30", "40", "50", "60", "70", "80", "90"], fontsize=font_size)
ax.set_xlim(-0.1, 1.1)
ax.set_xlabel('Number of Mutations', fontsize=font_size)
ax.set_ylabel('Age at Draw', fontsize=font_size)
ax.set_title('Number of mutations and age at blood draw', fontsize=font_size)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/figures/age_and_mutcounts.png")
#==================================================================

#==================================================================
# VAF and age
gene_colors = {
    'DNMT3A': 'dodgerblue',   # Blue
    'TET2': 'red',     # Red
    'ASXL1': 'limegreen',    # Green
    'Other': 'Black'
}

age_bins = [20, 40, 50, 60, 70, 80, 90]
age_labels = ['21-40', '41-50', '51-60', '61-70', '71-80', '81-90']

filtered_muts = muts[["Gene", "Age_at_draw", "VAF_n"]]
filtered_muts["VAF_n"] = filtered_muts["VAF_n"] * 100
filtered_muts['Mutation_Bins'] = pd.cut(filtered_muts['Age_at_draw'], age_bins, labels=age_labels)
filtered_muts["Bins_enumerated"] = pd.cut(filtered_muts['Age_at_draw'], age_bins, labels=False)

boxprops = dict(facecolor='lightgray', color='black', linewidth=1.5)
medianprops = dict(color='black', linewidth=2)  


# ==================================================
# Below 2 %
below_2_perc = filtered_muts[filtered_muts["VAF_n"] < 2]
boxplot_data = [group['VAF_n'] for _, group in below_2_perc.groupby('Bins_enumerated')]

fig, ax = plt.subplots(figsize=(7, 5))
ax.boxplot(boxplot_data, positions=np.sort(below_2_perc["Bins_enumerated"].unique()), patch_artist=True, medianprops=medianprops, boxprops=boxprops, showfliers=False)

for bin_value, group in below_2_perc.groupby('Bins_enumerated'):
    jitter_x = np.random.uniform(-0.03, 0.03, len(group))  # Add jitter within the range [-0.02, 0.02]
    ax.scatter(group["Bins_enumerated"] + jitter_x, group['VAF_n'], s = 14, alpha=1, color="black", zorder = 10)

ax.set_xticklabels(np.sort(list(below_2_perc['Mutation_Bins'].unique())))
ax.set_xlabel('Age at blood draw')
ax.set_ylabel('Variant allele frequency (%)')
ax.set_title('VAF and Age, for mutations with VAF < 2%')
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

group_values = below_2_perc.groupby('Bins_enumerated')['VAF_n'].apply(list)
_, p_value = kruskal(*group_values.values)
plt.text(0.05, 0.95, f'Kruskal-Wallis p-value = {p_value:.4f}', transform=plt.gca().transAxes)
plt.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/figures/age_and_vaf_below_2_perc.png")
# ==================================================

# ==================================================
# Above 2 %
above_2_perc = filtered_muts[filtered_muts["VAF_n"] >= 2]
boxplot_data = [group['VAF_n'] for _, group in above_2_perc.groupby('Bins_enumerated')]

fig, ax = plt.subplots(figsize=(7, 5))
ax.boxplot(boxplot_data, positions=np.sort(above_2_perc["Bins_enumerated"].unique()), patch_artist=True, medianprops=medianprops, boxprops=boxprops, showfliers=False)

for bin_value, group in above_2_perc.groupby('Bins_enumerated'):
    jitter_x = np.random.uniform(-0.03, 0.03, len(group))  # Add jitter within the range [-0.02, 0.02]
    ax.scatter(group["Bins_enumerated"] + jitter_x, group['VAF_n'], s = 14, alpha=1, color="black", zorder = 10)

ax.set_xticklabels(np.sort(list(above_2_perc['Mutation_Bins'].unique())))
ax.set_xlabel('Age at blood draw')
ax.set_ylabel('Variant allele frequency (%)')
ax.set_title('VAF and Age, for mutations with VAF \u2265 2%')
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
group_values = above_2_perc.groupby('Bins_enumerated')['VAF_n'].apply(list)
_, p_value = kruskal(*group_values.values)
plt.text(0.05, 0.95, f'Kruskal-Wallis p-value = {p_value:.4f}', transform=plt.gca().transAxes)

plt.tight_layout()
plt.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/figures/age_and_vaf_above_2_perc.png")
# ==================================================



# ==================================================
# FRACTION OF PTS WITH CHIP IN EACH AGE GROUP AND CHI SQUARE TEST
age_bins = [20, 40, 50, 60, 70, 80, 90]
age_labels = ['21-40', '41-50', '51-60', '61-70', '71-80', '81-90']

clinical = pd.read_csv(PATH_clinical, "\t")
muts = pd.read_csv(PATH_muts)
muts = muts[muts["n_callers"] != 1].reset_index(drop = True)
combined = clinical.merge(muts[["Patient_id", "Gene", "VAF_n"]], "left")

combined["VAF_n"] = combined["VAF_n"] * 100
combined['Mutation_Bins'] = pd.cut(combined['Age_at_draw'], age_bins, labels=age_labels)
combined["Bins_enumerated"] = pd.cut(combined['Age_at_draw'], age_bins, labels=False)
combined["CHIP"] = combined["VAF_n"] > 2
# combined["CHIP"] = ~pd.isnull(combined["Gene"])
combined = combined.groupby('Patient_id').apply(lambda group: group.iloc[group['CHIP'].values.argmax()]).reset_index(drop=True)

# Group the data by 'Mutation_Bins' and 'CHIP', and calculate the count for each combination
grouped_data = combined.groupby(['Mutation_Bins', 'CHIP']).size().unstack(fill_value=0)
grouped_data['Total'] = grouped_data[True] + grouped_data[False]
grouped_data['Fraction_True'] = grouped_data[True] / grouped_data['Total']
grouped_data['Fraction_False'] = grouped_data[False] / grouped_data['Total']

# Plot the stacked bar chart
fig, ax = plt.subplots(figsize=(7, 5))
grouped_data[['Fraction_True', 'Fraction_False']].plot(kind='bar', stacked=True, color=['black', 'lightgrey'], ax=ax)

# Set axis labels and title
font_size = 13
ax.set_xlabel('Number of patients in each age category', fontsize = font_size)
ax.set_xticklabels(grouped_data.index, rotation=0, fontsize = font_size)
ax.set_xticklabels([f'{label}\nn = {count}' for label, count in zip(age_labels, grouped_data["Total"])], rotation=0, fontsize = font_size)
ax.set_ylabel('Fraction of patients', fontsize = font_size)
ax.set_title('Fraction of patients with CH \u2265 2% in each age category', fontsize = font_size)
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_yticklabels(["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(title='CH with VAF \u2265 2%', labels=['Present', 'Absent'])

# chi square test
contingency_table = np.column_stack((grouped_data[True], grouped_data[False]))
chi2_stat, p_value, dof, expected = stats.chi2_contingency(contingency_table)
plt.text(0.6, 0.7, f'Chi-square test p-value = {p_value:.4f}', transform=plt.gca().transAxes)

plt.tight_layout()
plt.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/figures/above_2_perc_stacked.png")
# ==================================================

# ==================================================
# ASSOCIATING PRESENCE ABSENCE WITH SEX
clinical = pd.read_csv(PATH_clinical, "\t")
muts = pd.read_csv(PATH_muts)
muts = muts[muts["n_callers"] != 1].reset_index(drop = True)

combined = clinical.merge(muts[["Patient_id", "Gene", "VAF_n"]], "left")
combined["CH"] = ~pd.isnull(combined["Gene"])
combined = combined.groupby('Patient_id').apply(lambda group: group.iloc[group['CH'].values.argmax()]).reset_index(drop=True)

grouped_data = combined.groupby(['Sex', 'CH']).size().unstack(fill_value=0)
fractions = grouped_data.div(grouped_data.sum(axis=1), axis=0)
fractions.rename(columns={True: 'Present', False: 'Absent'}, inplace=True)

# Plot the stacked bar plot
colors = {"Present": 'Black', "Absent": 'lightgrey'}
fig, ax = plt.subplots(figsize=(2, 2))  # Adjust the figure size here
ax = fractions.plot(kind='bar', stacked=True, color=[colors[col] for col in fractions.columns], width=0.3)

# Customize the plot
font_size = 14
ax.set_xlabel('Biological Sex', fontsize=font_size)
ax.set_xticklabels(["Female\nn = 46", "Male\nn = 131"], rotation=0, fontsize=font_size)
ax.set_ylabel('Fraction', fontsize=font_size)
ax.set_title('CH in Females and Males', fontsize=font_size)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Show the plot
plt.tight_layout()
plt.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/figures/sex_CH_presence_absence.png")




# ==================================================
# Associating the presence/absence of genes with age
font_size = 12
age_bins = [20, 40, 50, 60, 70, 80, 90]
age_labels = ['21-40', '41-50', '51-60', '61-70', '71-80', '81-90']

filtered_muts = muts[muts["Diagnosis"] == "Kidney"]
clinical = pd.read_csv(PATH_clinical, "\t")

filtered_muts = filtered_muts.merge(clinical, "left")


filtered_muts = filtered_muts[["Patient_id", "Gene", "Age_at_draw", "VAF_n"]]
filtered_muts.loc[filtered_muts['Age_at_draw'] < 40, 'Age_at_draw'] = 40 # Replace ages less than 50 with 50 in the Age_at_draw column
selected_genes = ["DNMT3A", "ASXL1", "TET2", "TP53", "CHEK2", "PPM1D", "ATM"]
filtered_muts["selected_genes"] = filtered_muts["Gene"].isin(selected_genes)
filtered_muts["Gene"] = filtered_muts.apply(lambda row: row["Gene"] if row["selected_genes"] else "Other", axis=1)
filtered_muts.drop(columns=["selected_genes"], inplace=True)
filtered_muts = filtered_muts.groupby(["Patient_id", "Gene"]).first().reset_index()

filtered_muts["VAF_n"] = filtered_muts["VAF_n"] * 100
filtered_muts["Age_Bins"] = pd.cut(filtered_muts["Age_at_draw"], bins = age_bins, labels = age_labels)

# Calculate the proportion of samples with mutations in each gene within each age group
sample_counts = filtered_muts.groupby(['Age_Bins', 'Gene']).size().unstack(fill_value=0)
sample_proportions = sample_counts.div(sample_counts.sum(axis=1), axis=0)

color_map = {gene: 'C'+str(i) for i, gene in enumerate(["DNMT3A", "ASXL1", "TET2", "TP53", "CHEK2", "PPM1D", "ATM"])}
color_map.update({"Other": "Gray"})

# Reorder columns in DataFrame to move "Other" to the top
cols = selected_genes + ["Other"]
sample_proportions = sample_proportions[cols]

fig, ax = plt.subplots(figsize=(6, 5))
sample_proportions.plot(kind='bar', stacked=True, ax=ax, color=[color_map[gene] for gene in sample_proportions.columns], width=0.5)

# Customize the plot for sample proportions
age_bin_counts = filtered_muts["Age_Bins"].value_counts().sort_index(ascending=True)

x_ticks_mutation_counts = [f'{label}\nn = {count}' for label, count in zip(age_labels, age_bin_counts)]
x_ticks_patient_counts = ["\nn = 9", "\nn = 13", "\nn = 36", "\nn = 64", "\nn = 46", "\nn = 8"]
x_ticks_combined = [x_ticks_mutation_counts[i] + x_ticks_patient_counts[i] for i in range(len(x_ticks_mutation_counts))]

ax.set_xticklabels(x_ticks_combined, rotation=0, fontsize=font_size)
ax.set_yticklabels(["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=font_size)
ax.set_xlabel('Number of unique gene mutations in each age group', fontsize=font_size)
ax.set_ylabel('Fraction of Patients', fontsize=font_size)
ax.set_title('Fraction of Patients with CH in Select Genes by Age', fontsize=font_size)
ax.legend(title='Gene', bbox_to_anchor=(1.05, 1), loc='upper left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/figures/gene_presence_absence_by_age_samples.png")
plt.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/figures/gene_presence_absence_by_age_samples.pdf")


# ==================================================
# Scatter plot in front of box plot. X is categorical, M and F. Y axis is the VAFs.
filtered_muts = muts[["Patient_id", "Gene", "Sex", "VAF_n"]]
filtered_muts["VAF_n"] = filtered_muts["VAF_n"] * 100

# Group the data by 'Mutation_Bins' and 'CHIP', and calculate the count for each combination
grouped_data = filtered_muts.groupby("Sex")
male_vaf = grouped_data.get_group("M")["VAF_n"]
female_vaf = grouped_data.get_group("F")["VAF_n"]
data = [male_vaf, female_vaf]

fig, ax = plt.subplots(figsize=(5, 4))
ax.boxplot(female_vaf, positions=[0], patch_artist=True, medianprops=medianprops, boxprops=boxprops, showfliers=False)
ax.boxplot(male_vaf, positions=[1], patch_artist=True, medianprops=medianprops, boxprops=boxprops, showfliers=False)

for _, vaf in enumerate(female_vaf):
    ax.scatter(0 + np.random.uniform(-0.2, 0.2), vaf, s = 14, alpha=1, color="black", zorder = 10)

for _, vaf in enumerate(male_vaf):
    ax.scatter(1 + np.random.uniform(-0.2, 0.2), vaf, s = 14, alpha=1, color="black", zorder = 10)

font_size = 12
ax.set_xticklabels(["Female", "Male"], fontsize=font_size)
ax.set_xlabel('Biological Sex', fontsize=font_size)
ax.set_ylabel('Variant allele frequency (%)', fontsize=font_size)
ax.set_yticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
ax.set_yticklabels(["0", "10", "20", "30", "40", "50", "60", "70", "80", "90"], fontsize=font_size)
ax.set_title('Distribution of VAFs in Males and Females', fontsize=font_size)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

statistic, p_value = stats.mannwhitneyu(male_vaf, female_vaf, alternative='two-sided')
plt.text(0.05, 0.90, f'Mann-Whitney p-value = {p_value:.4f}', transform=plt.gca().transAxes)

# f_median = np.median(female_vaf)
# m_median = np.median(male_vaf)
# plt.text(0.15, 0.8, f'median = {f_median:.2f}', transform=plt.gca().transAxes)
# plt.text(0.75, 0.8, f'median = {m_median:.2f}', transform=plt.gca().transAxes)


plt.tight_layout()
plt.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/figures/sex_vaf.png")
# ==================================================


# ==================================================
# Linear model for all ages and all mutations
filtered_muts = muts[["Patient_id", "VAF_n", "Age_at_draw"]]
filtered_muts["VAF_n"] = filtered_muts["VAF_n"]*100
filtered_muts = filtered_muts[filtered_muts["VAF_n"] > 2]
fig, ax = plt.subplots(figsize=(6, 4))
ax.scatter(filtered_muts["Age_at_draw"], filtered_muts["VAF_n"], s = 14, alpha=1, color="black", zorder = 10)

# Calculate the linear regression
slope, intercept, r_value, p_value, std_err = linregress(filtered_muts["Age_at_draw"], filtered_muts["VAF_n"])
line_x = np.linspace(filtered_muts["Age_at_draw"].min(), filtered_muts["Age_at_draw"].max(), 100)
line_y = slope * line_x + intercept
ax.plot(line_x, line_y, color='blue', linestyle='dashed', label=f"Regression Line: y = {slope:.2f}x + {intercept:.2f}")

# Add equation and R-squared value to the plot
equation_text = f"y = {slope:.2f}x + {intercept:.2f}"
r2_text = f"R\u00b2 = {r_value**2:.2f}"
ax.text(0.05, 0.95, equation_text, transform=ax.transAxes, fontsize=12, verticalalignment='top')
ax.text(0.05, 0.9, r2_text, transform=ax.transAxes, fontsize=12, verticalalignment='top')

ax.set_xlabel('Age at Blood Draw')
ax.set_ylabel('Variant allele frequency (%)')
ax.set_title('Age at Blood Draw and VAF(%) for CH > 2%')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/figures/age_vaf.png")
