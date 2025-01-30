"""
Single supp figure. Sex ratios in each gene.
"""

import pandas as pd
import os
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact

def plot_mut_prevalence_and_sex_per_gene(muts_df_main, sex_df, ax0, subset_to_top10 = False):
    """
    muts_df: Can be CHIP or ctDNA.
    sex_df: Two column df showing patient ID and sex.
    ax0 = barchart
    ax1 = counts (broken ax)
    ax2 = counts (broken ax)
    """
    
    # Generate the plotting DF
    n_total_male = sex_df[sex_df["Sex"] == "Male"].shape[0]
    n_total_female = sex_df[sex_df["Sex"] == "Female"].shape[0]
    
    # add sex data to the genomic DFs
    muts_df = muts_df_main.merge(sex_df, how = "inner", on = "Patient_id")
    muts_df = muts_df[["Patient_id", "Gene", "Sex"]].drop_duplicates()
    muts_df = muts_df[["Gene", "Sex"]].value_counts().reset_index().sort_values(by = "Gene").rename(columns = {0: "Count"})
    muts_df["ratio"] = np.where(
        muts_df["Sex"] == "Male",
        muts_df["Count"] / n_total_male,
        muts_df["Count"] / n_total_female
    )
    muts_df["ratio"] = muts_df["ratio"]*100
    filtered_df = muts_df.groupby('Gene').filter(
        lambda x: (
            'Male' in x['Sex'].values and 
            'Female' in x['Sex'].values and
            x.loc[x['Sex'] == 'Male', 'ratio'].values[0] > 0 and
            x.loc[x['Sex'] == 'Female', 'ratio'].values[0] > 0
        )
    )    
    # genes not in this criteria will be other
    other_genes = muts_df[~muts_df['Gene'].isin(filtered_df['Gene'])].copy()["Gene"].unique().tolist()
    other_genes_counts = muts_df_main[muts_df_main["Gene"].isin(other_genes)].merge(sex_df)[["Patient_id", "Sex"]].drop_duplicates()["Sex"].value_counts().reset_index()
    other_genes_counts["ratio"] = np.where(
        other_genes_counts["index"] == "Male",
        other_genes_counts["Sex"] / n_total_male,
        other_genes_counts["Sex"] / n_total_female
    )
    other_genes_counts["ratio"] = other_genes_counts["ratio"]*100
    other_genes_counts["Gene"] = "Other"
    other_genes_counts.columns = ["Sex", "Count", "ratio", "Gene"]
    other_genes_counts = other_genes_counts[["Gene", "Sex", "Count", "ratio"]]
    filtered_df = pd.concat([filtered_df, other_genes_counts])
    filtered_df = filtered_df.reset_index(drop = True).reset_index()
    
    # Sorting
    grouped_counts = filtered_df.groupby('Gene')['Count'].sum().reset_index()
    gene_order = grouped_counts.sort_values(by=['Count'], ascending=False)
    gene_order['sort_key'] = gene_order['Gene'].apply(lambda x: 1 if x == 'Other' else 0)
    gene_order = gene_order.sort_values(by=['sort_key', 'Count'], ascending=[True, False]).drop(columns=['sort_key'])
    filtered_df = gene_order.merge(filtered_df, how="inner", on="Gene")
    
    male_color = "#ffff00"
    female_color = "#6a5acd"
    
    if subset_to_top10:
        filtered_df= filtered_df.iloc[0: 20]
    
    # Plot the prevalence
    p_values = []
    for idx, (gene, group) in enumerate(filtered_df.groupby("Gene", sort = False)):
        gene = group['Gene'].unique()  # Get unique genes
        male_ratio = group[group['Sex'] == 'Male'].set_index('Gene')['ratio'].reindex(gene, fill_value=0)[0]
        female_ratio = group[group['Sex'] == 'Female'].set_index('Gene')['ratio'].reindex(gene, fill_value=0)[0]
        
        # Determine where to plot the bar
        ax0.bar(idx-0.2, male_ratio, color=male_color, width = 0.4)
        ax0.bar(idx+0.2, female_ratio, color=female_color, width = 0.4)
        
        # Annotate the n for males and females
        male_count = group[group['Sex'] == 'Male'].set_index('Gene')['Count_y'].reindex(gene, fill_value=0)[0]
        female_count = group[group['Sex'] == 'Female'].set_index('Gene')['Count_y'].reindex(gene, fill_value=0)[0]
        male_neg = n_total_male - male_count
        female_neg = n_total_female - female_count
        contingency_table = [[male_count, female_count], [male_neg, female_neg]]
        odds_ratio, p_value = fisher_exact(contingency_table)
        
        p_values.append(p_value)
        
        ax0.text(idx - 0.2, male_ratio, str(male_count), ha='center', va='bottom', fontsize=7, color='black')
        ax0.text(idx + 0.2, female_ratio, str(female_count), ha='center', va='bottom', fontsize=7, color='black')
        
        if p_value < 0.01: 
            p_rounded = round(p_value, 4)
            ax0.text(idx, 15, f"p={p_rounded}", ha='center', va='bottom', fontsize=6, color='black')
    
    ax0.set_ylabel("CH mutation frequency %", rotation = 90)
    ax0.spines[["top", "right"]].set_visible(False)
    ax0.tick_params(axis='x', bottom=False)
    ax0.set_ylim((0, 50))
    
    # Determine x ticks and tick labels
    ax0.set_xticks(range(0, filtered_df["Gene"].unique().shape[0]))
    ax0.set_xticklabels(filtered_df["Gene"].unique(), rotation = 90, fontstyle = "italic")
    
    # add legend
    legend_colors = [male_color, female_color]
    legend_labels = ["Male", "Female"]
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax0.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.5)    
        
    return(ax0)

PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_treatment = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/treatment.csv"
PATH_date_surgery = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/operations.csv"
PATH_CHIP = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"
PATH_panel_genes = "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/chip_panel_gene_categories.tsv"

PATH_clinical_bladder = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv"
PATH_clinical_kidney = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Supplementary tables - Clinical data - mRCC.csv"

# LOAD CHIP DATASETS
all_vars_chip_main = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip_main = all_vars_chip_main[(all_vars_chip_main["Timepoint"] == "Baseline") & (all_vars_chip_main["Dependent"] == False)].reset_index(drop = True)
bladder_chip = all_vars_chip_main[all_vars_chip_main["Diagnosis"] == "Bladder"]
kidney_chip = all_vars_chip_main[all_vars_chip_main["Diagnosis"] == "Kidney"]

# LOAD SOMATIC DATASETS
all_vars_ctDNA_main  = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_ctDNA_main = all_vars_ctDNA_main[~all_vars_ctDNA_main["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
all_vars_ctDNA_main = all_vars_ctDNA_main[(all_vars_ctDNA_main["Timepoint"] == "Baseline") & (all_vars_ctDNA_main["Dependent"] == False)].reset_index(drop = True)

# pull sex data
sex_bladder = pd.read_csv(PATH_clinical_bladder)[["Patient_id", 'Sex']].drop_duplicates().reset_index(drop = True)
sex_kidney = pd.read_csv(PATH_clinical_kidney)[["Patient_id", 'Sex']].drop_duplicates().reset_index(drop = True)

# Plotting
fig = plt.figure(figsize=(6, 4))
outer_gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace = 0.8) # outer most gs with 3 rows

ax_bladder = fig.add_subplot(outer_gs[0])
ax_kidney = fig.add_subplot(outer_gs[1])
ax_bladder = plot_mut_prevalence_and_sex_per_gene(muts_df_main = bladder_chip, sex_df = sex_bladder, ax0 = ax_bladder, subset_to_top10 = True)
ax_kidney = plot_mut_prevalence_and_sex_per_gene(muts_df_main = kidney_chip, sex_df = sex_kidney, ax0 = ax_kidney, subset_to_top10 = True)
ax_bladder.set_title("mUC")
ax_kidney.set_title("mRCC")
outer_gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_ch_and_sex.png")
fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/SUPP_ch_and_sex.pdf")

# Are bladder men older than bladder women that have CH?
age_bladder = pd.read_csv(PATH_clinical_bladder)
age_bladder = age_bladder[age_bladder["First sample?"]][["Patient_id", "Age at blood draw", "Sex"]]

male_ages = age_bladder[age_bladder["Sex"] == "Male"]["Age at blood draw"]
female_ages = age_bladder[age_bladder["Sex"] == "Female"]["Age at blood draw"]

np.median(male_ages)
np.median(female_ages)
mannwhitneyu(male_ages, female_ages)

# Are CH+ bladder men older than CH+ bladder women that have CH?
mut_status = annotate_mutation_status(bladder_chip, "Bladder", PATH_sample_information, annotate_what = "CHIP")
mut_status = mut_status[mut_status["Timepoint"] == "Baseline"]
ch_pos_pts = mut_status[mut_status["CHIP status"] == "Positive"]

age_bladder = pd.read_csv(PATH_clinical_bladder)
age_bladder = age_bladder[age_bladder["First sample?"]][["Patient_id", "Age at blood draw", "Sex"]]
age_bladder = age_bladder[age_bladder["Patient_id"].isin(ch_pos_pts["Patient_id"])]

male_ages = age_bladder[age_bladder["Sex"] == "Male"]["Age at blood draw"]
female_ages = age_bladder[age_bladder["Sex"] == "Female"]["Age at blood draw"]

np.median(male_ages)
np.median(female_ages)
mannwhitneyu(male_ages, female_ages)