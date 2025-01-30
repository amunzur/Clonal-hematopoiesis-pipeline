import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from itertools import cycle
from scipy.stats import ttest_ind
import seaborn as sns
from lifelines import KaplanMeierFitter
from adjustText import adjust_text
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.cm import get_cmap
from datetime import datetime
import matplotlib.cm as cm
from lifelines.statistics import logrank_test
from scipy.stats import fisher_exact
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
import matplotlib.ticker as ticker
import upsetplot
from scipy.stats import mannwhitneyu
from lifelines import CoxPHFitter
import matplotlib.dates as mdates




mpl.rcParams['font.size'] = 10
mpl.rcParams['text.color'] = 'k'
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.handletextpad'] = '0.8'
mpl.rcParams['legend.labelspacing'] = '0.4'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['axes.labelsize'] = 10

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
PATH_gene_categories = os.path.join(DIR_working, "resources/panel/chip_panel_gene_categories.tsv")
PATH_kidney_clinical = os.path.join(DIR_working, "resources/clinical_data/Kidney/mRCC clinical Data_clean.csv")
PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/Bladder_enrollment.csv")
PATH_treatment_landscape = os.path.join(DIR_working, "resources/clinical_data/bladder/treatment.csv")
PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"
figure_dir = os.path.join(DIR_working, "results/figures/amazing_figures")
source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")

color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
all_vars_chip = sample_info.merge(all_vars_chip, how = "inner")

base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

###################################################################

for ext in ["pdf"]:
    for baseline_chip, prog_chip, baseline_somatic, prog_somatic, path_clinical, diagnosis in zip([base_bladder_chip, base_kidney_chip], [prog_bladder_chip, prog_kidney_chip], [base_bladder_somatic, base_kidney_somatic], [prog_bladder_somatic, prog_kidney_somatic], [PATH_clinical_bladder, PATH_kidney_clinical], ["bladder", "kidney"]): 
        figure_dir_diagnosis = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/bladder/Patient_profiles_line_graph_pdf"
        # figure_dir_diagnosis = os.path.join(figure_dir, diagnosis)
        # # patient profiles showing the change in VAF in CH mutations
        CALL_generate_patient_profiles_dot_plot(baseline_chip, prog_chip, baseline_somatic, prog_somatic, figure_dir_diagnosis, PATH_sample_information, PATH_treatment_landscape, ext = "pdf")
        
        # # This one is patient profiles too, but this time this is plotting a bar plot indicating the VAF of the mutations.
        # generate_patient_profiles(baseline_chip, prog_chip, figure_dir_diagnosis, ext)
        
        # plot_vaf(baseline_chip, figure_dir_diagnosis, f"{diagnosis}_vaf_baseline.{ext}") # What is the VAF and gene distribution of mutations?
        # plot_vaf(prog_chip, figure_dir_diagnosis, f"{diagnosis}_vaf_prog.{ext}") # What is the VAF and gene distribution of mutations that appear in the progression - they may also have appeared in the baseline before as well.
        
        # # Barchart showing the number of mutations at baseline and progression, genes on the x axis, number of times they appear in the cohort is on the y axis.
        # plot_gene_counts(baseline_chip, figure_dir_diagnosis, f"{diagnosis}_gene_counts_ALL_MUTATIONS_baseline.{ext}", f"{diagnosis} - all mutations at baseline", unique = False)
        # plot_gene_counts(baseline_chip, figure_dir_diagnosis, f"{diagnosis}_gene_counts_UNIQUE_MUTATIONS_baseline.{ext}", f"{diagnosis} - unique mutations at baseline", unique = True)
        
        # plot_gene_counts(prog_chip, figure_dir_diagnosis, f"{diagnosis}_gene_counts_ALL_MUTATIONS_progression.{ext}", f"{diagnosis} - all mutations at progression", unique = False)
        # plot_gene_counts(prog_chip, figure_dir_diagnosis, f"{diagnosis}_gene_counts_UNIQUE_MUTATIONS_progression.{ext}", f"{diagnosis} - unique mutations at progression", unique = True)
        
        # What percentage of the total cohort is CH positive at baseline and progression timepoints?
        plot_pie_for_ch_presence_in_df(baseline_chip, PATH_sample_information, figure_dir_diagnosis, f"Baseline {diagnosis} samples", f"{diagnosis}_pie_cohort_baseline.{ext}")
        plot_pie_for_ch_presence_in_df(prog_chip, PATH_sample_information, figure_dir_diagnosis, f"Progression {diagnosis} samples", f"{diagnosis}_pie_cohort_progression.{ext}")
        
        # Histogram showing the number of patients that have n mutations
        plot_per_patient_counts(baseline_chip, figure_dir_diagnosis, f"{diagnosis}_patient_mutation_counts_baseline_CH.{ext}")
        plot_per_patient_counts(baseline_chip, figure_dir_diagnosis, f"{diagnosis}_patient_mutation_counts_baseline_CH_both_cohorts.{ext}", colorby = "Diagnosis")
        plot_per_patient_counts(baseline_somatic, figure_dir_diagnosis, f"{diagnosis}_patient_mutation_counts_baseline_ctDNA.{ext}")
        
        # Age plots that show the distirbution of mutation + patients in each age bin
        plot_age_plots(baseline_chip, PATH_sample_information, figure_dir_diagnosis, path_clinical, fontsize = 10, ext = ext)
        plot_age_plots(baseline_somatic, PATH_sample_information, figure_dir_diagnosis, path_clinical, fontsize = 15, ext = ext)
        if diagnosis == "bladder":
            plot_age_plots(prog_chip, PATH_sample_information, figure_dir_diagnosis, path_clinical, fontsize = 15)
            plot_age_plots(prog_somatic, PATH_sample_information, figure_dir_diagnosis, path_clinical, fontsize = 15)
        
        # How does the VAF in cfDNA correlate with the WBC VAF?
        plot_vaf_scatter(baseline_chip, figure_dir_diagnosis, f"{diagnosis}_baseline_vaf_correlation.{ext}")
        plot_vaf_scatter(prog_chip, figure_dir_diagnosis, f"{diagnosis}_progression_vaf_correlation.{ext}")

##########################################################################
# PLOTS THAT INTERSECT GENOMIC DATA WITH CLINICAL DATA - BLADDER
###########################################################################
PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/Bladder_enrollment.csv")
figures_clinical_bladder = os.path.join(DIR_working, "results/figures/amazing_figures/bladder/clinical")

clin_bladder = pd.read_csv(PATH_clinical_bladder)[["Patient_id", "Sex", "Previous smoking history", "Age at baseline blood draw"]]
clin_bladder_mibc = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/nmibc_mibc.csv")
base_bladder_chip = base_bladder_chip.merge(clin_bladder, how = "inner", on = "Patient_id")

# CHIP AND AGE
age_vs_CH_presence(df_CH = base_bladder_chip, numerical_col_to_plot = "Age at baseline blood draw", PATH_sample_information = PATH_sample_information, clin = clin_bladder, figure_title = "", figure_dir = figures_clinical_bladder, fontsize = 10)
age_vs_CH_presence(df_CH = base_bladder_chip, numerical_col_to_plot = "Age at baseline blood draw", PATH_sample_information = PATH_sample_information, clin = clin_bladder, figure_title = "Age and presence of CH>2 percent in bladder", figure_dir = figures_clinical_bladder, min_WBC_VAF_threshold = 2, fontsize = 10)
age_vs_CH_presence(df_CH = base_bladder_chip, numerical_col_to_plot = "Age at baseline blood draw", PATH_sample_information = PATH_sample_information, clin = clin_bladder, figure_title = "Age and presence of CH>10 percent in bladder", figure_dir = figures_clinical_bladder, min_WBC_VAF_threshold = 10, fontsize = 10)

# SEX AND CHIP
categorical_vs_CH_vaf(df_CH = base_bladder_chip, categorical_col_to_plot = "Sex", clin =  clin_bladder, figure_title = "Sex and WBC VAF of any CH\nShowing only the mutation with the highest VAF per patient.", figure_dir = figures_clinical_bladder)

# SMOKING AND CHIP
clin_bladder_smoking = clin_bladder.copy()
clin_bladder_smoking["Previous smoking history"]= clin_bladder_smoking["Previous smoking history"].replace({"Current smoker": "Current or previous smoker", "Previous smoker": "Current or previous smoker"})
categorical_vs_CH_vaf(df_CH = base_bladder_chip, categorical_col_to_plot = "Previous smoking history", clin =  clin_bladder_smoking, figure_title = "Smoking history and WBC VAF of any CH\nShowing only the mutation with the highest VAF per patient.", figure_dir = figures_clinical_bladder)

# CH presence absence and age
figure_dir_age_plots = os.path.join(figure_dir, "clinical_bladder")
age_vs_CH_presence(base_bladder_chip, "Age at last follow-up or death", PATH_sample_information, clin_bladder, "CH status at baseline vs age", figure_dir_age_plots)
age_vs_CH_presence(base_bladder_chip, "Age at last follow-up or death", PATH_sample_information, clin_bladder, f"CH>2% status at baseline vs age", figure_dir_age_plots, min_WBC_VAF_threshold = 2)

# CH presence absence in a select gene vs age
for i, gene in enumerate(base["Gene"].unique()):
    age_vs_CH_presence(base_bladder_chip, "Age at last follow-up or death", PATH_sample_information, clin_bladder, f"{gene} CH status at baseline vs age", figure_dir_age_plots, gene)
    age_vs_CH_presence(prog, "Age", PATH_sample_info, clin, f"{gene} CH status at progression vs age", figure_dir_age_plots, gene)

# CH presence absence and PSA
age_vs_CH_presence(base, "PSA", PATH_sample_info, clin, "CH status at baseline vs PSA", figure_dir)
age_vs_CH_presence(prog, "PSA", PATH_sample_info, clin, "CH status at progression vs PSA", figure_dir)

# CH presence absence and Haemoglobin
age_vs_CH_presence(base, "Haemoglobin", PATH_sample_info, clin, "CH status at baseline vs Haemoglobin", figure_dir)
age_vs_CH_presence(prog, "Haemoglobin", PATH_sample_info, clin, "CH status at progression vs Haemoglobin", figure_dir)

# CH presence absence and LDH
age_vs_CH_presence(base, "LDH", PATH_sample_info, clin, "CH status at baseline vs LDH", figure_dir)
age_vs_CH_presence(prog, "LDH", PATH_sample_info, clin, "CH status at progression vs LDH", figure_dir)

# SURVIVAL ANALYSISf
figure_dir_bladder = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/bladder"
dir_bladder_survival = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/bladder/survival"
do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, base_bladder_chip, "CHIP", "Bladder_KM_CHIP.png", plot_title = "KM stratified by presence of any CH event", annotate_gene = False, figure_dir = dir_bladder_survival)
do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, base_bladder_somatic, "ctDNA", "Bladder_KM_ctDNA.png", plot_title = "KM stratified by presence of ctDNA", annotate_gene = False, figure_dir = dir_bladder_survival)



# 1% and 2% cutoff
base_bladder_chip_1_perc = base_bladder_chip[base_bladder_chip["VAF_n"] >= 1].reset_index(drop = True)
do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, base_bladder_chip_1_perc, "CHIP", "Bladder_KM_CHIP_1_percent.pdf", plot_title = "KM stratified by presence of CH>1%", annotate_gene = False, figure_dir = dir_bladder_survival)

base_bladder_chip_2_perc = base_bladder_chip[base_bladder_chip["VAF_n"] >= 2].reset_index(drop = True)
do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, base_bladder_chip_2_perc, "CHIP", "Bladder_KM_CHIP_2_percent.pdf", plot_title = "KM stratified by presence of CH>2%", annotate_gene = False, figure_dir = dir_bladder_survival)

# for gene in ["DNMT3A", "ASXL1", "TET2", "TP53", "ATM", "PPM1D", "CHEK2"]:
panel_genes = pd.read_csv(PATH_gene_categories, sep = "\t")["Panel genes"].tolist()
for gene in panel_genes:
    print(f"Working on {gene}")
    do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, base_bladder_chip, "CHIP", f"{gene}_CHIP_Bladder_KM.pdf", plot_title = f"KM stratified by presence of any CH in {gene}", annotate_gene = gene, figure_dir = dir_bladder_survival)
    do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, base_bladder_chip_2_perc, "CHIP", f"{gene}_CHIP_2_percent_Bladder_KM.pdf", plot_title = f"KM stratified by presence CH>2% in {gene}", annotate_gene = gene, figure_dir = dir_bladder_survival)
    do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, base_bladder_somatic, "ctDNA", f"{gene}_ctDNA_Bladder_KM.pdf", f"KM stratified by presence of ctDNA in {gene}", annotate_gene = gene, figure_dir = dir_bladder_survival)

# Now do a gene category specific
cats = {"DTA genes": ["DNMT3A", "TET2", "ASXL1"], "DDR genes": ["BRCA1", "BRCA2", "TP53", "CHEK2", "ATM"], "genes involved in splicing": ['SF3B1', 'SRSF2', 'U2AF1', 'ZRSR2']}
for cat in cats.keys(): 
    genes = cats[cat]
    print(f"Working on {cat}")
    do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, base_bladder_chip, "CHIP", f"{cat}_CHIP_Bladder_KM.png", plot_title = f"KM stratified by presence of any CH in {cat}", annotate_gene = genes, figure_dir = dir_bladder_survival)
    do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, base_bladder_chip_2_perc, "CHIP", f"{cat}_CHIP_2_percent_Bladder_KM.png", plot_title = f"KM stratified by presence CH>2% in {cat}", annotate_gene = genes, figure_dir = dir_bladder_survival)
    do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, base_bladder_somatic, "ctDNA", f"{cat}_ctDNA_Bladder_KM.png", f"KM stratified by presence of ctDNA in {cat}", annotate_gene = genes, figure_dir = dir_bladder_survival)

# Survival in non-DNMT3A CHIP
non_dnmt_chip = base_bladder_chip[base_bladder_chip["Gene"] != "DNMT3A"].reset_index(drop = True)
do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, non_dnmt_chip, "CHIP", f"non-DNMT3A_CHIP_Bladder_KM.pdf", plot_title = f"KM stratified by presence of any non-DNMT3A CH", annotate_gene = False, figure_dir = dir_bladder_survival)

non_dnmt_chip_2_percent = base_bladder_chip[(base_bladder_chip["Gene"] != "DNMT3A") & (base_bladder_chip["VAF_n"] >= 2)].reset_index(drop = True)
do_survival_analysis(PATH_clinical_bladder, PATH_sample_information, non_dnmt_chip_2_percent, "CHIP", f"non-DNMT3A_CHIP_2_percent_Bladder_KM.pdf", plot_title = f"KM stratified by presence of non-DNMT3A CH>2%", annotate_gene = False, figure_dir = dir_bladder_survival)

###################################################################################
# Now the next part runs Fisher's exact test in two groups, DOESN'T MAKE PLOTS
bladder_base_ctDNA_status = annotate_mutation_status(base_bladder_somatic, "Bladder", PATH_sample_information, annotate_what = "ctDNA", annotate_gene = False)
bladder_base_CHIP_status = annotate_mutation_status(base_bladder_chip, "Bladder", PATH_sample_information, annotate_what = "CHIP", annotate_gene = False)
base_bladder_status = bladder_base_ctDNA_status.merge(bladder_base_CHIP_status, how = "inner") # annotate patients based on ctDNA and CHIP status
clin_df = pd.read_csv(PATH_clinical_bladder)[["Patient_id", "Sex", "History of non-urothelial malignancy", "Previous smoking history"]]
clin_df["Smoking_general"] = clin_df["Previous smoking history"].isin(["Current smoker", "Previous smoker"]) # if they ever smoked or not

df = clin_df.merge(base_bladder_status, how = "inner")

for variable in ["Sex", "History of non-urothelial malignancy", "Smoking_general"]:
    contingency_table = pd.crosstab(df[variable], df['CHIP status'])
    odds_ratio, p_value = fisher_exact(contingency_table)
    print(variable)
    print(f"Odds Ratio: {odds_ratio}")
    print(f"P-value: {p_value}")
    print(" ")

# Same thing, now in a gene specific manner
for gene in ["DNMT3A", "ASXL1", "TET2", "TP53", "ATM", "CHEK2", "PPM1D"]:
    bladder_base_ctDNA_status = annotate_mutation_status(base_bladder_somatic, "Bladder", PATH_sample_information, annotate_what = "ctDNA", annotate_gene = gene)
    bladder_base_CHIP_status = annotate_mutation_status(base_bladder_chip, "Bladder", PATH_sample_information, annotate_what = "CHIP", annotate_gene = gene)
    base_bladder_status = bladder_base_ctDNA_status.merge(bladder_base_CHIP_status, how = "inner") # annotate patients based on ctDNA and CHIP status
    clin_df = pd.read_csv(PATH_clinical_bladder)[["Patient_id", "Sex", "History of non-urothelial malignancy", "Previous smoking history"]]
    clin_df["Smoking_general"] = clin_df["Previous smoking history"].isin(["Current smoker", "Previous smoker"]) # if they ever smoked or not
    df = clin_df.merge(base_bladder_status, how = "inner")
    for variable in ["Sex", "History of non-urothelial malignancy", "Smoking_general"]:
        contingency_table = pd.crosstab(df[variable], df['CHIP status'])
        odds_ratio, p_value = fisher_exact(contingency_table)
        print(variable, gene)
        print(f"Odds Ratio: {odds_ratio}")
        print(f"P-value: {p_value}")
        print(" ")

#############################################
# Next, we plot the patient profiles (dot plot) together in the same fig for pts that have a TP53 CH event.
tp53_pts = base_bladder_chip[base_bladder_chip["Gene"] == "TP53"]["Patient_id"].tolist()
fig = plt.figure(figsize=(12, 5))  # Increase the height to accommodate two rows
fig.suptitle(f'Patients with bladder cancer who have TP53 CH mutations')
num_cols = int(len(tp53_pts)/2)
gs = fig.add_gridspec(nrows=2, ncols=num_cols, width_ratios=[1] * num_cols, hspace = 0.4)
diagnosis = "Bladder"
show_dates_on_x_axis = False

for i, patient_id in enumerate(tp53_pts): 
    if i < 5: 
        row_index = 0
        col_index = i
    else: 
        row_index = 1
        col_index = i - 5
    ax = fig.add_subplot(gs[row_index, col_index])
    patient_chip_df = all_vars_chip[all_vars_chip["Patient_id"] == patient_id]
    ax = generate_patient_profiles_dot_plot(patient_chip_df, figure_dir, PATH_sample_information, show_dates_on_x_axis, ax, "tab10", patient_id, diagnosis, "Clonal hematopoiesis", title=patient_id, gene = "TP53", show_legend = False)
    if col_index != 0:
        ax.set_ylabel("")

fig.savefig(os.path.join(figure_dir, "TP53_positive_pts_profile_bladder.pdf"))
#############################################

make_mutation_count_histograms(base_bladder_chip, figure_dir)











##########################################################################
# PLOTS THAT INTERSECT GENOMIC DATA WITH CLINICAL DATA - KIDNEY DATASET
###########################################################################
kidney_figure_dir = os.path.join(DIR_working, "results/figures/amazing_figures/kidney")
figures_survival_kidney = os.path.join(DIR_working, "results/figures/amazing_figures/kidney/survival")
PATH_kidney_clean = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/kidney_clin_clean.csv"

diagnosis = "Kidney"

clin_df_kidney = pd.read_csv(PATH_kidney_clinical)
clin_df_kidney = clin_df_kidney[clin_df_kidney["Patient_id"] != "23-083"]
clin_df = clin_df_kidney[["Patient_id", "Sex", "DOB", "DOD", "Date of GUBB draw", "irAE"]].rename(columns = {"GUBB ID": "Patient_id", "DOD": "Date of last follow-up or death"})
clin_df["Death at last follow up"] = ~clin_df["Date of last follow-up or death"].isin(["X", "XX"])
clin_df.loc[clin_df["Date of last follow-up or death"].isin(["X", "XX"]) , "Date of last follow-up or death"] = pd.to_datetime("2023-10-17") # thats when we had the clinical data sent to us
clin_df[['DOB', 'Date of last follow-up or death', 'Date of GUBB draw']] = clin_df[['DOB', 'Date of last follow-up or death', 'Date of GUBB draw']].apply(pd.to_datetime)
clin_df["Age at last follow-up or death"] = (clin_df["Date of last follow-up or death"] - clin_df["DOB"]).astype('<m8[Y]').astype(int)

PATH_kidney_clean = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/kidney_clin_clean.csv"
clin_df.to_csv(PATH_kidney_clean)

# CHIP AND AGE
age_vs_CH_presence(df_CH = base_kidney_chip, numerical_col_to_plot = "Age at last follow-up or death", PATH_sample_information = PATH_sample_information, clin = clin_df, figure_title = "Age and any CH in kidney", figure_dir = kidney_figure_dir)
age_vs_CH_presence(df_CH = base_kidney_chip, numerical_col_to_plot = "Age at last follow-up or death", PATH_sample_information = PATH_sample_information, clin = clin_df, figure_title = "Age and CH>2 in kidney", figure_dir = kidney_figure_dir, min_WBC_VAF_threshold = 2)
age_vs_CH_presence(df_CH = base_kidney_chip, numerical_col_to_plot = "Age at last follow-up or death", PATH_sample_information = PATH_sample_information, clin = clin_df, figure_title = "Age and CH>10 in kidney", figure_dir = kidney_figure_dir, min_WBC_VAF_threshold = 10)

# SEX AND CHIP
categorical_vs_CH_vaf(df_CH = base_kidney_chip, categorical_col_to_plot = "Sex", clin =  clin_df, figure_title = "Kidney - Sex and WBC VAF of any CH\nShowing only the mutation with the highest VAF per patient.", figure_dir = kidney_figure_dir)

# CH presence absence in a select gene vs age
# for gene in panel_genes:
#     age_vs_CH_presence(df_CH = base_kidney_chip, numerical_col_to_plot = "Age at last follow-up or death", PATH_sample_information = PATH_sample_information, clin = clin_df, figure_title = f"Age and presence of CH in {gene} in kidney", figure_dir = figures_clinical, gene = gene)

# SURVIVAL ANALYSIS
do_survival_analysis(PATH_kidney_clean, PATH_sample_information, base_kidney_chip, "CHIP", f"{diagnosis}_KM_CHIP.pdf", plot_title = "KM stratified by presence of any CH event", annotate_gene = False, diagnosis = diagnosis, figure_dir = figures_survival_kidney)
do_survival_analysis(PATH_kidney_clean, PATH_sample_information, base_kidney_somatic, "ctDNA", f"{diagnosis}_KM_ctDNA.pdf", plot_title = "KM stratified by presence of ctDNA", annotate_gene = False, diagnosis = diagnosis, figure_dir = figures_survival_kidney)

# 1% and 2% cutoff
df_1_perc = base_kidney_chip[base_kidney_chip["VAF_n"] >= 1].reset_index(drop = True)
do_survival_analysis(PATH_kidney_clean, PATH_sample_information, df_1_perc, "CHIP", f"{diagnosis}_KM_CHIP_1_percent.pdf", plot_title = "KM stratified by presence of CH>1%", annotate_gene = False, diagnosis = diagnosis, figure_dir = figures_survival_kidney)

df_2_perc = base_kidney_chip[base_kidney_chip["VAF_n"] >= 2].reset_index(drop = True)
do_survival_analysis(PATH_kidney_clean, PATH_sample_information, df_2_perc, "CHIP", f"{diagnosis}_KM_CHIP_2_percent.pdf", plot_title = "KM stratified by presence of CH>2%", annotate_gene = False, diagnosis = diagnosis, figure_dir = figures_survival_kidney)

# for gene in ["DNMT3A", "ASXL1", "TET2", "TP53", "ATM", "PPM1D", "CHEK2"]:
panel_genes = pd.read_csv(PATH_gene_categories, sep = "\t")["Panel genes"].tolist()
for gene in panel_genes:
    print(f"Working on {gene}")
    do_survival_analysis(PATH_kidney_clean, PATH_sample_information, base_kidney_chip, "CHIP", f"{gene}_CHIP_{diagnosis}_KM.pdf", plot_title = f"KM stratified by presence of any CH in {gene}", annotate_gene = gene, diagnosis = diagnosis, figure_dir = figures_survival_kidney)
    do_survival_analysis(PATH_kidney_clean, PATH_sample_information, df_2_perc, "CHIP", f"{gene}_CHIP_2_percent_{diagnosis}_KM.pdf", plot_title = f"KM stratified by presence CH>2% in {gene}", annotate_gene = gene, diagnosis = diagnosis, figure_dir = figures_survival_kidney)
    do_survival_analysis(PATH_kidney_clean, PATH_sample_information, base_kidney_somatic, "ctDNA", f"{gene}_ctDNA_{diagnosis}_KM.pdf", f"KM stratified by presence of ctDNA in {gene}", annotate_gene = gene, diagnosis = diagnosis, figure_dir = figures_survival_kidney)

# Now do a gene category specific
cats = {"DTA genes": ["DNMT3A", "TET2", "ASXL1"], "DDR genes": ["BRCA1", "BRCA2", "TP53", "CHEK2", "ATM"], "genes involved in splicing": ['SF3B1', 'SRSF2', 'U2AF1', 'ZRSR2']}
for cat in cats.keys(): 
    genes = cats[cat]
    print(f"Working on {cat}")
    do_survival_analysis(PATH_kidney_clean, PATH_sample_information, base_kidney_chip, "CHIP", f"{cat}_CHIP_{diagnosis}_KM.pdf", plot_title = f"KM stratified by presence of any CH in {cat}", annotate_gene = genes, diagnosis = diagnosis, figure_dir = figures_survival_kidney)
    do_survival_analysis(PATH_kidney_clean, PATH_sample_information, df_2_perc, "CHIP", f"{cat}_CHIP_2_percent_{diagnosis}_KM.pdf", plot_title = f"KM stratified by presence CH>2% in {cat}", annotate_gene = genes, diagnosis = diagnosis, figure_dir = figures_survival_kidney)
    do_survival_analysis(PATH_kidney_clean, PATH_sample_information, base_kidney_somatic, "ctDNA", f"{cat}_ctDNA_{diagnosis}_KM.pdf", f"KM stratified by presence of ctDNA in {cat}", annotate_gene = genes, diagnosis = diagnosis, figure_dir = figures_survival_kidney)

# Survival in non-DNMT3A CHIP
non_dnmt_chip = base_kidney_chip[base_kidney_chip["Gene"] != "DNMT3A"].reset_index(drop = True)
do_survival_analysis(PATH_kidney_clean, PATH_sample_information, non_dnmt_chip, "CHIP", f"non-DNMT3A_CHIP_{diagnosis}_KM.pdf", plot_title = f"KM stratified by presence of any non-DNMT3A CH", annotate_gene = False, diagnosis = diagnosis, figure_dir = figures_survival_kidney)

non_dnmt_chip_2_percent = base_kidney_chip[(base_kidney_chip["Gene"] != "DNMT3A") & (base_kidney_chip["VAF_n"] >= 2)].reset_index(drop = True)
do_survival_analysis(PATH_kidney_clean, PATH_sample_information, non_dnmt_chip_2_percent, "CHIP", f"non-DNMT3A_CHIP_2_percent_{diagnosis}_KM.pdf", plot_title = f"KM stratified by presence of non-DNMT3A CH>2%", annotate_gene = False, diagnosis = diagnosis, figure_dir = figures_survival_kidney)

##########################################################################
# PLOTS THAT INTERSECT GENOMIC DATA WITH CLINICAL DATA - MEGACOHORT! 
# Plotting survival curves for kidney and bladder together
###########################################################################

base_combined = pd.concat([base_kidney_chip, base_bladder_chip]).reset_index(drop = True)
base_combined_somatic = pd.concat([base_kidney_somatic, base_bladder_somatic]).reset_index(drop = True)
diagnosis = "Both"

c1 = pd.read_csv(PATH_clinical_bladder)[["Patient_id", 'Date of last follow-up or death', 'Death at last follow up', 'Age at last follow-up or death']]
c2 = pd.read_csv(PATH_kidney_clean)[["Patient_id", 'Date of last follow-up or death', 'Death at last follow up', 'Age at last follow-up or death']]
c3 = pd.concat([c1, c2]).reset_index(drop = True)
c3.to_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder_kidney_combined_survival.csv", index = False)

PATH_clinical = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder_kidney_combined_survival.csv"
dir_combined_survival = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/survival_combined"
do_survival_analysis(PATH_clinical, PATH_sample_information, base_combined, "CHIP", "Any_ch_megacohort.pdf", plot_title = f"KM stratified by presence of CH", annotate_gene = False, diagnosis = "Both", figure_dir = dir_combined_survival)
do_survival_analysis(PATH_clinical, PATH_sample_information, base_combined_somatic, "ctDNA", "Any_ctdna_megacohort.pdf", plot_title = f"KM stratified by presence of ctDNA", annotate_gene = False, diagnosis = "Both", figure_dir = dir_combined_survival)

# 1% and 2% cutoff
df_1_perc = base_combined[base_combined["VAF_n"] >= 1].reset_index(drop = True)
do_survival_analysis(PATH_clinical, PATH_sample_information, df_1_perc, "CHIP", f"{diagnosis}_KM_CHIP_1_percent.png", plot_title = "KM stratified by presence of CH>1%", annotate_gene = False, diagnosis = diagnosis, figure_dir = dir_combined_survival)

df_2_perc = base_combined[base_combined["VAF_n"] >= 2].reset_index(drop = True)
do_survival_analysis(PATH_clinical, PATH_sample_information, df_2_perc, "CHIP", f"{diagnosis}_KM_CHIP_2_percent.png", plot_title = "KM stratified by presence of CH>2%", annotate_gene = False, diagnosis = diagnosis, figure_dir = dir_combined_survival)

# for gene in ["DNMT3A", "ASXL1", "TET2", "TP53", "ATM", "PPM1D", "CHEK2"]:
panel_genes = pd.read_csv(PATH_gene_categories, sep = "\t")["Panel genes"].tolist()
for gene in panel_genes:
    print(f"Working on {gene}")
    do_survival_analysis(PATH_clinical, PATH_sample_information, base_combined, "CHIP", f"{gene}_CHIP_{diagnosis}_KM.png", plot_title = f"KM stratified by presence of any CH in {gene}", annotate_gene = gene, diagnosis = diagnosis, figure_dir = dir_combined_survival)
    do_survival_analysis(PATH_clinical, PATH_sample_information, df_2_perc, "CHIP", f"{gene}_CHIP_2_percent_{diagnosis}_KM.png", plot_title = f"KM stratified by presence CH>2% in {gene}", annotate_gene = gene, diagnosis = diagnosis, figure_dir = dir_combined_survival)
    do_survival_analysis(PATH_clinical, PATH_sample_information, base_combined_somatic, "ctDNA", f"{gene}_ctDNA_{diagnosis}_KM.png", f"KM stratified by presence of ctDNA in {gene}", annotate_gene = gene, diagnosis = diagnosis, figure_dir = dir_combined_survival)

# Now do a gene category specific
cats = {"DTA genes": ["DNMT3A", "TET2", "ASXL1"], "DDR genes": ["BRCA1", "BRCA2", "TP53", "CHEK2", "ATM"], "genes involved in splicing": ['SF3B1', 'SRSF2', 'U2AF1', 'ZRSR2']}
for cat in cats.keys(): 
    genes = cats[cat]
    print(f"Working on {cat}")
    do_survival_analysis(PATH_clinical, PATH_sample_information, base_combined, "CHIP", f"{cat}_CHIP_{diagnosis}_KM.pdf", plot_title = f"KM stratified by presence of any CH in {cat}", annotate_gene = genes, diagnosis = diagnosis, figure_dir = dir_combined_survival)
    do_survival_analysis(PATH_clinical, PA
    TH_sample_information, base_combined_somatic, "ctDNA", f"{cat}_ctDNA_{diagnosis}_KM.pdf", f"KM stratified by presence of ctDNA in {cat}", annotate_gene = genes, diagnosis = diagnosis, figure_dir = dir_combined_survival)

# non DNMT3A chip
base_combined_non_dnmt3a = base_combined[base_combined["Gene"] != "DNMT3A"]
do_survival_analysis(PATH_clinical, PATH_sample_information, base_combined_non_dnmt3a, "CHIP", f"nonDNMT3A_chip_megacohort.png", plot_title = f"KM stratified by presence of non-DNMT3A CH", annotate_gene = False, diagnosis = "Both", figure_dir = dir_combined_survival)

