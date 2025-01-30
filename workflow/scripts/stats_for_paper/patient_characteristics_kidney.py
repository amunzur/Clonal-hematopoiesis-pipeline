import pandas as pd
import numpy as np
import os 
import re
from scipy.stats import fisher_exact

# Assuming df is already defined

source_functions = "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/visualization/UTILITIES_make_chip_plots.py"
with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

PATH_clinical_bladder = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Bladder_enrollment.csv"
PATH_clinical_kidney = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Supplementary tables - Clinical data - mRCC.csv"
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

all_vars_chip = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv")
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False].reset_index(drop = True)

baseline_chip = all_vars_chip[all_vars_chip["Timepoint"] == "Baseline"]
ot_chip = all_vars_chip[all_vars_chip["Timepoint"] == "During treatment"]

base_kidney_chip = baseline_chip[baseline_chip["Diagnosis"] == "Kidney"].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[all_vars_chip["Diagnosis"] == "Kidney"].reset_index(drop = True)
base_bladder_chip = baseline_chip[all_vars_chip["Diagnosis"] == "Bladder"].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[all_vars_chip["Diagnosis"] == "Bladder"].reset_index(drop = True)

# Load and rename columns
kidney_clin = pd.read_csv(PATH_clinical_kidney).rename(columns={"GUBB ID": "Patient_id", "Age at GUBB draw": "Age"})
bladder_clin = pd.read_csv(PATH_clinical_bladder).rename(columns={"GUBB ID": "Patient_id", "Age at baseline blood draw": "Age"})

# Subset to samples we are including in the study
kidney_pts = sample_info[(sample_info["Diagnosis"] == "Kidney") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()
bladder_pts = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()

# Filter the clinical data based on the patients included in the study
kidney_clin = kidney_clin[kidney_clin["Patient_id"].isin(kidney_pts)]
bladder_clin = bladder_clin[bladder_clin["Patient_id"].isin(bladder_pts)]

n_denom = kidney_clin.shape[0] # number of patient for whom clinical data is available
############################### Kidney
# Number of patients exposed/not exposed to treatment prior to baseline sample collection
treatment_dates_df = kidney_clin[["Patient_id", "Date of start of 1L systemic therapy", "Date of end of 1L systemic therapy", "Type of treatment"]]
treatment_dates_df = treatment_dates_df.merge(sample_info[sample_info["Timepoint"] == "Baseline"].reset_index(drop = True), on = "Patient_id")

# Convert date columns from strings to datetime objects
treatment_dates_df["Date of start of 1L systemic therapy"] = pd.to_datetime(treatment_dates_df["Date of start of 1L systemic therapy"],  format='%d-%b-%y', errors='coerce')
treatment_dates_df["Date of end of 1L systemic therapy"] = pd.to_datetime(treatment_dates_df["Date of end of 1L systemic therapy"], format='%d-%b-%y', errors='coerce')
treatment_dates_df["Date_collected"] = pd.to_datetime(treatment_dates_df["Date_collected"], format='%Y%b%d')

treatment_dates_df["Draw before tx"] = treatment_dates_df["Date_collected"] <= treatment_dates_df["Date of start of 1L systemic therapy"] 
treatment_dates_df["Draw minus DOSxStart"] = treatment_dates_df["Date_collected"] - treatment_dates_df["Date of start of 1L systemic therapy"] 
treatment_dates_df["Tx duration"] = treatment_dates_df["Date of end of 1L systemic therapy"] - treatment_dates_df["Date of start of 1L systemic therapy"] 

# GET THE NUMBER OF TRUE BASELINE SAMPLES
no_treatment_patients_df = treatment_dates_df[treatment_dates_df["Date of start of 1L systemic therapy"].isna()]
n_no_systemic = no_treatment_patients_df.shape[0] # number of patients that didn't receive systemic treatment at all

draw_before_tx_df = treatment_dates_df[treatment_dates_df["Draw minus DOSxStart"] <= pd.Timedelta(days=8)]
n_draw_before_tx = draw_before_tx_df.shape[0]
round((n_no_systemic+n_draw_before_tx)/n_denom*100) # percent and number of true baseline

# Patients with prior exposure to treatment when their blood was drawn
draw_after_tx_df = treatment_dates_df[treatment_dates_df["Draw minus DOSxStart"] > pd.Timedelta(days=8)]

treatment_dict = {1: "Ipilimumab-Nivolumab", 
                  2: "Pembrolizumab-Axitinib", 
                  3: "Avelumab-Axitinib", 
                  4: "Nivolumab-Cabozantinib", 
                  5: "Pembrolizumab-Lenvatinib", 
                  6: "Sunitinib", 
                  7: "Pazopanib", 
                  8: "Quavonlimab + Pembrolizumab + Lenvatinib", 
                  9: "NKTR214 + Nivolumab + Axitinib"}

draw_after_tx_df["Treatment status"] = "Exposed"
draw_after_tx_df.shape[0]

draw_after_tx_df["Type of treatment"].value_counts()
draw_after_tx_df["Type of treatment"].value_counts(normalize=True).round(2) * 100


# Is CH more common in patients with prior exposure?
# Now combine these into one df
no_exposure_df = pd.concat([no_treatment_patients_df, draw_before_tx_df]).reset_index(drop = True)
no_exposure_df["Treatment status"] = "Unexposed"
no_exposure_df["Regimen"] = np.nan

df = pd.concat([no_exposure_df, draw_after_tx_df]).reset_index(drop = True)[["Patient_id", "Treatment status", "Regimen"]]
chip_status = annotate_mutation_status(base_kidney_chip, "Kidney", PATH_sample_information, annotate_what = "CHIP")
chip_status = chip_status[chip_status["Timepoint"] == "Baseline"]
df = df.merge(chip_status[["Patient_id", "CHIP status"]].drop_duplicates())

# Perform Fisher's Exact Test
contingency_table = pd.crosstab(df['Treatment status'], df['CHIP status'])
odds_ratio, p_value = fisher_exact(contingency_table)

# PART 2. METS INFORMATION
# percent with bone mets
kidney_clin["Lung_Met"].value_counts(normalize = True)
kidney_clin["Lung_Met"].value_counts()

kidney_clin["Bone_Met"].value_counts(normalize = True)
kidney_clin["Bone_Met"].value_counts()

kidney_clin["LN_Met"].value_counts(normalize = True)
kidney_clin["LN_Met"].value_counts()

# PART 3 
# Pathological subtype
subtype_dict = {"0": "Not available", 
                "1": "Clear cell", 
                "2": "Papillary", 
                "2a": "Papillary", 
                "2b": "Papillary", 
                "3": "Chromophobe", 
                "4": "Oncocytic / oncocytoma", 
                "5": "Collecting duct / Bellini tumour", 
                "6": "Undifferenciated or unclassifiable", 
                "7": "Mixed"}

kidney_clin["subtype_coded"] = kidney_clin["subtype"].map(subtype_dict)
kidney_clin["subtype_coded"].value_counts()
kidney_clin["subtype_coded"].value_counts(normalize=True)

