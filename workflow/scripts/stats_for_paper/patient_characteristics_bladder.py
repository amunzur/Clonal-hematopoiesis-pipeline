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
PATH_treatment_bladder = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/treatment.csv"

PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

all_vars_chip = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv")
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False].reset_index(drop = True)

baseline_chip = all_vars_chip[all_vars_chip["Timepoint"] == "Baseline"]
ot_chip = all_vars_chip[all_vars_chip["Timepoint"] == "During treatment"]

base_bladder_chip = baseline_chip[all_vars_chip["Diagnosis"] == "Bladder"].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[all_vars_chip["Diagnosis"] == "Bladder"].reset_index(drop = True)

bladder_pts = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "Baseline")].reset_index(drop = True)
n_denom = bladder_pts["Patient_id"].shape[0]
# Load and rename columns
bladder_clin = pd.read_csv(PATH_clinical_bladder).rename(columns={"GUBB ID": "Patient_id", "Age at baseline blood draw": "Age"})
treatment = pd.read_csv(PATH_treatment_bladder)
treatment = treatment[treatment["Treatment"] != "Palliative / best supportive care"].reset_index(drop = True)

# Filter the clinical data based on the patients included in the study
bladder_clin = bladder_clin.merge(bladder_pts, on = "Patient_id", how = "inner")
treatment = treatment.merge(bladder_pts, on = "Patient_id", how = "inner").drop("Comments", axis = 1)

treatment["Date start"] = pd.to_datetime(treatment["Date start"], errors='coerce')
treatment["Date discontinuation"] = pd.to_datetime(treatment["Date discontinuation"], errors='coerce')
treatment["Date_collected"] = pd.to_datetime(treatment["Date_collected"], format='%Y%b%d')

# treatment["Draw before tx"] = treatment["Date_collected"] <= treatment["Date start"] 
treatment["Draw minus start"] = treatment["Date_collected"] - treatment["Date start"] 
treatment["Tx duration"] = treatment["Date discontinuation"] - treatment["Date start"] 

# Identify patients whose baseline was taken after at least one line of treatment
prior_exposure_dict = {}
for i, group in treatment.groupby("Patient_id"): 
    # For each patient determine their prior therapy exposure.
    patient = group["Patient_id"].iloc[0]
    earliest_line = # start date of the earliest line of treatment
    for i, row in group.iterrows(): # iterate through the lines of treatment patient received
        date_difference = row["Date_collected"] - row["Date start"]
        if date_difference <= pd.Timedelta(days=-30): # Check if the difference is more than 30 days
            prior_exposure_dict[patient] = row["Treatment"]

prior_exposure_df = pd.DataFrame(list(prior_exposure_dict.items()), columns=['Patient_id', 'Treatment'])
prior_exposure_df.shape[0] # number of patients
round(prior_exposure_df.shape[0]/n_denom, 2)*100

prior_exposure_df["Treatment"].value_counts()
prior_exposure_df["Treatment"].value_counts(normalize = True)