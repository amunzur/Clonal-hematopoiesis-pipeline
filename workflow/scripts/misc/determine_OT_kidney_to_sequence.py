import pandas as pd 
import numpy as np
import os
from datetime import datetime


DIR_bams = "/groups/wyattgrp/users/amunzur/pipeline/results/data/bam/SSCS1_filtered"

# The purpose of this script is to identify which OT kidney samples to do. We choose samples that are at least 3 months apart from the baseline
need_library = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/sequencing_sheets/Copy of Bladder & Kidney sequencing - June 5, 10_01 a.m. - Just OT - Library_needs_to_be_made.csv")
need_hyb = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/sequencing_sheets/Copy of Bladder & Kidney sequencing - June 5, 10_01 a.m. - Just OT - Library_already_made.csv")
need_exraction = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/sequencing_sheets/Copy of Bladder & Kidney sequencing - June 5, 10_01 a.m. - Just OT - DNA_needs_to_be_extracted .csv")

# On treatment
OT_df_list = []
for file in [need_library, need_hyb, need_exraction]: 
    file = file[file["Timepoint"] == "During treatment"].reset_index()
    if not file.empty:     
        pt_col = file.columns[file.columns.str.contains("Patient")][0] # the col that has patient name
        file.rename(columns={pt_col: "Patient"}, inplace=True)
        file["Timepoint"] = "During treatment"
        OT_df_list.append(file[[pt_col, "Date collected", "Status", "Timepoint"]])

OT_df = pd.concat(OT_df_list) # the patient IDs that have at least one OT sample

# Baseline samples we already profiled
bams = [file for file in os.listdir(DIR_bams) if file.endswith(".bam") and "cfDNA" in file]
baseline_kidney = [item for item in bams if any(string in item for string in set(OT_df["Patient"]))] # these are kidney baseline samples

baseline_dict = {}
for path in baseline_kidney:
    patient = path.split("_")[0].replace("GUBB-", "")
    baseline_date_str = path.split("_")[-1].replace(".bam", "").replace(".fixed", "")
    # Convert baseline date format from "2020Aug05" to "5-Aug-2020"
    baseline_date = datetime.strptime(baseline_date_str, "%Y%b%d")
    formatted_baseline_date = baseline_date.strftime("%d-%b-%Y").lstrip("0")
    baseline_dict.update({patient: formatted_baseline_date})

baseline_df = pd.DataFrame(list(baseline_dict.items()), columns=['Patient', 'Date collected'])
baseline_df["Timepoint"] = "Baseline"
baseline_df["Status"] = "Done"
baseline_df = baseline_df.sort_values(by = "Patient").reset_index(drop = True)

combined = pd.concat((OT_df, baseline_df)).reset_index(drop = True)
combined.to_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/sequencing_sheets/all_kidney_OT_and_their_baseline.csv", index = False)



#======================================================================
data = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv", "\t")

data = data[data["Diagnosis"] == "Kidney"]

baseline_data = {}
during_data = {}

for index, row in data.iterrows():
    patient_id = row["Patient_id"]
    date_collected = row["Date_collected"]
    timepoint = row["Timepoint"]
    
    if timepoint == "Baseline":
        baseline_data[patient_id] = datetime.strptime(date_collected, "%Y%b%d")
    elif timepoint == "During treatment":
        if patient_id not in during_data:
            during_data[patient_id] = []
        during_data[patient_id].append(datetime.strptime(date_collected, "%Y%b%d"))

with open("/groups/wyattgrp/users/amunzur/pipeline/resources/sequencing_sheets/time_differences.txt", "w") as f:
    for patient_id in during_data:
        if patient_id in baseline_data:
            f.write(f"Patient ID: {patient_id}\n")
            f.write("Time differences:\n")
            
            baseline_date = baseline_data[patient_id]
            during_dates = during_data[patient_id]
            
            for i, during_date in enumerate(during_dates):
                time_diff_to_during = (during_date - baseline_date).days / 7  # Calculate in weeks
                f.write(f"Baseline sample taken on: {baseline_date.strftime('%Y-%b-%d')}\n")
                f.write(f"During treatment sample {i+1} taken on: {during_date.strftime('%Y-%b-%d')}\n")
                f.write(f"Time difference to during treatment {i+1}: {time_diff_to_during:.2f} weeks\n")
                
                if time_diff_to_during >= 24:
                    f.write("WARNING: TIME DIFFERENCE IS MORE THAN OR EQUAL TO 24 WEEKS\n")
                
                if i > 0:
                    time_diff_between_during = (during_date - during_dates[i - 1]).days / 7  # Calculate in weeks
                    f.write(f"Time difference between treatment {i} and treatment {i + 1}: {time_diff_between_during:.2f} weeks\n")
            
            f.write("------------------------\n")
