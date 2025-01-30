"""
Removes dates from supp tables.
"""
import pandas as pd
import os

def enumerate_bladder_OTs(PATH_sample_information):
    """
    Instead of dates enumerate them. 
    """
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    bladder_samples = sample_info[sample_info["Diagnosis"] == "Bladder"].reset_index(drop = True)
    bladder_samples["Date_collected"] = pd.to_datetime(bladder_samples["Date_collected"], format='%Y%b%d')
    
    # Create a helper function to label the timepoints
    def label_timepoints(dates):
        labels = ['Baseline']
        if len(dates) > 1:
            labels += [f"On treatment {i}" for i in range(1, len(dates))]
        return labels
    
    # Apply the logic for labeling timepoints based on the dates for each patient
    for i, group in bladder_samples.groupby("Patient_id"):
        if group.shape[0] > 1:
            # Sort the group by date (assumes there's a date column named "Date")
            group = group.sort_values("Date_collected")
            
            # Label the timepoints
            timepoints = label_timepoints(group["Date_collected"])
            
            # Assign the labels to the Timepoint column
            bladder_samples.loc[group.index, "Timepoint"] = timepoints
        else:
            bladder_samples.loc[group.index, "Timepoint"] = "Baseline"
    
    return(bladder_samples)
    
    


kidney_clinical = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Supplementary tables - Clinical data - mRCC.csv")
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
dir_output = "/groups/wyattgrp/users/amunzur/pipeline/results/submission_materials"
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

bladder_samples = enumerate_bladder_OTs(PATH_sample_information)

###################################################
# 1. CHIP DATASET
all_vars_chip = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv")
all_vars_chip["Date_collected"] = pd.to_datetime(all_vars_chip["Date_collected"], format='%Y%b%d')
idx_to_remove = all_vars_chip[(all_vars_chip["Diagnosis"] == "Kidney") & (all_vars_chip["Timepoint"] == "During treatment")].index # removing OT kidney samples
all_vars_chip = all_vars_chip.drop(idx_to_remove)

chip_cols = [
    'Patient_id', 'Diagnosis', 'Date_collected', 'Chrom', 'Position', 'Ref', 'Alt', 'Gene', 'Type', 'VAF_t', 'Depth_t', 'Alt_forward_t', 'Alt_reverse_t',
    'VAF_n', 'Depth_n', 'Alt_forward_n', 'Alt_reverse_n', 'Effects', 'Protein_annotation', 'error_rate_n', 'cosmic97_coding', 'avsnp150', 'gnomad40_exome_AF_t', 'CLNALLELEID', 'CLNSIG', 'Dependent']

all_vars_chip = all_vars_chip[chip_cols]
all_vars_chip = all_vars_chip.merge(bladder_samples, how = "left")
all_vars_chip["Timepoint"]= all_vars_chip["Timepoint"].fillna("Baseline")
all_vars_chip["cfDNA alt reads"] = all_vars_chip["Alt_forward_t"] + all_vars_chip["Alt_reverse_t"]
all_vars_chip["WBC alt reads"] = all_vars_chip["Alt_forward_n"] + all_vars_chip["Alt_reverse_n"]

# Reorder cols
ordered_cols = [
    'Patient_id', 'Diagnosis', 'Timepoint', 'Chrom', 'Position', 'Ref', 'Alt', 'Gene', 'Type', 'VAF_t', 'Depth_t', 'cfDNA alt reads', 'VAF_n', 'Depth_n', 'WBC alt reads', 'Effects', 
    'Protein_annotation', 'error_rate_n', 'cosmic97_coding', 'avsnp150', 'gnomad40_exome_AF_t', 'CLNALLELEID', 'CLNSIG', 'Dependent']
all_vars_chip = all_vars_chip[ordered_cols]

# Rename cols
chip_cols_renamed = [
    'Patient_id', 'Diagnosis', 'Timepoint', 'Chromosome', 'Position', 'Ref', 'Alt', 'Gene', 'Type', 'cfDNA VAF', 'cfDNA depth', 'cfDNA alt reads',
    'WBC VAF', 'WBC depth', 'WBC alt reads', 'Effects', 'Protein annotation', 'Error rate', 'cosmic97_coding', 'avsnp150', 'gnomad40_exome_AF', 'CLNALLELEID', 'CLNSIG', 'Dependent mutation calling?']
all_vars_chip.columns = chip_cols_renamed
all_vars_chip["Diagnosis"] = all_vars_chip["Diagnosis"].replace({"Kidney": "mRCC", "Bladder": "mUC"})
all_vars_chip = all_vars_chip.sort_values(by = ["Diagnosis", "Patient_id"]).reset_index(drop = True)
chip_path_output = os.path.join(dir_output, "Supplementary table 3 - CH variants.csv")
all_vars_chip.to_csv(chip_path_output, index = False)

###################################################
# 2. CTDNA MUTATION DATASET
all_vars_somatic = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv")
all_vars_somatic["Date_collected"] = pd.to_datetime(all_vars_somatic["Date_collected"], format='%Y%b%d')
idx_to_remove = all_vars_somatic[(all_vars_somatic["Diagnosis"] == "Kidney") & (all_vars_somatic["Timepoint"] == "During treatment")].index # removing OT kidney samples
all_vars_somatic = all_vars_somatic.drop(idx_to_remove)

ctdna_cols = [
    'Patient_id', 'Diagnosis', 'Date_collected', 'Chrom', 'Position', 'Ref', 'Alt', 'Gene', 'Type', 'VAF_t', 'Depth_t', 'Alt_forward_t', 'Alt_reverse_t', 'Effects', 
    'Protein_annotation', 'error_rate', 'cosmic97_coding', 'avsnp150', 'gnomad40_exome_AF', 'CLNALLELEID', 'CLNSIG', 'Dependent']

all_vars_somatic = all_vars_somatic[ctdna_cols]
all_vars_somatic = all_vars_somatic.merge(bladder_samples, how = "left")
all_vars_somatic["Timepoint"] = all_vars_somatic["Timepoint"].fillna("Baseline")
all_vars_somatic["cfDNA alt reads"] = all_vars_somatic["Alt_forward_t"] + all_vars_somatic["Alt_reverse_t"]

# Reorder cols
ordered_cols = [
    'Patient_id', 'Diagnosis', 'Timepoint', 'Chrom', 'Position', 'Ref', 'Alt', 'Gene', 'Type', 'VAF_t', 'Depth_t', 'cfDNA alt reads', 'Effects', 
    'Protein_annotation', 'error_rate', 'cosmic97_coding', 'avsnp150', 'gnomad40_exome_AF', 'CLNALLELEID', 'CLNSIG', 'Dependent']
all_vars_somatic = all_vars_somatic[ordered_cols]

# Rename cols
ctdna_cols_renamed = [
    'Patient_id', 'Diagnosis', 'Timepoint', 'Chromosome', 'Position', 'Ref', 'Alt', 'Gene', 'Type', 'cfDNA VAF', 'cfDNA depth', 'cfDNA alt reads', 'Effects', 
    'Protein annotation', 'Error rate', 'cosmic97_coding', 'avsnp150', 'gnomad40_exome_AF', 'CLNALLELEID', 'CLNSIG', 'Dependent mutation calling?']
all_vars_somatic.columns = ctdna_cols_renamed
all_vars_somatic["Diagnosis"] = all_vars_somatic["Diagnosis"].replace({"Kidney": "mRCC", "Bladder": "mUC"})

ctdna_path_output = os.path.join(dir_output, "Supplementary table 4 - ctDNA variants.csv")
all_vars_somatic.to_csv(ctdna_path_output, index = False)

# 3. Sequencing stats
path_input = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/seq_quality_metrics/CHIP Supplementary tables - Supplementary table 2 - Sequencing quality metrics.csv"
df = pd.read_csv(path_input).rename(columns = {"Patient": "Patient_id"})
df["Sample type"] = df["Sample_name"].str.extract(r"(WBC|cfDNA)")
df = df[df["Patient_id"].isin(sample_info["Patient_id"].unique())].reset_index(drop = True)
df["Date_collected"] = df["Sample_name"].str.replace(r"[_\-]", " ", regex=True).str.split().str[-1]
df["Date_collected"] = pd.to_datetime(df["Date_collected"], format='%Y%b%d')
idx_to_remove = df[(df["Diagnosis"] == "Kidney") & (df["Timepoint"] == "During treatment")].index # removing OT kidney samples
df = df.drop(idx_to_remove).reset_index(drop = True)
del df["Timepoint"]
df = df.merge(bladder_samples, how = "left")
df["Timepoint"]= df["Timepoint"].fillna("Baseline")
df = df[["Patient_id", "Diagnosis", "Sample type", "Timepoint", "Median depth", "On target rate", "Duplicate fraction", "Number of reads (millions)"]]
df = df.sort_values(by = ["Diagnosis", "Patient_id"]).reset_index(drop = True)
df["Diagnosis"] = df["Diagnosis"].replace({"Kidney": "mRCC", "Bladder": "mUC"})
seq_stats = os.path.join(dir_output, "Supplementary table 2 - Sequencing quality metrics.csv")
df.to_csv(seq_stats, index = False)

# 4. Bladder clinical data
bladder_clinical = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv")

# Function to calculate time difference from baseline
def calculate_time_diff(group):
    # Get the baseline date
    baseline_date = group.loc[group["Timepoint"] == "Baseline", "Date_collected"].values
    if len(baseline_date) > 0:
        baseline_date = pd.to_datetime(baseline_date[0])  # Extract the baseline date
        # Calculate the time difference for OT samples
        group["If OT sample, time difference from baseline (weeks)"] = group.apply(
            lambda row: round((row["Date_collected"] - baseline_date).days / 7) 
            if "On treatment" in row["Timepoint"] else pd.NA, axis=1
        )
    return group

bladder_clinical["Date_collected"] = pd.to_datetime(bladder_clinical["Date_collected"])
bladder_clinical = bladder_clinical.merge(bladder_samples)
bladder_clinical["If OT sample, time difference from baseline (weeks)"] = np.nan
bladder_clinical = bladder_clinical.groupby("Patient_id").apply(calculate_time_diff) # Apply the function for each patient group
bladder_clinical = bladder_clinical.reset_index(drop=True)

bladder_clinical = bladder_clinical.merge(bladder_samples)

cols_to_retain = [
    "Patient_id", "Diagnosis", "Timepoint", "Sex", "Previous smoking history", 
    "Total pack years", "Age at blood draw", "Radical surgery performed", "Histology of pathology specimen", "Death", 
    "OS from cfDNA collection (mo)", "Event", "Treatment", "Drug", "Reason discontinuation", "Progression", "PFS (mo)", 
    "Carboplatin", "Cisplatin", "Gemcitabine", "Docetaxel", "Pembrolizumab", "Durvalumab", "Avelumab", "Enfortumab Vedotin", "Erdafitinib"]

bladder_clinical = bladder_clinical[cols_to_retain]
bladder_clinical["Diagnosis"] = "mUC"
bladder_clinical.to_csv(os.path.join(dir_output, "Supplementary table 5 - mUC clinical data.csv"), index = False)