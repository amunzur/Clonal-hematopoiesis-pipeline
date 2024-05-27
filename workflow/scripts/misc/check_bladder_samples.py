"""
This script goes through all the FASTQs from bladder samples and ensures that we have both cfDNA and WBC. 
Returns stats like number of samples, number of baseline and OT etc.
"""

import os 
import pandas as pd
import numpy as np
import re
from collections import Counter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def return_OT_samples(file_list, patient_name):
    """
    Given a patient name, return number of time points we have samples available.
    """
    files = [item for item in file_list if patient_name in item]
    timepoints = set([re.split(r'[-_]', file)[-2] for file in files])
    n_timepoints = len(timepoints)
    return(n_timepoints)

def filter_strings(input_strings, cfDNA=True, baseline=True):
    """
    Based on the specified arguments, return the relevant samples. This is used to count how many cfDNA and WBC samples we have from baseline and OT.
    """
    if baseline and cfDNA:
        pattern = re.compile(r'(?=.*aseline)(?=.*cfDNA)(?=.*IDT)(?=.*_2\.fq\.gz)')
    elif baseline and not cfDNA:
        pattern = re.compile(r'(?=.*aseline)(?=.*WBC)(?=.*IDT)(?=.*_2\.fq\.gz)')
    elif not baseline and cfDNA:
        pattern = re.compile(r'(?=.*reatment)(?=.*cfDNA)(?=.*IDT)(?=.*_2\.fq\.gz)')
    elif not baseline and not cfDNA:
        pattern = re.compile(r'(?=.*reatment)(?=.*WBC)(?=.*IDT)(?=.*_2\.fq\.gz)')
    return [s for s in input_strings if re.match(pattern, s)]

def return_samples_from_patient(patient_id, file_list):
    """
    Given a single PATIENT ID return all of their fastq files across all timepoints.
    """
    return [s for s in file_list if re.match(patient_id, s)]

def check_cfDNA_and_WBC_concordance(sample_list):
    """
    Given a list of FastQ files from a single patient, check to see if they have two cfDNA files and two WBC files per timepoint. Two because read1 and read2 from FastQ. 
    """
    timepoints = set([re.split(r'[-_]', sample)[-2] for sample in sample_list])
    for timepoint in timepoints:
        timepoint_samples = [s for s in sample_list if re.match(re.compile(fr'.*{timepoint}.*'), s)]
        cfDNA_samples = [s for s in timepoint_samples if re.match(re.compile(fr'.*cfDNA.*'), s)]
        WBC_samples = [s for s in timepoint_samples if re.match(re.compile(fr'.*WBC.*'), s)]
        if len(cfDNA_samples) == len(WBC_samples) == 2: 
            print("Sample list ok.")
        else: 
            print("ERROR.")

#=========================================================================================
DIR_bladder_fastq = "/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/merged/bladder"
DIR_figures = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/misc"

file_list = os.listdir(DIR_bladder_fastq)
patient_list = set([x.split("_")[0] for x in file_list])

# number of patients
len(patient_list)

# FIGURE 1
# number of patients across timepoints
timepoints = [return_OT_samples(file_list, patient_name) for patient_name in patient_list]
element_counts = Counter(timepoints) # how many times each element appears
df = pd.DataFrame(list(element_counts.items()), columns=['Element', 'Count'])

fig, ax = plt.subplots()
ax.bar(df['Element'], df['Count'], color='gray', edgecolor='black')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('Number of timepoints')
ax.set_ylabel('Number of patients')
ax.set_title('Timepoint distribution of bladder samples')
ax.set_xticks([1, 2, 3, 4])

fig.savefig(os.path.join(DIR_figures, "bladder_sample_distribution.png"))

# Check for wbc and cfDNA concordance, does each patient have 2 cfDNA and 2 WBC FASTQ files for each time point?
for patient_id in set([re.split(r'[_]', file)[0] for file in file_list]): 
    print(patient_id)
    sample_list = return_samples_from_patient(patient_id, file_list)
    check_cfDNA_and_WBC_concordance(sample_list)
    print(" ")

# FIGURE 2
# Showing the number of baseline and OT cfDNA samples in a bar chart
input_strings = ["sample1_baseline_fDNA_2.fq.gz", "sample2_fDNA_2.fq.gz", "sample3_baseline_2.fq.gz", "sample4_baseline_fDNA_2.fq.gz"]
n_baseline_cfDNA = len(filter_strings(file_list, cfDNA=True, baseline=True)) # number of baseline cfDNA samples
n_OT_cfDNA = len(filter_strings(file_list, cfDNA=True, baseline=False)) # number of baseline cfDNA samples
n_baseline_WBC = len(filter_strings(file_list, cfDNA=False, baseline=True)) # number of baseline cfDNA samples
n_OT_WBC = len(filter_strings(file_list, cfDNA=False, baseline=False)) # number of baseline cfDNA samples

cfDNA_dict = {"Baseline cfDNA samples": n_baseline_cfDNA, "On treatment cfDNA samples": n_OT_cfDNA}
fig, ax = plt.subplots(figsize=(3, 5))
ax.bar(cfDNA_dict.keys(), cfDNA_dict.values(), color='gray', edgecolor='black', width = 0.6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('')
ax.set_ylabel('Number of samples')
ax.set_xticklabels(["Baseline cfDNA\nsamples", "On treatment cfDNA\nsamples"])
ax.set_title('Bladder cfDNA samples')

fig.tight_layout()
fig.savefig(os.path.join(DIR_figures, "baseline_and_OT_bladder_samples.png"))