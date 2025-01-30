"""
For each ctDNA and CH mutation this script asks the question, just based on fragment length is it possible to distingush CH 
from normal DNA?
"""

import pandas as pd
import numpy as np
import os
import pysam
import matplotlib.gridspec as gridspec
from itertools import chain
import ast
import random
from scipy.stats import gaussian_kde

# conda activate pysam_environment
with open("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/fragmentomics/source_functions.py", 'r') as file:
    script_code = file.read()

exec(script_code)

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
DIR_output = "/groups/wyattgrp/users/amunzur/pipeline/results/fragmentomics"

sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
sample_info = sample_info[sample_info["Diagnosis"] != "Healthy"]

all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = sample_info.merge(all_vars_chip, how = "inner")
all_vars_chip = all_vars_chip[all_vars_chip["Timepoint"] == "Baseline"].reset_index(drop = True)
all_vars_chip["Var_reads"] = all_vars_chip["Alt_forward_t"] + all_vars_chip["Alt_reverse_t"]
all_vars_chip = all_vars_chip[all_vars_chip["Var_reads"] > 15].reset_index(drop = True)

all_vars_ctDNA = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
all_vars_ctDNA = sample_info.merge(all_vars_ctDNA, how = "inner")
all_vars_ctDNA = all_vars_ctDNA[all_vars_ctDNA["Timepoint"] == "Baseline"].reset_index(drop = True)
all_vars_ctDNA["Var_reads"] = all_vars_ctDNA["Alt_forward_t"] + all_vars_ctDNA["Alt_reverse_t"]
all_vars_ctDNA = all_vars_ctDNA[all_vars_ctDNA["Var_reads"] > 15].reset_index(drop = True)

#######################################################################
# CH mutations
#######################################################################
# Get the fragment lengths from CH mutations
# these lists keep info from the entire cohort. In the loop we also collect info on a per sample and per position basis.
deletion_fragments = []
insertion_fragments = []
snv_fragments = []
row_list = []
for i, row in all_vars_chip.iterrows():
    if not row["Dependent"]:
        print(i)
        chromosome = row["Chrom"]
        position = row["Position"]
        event_type = row["Type"]
        sample_name = row["Sample_name_t"]
        bam_file = os.path.join(DIR_working, "results/data/bam/SSCS2_final", sample_name + ".bam")
        if event_type == "SNV":
            ref_base = row["Ref"]
            alt_base = row["Alt"]
            ch_fragments = get_fragments_with_snv(bam_file, chromosome, position, ref_base, alt_base)
            snv_fragments.extend(ch_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position, variant_type=event_type, ref_base=ref_base, alt_base=alt_base)
        elif event_type == "Deletion":
            ch_fragments = get_fragments_with_deletion(bam_file, chromosome, position)
            deletion_fragments.extend(ch_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position, variant_type=event_type)
        else:
            ch_fragments = get_fragments_with_insertion(bam_file, chromosome, position)
            insertion_fragments.extend(ch_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position+1, variant_type=event_type)
        # Add the fragment lengths to the row
        row["CH fragments"] = ch_fragments
        row["WT fragments"] = wt_fragments
        row_list.append(row)      

df = pd.DataFrame(row_list)
df.to_csv(os.path.join(DIR_output, "CH.csv"))

#######################################################################
# ctDNA mutations
#######################################################################
# Get the fragment lengths from ctDNA mutations
ctDNA_deletion_fragments = []
ctDNA_insertion_fragments = []
ctDNA_snv_fragments = []
row_list = []
for i, row in all_vars_ctDNA.iterrows():
    if not row["Dependent"]:
        print(i)
        chromosome = row["Chrom"]
        position = row["Position"]
        event_type = row["Type"]
        sample_name = row["Sample_name_t"]
        bam_file = os.path.join(DIR_working, "results/data/bam/SSCS2_final", sample_name + ".bam")
        if row["Sample_name_t"] == "GU-23-022_cfDNA-Baseline-IDT-2023Jan11" and position == 26696588:
            position += 1
        if event_type == "SNV":
            ref_base = row["Ref"]
            alt_base = row["Alt"]
            ctDNA_fragments = get_fragments_with_snv(bam_file, chromosome, position, ref_base, alt_base)
            ctDNA_snv_fragments.extend(ctDNA_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position, variant_type=event_type, ref_base=ref_base, alt_base=alt_base)
        elif event_type == "Deletion":
            ctDNA_fragments = get_fragments_with_deletion(bam_file, chromosome, position)
            ctDNA_deletion_fragments.extend(ctDNA_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position, variant_type=event_type)
        else:
            ctDNA_fragments = get_fragments_with_insertion(bam_file, chromosome, position)
            ctDNA_insertion_fragments.extend(ctDNA_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position+1, variant_type=event_type)
        # Add the fragment lengths to the row
        row["ctDNA fragments"] = ctDNA_fragments
        row["WT fragments"] = wt_fragments
        row_list.append(row)      

df = pd.DataFrame(row_list)
df.to_csv(os.path.join(DIR_output, "ctDNA.csv"))


