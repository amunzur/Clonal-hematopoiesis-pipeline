import numpy as np
import pandas as pd 
import os
import re

"""
After curating mutations by going through IGV snapshots, this script detects which mutations have been removed, and then adds them to the blacklist.
"""

consensus = "SSCS2"
DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")

blacklist = "/groups/wyattgrp/users/amunzur/pipeline/resources/validated_variants/blacklisted_variants.csv"

somatic_PATH_muts = f"/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic_SSCS2.csv"
somatic_DIR_curated_muts = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/IGV_snapshots/kidney_bladder_somatic_curated"
somatic_muts_to_exclude = os.path.join(DIR_working, "resources/validated_variants/bladder_kidney_somatic_to_exclude_IGV.csv")
somatic_muts_to_keep = os.path.join(DIR_working, "results/variant_calling/somatic_curated_SSCS2.csv")

chip_PATH_muts = f"/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip_SSCS2.csv"
chip_DIR_curated_muts = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/IGV_snapshots/chip_sscs2_newbg/"
chip_muts_to_exclude = os.path.join(DIR_working, "resources/validated_variants/chip_to_exclude_IGV.csv")
chip_muts_to_keep = os.path.join(DIR_working, "results/variant_calling/chip_SSCS2_curated.csv")

def add_timepoint(all_muts, PATH_sample_information):
    # Add the correct time point, some are problematic
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    del all_muts["Timepoint"]
    all_muts = all_muts.merge(sample_info, how = "inner")
    # reorder cols
    columns = list(all_muts.columns)
    timepoint_index = columns.index('Timepoint')
    date_collected_index = columns.index('Date_collected')
    columns.pop(timepoint_index)
    columns.insert(date_collected_index + 1, 'Timepoint')
    all_muts = all_muts[columns]
    return(all_muts)

def curate(PATH_muts, DIR_curated_muts, path_to_keep, path_to_exclude, PATH_sample_information): 
    
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])[["Patient_id", "Diagnosis"]].drop_duplicates()
    all_muts = pd.read_csv(PATH_muts)
    all_muts.loc[all_muts["Diagnosis"] == "Kidney During treatment", "Diagnosis"] = "Kidney"
    all_muts = add_timepoint(all_muts, PATH_sample_information)
    
    all_muts["IGV_screenshot_name"] = all_muts.apply(lambda row: "_".join(map(str, row[['Gene', 'Protein_annotation', 'Chrom', 'Position', 'Sample_name_t']])), axis=1) + ".png" # all called muts
    
    curated_muts = pd.DataFrame({"IGV_status": "keep", "IGV_screenshot_name": os.listdir(DIR_curated_muts)}) # muts we are keeping
    merged = all_muts.merge(curated_muts, how = "left")
    
    muts_to_keep = merged[merged["IGV_status"] == "keep"]
    muts_to_exclude = merged[pd.isnull(merged["IGV_status"])]
    
    # merged.loc[(merged["Gene"] == "FLT3") & (merged["VAF_n"] > 40)] = ""
    
    del muts_to_keep["IGV_status"]
    del muts_to_exclude["IGV_status"]
    
    # Only keep muts present in the patient id list
    muts_to_keep = sample_info.merge(muts_to_keep, how = "inner")
    
    muts_to_exclude.to_csv(path_to_exclude, index = False)
    muts_to_keep.to_csv(path_to_keep, index = False)

curate(PATH_muts = chip_PATH_muts, DIR_curated_muts = chip_DIR_curated_muts, path_to_keep = chip_muts_to_keep, path_to_exclude = chip_muts_to_exclude, PATH_sample_information = PATH_sample_information)
curate(PATH_muts = somatic_PATH_muts, DIR_curated_muts = somatic_DIR_curated_muts, path_to_keep = somatic_muts_to_keep, path_to_exclude = somatic_muts_to_exclude, PATH_sample_information = PATH_sample_information)

# # compare with jack's list
# curated_muts_kidney = muts_to_keep[muts_to_keep["Diagnosis"] == "Kidney"]
# jack = pd.read_csv("/groups/wyattgrp/users/jbacon/projects/CHIP/variant_calling/chip/new_intersected_list_jack.tsv", sep = "\t")
# jack["Sample_ID"] = jack["Sample_ID"].str.replace(".table", "")
# jack["Sample_ID"] = jack["Sample_ID"].str.replace("GUBB-", "")
# jack = jack.rename(columns = {"Sample_ID": "Patient_id", "AA_Change":"Protein_annotation"})
# # jack = jack.iloc[:-3]

# # merge
# combined = jack[["CHROM", "POS", "full_ID", "Patient_id", "Gene", "Protein_annotation", "VAF", "GNOMAD_AF"]].merge(curated_muts_kidney[["Patient_id", "Protein_annotation", "Gene", "VAF_t"]], indicator = True, how = "outer")
# jack_only = combined[combined["_merge"] == "left_only"].reset_index(drop = True)
# asli_only = combined[combined["_merge"] == "right_only"].reset_index(drop = True)
# jack_only.to_csv("/groups/wyattgrp/users/amunzur/COMPOST_BIN/jack_only.csv")
