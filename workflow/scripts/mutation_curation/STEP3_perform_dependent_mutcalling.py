import pandas as pd
import numpy as np
import os
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import concurrent.futures

def clean_up_dataset2(path, patient_to_exclude): 
    """
    Minor aesthetic changes to make it easier to plot.
    """
    all_vars = pd.read_csv(path)
    all_vars["VAF_t"] = all_vars["VAF_t"]*100
    all_vars["VAF_n"] = all_vars["VAF_n"]*100
    all_vars.loc[all_vars["Gene"] == "TERT", "Protein_annotation"] = "Promoter"
    all_vars.loc[(all_vars["Gene"] != "TERT") & pd.isnull(all_vars["Protein_annotation"]), "Protein_annotation"] = "Splicing"
    all_vars["Consequence"] = all_vars["Consequence"].replace({'stopgain': 'Stopgain', 'nonsynonymous_SNV': 'Missense', '.':'Splicing', 'frameshift_insertion': 'Frameshift insertion', 'frameshift_deletion': 'Frameshift Deletion', 'nonframeshift_insertion': 'Nonframeshift insertion', 'frameshift_substitution': 'Frameshift deletion', 'nonframeshift_deletion': 'Nonframeshift deletion', 'startloss': 'Startloss'})
    if len(patient_to_exclude) > 0:
        for pt in patient_to_exclude: 
            all_vars = all_vars[all_vars["Patient_id"] != pt].reset_index(drop = True)
    return(all_vars)

def run_system_grep(keyword, input_path, output_path): 
    """
    From python runs grep on a file saves the output.
    """
    os.system(f"grep {keyword} {input_path} > {output_path} &")     

def return_dependent_mut_calling_details(chrom, position, ref, alt, path_pileup):
    """
    Return details for dependent mutation calling based on the pileup file.
    """
    df = pd.read_csv(path_pileup, sep="\t", names=["chrom", "position", "ref", "Coverage", "Consensus_read_bases", "MapQ"])
    df = df[df["chrom"] == chrom]
    position = int(df["position"])
    depth = int(df["Coverage"])
    n_alt = int(df["Consensus_read_bases"].str.count(alt))
    vaf = float(n_alt / depth)
    return [chrom, position, ref, alt, n_alt, vaf, depth]

def generate_mutation_samples_list(all_vars_chip):
    df_list = []
    for patient_id, group in all_vars_chip.groupby("Patient_id"):
        unique_mutations = group[['Chrom', 'Position', 'Ref', 'Alt', 'Gene', 'Protein_annotation', 'Consequence', 'Patient_id']].drop_duplicates()
        unique_mutations["Mutation_ID"] = range(len(unique_mutations))  # Create a dictionary to store sample information for each mutation
        mutation_samples = {}
        # Iterate over unique mutations
        for index, mutation in unique_mutations.iterrows():
            mutation_id = mutation['Mutation_ID']
            patient_id = mutation['Patient_id']
            # Filter 'group' DataFrame for the current mutation and patient
            mutation_rows = group[(group['Chrom'] == mutation['Chrom']) & 
                                (group['Position'] == mutation['Position']) & 
                                (group['Ref'] == mutation['Ref']) & 
                                (group['Alt'] == mutation['Alt']) & 
                                (group['Gene'] == mutation['Gene']) & 
                                (group['Protein_annotation'] == mutation['Protein_annotation']) & 
                                (group['Consequence'] == mutation['Consequence']) & 
                                (group['Patient_id'] == patient_id)]
            # Store the sample information in the dictionary
            result_df = pd.DataFrame({
                'Patient_id': patient_id,
                'Chrom': mutation['Chrom'],
                'Position': mutation['Position'],
                'Ref': mutation['Ref'],
                'Alt': mutation['Alt'],
                'Gene': mutation['Gene'],
                'Protein_annotation': mutation['Protein_annotation'],
                'Consequence': mutation['Consequence'],
                'Samples_Present': ", ".join(mutation_rows['Date_collected'].tolist()),
            }, index=[0])
            df_list.append(result_df)
    #          
    mutation_samples_df = pd.concat(df_list).reset_index(drop=True)
    return mutation_samples_df

def get_missing_samples_df(mutation_samples_df, sample_info_main):
    missing_samples_list = []
    for i, row in mutation_samples_df.iterrows():
        patient_id = row["Patient_id"]
        if patient_id == "21-408":
            pass
        all_timepoints = sample_info_main[sample_info_main["Patient_id"] == patient_id]["Date_collected"].tolist()
        if not set(row["Samples_Present"].split(", ")) == set(all_timepoints):
            row["Missing_samples"] = ", ".join(set(all_timepoints) - set(row["Samples_Present"].split(", ")))
            missing_samples_list.append(row)
    return pd.DataFrame(missing_samples_list).reset_index(drop=True)

def get_dependent_mut_calls(missing_samples_df, DIR_mpileup_SSCS2, dir_filtered_pileups, run_grep = False):
    dependent_mut_calls_list = []
    for i, row in missing_samples_df.iterrows():
        patient_id = row["Patient_id"]
        timepoints = row["Missing_samples"].split(", ")
        chrom = row["Chrom"]
        position = row["Position"]
        ref = row["Ref"]
        alt = row["Alt"]
        protein_annot = row["Protein_annotation"]
        for timepoint in timepoints:
            patterns = {
                'WBC': re.compile(f'GU-{re.escape(patient_id)}.*WBC.*{re.escape(timepoint)}.mpileup'),
                'cfDNA': re.compile(f'GU-{re.escape(patient_id)}.*cfDNA.*{re.escape(timepoint)}.mpileup')
            }
            matching_files = {pattern_name: [file for file in os.listdir(DIR_mpileup_SSCS2) if pattern.match(file)] for pattern_name, pattern in patterns.items()}
            for key in matching_files.keys():
                input_path_pileup = os.path.join(DIR_mpileup_SSCS2, matching_files[key][0]) # not filtered all there is to it
                output_path_grep = os.path.join(dir_filtered_pileups, f"{patient_id}_{timepoint}_{key}_{ref}_{position}_{protein_annot}.tsv")
                if run_grep:
                    run_system_grep(keyword = position, input_path = input_path_pileup, output_path = output_path_grep) # Runs grep in the background
                else:
                    output_list = return_dependent_mut_calling_details(chrom, position, ref, alt, output_path_grep) # dependent mut calling results
                    output_list = output_list + [row["Protein_annotation"], row["Patient_id"], timepoint, key] # add a few more details like which sample, which timepoint etc.
                    dependent_mut_calls_list.append(output_list)
    return pd.DataFrame(dependent_mut_calls_list)

def create_df_wide(dependent_mut_calls_df, all_vars_chip):
    dependent_mut_calls_df.columns = ["Chrom", "Position", "Ref", "Alt", "Alt_forward", "VAF", "Depth", "Protein_annotation", "Patient_id", "Date_collected", "Sample_type"]
    dependent_mut_calls_df["Alt_forward"] = dependent_mut_calls_df["Alt_forward"].astype(int)
    dependent_mut_calls_df = dependent_mut_calls_df.pivot_table(index=["Patient_id", "Chrom", "Position", "Ref", "Alt", "Protein_annotation", "Date_collected"], columns='Sample_type', values = ["Depth", "VAF"]).reset_index()
    dependent_mut_calls_df.columns = [f'{col}_{suffix}' if suffix else col for col, suffix in dependent_mut_calls_df.columns]
    dependent_mut_calls_df.columns = [word.replace("_cfDNA", "_t") for word in dependent_mut_calls_df.columns]
    dependent_mut_calls_df.columns = [word.replace("_WBC", "_n") for word in dependent_mut_calls_df.columns]
    dependent_mut_calls_df["Dependent"] = True
    return dependent_mut_calls_df

def add_sample_name_from_patientID_and_date(row):
    patient_id = row["Patient_id"]
    date = row["Date_collected"]
    dir_bams = "/groups/wyattgrp/users/amunzur/pipeline/results/data/bam/SSCS2_final"
    # cfDNA
    filename_pattern = re.compile(f"GU-.*{patient_id}.*cfDNA.*{date}.*.bam")
    cfDNA_file = [file for file in os.listdir(dir_bams) if filename_pattern.match(file) and file.endswith(".bam")][0].replace(".bam", "")
    # WBC
    filename_pattern = re.compile(f"GU-.*{patient_id}.*WBC.*{date}.*.bam")
    WBC_file = [file for file in os.listdir(dir_bams) if filename_pattern.match(file) and file.endswith(".bam")][0].replace(".bam", "")
    return pd.Series([cfDNA_file, WBC_file], index=["Sample_name_t", "Sample_name_n"])


def concat_and_save_to_csv(variants, df_wide, sample_info_main, output_path):
    # now we need to use the all_chip_vars df to add a few cols to the df_wide
    variants["Dependent"] = False
    # if somatic
    if "Status" in variants.columns: 
        cols_to_subset = ['Patient_id', 'Path_bam_t', 'Chrom', 'Position', 'Ref', 'Alt', 'Type',
        'Function', 'Gene', 'Consequence', 'AAchange', 'cosmic97_coding', 'avsnp150', 'CLNALLELEID',
        'CLNSIG', 'error_type', 'error_rate', 'Protein_annotation', 'Effects', 'Status', 'Path_bam_n']
        
        variants_subset = variants[cols_to_subset]
        variants_subset = variants_subset.drop_duplicates(["Patient_id", "Position", "Chrom", "Ref", "Alt", "Gene", "Protein_annotation"])
        df_wide[["Sample_name_t", "Sample_name_n"]] = df_wide.apply(add_sample_name_from_patientID_and_date, axis=1)
        df_wide["Timepoint"] = df_wide["Sample_name_t"].str.extract("(?i)(Baseline|OnTreatment)").replace("(?i)OnTreatment", "During treatment", regex=True)
        df_wide = df_wide.merge(variants_subset, how = "left", on = ["Patient_id", "Position", "Chrom", "Ref", "Alt", "Protein_annotation"])
        df_wide[["Ref_forward_t", "Alt_reverse_t", "Ref_reverse_t", "Ref_forward_n", "Alt_reverse_n", "Ref_reverse_n", "Strand_bias_fishers"]] = None
        df_wide["VAF_bg_ratio"] = df_wide["VAF_t"]/df_wide["error_rate"]
        df_wide["tumor_to_normal_VAF_ratio"] = df_wide["VAF_t"]/df_wide["VAF_n"]
        df_wide["IGV_screenshot_name"] = None
    else: 
        # if CHIP
        cols_to_subset = ['Patient_id', 'Path_bam_t', 'Chrom', 'Position', 'Ref', 'Alt', 'Type',
        'Function', 'Gene', 'Consequence', 'AAchange', 'cosmic97_coding', 'avsnp150', 'gnomad40_exome_AF_t', 'CLNALLELEID',
        'CLNSIG', 'error_type_t', 'error_rate_t', 'Protein_annotation', 'Effects', 'Status_t', 'Path_bam_n',
        'gnomad40_exome_AF_n', 'error_type_n', 'error_rate_n', 'Status_n', 
        'Consensus']
        
        variants_subset = variants[cols_to_subset].drop_duplicates(["Patient_id", "Position", "Chrom", "Ref", "Alt", "Protein_annotation"]).reset_index(drop = True)
        variants_subset = variants_subset.drop_duplicates(["Patient_id", "Position", "Chrom", "Ref", "Alt", "Protein_annotation"])
        df_wide[["Sample_name_t", "Sample_name_n"]] = df_wide.apply(add_sample_name_from_patientID_and_date, axis=1)
        df_wide["Timepoint"] = df_wide["Sample_name_t"].str.extract("(?i)(Baseline|OnTreatment)").replace("(?i)OnTreatment", "During treatment", regex=True)
        df_wide = df_wide.merge(variants_subset, how = "left", on = ["Patient_id", "Position", "Chrom", "Ref", "Alt", "Protein_annotation"])
        df_wide[["Ref_forward_t", "Alt_reverse_t", "Ref_reverse_t", "Ref_forward_n", "Alt_reverse_n", "Ref_reverse_n", "Strand_bias_fishers_n", "Strand_bias_fishers_t", 'Variant_caller', 'Mutect2', 'Vardict', 'freebayes', 'n_callers', 'IGV_screenshot_name']] = None
        df_wide["VAF_bg_ratio_t"] = df_wide["VAF_t"]/df_wide["error_rate_t"]
        df_wide["VAF_bg_ratio_n"] = df_wide["VAF_n"]/df_wide["error_rate_n"]
        df_wide["tumor_wbc_vaf_ratio"] = df_wide["VAF_t"]/df_wide["VAF_n"]
        df_wide["tumor_wbc_depth_ratio"] = df_wide["Depth_t"]/df_wide["Depth_n"]
    #    
    output = pd.concat([variants, df_wide]).reset_index(drop = True)
    del output["Diagnosis"]
    output = output.merge(sample_info_main[["Patient_id", "Diagnosis"]].drop_duplicates(), how = "left")
    output = output[output["Patient_id"] != "21-408"].reset_index(drop = True) # this one is excluded 
    output.to_csv(output_path, index = False)

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
DIR_mpileup_SSCS2 = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/mpileup/SSCS2"

chip_path = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip_SSCS2_curated.csv"
chip_output = os.path.join(DIR_working, "results/variant_calling/chip.csv")
DIR_filtered_pileups_chip = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/mpileup_filtered_SSCS2_chip_dependent_mut_calling"

somatic_path = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic_curated_SSCS2.csv"
somatic_output = os.path.join(DIR_working, "results/variant_calling/somatic.csv")
DIR_filtered_pileups_somatic = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/mpileup_filtered_SSCS2_somatic_dependent_mut_calling"

# patient_to_exclude = ["21-408", "20-265"]

# path = chip_path
# dir_filtered_pileups = DIR_filtered_pileups_chip
# path_output = chip_output

path = somatic_path
dir_filtered_pileups = DIR_filtered_pileups_somatic
path_output = somatic_output

variants = clean_up_dataset2(path, patient_to_exclude = [])
sample_info_main = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
mutation_samples_df = generate_mutation_samples_list(variants) # Single DataFrame that tells us in which samples a mutation is called
missing_samples_df = get_missing_samples_df(mutation_samples_df, sample_info_main)
# get_dependent_mut_calls(missing_samples_df, DIR_mpileup_SSCS2, dir_filtered_pileups, run_grep = False) # if you need to run grep run it like this.
dependent_mut_calls_df = get_dependent_mut_calls(missing_samples_df, DIR_mpileup_SSCS2, dir_filtered_pileups, run_grep = False) # Get the dependent calls at those positions, run grep if needed. Wait before proceeding to the next step that requires grep outputs.
df_wide = create_df_wide(dependent_mut_calls_df, variants) # slight reformatting and save to file
concat_and_save_to_csv(variants, df_wide, sample_info_main, path_output)