"""
Anonymizes sample names for submitting the bams to EGA
"""
import pandas as pd
import os
import re
import numpy as np

PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["GU_ID", "Date_collected_str", "Diagnosis", "Timepoint"])
sample_info["Date_collected"] =  pd.to_datetime(sample_info["Date_collected_str"], format = '%Y%b%d')

# Step 1: Enumerate unique patients
sample_info['Patient_num'] = sample_info['GU_ID'].rank(method='dense').astype(int)
sample_info['Patient_id'] = 'P' + sample_info['Patient_num'].astype(str)

# Step 2: Sort by `Patient_id` and `Date_collected` to get chronological order
sample_info = sample_info.sort_values(by=['Patient_id', 'Date_collected'])

# Step 3: Generate sample identifiers
sample_info['Sample_num'] = sample_info.groupby('Patient_id').cumcount() + 1
sample_info['Sample_id'] = sample_info['Patient_id'] + '_S' + sample_info['Sample_num'].astype(str)

sample_info.to_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/anonymized.tsv", sep = "\t", index = False)

# STEP 1
# Next section writes a script to move and rename the files.
dir_all_fastqs = "/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/merged/all"
all_files = [os.path.join(dir_all_fastqs, f) for f in os.listdir(dir_all_fastqs) if f.endswith(".gz")]
dir_ega = "/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/ega_submission"
script_path = "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/ega_submission/ega_copy_commands.bash"
os.remove(script_path)

for i, row in sample_info.iterrows():
    gu_id = row["GU_ID"]
    date_collected_str = row["Date_collected_str"]
    sample_id = row["Sample_id"]
    files = [f for f in all_files if gu_id in f and ("cfDNA" in f or "WBC" in f) and date_collected_str in f]
    if len(files) > 4: 
        ValueError("More than 4 files detected.")
    for f in files:
        # Replace the string before "cfDNA" or "WBC" with gu_id and after with sample_id\
        sample_type = re.search(r"(cfDNA|WBC)", f).group(0) # cfDNA or WBC
        which_read = re.search(r"(1.fq.gz|2.fq.gz)", f).group(0)
        new_name = "_".join([sample_id, sample_type, which_read])
        
        # Write to file for renaming
        path_in = os.path.join(dir_all_fastqs, f)
        path_out = os.path.join(dir_ega, new_name)
        command = f"cp -n {path_in} {path_out} &\n"
        
        with open(script_path, 'a') as file:
            file.write(command)

# STEP 2
# Ega requires the sample names to be provided in a specific format. 
ega_samples = pd.DataFrame([f.replace("_2.fq.gz", "").replace("_1.fq.gz", "") for f in os.listdir(dir_ega)])
ega_samples.columns = ["alias"]
ega_samples= ega_samples.drop_duplicates().reset_index(drop = True)

ega_samples["title"] = "Targeted DNA sequencing"
ega_samples["phenotype"] = "Cancer"
ega_samples["subject_id"] = ega_samples["alias"].str.split("_").str[0]

# We need sex data
kidney_sex = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Supplementary tables - Clinical data - mRCC.csv")[["Patient_id", "Sex"]]
bladder_sex = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv")[["Patient_id", "Sex"]].drop_duplicates()
sex_df = pd.concat([kidney_sex, bladder_sex]).reset_index(drop = True).rename(columns = {"Patient_id": "GU_ID"})
sex_df = sex_df.merge(sample_info[["GU_ID", "Patient_id"]].drop_duplicates(), how = "inner")
ega_samples = ega_samples.merge(sex_df, left_on = "subject_id", right_on = "Patient_id").drop(["GU_ID", "Patient_id"], axis = 1).rename(columns = {"Sex": "biological_sex"})
ega_samples["organism_part"] = ega_samples["alias"].str.split("_").str[-1]
ega_samples["biological_sex"] = ega_samples["biological_sex"].str.lower()

# fill out empty columns
ega_samples["description"] = np.nan
ega_samples["biosample_id"] = np.nan
ega_samples["case_control"] = np.nan
ega_samples["cell_line"] = np.nan

# Minor fix
ega_samples["organism_part"] = ega_samples["organism_part"].replace("cfDNA", "cell-free DNA")
ega_samples["organism_part"] = ega_samples["organism_part"].replace("WBC", "White blood cells")

# ega_samples["description"] = "FastQ with " + "READ" + ega_samples["alias"].str.extract(r'(_1|_2)')

# Put them into order
ega_samples = ega_samples[["alias", "title", "description", "biological_sex", "subject_id", "phenotype", "biosample_id", "case_control", "organism_part", "cell_line"]]

# Sort
ega_samples["Patient_id_numerical"] = ega_samples["subject_id"].str.replace("P", "").astype(int)
ega_samples.sort_values(by = "Patient_id_numerical", inplace = True)
del ega_samples["Patient_id_numerical"]

ega_samples.to_csv("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/ega_submission/samples.csv", index = False)

# STEP 3. Add runs information
# Ega requires the sample names to be provided in a specific format. 
ega_runs = pd.DataFrame(ega_samples["alias"].reset_index(drop=True))
ega_runs = ega_runs.rename(columns={"alias": "sample"})
ega_runs["file1"] = ega_runs["sample"] + "_1.fq.gz"
ega_runs["file2"] = ega_runs["sample"] + "_2.fq.gz"

ega_runs.to_csv("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/ega_submission/runs.csv", index = False)



















ega_samples = pd.DataFrame([f.replace("_2.fq.gz", "_2").replace("_1.fq.gz", "_1") for f in os.listdir(dir_ega)])
ega_samples.columns = ["alias"]
ega_samples["title"] = "Targeted DNA sequencing"
ega_samples["phenotype"] = "Cancer"
ega_samples["subject_id"] = ega_samples["alias"].str.split("_").str[0]

# We need sex data
kidney_sex = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Supplementary tables - Clinical data - mRCC.csv")[["Patient_id", "Sex"]]
bladder_sex = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv")[["Patient_id", "Sex"]].drop_duplicates()
sex_df = pd.concat([kidney_sex, bladder_sex]).reset_index(drop = True).rename(columns = {"Patient_id": "GU_ID"})
sex_df = sex_df.merge(sample_info[["GU_ID", "Patient_id"]].drop_duplicates(), how = "inner")
ega_samples = ega_samples.merge(sex_df, left_on = "subject_id", right_on = "Patient_id").drop(["GU_ID", "Patient_id"], axis = 1).rename(columns = {"Sex": "biological_sex"})
ega_samples["organism_part"] = ega_samples["alias"].str.split("_").str[-2]
ega_samples["biological_sex"] = ega_samples["biological_sex"].str.lower()

# fill out empty columns
ega_samples["description"] = np.nan
ega_samples["biosample_id"] = np.nan
ega_samples["case_control"] = np.nan
ega_samples["cell_line"] = np.nan

# Minor fix
ega_samples["organism_part"] = ega_samples["organism_part"].replace("cfDNA", "cell-free DNA")
ega_samples["organism_part"] = ega_samples["organism_part"].replace("WBC", "White blood cells")

ega_samples["description"] = "FastQ with " + "READ" + ega_samples["alias"].str.extract(r'(_1|_2)')

# Put them into order
ega_samples = ega_samples[["alias", "title", "description", "biological_sex", "subject_id", "phenotype", "biosample_id", "case_control", "organism_part", "cell_line"]]

# Sort
ega_samples["Patient_id_numerical"] = ega_samples["subject_id"].str.replace("P", "").astype(int)
ega_samples.sort_values(by = "Patient_id_numerical", inplace = True)
del ega_samples["Patient_id_numerical"]

ega_samples.to_csv("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/ega_submission/samples.csv", index = False)

