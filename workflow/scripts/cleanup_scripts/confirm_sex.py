import os
import pandas as pd

DIR_cnvkit = "/groups/wyattgrp/users/amunzur/pipeline/results/cnvkit/final_result"
path_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
path_clean = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Kidney/mRCC clinical Data_clean.csv" # clinical data



sample_info = pd.read_csv(path_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
sample_info = sample_info[sample_info["Diagnosis"] == "Kidney"][["Patient_id", "Diagnosis"]].drop_duplicates()
sample_info["Patient_id"] = "GU-" + sample_info["Patient_id"] 

cnvkit_files = [os.path.join(DIR_cnvkit, file) for file in os.listdir(DIR_cnvkit) if "WBC" in file]

def get_chrX_CN(file_path):
    df = pd.read_csv(file_path, sep = "\t")
    df = df[["chromosome", "gene", "probe_start", "log2"]]
    df = df[df["chromosome"] == "chrX"]
    median_value = np.median(df["log2"])
    return(median_value)

results = {}
for file in cnvkit_files: 
    patient_id = os.path.basename(file).split("_WBC")[0]
    value = get_chrX_CN(file)
    results[patient_id] = value
    
results_df = pd.DataFrame(list(results.items()), columns=['Patient_id', 'Value'])
clinical = pd.read_csv(path_clean)[["Patient_id", "Sex"]]
clinical["Patient_id"] = "GU-" + clinical["Patient_id"]
merged = sample_info.merge(results_df, how = "left", on = "Patient_id").merge(clinical, how = "left")

merged.to_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/check_sex.csv", index = False)