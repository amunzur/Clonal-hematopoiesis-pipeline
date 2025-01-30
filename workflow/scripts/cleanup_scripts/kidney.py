"""
There are some missing columns in the KIDNEY clinical data. This script cleans it adds the following relevant columns: 
- Date of GUBB blood draw
- Age at GUBB blood draw
"""

path_raw = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Kidney/mRCC clinical Data.csv"
path_clean = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Kidney/mRCC clinical Data_clean.csv"
path_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"

sample_info = pd.read_csv(path_sample_information, sep = "\t", names = ["GUBB ID", "Date_collected", "Diagnosis", "Timepoint"])
sample_info = sample_info[sample_info["Timepoint"] == "Baseline"]
df = pd.read_csv(path_raw)
df = df.merge(sample_info[["GUBB ID", "Date_collected", "Timepoint"]], how = "inner", on = "GUBB ID")
df["DOB"] = pd.to_datetime(df["DOB"])
df["Date_collected"] = pd.to_datetime(df["Date_collected"], format = '%Y%b%d')
df["Age at GUBB draw"] = (df["Date_collected"] - df["DOB"]).dt.days//365.25

del df["Date of GUBB draw"]
del df["Timepoint"]

# df = df[df["GUBB ID"] != "22-569"]
df.loc[df["GUBB ID"] == "23-397", "irAE"] = 1

df = df.rename(columns = {"GUBB ID": "Patient_id", "Date_collected": "Date of GUBB draw"})

# add OS and PFS data 
clin_df = clin_df_kidney[["Patient_id", "Sex", "DOB", "DOD", "Date of GUBB draw", "irAE"]].rename(columns = {"GUBB ID": "Patient_id", "DOD": "Date of last follow-up or death"})
clin_df["Death at last follow up"] = ~clin_df["Date of last follow-up or death"].isin(["X", "XX"])
clin_df.loc[clin_df["Date of last follow-up or death"].isin(["X", "XX"]) , "Date of last follow-up or death"] = pd.to_datetime("2023-10-17") # thats when we had the clinical data sent to us
clin_df[['DOB', 'Date of last follow-up or death', 'Date of GUBB draw']] = clin_df[['DOB', 'Date of last follow-up or death', 'Date of GUBB draw']].apply(pd.to_datetime)
clin_df["Age at last follow-up or death"] = (clin_df["Date of last follow-up or death"] - clin_df["DOB"]).astype('<m8[Y]').astype(int)

PATH_kidney_clean = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/kidney_clin_clean.csv"
clin_df.to_csv(PATH_kidney_clean)

df.to_csv(path_clean, index = False)