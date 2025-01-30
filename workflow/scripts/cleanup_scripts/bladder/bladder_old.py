"""
This script puts the bladder clinical tables into a usable format.
"""
import pandas as pd
import numpy as np
import os

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
PATH_raw_table = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/UrothelialCancerData_DATA_LABELS_2024-05-24_0957.csv"
# output tables
PATH_enrollment = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Bladder_enrollment.csv"
PATH_nmibc_mibc = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/nmibc_mibc.csv"
DIR_output = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder"

sample_info_main = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])

df_main = pd.read_csv(PATH_raw_table)
df_main = df_main[~((df_main["Repeat Instrument"] == "Pathology") | (df_main["Repeat Instrument"] == "Form 1 Urothelial Cancer And Treatment (old, read-only)"))] # drop old forms no longer used. This one 
# cols_to_keep = ["Removed" not in col for col in df.columns]
# df_main = df_main.loc[:, cols_to_keep]

######################################################################################################
# Confirm mUC status
# x = df_main[["STUDY ID for Urothelial Cancer ", "Date of diagnosis of mUC"]].dropna(subset = "Date of diagnosis of mUC").reset_index(drop = True).rename(columns = {"STUDY ID for Urothelial Cancer ": "Patient_id"})
# bladder_info = sample_info_main[(sample_info_main["Diagnosis"] == "Bladder") & (sample_info_main["Timepoint"] == "Baseline")].reset_index(drop = True)
# bladder_info = bladder_info.merge(x, how = "left")
# bladder_info['Date of diagnosis of mUC'] = pd.to_datetime(bladder_info['Date of diagnosis of mUC'])
# bladder_info['Date_collected'] = pd.to_datetime(bladder_info['Date_collected'], format='%Y%b%d')

# # Check if 'Date_collected' is after 'Date of diagnosis of mUC'
# bladder_info['Date_Collected_After_Diagnosis'] = bladder_info['Date_collected'] > bladder_info['Date of diagnosis of mUC']
# before_muc = bladder_info[bladder_info['Date_Collected_After_Diagnosis'] == False]
# These will be dropped: 20-137|20-299|20-369|20-371|21-010|21-453|20-384.
# 
# ######################################################################################################
# Generate a clean enrollment form
enrollment = df_main[df_main["Event Name"] == "Enrollment"].reset_index(drop = True)
enrollment = enrollment[["STUDY ID for Urothelial Cancer ",
              "Date of birth",
              "Gender",
              "Previous smoking history",
              " => If smoking history is documented, total pack-years?",
              "Date of last follow-up or death",
              "Alive or dead at last follow-up?"]]

dict_rename = {
    "STUDY ID for Urothelial Cancer ": "Patient_id",
    "Gender": "Sex",
    " => If smoking history is documented, total pack-years?": "Total pack years",
    "Alive or dead at last follow-up?": "Death at last follow up",
    "=> Type of intravesical therapy administered for NMIBC?": "Type of intravesical therapy administered for NMIBC"
}

enrollment = enrollment.rename(columns = dict_rename)

enrollment["Death"] = enrollment["Death at last follow up"].replace({"Dead": True, "Alive": False})
enrollment.loc[enrollment["Previous smoking history"] == "Never smoker", "Total pack years"] = 0
# Calculate age at last follow up
enrollment['Date of birth'] = pd.to_datetime(enrollment['Date of birth'])
enrollment['Date of last follow-up or death'] = pd.to_datetime(enrollment['Date of last follow-up or death'])
enrollment['Age at last follow-up or death'] = (enrollment['Date of last follow-up or death'] - enrollment['Date of birth']).astype('<m8[Y]')
mask = ~enrollment['Patient_id'].astype(str).str.startswith("BC")
enrollment = enrollment[mask].reset_index(drop = True)
enrollment = enrollment.merge(sample_info_main["Patient_id"], how = "inner") # subset tp pts in the cohort
enrollment["Total pack years"] = enrollment["Total pack years"].fillna(value=pd.NA)

# Add the age at gubb blood draw
sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
sample_info = sample_info[sample_info["Timepoint"] == "Baseline"]
enrollment = enrollment.merge(sample_info, how = "left")
enrollment["Date_collected"] = pd.to_datetime(enrollment['Date_collected'], format='%Y%b%d')
enrollment["Age at baseline blood draw"] = (round((enrollment["Date_collected"] - enrollment["Date of birth"]).dt.days / 365.25)).astype(int)
enrollment = enrollment.drop_duplicates()

# add some PFS information
pfs_df = df_main[df_main["Event Name"] == "first-line"][["STUDY ID for Urothelial Cancer ", "Treatment administered in [event-label] for metastatic urothelial cancer", " => Date of [event-label] treatment initiation", " => Discontinuation date of [event-label] treatment", " => Reason for discontinuation of [event-label] treatment"]]
pfs_df.rename(columns = {
    "STUDY ID for Urothelial Cancer ": "Patient_id",
    " => Date of [event-label] treatment initiation": "DOSxStart",
    " => Discontinuation date of [event-label] treatment": "DOSxStop",
    " => Reason for discontinuation of [event-label] treatment": "Reason discontinuation", 
    "Treatment administered in [event-label] for metastatic urothelial cancer": "Type of treatment"
}, inplace = True)
pfs_df["Progression"] = pfs_df["Reason discontinuation"].isin(["Clinical progression", "Radiographic progression", "Death"])
pfs_df = pfs_df[pfs_df["Type of treatment"] != "Palliative / best supportive care"].reset_index(drop = True)

pfs_df["DOSxStop"] = pd.to_datetime(pfs_df['DOSxStop'])
pfs_df["DOSxStart"] = pd.to_datetime(pfs_df['DOSxStart'])
pfs_df["Treatment line"] = "1L"

enrollment = enrollment.merge(pfs_df, how = "left")

enrollment["OS from cfDNA collection (mo)"] = round((enrollment["Date of last follow-up or death"] - enrollment["Date_collected"]).dt.days/30, 2)
enrollment["PFS (mo)"] = round((enrollment["DOSxStop"] - enrollment["DOSxStart"]).dt.days/30, 2)

# add mets information
mets_cols = [
    "STUDY ID for Urothelial Cancer ",
    "Event Name",
    "Was baseline CT scan performed prior to [event-label] treatment initiation?",
    " => Date of baseline CT scan",
    "Bone metastasis on baseline bone scan or CT prior to [event-label] treatment initiation?",
    "Lymph node metastasis on baseline CT scan prior to [event-label] treatment initiation?",
    "Lung metastasis on baseline chest CT scan prior to [event-label] treatment initiation?",
    "Liver metastasis on baseline CT scan prior to [event-label] treatment initiation?",
    "Other sites of metastasis on baseline CT prior to [event-label] treatment initiation?",
    " => List other sites of metastasis",
    "ECOG performance status prior to [event-label] treatment initiation?"]

mets_cols_rename_dict = {
    "STUDY ID for Urothelial Cancer ": "Patient_id",
    "Event Name": "Line",
    "Was baseline CT scan performed prior to [event-label] treatment initiation?": "Baseline CT performed",
    " => Date of baseline CT scan": "Baseline CT date",
    "Bone metastasis on baseline bone scan or CT prior to [event-label] treatment initiation?": "Bone mets",
    "Lymph node metastasis on baseline CT scan prior to [event-label] treatment initiation?": "Lymph node mets",
    "Lung metastasis on baseline chest CT scan prior to [event-label] treatment initiation?": "Lung mets",
    "Liver metastasis on baseline CT scan prior to [event-label] treatment initiation?": "Liver mets",
    "Other sites of metastasis on baseline CT prior to [event-label] treatment initiation?": "Other mets",
    " => List other sites of metastasis": "Other mets site",
    "ECOG performance status prior to [event-label] treatment initiation?": "ECOG"}


mets_df = df_main[df_main["Event Name"].isin(["first-line", "second-line", "third-line", "fourth-line", "fifth-line"])][mets_cols].reset_index(drop = True)
mets_df.rename(columns = mets_cols_rename_dict, inplace = True)

mets_df["Baseline CT performed"] = mets_df["Baseline CT performed"].fillna("No")

# Now find the closest ctscan date to the date of blood draw
bladder_pts_df = sample_info_main[(sample_info_main["Diagnosis"] == "Bladder") & (sample_info_main["Timepoint"] == "Baseline")][["Patient_id", "Date_collected"]].reset_index(drop = True)
mets_df = mets_df.merge(bladder_pts_df, how = "inner")
mets_df["Baseline CT date"] = pd.to_datetime(mets_df['Baseline CT date'])
mets_df['Date_collected'] = pd.to_datetime(mets_df['Date_collected'], format='%Y%b%d')
mets_df["Draw minus CT date"] = (mets_df["Date_collected"] - mets_df["Baseline CT date"]).dt.days

# Define a custom function to get the minimum absolute value and corresponding event name
def get_min_abs_and_event(group):
    min_abs_value = group["Draw minus CT date"].abs().min()
    matching_rows = group.loc[group["Draw minus CT date"].abs() == min_abs_value]
    if not matching_rows.empty:
        event_name = matching_rows["Line"].values[0]
    else:
        event_name = None
    return pd.Series({"Min Abs Draw minus CT date": min_abs_value, "Line": event_name})

# Apply the custom function to each group
ct_dates = mets_df.groupby("Patient_id").apply(get_min_abs_and_event).reset_index()
ct_dates = ct_dates[ct_dates["Min Abs Draw minus CT date"] < 60].reset_index(drop = True) # limit to CT done within 60 days
mets_df = mets_df.merge(ct_dates, how = "inner")
mets_df = mets_df[['Patient_id', 'Baseline CT date', 'Draw minus CT date', 'Bone mets', 'Lymph node mets', 'Lung mets', 'Liver mets']]
# enrollment = enrollment.merge(mets_df, how = "left")

# organize the order of cols

######################################################################################################

######################################################################################################
# Generate a NMIBC\MIBC form
df = df_main[df_main["Event Name"] == "NMIBC/MIBC"].reset_index(drop = True)
df = df[["STUDY ID for Urothelial Cancer ", 
         "Date of initial diagnosis",
         "Disease stage at initial diagnosis",
         "=> Type of intravesical therapy administered for NMIBC?",
         " => If other intravesical therapy, specify",
         "NMIBC recurrence",
         "Radical radiotherapy (+/- chemotherapy) administered? ",
         "Radical surgery performed",
         "Neoadjuvant chemotherapy?",
         " => If other neoadjuvant chemotherapy regimen, describe here",
         " => Number of neoadjuvant chemotherapy cycles administered",
         " => Response to neoadjuvant chemotherapy noted on pathology    ",
         "Adjuvant therapy?",
         " => If other adjuvant therapy regimen, describe here",
         " => Date of initiation of adjuvant therapy",
         " => Number of adjuvant therapy cycles administered",
         "Date of diagnosis of mUC", 
         "Comments"]]

condition = df["=> Type of intravesical therapy administered for NMIBC?"] == "Other (describe)"
df.loc[condition, "=> Type of intravesical therapy administered for NMIBC?"] = df.loc[condition, " => If other intravesical therapy, specify"]
del df[" => If other intravesical therapy, specify"]
df["Type of intravesical therapy administered for NMIBC"] = df["Type of intravesical therapy administered for NMIBC"].replace({"Bacillus Calmette-Guérin (BCG)": "BCG"})

condition = df["Neoadjuvant chemotherapy?"] == "Other (describe)"
df.loc[condition, "Neoadjuvant chemotherapy?"] = df.loc[condition, " => If other neoadjuvant chemotherapy regimen, describe here"]
del df[" => If other neoadjuvant chemotherapy regimen, describe here"]

condition = df["Adjuvant therapy?"] == "Other (describe)"
df.loc[condition, "Adjuvant therapy?"] = df.loc[condition, " => If other adjuvant therapy regimen, describe here"]
del df[" => If other adjuvant therapy regimen, describe here"]

df["NMIBC recurrence"] = df["NMIBC recurrence"].fillna(pd.NA).replace({"Yes": True, "No": False})
df["Radical radiotherapy (+/- chemotherapy) administered? "] = df["Radical radiotherapy (+/- chemotherapy) administered? "].fillna(pd.NA).replace({"Yes": True, "No": False})
df["Radical surgery performed"] = df["Radical surgery performed"].fillna(pd.NA).replace({"Yes": True, "No": False})
df["Neoadjuvant chemotherapy?"] = df["Neoadjuvant chemotherapy?"].fillna(pd.NA).replace({"None": pd.NA})
df["Adjuvant therapy?"] = df["Adjuvant therapy?"].fillna(pd.NA).replace({"None": pd.NA})

df = df.rename(columns = {"STUDY ID for Urothelial Cancer ": "Patient_id",
         "=> Type of intravesical therapy administered for NMIBC?": "Type of intravesical therapy administered for NMIBC",
         "Radical radiotherapy (+/- chemotherapy) administered? ": "Radical radiotherapy (+/- chemotherapy)",
         "Radical surgery performed": "Radical surgery",
         "Neoadjuvant chemotherapy?": "Neoadjuvant chemotherapy",
         " => Number of neoadjuvant chemotherapy cycles administered": "Number of neoadjuvant chemotherapy cycles",
         " => Response to neoadjuvant chemotherapy noted on pathology    ": "Response to neoadjuvant chemotherapy",
         "Adjuvant therapy?": "Adjuvant therapy",
         " => Date of initiation of adjuvant therapy": "Date of initiation of adjuvant therapy",
         " => Number of adjuvant therapy cycles administered": "Number of adjuvant therapy cycles"})


# subset to pts of interest
mask = ~df['Patient_id'].astype(str).str.startswith("BC")
df = df[mask].reset_index(drop = True)
df = df.merge(sample_info_main["Patient_id"], how = "inner") # subset tp pts in the cohort

df = df.drop_duplicates().reset_index(drop = True)
df = df[~pd.isnull(df["Date of diagnosis of mUC"])]
df.to_csv(PATH_nmibc_mibc, index = False)

df_subset = df[['Patient_id', 'Date of initial diagnosis',
    'Disease stage at initial diagnosis', 
    'Type of intravesical therapy administered for NMIBC',
    'NMIBC recurrence',
    'Radical surgery',
    "Date of diagnosis of mUC",
]]

enrollment = enrollment.merge(df_subset, how = "left", on = "Patient_id")

# Now we tackle the variant histology information from the pathology reports.
df_main = pd.read_csv(PATH_raw_table)
path_cols = [
    "STUDY ID for Urothelial Cancer ",
    "Histology of pathology specimen",
    "Date of pathology specimen collection",
    " => Histology of pathology specimen variant component  (choice=Squamous)",
    " => Histology of pathology specimen variant component  (choice=Micropapillary)",
    " => Histology of pathology specimen variant component  (choice=Sarcomatoid)",
    " => Histology of pathology specimen variant component  (choice=Small cell)",
    " => Histology of pathology specimen variant component  (choice=Other (describe))",
    " => Histology of pathology specimen variant component  (choice=Unavailable/not applicable)",
    ]

cols_rename = {
    " => Histology of pathology specimen variant component  (choice=Squamous)" : "Squamous",
    " => Histology of pathology specimen variant component  (choice=Micropapillary)" : "Micropapillary",
    " => Histology of pathology specimen variant component  (choice=Sarcomatoid)" : "Sarcomatoid",
    " => Histology of pathology specimen variant component  (choice=Small cell)" : "Small cell",
    " => Histology of pathology specimen variant component  (choice=Other (describe))" : "Other",
    " => Histology of pathology specimen variant component  (choice=Unavailable/not applicable)" : "NA",
}

path = df_main[df_main["Repeat Instrument"] == "Pathology"][path_cols]

# Choose the closest path report to diagnosis of mets disease
path = path.merge(enrollment[["Patient_id", "Date of initial diagnosis"]].drop_duplicates(), left_on = "STUDY ID for Urothelial Cancer ", right_on = "Patient_id")
del path["STUDY ID for Urothelial Cancer "]

path['Date of initial diagnosis'] = pd.to_datetime(path['Date of initial diagnosis'])
path['Date of pathology specimen collection'] = pd.to_datetime(path['Date of pathology specimen collection'])
path["dx date minus path date"] = abs(path["Date of initial diagnosis"] - path["Date of pathology specimen collection"]).dt.days

# find the path report closest to the date of mets diagnosis
idx_list = []
for i, group in path.groupby("Patient_id"):
    idx = group["dx date minus path date"].idxmin()
    idx_list.append(idx)

path = path.iloc[idx_list, ].reset_index(drop = True)
path = path.rename(columns = cols_rename)
path = path.replace({"Checked": True, "Unchecked": False})

enrollment = enrollment.merge(path[["Patient_id", "Histology of pathology specimen"]])

cols_order = [
    'Patient_id', 
    'Date_collected',
    'Diagnosis',
    'Timepoint',
    'Date of birth', 
    'Age at baseline blood draw',
    'Sex',
    'Previous smoking history',
    'Total pack years',
    'Date of initial diagnosis',
    'Disease stage at initial diagnosis', 
    'Type of intravesical therapy administered for NMIBC',
    'NMIBC recurrence',
    'Radical surgery',
    'Date of diagnosis of mUC',
    'Type of treatment',
    'DOSxStart',
    'DOSxStop',
    'Reason discontinuation',
    'Date of last follow-up or death',
    'Death',
    'Age at last follow-up or death',
    'Progression',
    'PFS (mo)',
    'OS from cfDNA collection (mo)',
    'Bone mets',
    'Lymph node mets',
    'Lung mets',
    'Liver mets', 
    'Histology of pathology specimen']

enrollment = enrollment[cols_order]
enrollment.to_csv(PATH_enrollment, index = False)

# ######################################################################################################
# # EXTRACTING BLOOD COUNTS.
# blood_cols = ["STUDY ID for Urothelial Cancer ", 
#               "Event Name",
#               "What was the hemoglobin (Hb) value prior to [event-label] treatment initiation?",
#               "What was the absolute neutrophil count prior to [event-label] treatment initiation?",
#               "Date of baseline neutrophil value",
#               "What was the absolute lymphocyte count prior to [event-label] treatment initiation?",
#               "Date of baseline lymphocyte value",
#               "What was the serum lactate dehydrogenase (LDH) value prior to [event-label] treatment initiation?",
#               "What was the serum lactate dehydrogenase (LDH) upper limit of normal (ULN) at the laboratory where the LDH was recorded prior to [event-label] treatment initiation?",
#               "Date of baseline serum lactate dehydrogenase (LDH) value",
#               "What was the serum albumin value prior to [event-label] treatment initiation?",
#               "What was the serum albumin lower limit of normal (LLN) at the laboratory where the albumin was recorded prior to [event-label] treatment initiation?",
#               "Date of baseline serum albumin value"]

# dict_rename = {"STUDY ID for Urothelial Cancer ": "Patient_id",              
#                "What was the absolute neutrophil count prior to [event-label] treatment initiation?": "neutrophil count",
#                "Date of baseline neutrophil value": "neutrophil date",
#                "What was the absolute lymphocyte count prior to [event-label] treatment initiation?": "lymphocyte count",
#                "Date of baseline lymphocyte value": "lymphocyte date",
#                "What was the serum lactate dehydrogenase (LDH) value prior to [event-label] treatment initiation?": "LDH",
#                "What was the serum lactate dehydrogenase (LDH) upper limit of normal (ULN) at the laboratory where the LDH was recorded prior to [event-label] treatment initiation?": "LDH ULN",
#                "Date of baseline serum lactate dehydrogenase (LDH) value": "LDH date",
#                "What was the serum albumin value prior to [event-label] treatment initiation?": "albumin",
#                "What was the serum albumin lower limit of normal (LLN) at the laboratory where the albumin was recorded prior to [event-label] treatment initiation?": "albumin LLN",
#                "Date of baseline serum albumin value": "albumin date"}

# blood_df = df_main[blood_cols].rename(columns = dict_rename)
# blood_df = blood_df.merge(sample_info_main[["Patient_id", "Date_collected", "Timepoint"]], on = "Patient_id", how = "inner") # subset tp pts in the cohort

# blood_df_dates = blood_df[['Patient_id', 'Event Name', 'Timepoint', 'neutrophil date', 'lymphocyte date', 'LDH date', 'albumin date', 'Date_collected']]
# blood_df_dates = blood_df_dates.dropna(subset=['neutrophil date', 'lymphocyte date', 'LDH date', 'albumin date']).reset_index(drop = True) # Drop rows with NaN values in specified columns

# # Convert to date time format
# date_cols = ['neutrophil date', 'lymphocyte date', 'LDH date', 'albumin date']
# blood_df_dates[date_cols] = blood_df_dates[date_cols].apply(pd.to_datetime)
# blood_df_dates['Date_collected'] = pd.to_datetime(blood_df_dates['Date_collected'], format='%Y%b%d')

# # Find the time difference between the date of blood draw for GUBB and the blood counts dates
# blood_df_dates["neutrophil date difference"] = (blood_df_dates["Date_collected"] - blood_df_dates["neutrophil date"]).abs()
# blood_df_dates["lymphocyte date difference"] = (blood_df_dates["Date_collected"] - blood_df_dates["lymphocyte date"]).abs()
# blood_df_dates["LDH date difference"] = (blood_df_dates["Date_collected"] - blood_df_dates["LDH date"]).abs()
# blood_df_dates["albumin date difference"] = (blood_df_dates["Date_collected"] - blood_df_dates["albumin date"]).abs()

# # Generate a df that df contains the blood count dates closest to the GUBB blood draw for each patient / timepoint combination
# results_dict = {}
# for (patient_id, Timepoint), group in blood_df_dates.groupby(["Patient_id", "Timepoint"]): 
#     counts_dict = {}
#     for count in ["neutrophil", "lymphocyte", "LDH", "albumin"]: 
#         group = group.reset_index(drop = True)
#         idx_min = group[f"{count} date difference"].idxmin()
#         count_date = group.iloc[idx_min][f"{count} date"] # date of the blood draw that is as close to the GUBB draw as possible
#         count_date_diff = group.iloc[idx_min][f"{count} date difference"] # time difference between the blood draw that is as close to the GUBB draw as possible, and the GUBB draw
#         counts_dict[f"{count} date"] = count_date
#         counts_dict[f"{count} date difference"] = count_date_diff
#     results_dict[(patient_id, Timepoint)] = counts_dict

# blood_counts_closest = pd.DataFrame(results_dict).T.reset_index().rename(columns = {"level_0" : "Patient_id", "level_1": "Timepoint"}) # this df contains the blood count dates closest to the GUBB blood draw for each patient / timepoint combination

# blood_counts_closest["neutrophil date difference"] = blood_counts_closest["neutrophil date difference"].dt.days
# blood_counts_closest["lymphocyte date difference"] = blood_counts_closest["lymphocyte date difference"].dt.days
# blood_counts_closest["albumin date difference"] = blood_counts_closest["albumin date difference"].dt.days
# blood_counts_closest["LDH date difference"] = blood_counts_closest["LDH date difference"].dt.days



# blood_counts_closest.to_csv("/groups/wyattgrp/users/amunzur/COMPOST_BIN/blood.csv", index = False)

# blood_counts_closest["neutrophil date difference"].replace("days", "")
# blood_counts_closest["lymphocyte date difference"].replace(" days", "")
# blood_counts_closest["albumin date difference"].replace(" days", "")


