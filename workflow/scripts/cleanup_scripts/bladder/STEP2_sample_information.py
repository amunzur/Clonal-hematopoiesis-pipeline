"""
Following coloumns:

Patient id
First sample (True or False)
Date collected
Number of samples this patient has
Disease stage at sample collection
"""

import pandas as pd
import numpy as np
import os
import re


def generate_sample_info_df(path_sample_info):
    """
    Patient_id
    Date_collected
    Number of samples
    First sample
    Disease state at sample collection
    """
    df = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    df = df[df["Diagnosis"] == "Bladder"]
    df["First sample?"] = False
    
    df.loc[df["Timepoint"] == "Baseline", "First sample?"] = True
    df["Disease state at sample collection"] = "mUC"
    
    del df["Timepoint"]
    del df["Diagnosis"]
    
    # number of samples for each pt
    df_n_samples = df["Patient_id"].value_counts().reset_index()
    df_n_samples.columns = ["Patient_id", "Number of samples"]
    
    df = df.merge(df_n_samples)
    df["Date_collected"] = pd.to_datetime(df['Date_collected'], format='%Y%b%d')
    
    return(df)

def generate_patient_info_df(PATH_sample_information, PATH_raw_table):
    """
    Sex
    Date of birth
    Age at GUBB draw
    Previous smoking history
    Total pack years
    Disease stage at initial diagnosis
    Date of initial diagnosis
    Date of diagnosis of mUC
    """
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    bladder_pts = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()
    
    df = pd.read_csv(PATH_raw_table)
    df = df[[
        "STUDY ID for Urothelial Cancer ", 
        "Date of birth",
        "Gender", 
        "Previous smoking history",
        " => If smoking history is documented, total pack-years?"]]
    
    df.columns = ["Patient_id", "Date of Birth", "Sex", "Previous smoking history", "Total pack years"]
    df = df.dropna(subset=["Date of Birth", "Sex", "Previous smoking history", "Total pack years"], how='all')
    
    df = df[df["Patient_id"].isin(bladder_pts)].reset_index(drop = True)
    
    # Add age at blood draw
    df["Date of Birth"] = pd.to_datetime(df["Date of Birth"])
    sample_info = sample_info[sample_info["Diagnosis"] == "Bladder"]
    sample_info["Date_collected"] = pd.to_datetime(sample_info["Date_collected"], format='%Y%b%d')
    sample_info = sample_info[["Patient_id", "Date_collected"]]
    df = df.merge(sample_info)
    df["Age at blood draw"] = ((df["Date_collected"] - df["Date of Birth"]).dt.days/365.25).astype(int)
    
    # add disease stage at initial diagnosis
    xx = pd.read_csv(PATH_raw_table)[["STUDY ID for Urothelial Cancer ", "Disease stage at initial diagnosis"]].rename(columns = {"STUDY ID for Urothelial Cancer ": "Patient_id"})
    xx = xx.dropna(subset = "Disease stage at initial diagnosis").reset_index(drop = True)
    df = df.merge(xx, how = "left")
    
    # add date of diagnosis of initial disease
    yy = pd.read_csv(PATH_raw_table)[["STUDY ID for Urothelial Cancer ", "Date of initial diagnosis"]].rename(columns = {"STUDY ID for Urothelial Cancer ": "Patient_id"})
    yy = yy.dropna(subset = "Date of initial diagnosis").reset_index(drop = True)
    df = df.merge(yy, how = "left")
    
    # add date of diagnosis of mUC
    zz = pd.read_csv(PATH_raw_table)[["STUDY ID for Urothelial Cancer ", "Date of diagnosis of mUC"]].rename(columns = {"STUDY ID for Urothelial Cancer ": "Patient_id"})
    zz = zz.dropna(subset = "Date of diagnosis of mUC").reset_index(drop = True)
    df = df.merge(zz, how = "left")
    
    return(df)

def get_OS_info(PATH_sample_information, PATH_raw_table):
    """
    Alive or dead
    Date of last follow up
    """
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    sample_info = sample_info[sample_info["Timepoint"] == "Baseline"]
    bladder_pts = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()
    
    df = pd.read_csv(PATH_raw_table)
    df = df[[
        "STUDY ID for Urothelial Cancer ", 
        "Alive or dead at last follow-up?",
        "Date of last follow-up or death"]]
    df = df.dropna(subset=["Alive or dead at last follow-up?", "Date of last follow-up or death"], how='all')
    df = df[df["STUDY ID for Urothelial Cancer "].isin(bladder_pts)].reset_index(drop = True)
    df.columns = ["Patient_id", "Death", "OS date"]
    df["Death"] = df["Death"].replace({"Alive": False, "Dead": True, "Lost to follow-up": False})
    df["OS date"] = pd.to_datetime(df["OS date"])
    
    # Add OS time (months)
    df = df.merge(sample_info[["Patient_id", "Date_collected"]])
    df["Date_collected"] = pd.to_datetime(df["Date_collected"], format='%Y%b%d')
    df["OS from cfDNA collection (mo)"] = round((df["OS date"] - df["Date_collected"]).dt.days/30, 2)
    del df["Date_collected"]
    return(df)

def get_at_diagnosis(PATH_sample_information, PATH_raw_table):
    """
    Date of initial diagnosis
    Disease stage at initial diagnosis
    Date of diagnosis of mUC
    """
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    bladder_pts = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()
    
    df = pd.read_csv(PATH_raw_table)
    
    # Date of initial diagnosis 
    initial_diagnosis_df = df[['STUDY ID for Urothelial Cancer ', 'Date of initial diagnosis']]
    initial_diagnosis_df = initial_diagnosis_df.dropna(subset=["Date of initial diagnosis"], how='all').reset_index(drop = True)
    initial_diagnosis_df["Date of initial diagnosis"] = pd.to_datetime(initial_diagnosis_df["Date of initial diagnosis"])
    initial_diagnosis_df = initial_diagnosis_df.rename(columns = {"STUDY ID for Urothelial Cancer ": "Patient_id"})
    
    # Disease stage at initial diagnosis
    disease_stage_df = df[['STUDY ID for Urothelial Cancer ', 'Disease stage at initial diagnosis']]
    disease_stage_df = disease_stage_df.dropna(subset=["Disease stage at initial diagnosis"], how='all').reset_index(drop = True)
    disease_stage_df = disease_stage_df.rename(columns = {"STUDY ID for Urothelial Cancer ": "Patient_id"})
    
    # Get the date of diagnosis of muc
    muc_df = df[["STUDY ID for Urothelial Cancer ", "Date of diagnosis of mUC"]]
    muc_df = muc_df.dropna(subset=["Date of diagnosis of mUC"], how='all').reset_index(drop = True)
    muc_df["Date of diagnosis of mUC"] = pd.to_datetime(muc_df["Date of diagnosis of mUC"])
    muc_df = muc_df.rename(columns = {"STUDY ID for Urothelial Cancer ": "Patient_id"})
    
    combined_df = initial_diagnosis_df.merge(disease_stage_df).merge(muc_df)
    
    return(combined_df)

def get_pathology_histology(PATH_sample_information, PATH_raw_table):
    """
    Histology of pathology specimen
    Date of specimen collection
    Delta initial diagnosis and specimen collection
    """
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    bladder_pts = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()
    
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
        "STUDY ID for Urothelial Cancer ": "Patient_id",
        " => Histology of pathology specimen variant component  (choice=Squamous)" : "Squamous",
        " => Histology of pathology specimen variant component  (choice=Micropapillary)" : "Micropapillary",
        " => Histology of pathology specimen variant component  (choice=Sarcomatoid)" : "Sarcomatoid",
        " => Histology of pathology specimen variant component  (choice=Small cell)" : "Small cell",
        " => Histology of pathology specimen variant component  (choice=Other (describe))" : "Other",
        " => Histology of pathology specimen variant component  (choice=Unavailable/not applicable)" : "NA",
    }
    
    path = df_main[df_main["Repeat Instrument"] == "Pathology"][path_cols]
    path.rename(columns = cols_rename, inplace = True)
    path = path[path["Patient_id"].isin(bladder_pts)]
    
    # Choose the closest path report to diagnosis of initial disease
    initial_diagnosis_df = df_main[['STUDY ID for Urothelial Cancer ', 'Date of initial diagnosis']]
    initial_diagnosis_df = initial_diagnosis_df.dropna(subset=["Date of initial diagnosis"], how='all').reset_index(drop = True)
    initial_diagnosis_df["Date of initial diagnosis"] = pd.to_datetime(initial_diagnosis_df["Date of initial diagnosis"])
    initial_diagnosis_df = initial_diagnosis_df.rename(columns = {"STUDY ID for Urothelial Cancer ": "Patient_id"})
    
    path = path.merge(initial_diagnosis_df)
    path['Date of initial diagnosis'] = pd.to_datetime(path['Date of initial diagnosis'])
    path['Date of pathology specimen collection'] = pd.to_datetime(path['Date of pathology specimen collection'])
    
    # Calculate the difference in days and convert to numeric, coercing errors to NaN
    path["Delta initial diagnosis and specimen collection (days)"] = ((path["Date of initial diagnosis"] - path["Date of pathology specimen collection"]).dt.days)
    path["Delta initial diagnosis and specimen collection (days)"] = pd.to_numeric(path["Delta initial diagnosis and specimen collection (days)"], errors='coerce')
    
    # find the path report closest to the date of mets diagnosis
    idx_list = []
    for i, group in path.groupby("Patient_id"):
        idx = abs(group["Delta initial diagnosis and specimen collection (days)"]).idxmin()
        idx_list.append(idx)    
    
    path = path.iloc[idx_list, ].reset_index(drop = True)
    path = path.rename(columns = cols_rename)
    path = path.replace({"Checked": True, "Unchecked": False})
    
    path = path[[
        'Patient_id',
        'Histology of pathology specimen',
        'Date of pathology specimen collection',
        'Delta initial diagnosis and specimen collection (days)']]
    
    return(path)

def localized_setting_info(PATH_sample_information, PATH_raw_table):
    """
    Type of intravesical therapy administered for NMIBC
    NMIBC recurrence
    Radical surgery
    """
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    bladder_pts = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "Baseline")]["Patient_id"].tolist()
    
    df_main = pd.read_csv(PATH_raw_table)
    df_main = df_main[df_main["Event Name"] == "NMIBC/MIBC"]
    
    df = df_main[[
        "STUDY ID for Urothelial Cancer ",
        "=> Type of intravesical therapy administered for NMIBC?", 
        "NMIBC recurrence", 
        "Radical surgery performed"
    ]]
    
    df = df.rename(columns = {
        "STUDY ID for Urothelial Cancer ": "Patient_id",
        "=> Type of intravesical therapy administered for NMIBC?": "Type of intravesical therapy administered for NMIBC"
    })
    
    df = df.dropna(subset=["Type of intravesical therapy administered for NMIBC", "NMIBC recurrence", "Radical surgery performed"], how='all')
    df = df[df["Patient_id"].isin(bladder_pts)].reset_index(drop = True)
    
    df["NMIBC recurrence"] = df["NMIBC recurrence"].replace({"No": False, "Yes": True})
    df["Radical surgery performed"] = df["Radical surgery performed"].replace({"No": False, "Yes": True})
    
    return(df)

def get_prior_treatment_history(PATH_sample_information, PATH_treatment_history):
    
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    sample_info = sample_info[sample_info["Diagnosis"] == "Bladder"]
    sample_info["Date_collected"] = pd.to_datetime(sample_info["Date_collected"], format='%Y%b%d')
    
    df = pd.read_csv(PATH_treatment_history)
    
    # ADD PRIOR TREATMENT HISTORY
    df = df.merge(sample_info[["Patient_id", "Date_collected"]])
    df["Date start"] = pd.to_datetime(df["Date start"])
    df["Date start"] = pd.to_datetime(df["Date start"])
    df["Date discontinuation"] = pd.to_datetime(df["Date discontinuation"])
    df["Date_collected"] = pd.to_datetime(df["Date_collected"], format='%Y%b%d')
    
    df["Treatment duration (days)"] = (df["Date discontinuation"] - df["Date start"]).dt.days
    df.loc[df["Reason discontinuation"] == "Treatment on-going", "Treatment duration (days)"] = 99999
    
    df["Delta cfDNA collection and treatment start date"] = (df["Date_collected"] - df["Date start"]).dt.days
    
    df["Exposed"] = False
    # df.loc[(df["Delta cfDNA collection and treatment start date"] > 14) & (df["Treatment duration (days)"] > 14), "Exposed"] =  True
    df.loc[df["Delta cfDNA collection and treatment start date"] > 0, "Exposed"] =  True
    
    df = df[["Patient_id", "Date_collected", "Drug", "Exposed"]]
    
    # List of drugs of interest
    drugs_of_interest = [
        "Pembrolizumab", 
        "Durvalumab", 
        "Avelumab", 
        "Docetaxel", 
        "Cisplatin", 
        "Carboplatin",
        "Gemcitabine", 
        "Erdafitinib", 
        "Enfortumab vedotin"
    ]
    
    # Create a case-insensitive regex pattern for the drugs of interest
    pattern = re.compile(r'(?i)(' + '|'.join(re.escape(drug) for drug in drugs_of_interest) + ')')
    
    # Function to extract all mentioned drugs in an entry
    def extract_drugs(drug_entry):
        return pattern.findall(drug_entry)
        
        # Apply the extraction function to each entry in the Drug column
    df['Extracted_Drug'] = df['Drug'].apply(lambda x: extract_drugs(str(x)))
    
    df = df.explode('Extracted_Drug') # Explode the DataFrame to have one drug per row
    df = df[df['Extracted_Drug'].notna()] # Remove rows with no extracted drug
    df['Extracted_Drug'] = df['Extracted_Drug'].str.lower().str.title()
    aggregated_df = df.groupby(['Patient_id', 'Date_collected', 'Extracted_Drug'], as_index=False).agg({'Exposed': 'max'}) # Aggregate duplicates by taking logical OR
    
    # Pivot the DataFrame
    reshaped_df = aggregated_df.pivot(index=['Patient_id', 'Date_collected'], columns='Extracted_Drug', values='Exposed')
    reshaped_df = reshaped_df.fillna(False)
    reshaped_df = reshaped_df.astype(bool)
    
    # For everyone else that isn't in this df add False for all drugs.
    reshaped_df = reshaped_df.reset_index()
    merged_df = sample_info[['Patient_id', 'Date_collected']].merge(reshaped_df, on=['Patient_id', 'Date_collected'], how='left').fillna(False)
    
    cols = ["Patient_id", "Date_collected", "Carboplatin", "Cisplatin", "Gemcitabine", "Docetaxel", "Pembrolizumab", "Durvalumab", "Avelumab", "Enfortumab Vedotin", "Erdafitinib"]
    merged_df = merged_df[cols]
    
    return(merged_df)

def add_treatment_history(PATH_treatment_history, PATH_sample_classification): 
    """
    Adds information about when the sample was collected during the treatment journey, such as 1L, 2L etc.
    This information is encoded in a text file. 
    Also adds the following information about whichever treatment line the sample is associated with:
    1. Event (1L, 2L...)
    2. Treatment (Immunotherapy etc)
    3. Drug (Name of the drugs given to the patient)
    4. Date start
    5. Date discontinuation
    6. Reason discontinuation
    """
    treatment = pd.read_csv(PATH_treatment_history)
    df = pd.read_csv(PATH_sample_classification, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint", "Event"])
    df["Date_collected"] = pd.to_datetime(df["Date_collected"], format='%Y%b%d')
    
    # Palliative setting will have no PFS information.
    df_palliative = df[df["Event"] == "Palliative"].reset_index(drop = True)[["Patient_id", "Date_collected", "Event"]]
    df_palliative["Event"] = "Palliative / best supportive care"
    df = df[df["Event"] != "Palliative"].reset_index(drop = True)
    
    # Helps map the treatment information from one df to the other
    my_dict = {
        '1L': 'first-line',
        '1L OT': 'first-line',
        '1L EOT': 'first-line',
        '2L': 'second-line',
        '2L OT': 'second-line',
        '2L EOT': 'second-line',
        '3L': 'third-line',
        '3L OT': 'third-line',
        '3L EOT': 'third-line',
        '4L': 'fourth-line',
        '4L OT': 'fourth-line',
        '4L EOT': 'fourth-line',
        'Maintenance EOT': 'Maintenance',
        'Maintenance': 'Maintenance',
        'Maintenance OT': 'Maintenance'
        }
    
    df["Treatment line"] = df["Event"].map(my_dict)
    df = df.merge(treatment, on = ["Patient_id", "Treatment line"])
    
    # For patients whose treatment is ongoing, add the date I pulled the data as date discontinuation, which is 2024 May 24th.
    df.loc[df["Reason discontinuation"] == "Treatment on-going", "Date discontinuation"] = "2024-05-24"
    df["Date discontinuation"] = pd.to_datetime(df["Date discontinuation"])
    df["Date start"] = pd.to_datetime(df["Date start"])
    
    # To clean up treatment names
    treatment_rename_dict = {
        'GUTNILE Arm 2, durvalumab + tremelimumab + Chemotherapy (gemcitabine and cisplatin)': 'Durvalumab + Tremelimumab + Cisplatin + Gemcitabine',
        'oral etoposide': 'Oral etoposide',
        'Platin-based + etoposide': 'Cisplatin + Etoposide',
        'Enfortumab vedotin + pembrolizumab': 'Enfortumab vedotin + Pembrolizumab',
        'Oxaliplatin, 5-Fluorouracil and Folinic  Acid (Leucovorin) (FOLFOX)': 'Oxaliplatin + 5-Fluorouracil + Folinic acid'
        }
    
    df["Drug"] = df["Drug"].replace(treatment_rename_dict)
    
    # To classify platinum - non platinum chemo
    chemo_classification_dict = {
        "Cisplatin + Gemcitabine": "Platinum chemotherapy",
        "Oral etoposide": "Non-platinum chemotherapy",
        "Carboplatin + Gemcitabine": "Platinum chemotherapy",
        "Cisplatin + Etoposide": "Platinum chemotherapy",
        "Docetaxel": "Non-platinum chemotherapy",
        "Carboplatin": "Platinum chemotherapy",
        "Oxaliplatin + 5-Fluorouracil + Folinic acid": "Platinum chemotherapy"
        }
    
    for key, value in chemo_classification_dict.items():
        df.loc[df["Drug"] == key, "Treatment"] = value
    
    df = df[['Patient_id', 'Date_collected', 'Event', 'Treatment', 'Drug', 'Date start', 'Date discontinuation', 'Reason discontinuation']]
    result = pd.concat([df, df_palliative], ignore_index=True, sort=False)
    result = result.sort_values(by = "Patient_id")
    return(result)

def get_PFS_info(treatment_df):
    """
    Assumes treatment history has been added to the df.
    Adds the following two cols:
    1. Progression
    2. PFS (mo)
    """
    treatment_df["Progression"] = treatment_df["Reason discontinuation"].isin(["Clinical progression", "Radiographic progression", "Death"])
    treatment_df["PFS (mo)"] = round((treatment_df["Date discontinuation"] - treatment_df["Date start"]).dt.days/30, 2)
        
    treatment_df.loc[treatment_df["Treatment"] == "Palliative / best supportive care", "PFS (mo)"] = np.nan 
    treatment_df.loc[treatment_df["Treatment"] == "Palliative / best supportive care", "Progression"] = np.nan
    
    treatment_df = treatment_df.reset_index(drop = True)
    
    return(treatment_df)

PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_raw_table = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/UrothelialCancerData_DATA_LABELS_2024-05-24_0957.csv"
PATH_treatment_history = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/treatment.csv"
PATH_sample_classification = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/bladder_sample_classification.tsv"
# Sample info
sample_info_df = generate_sample_info_df(PATH_sample_information)

# Patient info
patient_info_df = generate_patient_info_df(PATH_sample_information, PATH_raw_table)

# OS info
os_df = get_OS_info(PATH_sample_information, PATH_raw_table)

# Treatment history
treatment_df = add_treatment_history(PATH_treatment_history, PATH_sample_classification)

# PFS info
treatment_df = get_PFS_info(treatment_df)

# pathology
path_df = get_pathology_histology(PATH_sample_information, PATH_raw_table)

# localized setting
localized_setting_df = localized_setting_info(PATH_sample_information, PATH_raw_table)

# prior treatment
treatment_history_df = get_prior_treatment_history(PATH_sample_information, PATH_treatment_history)

# Now combine all
clinical_data = sample_info_df.merge(patient_info_df).merge(localized_setting_df).merge(path_df).merge(os_df).merge(treatment_df).merge(treatment_history_df)
clinical_data.to_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv", index = False)
