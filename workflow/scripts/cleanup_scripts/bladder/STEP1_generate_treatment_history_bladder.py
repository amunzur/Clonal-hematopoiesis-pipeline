"""
This script puts the bladder clinical tables into a usable format.
"""
import pandas as pd
import os
import numpy as np

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

######################################################################################################
# TREATMENT HISTORY IN THE METASTATIC SETTING
mibc_treatmnent_cols = ["STUDY ID for Urothelial Cancer ",
                            "Event Name", 
                            "Treatment administered in [event-label] for metastatic urothelial cancer",
                            " => Chemotherapy administered in [event-label]",
                            " => Immunotherapy administered in [event-label]",
                            " => Targeted therapy administered in [event-label]",
                            " => Antibody-drug conjugate administered in [event-label]",
                            " => Combination therapy administered in [event-label]",
                            " => If other treatment, describe here",
                            " => Date of [event-label] treatment initiation",
                            " => Discontinuation date of [event-label] treatment",
                            " => Reason for discontinuation of [event-label] treatment",
                            " => Additional comments on [event-label] treatment"]

mibc_rename_cols_dict = {"STUDY ID for Urothelial Cancer ": "Patient_id",
                    "Event Name": "Treatment line",
                    "Treatment administered in [event-label] for metastatic urothelial cancer": "Treatment",
                    " => Chemotherapy administered in [event-label]": "Treatment chemo",
                    " => Immunotherapy administered in [event-label]": "Treatment IO",
                    " => Targeted therapy administered in [event-label]": "Treatment targeted",
                    " => Antibody-drug conjugate administered in [event-label]": "Treatment ADC",
                    " => Combination therapy administered in [event-label]" : "Treatment combination",
                    " => If other treatment, describe here": "Treatment other",
                    " => Date of [event-label] treatment initiation" : "Date start",
                    " => Discontinuation date of [event-label] treatment" : "Date discontinuation",
                    " => Reason for discontinuation of [event-label] treatment": "Reason discontinuation",
                    " => Additional comments on [event-label] treatment": "Comments"}

mibc_drug_rename_dict = {
    'Docetaxel alone': "Docetaxel",
    'GUTNILE Arm 2, durvalumab + tremelimumab + Chemotherapy ': 'Durvalumab + Tremelimumab + Chemotherapy',
    'Cisplatin and gemcitabine': 'Cisplatin + Gemcitabine',
    'etoposide alone': 'Etoposide',
    'Carboplatin and gemcitabine': 'Carboplatin + Gemcitabine',
    'Carboplatin, gemcitabine and Atezolizumab': "Carboplatin + Gemcitabine + Atezolizumab",
    'Paclitaxel chemotherapy - intended for breast cancer': 'Paclitaxel chemotherapy',
    'Durvalumab in Combination with Cisplatin and Gemcitabine Chemotherapy': 'Durvalumab + Cisplatin + Gemcitabine',
    'Disitamab Vedotin (RC48-ADC)': 'Disitamab vedotin',
    'Durvalumab and tremelimumab with standard-of-care chemotherapy (cisplatin/gemcitabine)': 'Durvalumab + Tremelimumab + Chemotherapy',
    'Enfortumab vedotin and pembrolizumab': 'Enfortumab vedotin + pembrolizumab',
    'Platin-based and etoposide': 'Platin-based + etoposide',
    'Carboplatin alone': 'Carboplatin',
    "Cisplatin, gemcitabine": 'Cisplatin + Gemcitabine',
    'Cisplatin, gemcitabine (GC)': 'Cisplatin + Gemcitabine',
    'Methotrexate, vinblastine, doxorubicin, cisplatin (MVAC)': 'Methotrexate + vinblastine + doxorubicin + cisplatin (MVAC)',
    "Cisplatin, gemcitabine": 'Cisplatin + Gemcitabine',
    'Carboplatin, gemcitabine': "Carboplatin + Gemcitabine",
    "Palliative / best supportive care": "Palliative / best supportive care"
}

mibc_treatments = df_main[mibc_treatmnent_cols]
mibc_treatments = mibc_treatments[mibc_treatments["Event Name"].isin(["first-line", "second-line", "third-line", "fourth-line", "fifth-line"])].reset_index(drop = True)
mibc_treatments = mibc_treatments.rename(columns = mibc_rename_cols_dict)
mibc_treatments = mibc_treatments.merge(sample_info_main["Patient_id"], on = "Patient_id", how = "inner") # subset tp pts in the cohort

# combine treatment received in one col
mibc_treatments['Combined Treatment'] = mibc_treatments[['Treatment chemo', 'Treatment IO', 'Treatment targeted', 'Treatment ADC', 'Treatment combination', 'Treatment other']].fillna('').sum(axis=1)
mask = mibc_treatments['Combined Treatment'] == ''
mibc_treatments.loc[mask, 'Combined Treatment'] = mibc_treatments.loc[mask, 'Treatment']
mibc_treatments.drop(columns=['Treatment chemo', 'Treatment IO', 'Treatment targeted', 'Treatment ADC', 'Treatment combination', 'Treatment other'], inplace=True)
mibc_treatments = mibc_treatments.drop_duplicates().reset_index(drop = True)
mibc_treatments = mibc_treatments.rename(columns = {"Combined Treatment": "Drug"})
mibc_treatments['Drug'] = mibc_treatments['Drug'].str.replace(r'Other \(describe\)', '', regex=True)
mibc_treatments["Drug_new"] = mibc_treatments["Drug"].map(mibc_drug_rename_dict)
mibc_treatments["Drug_new"].fillna(mibc_treatments["Drug"], inplace = True)

# Some detailed touches
mibc_treatments["Drug_new"] = mibc_treatments["Drug_new"].mask(mibc_treatments["Drug_new"] == "None - pt passed away after 1st infusion of EV", pd.NA)
mibc_treatments["Drug_new"] = mibc_treatments["Drug_new"].mask(mibc_treatments["Drug_new"] == "Pt passed away in the middle of tx due to an MI", pd.NA)
mibc_treatments["Drug_new"] = mibc_treatments["Drug_new"].mask(mibc_treatments["Drug_new"] == "pt currently being monitored for progression in consideration for pembro or further chemo.", pd.NA)
mibc_treatments["Drug_new"] = mibc_treatments["Drug_new"].mask(mibc_treatments["Drug_new"] == "Pt moved to Edmonton and lost to follow-up", pd.NA)
mibc_treatments["Drug_new"] = mibc_treatments["Drug_new"].mask(mibc_treatments["Drug_new"] == "Pt passed away 20-09-2021", pd.NA)

del mibc_treatments["Drug"]
mibc_treatments = mibc_treatments.rename(columns = {"Drug_new": "Drug"})
# mibc_treatments.to_csv("/groups/wyattgrp/users/amunzur/COMPOST_BIN/test.csv", index = False)

######################################################################################################
# MAINTENANCE THERAPY
maintenance_cols = ["STUDY ID for Urothelial Cancer ",
                            " => Maintenance immunotherapy administered after [event-label] treatment for metastatic urothelial cancer?",
                            " => If other maintenance treatment, describe here",
                            " => Date of initiation of maintenance immunotherapy",
                            " => Discontinuation date of maintenance immunotherapy",
                            " => Reason for discontinuation of maintenance immunotherapy",
                            " => Additional comments on maintenance immunotherapy"]

rename_cols_dict = {"STUDY ID for Urothelial Cancer ": "Patient_id",
                    " => Maintenance immunotherapy administered after [event-label] treatment for metastatic urothelial cancer?": "Maintenance IO",
                    " => If other maintenance treatment, describe here": "Other maintenance",
                    " => Date of initiation of maintenance immunotherapy": "Date start",
                    " => Discontinuation date of maintenance immunotherapy": "Date discontinuation",
                    " => Reason for discontinuation of maintenance immunotherapy": "Reason discontinuation",
                    " => Additional comments on maintenance immunotherapy": "Comments"}

maintenance_df = df_main[maintenance_cols].rename(columns = rename_cols_dict).merge(sample_info_main["Patient_id"], on = "Patient_id", how = "inner") # subset
maintenance_df = maintenance_df.drop_duplicates().reset_index(drop = True) # To avoid duplicate rows
# we further filter the pts who didn't receive maintenance therapy
maintenance_df = maintenance_df.dropna(subset=["Maintenance IO"], inplace=False)
maintenance_df = maintenance_df[maintenance_df["Maintenance IO"] != "None"].reset_index(drop=True)
maintenance_df["Other maintenance"] = maintenance_df["Other maintenance"].str.extract(r'(Atezolizumab|Durvalumab|Nivolumab|Pembrolizumab)')
maintenance_df['Maintenance IO'] = maintenance_df[['Maintenance IO', 'Other maintenance']].fillna("").sum(axis=1)
maintenance_df['Maintenance IO'] = maintenance_df['Maintenance IO'].str.replace(r'Other \(describe\)', '', regex=True)
del maintenance_df["Other maintenance"]
maintenance_df["Treatment"] = "Maintenance"
maintenance_df = maintenance_df.rename(columns = {"Maintenance IO": "Drug"})

# Classify the treatments
del maintenance_df["Treatment"]
maintenance_df["Treatment line"] = "Maintenance"
maintenance_df["Treatment"] = maintenance_df["Drug"].map({"Atezolizumab": "Immunotherapy", "Durvalumab": "Immunotherapy", "Pembrolizumab": "Immunotherapy", "Avelumab": "Immunotherapy", "Nivolumab": "Immunotherapy"})

# NEOADJUVANT AND ADJUVANT CHEMO
chemo_cols = ["STUDY ID for Urothelial Cancer ",
              "Neoadjuvant chemotherapy?",
              " => If other neoadjuvant chemotherapy regimen, describe here",
              " => Number of neoadjuvant chemotherapy cycles administered",
              "Adjuvant therapy?",
              " => If other adjuvant therapy regimen, describe here",
              " => Date of initiation of adjuvant therapy",
              " => Number of adjuvant therapy cycles administered",
              " => Date of last cycle of adjuvant therapy",
              "Comments"]

chemo_cols_rename = {"STUDY ID for Urothelial Cancer ": "Patient_id",
                     "Neoadjuvant chemotherapy?": "NAC",
                     " => If other neoadjuvant chemotherapy regimen, describe here": "Other NAC",
                     " => Number of neoadjuvant chemotherapy cycles administered": "NAC cycles",
                     "Adjuvant therapy?": "Adjuvant",
                     " => If other adjuvant therapy regimen, describe here": "Other adjuvant",
                     " => Date of initiation of adjuvant therapy": "Date start adjuvant",
                     " => Number of adjuvant therapy cycles administered": "Adjuvant cycles",
                     " => Date of last cycle of adjuvant therapy": "Date last adjuvant cycle"}

chemo = df_main[chemo_cols].rename(columns = chemo_cols_rename).drop_duplicates().reset_index(drop = True)
chemo = chemo.dropna(subset=chemo.columns.difference(['Patient_id']), how='all')
chemo = chemo.merge(sample_info_main["Patient_id"].drop_duplicates(), how = "inner", on = "Patient_id")
chemo["NAC"] = chemo["NAC"].replace("None", np.nan)
chemo["NAC"] = chemo["NAC"].replace("", np.nan)

chemo["Adjuvant"] = chemo["Adjuvant"].replace("None", np.nan)
chemo["Adjuvant"] = chemo["Adjuvant"].replace("", np.nan)

# We need the date of the surgery for patients that have a "YES" in the radical surgery column. For that we go back to the pathology report.
path_cols = ["STUDY ID for Urothelial Cancer ",
             "Date of pathology specimen collection",
             "Pathology specimen collection method",
             " => If other collection method, describe here"]

path_cols_rename = {"STUDY ID for Urothelial Cancer ": "Patient_id",
                    "Date of pathology specimen collection": "Date of path",
                    "Pathology specimen collection method": "Operation type",
                    " => If other collection method, describe here": "Operation type other"}

operations_of_interest = [
    "Cystectomy (includes cystoprostatectomy)", 
    "Metastasectomy", 
    "Nephroureterectomy", 
    "Cystectomy + left Nephroureterectomy", 
    "Partial Cystectomy and pelvic lymph node dissection"
    "partial urethrectomy",
    "Transurethral resection of prostate",
    "TURBT and TUR-prostate", 
    "TURP", 
    "TUR-P", 
    "TURP, TUR- prostate fossa tumor", 
    "Ureterectomy", 
    "Urethrectomy", 
    "Transurethral resection of bladder tumor (TURBT)", 
    "Transurethral resection of bladder tumor (TURBT)Transurethral resection, prostate", 
    "Transurethral resection of bladder tumor (TURBT)TURBT & TUR-urethra", 
    "Transurethral resection of bladder tumor (TURBT)TURBT & URS (right ureter)",
    "resection of peritoneal deposit and lymph node dissection"
    ]

renaming_operations_of_interest = {
    "Cystectomy (includes cystoprostatectomy)": "Cystoprostatectomy",
    "Cystectomy + left Nephroureterectomy": "Cystectomy, Nephroureterectomy",
    "Partial Cystectomy and pelvic lymph node dissection": "Partial cystectomy, PLND",
    "partial urethrectomy": "Partial urethrectomy",
    "Transurethral resection of prostate": "TURP",
    "TURBT and TUR-prostate": "TURBT, TURP",
    "TUR-P": "TURP",
    "TURP, TUR- prostate fossa tumor": "TURP",
    "Transurethral resection of bladder tumor (TURBT)": "TURBT",
    "Transurethral resection of bladder tumor (TURBT)Transurethral resection, prostate": "TURBT, TURP",
    "Transurethral resection of bladder tumor (TURBT)TURBT & TUR-urethra": "TURBT, TUR_urethra",
    "Transurethral resection of bladder tumor (TURBT)TURBT & URS (right ureter)": "TURBT, URS",
    "resection of peritoneal deposit and lymph node dissection": "LND"
}

path = pd.read_csv(PATH_raw_table)
path = path[path_cols].rename(columns = path_cols_rename)
path['Operation type'] = path[['Operation type', 'Operation type other']].fillna('').sum(axis=1)
path['Operation type'] = path['Operation type'].str.replace(r'Other \(describe\)', '', regex=True)
path["Operation type"] = path["Operation type"].str.replace("Excisional biopsy", "")
del path["Operation type other"]
path.replace('', np.nan, inplace=True) # replace empty strings with NA
path = path.dropna(subset=path.columns.difference(['Patient_id']), how='all')
path = path.merge(sample_info_main["Patient_id"].drop_duplicates())
path = path[path["Operation type"].isin(operations_of_interest)].reset_index(drop = True) # Only interested in these surgeries, not the pathology biopsies etc
path['Operation type']= path['Operation type'].replace(renaming_operations_of_interest) # Replace values in the "Operation type" column using the mapping
path = path[path["Operation type"].str.contains(r'cyst|resection of peritoneal deposit and lymph node dissection', case=False, regex=True)].reset_index(drop=True)
path["Surgery performed"] = True


# now combining data frames
chemo = chemo.merge(path, on = "Patient_id", how = "left")
chemo["Surgery performed"] = chemo["Surgery performed"].fillna(False)

# Now to see when the patients started having NAC we need to: 
# subtract 4 months from the date of surgery, plus 3 weeks * number of cycles
chemo["Date NAC start"] = np.nan # first all nan, then we populate
chemo["Date NAC discontinuation"] = np.nan
chemo["NAC_cycles_weeks"] = chemo["NAC cycles"] * 3

chemo['Date of path'] = pd.to_datetime(chemo['Date of path'])
chemo["Date NAC start"] = chemo.apply(lambda row: pd.NaT if pd.isna(row["Date of path"]) or pd.isna(row["NAC_cycles_weeks"]) else row["Date of path"] - pd.DateOffset(weeks=4) - pd.DateOffset(weeks=row["NAC_cycles_weeks"]), axis=1)
chemo["Date NAC discontinuation"] = chemo["Date of path"] - pd.DateOffset(weeks=4)

# Now we fill in some manually based on the comments
# 1
pt = chemo.loc[chemo["Comments"] == "Neoadjuvant chemo administered in South Korea.  note patient initial diagnosis/first line chemo done in South Korea. ", "Patient_id"].values[0]
chemo.loc[chemo["Patient_id"] == pt, "Date NAC start"] = chemo.loc[chemo["Patient_id"] == pt, "Date of path"] - pd.DateOffset(weeks=16) # 4 weeks of recovery, and 3 cycles of chemo
chemo.loc[chemo["Patient_id"] == pt, "Date NAC discontinuation"] = chemo.loc[chemo["Patient_id"] == pt, "Date of path"] - pd.DateOffset(weeks=4)

# 2
pt = chemo.loc[chemo["Comments"] == "Patient only had Day 1 of cycle 1 and did not want the rest of the treatment due to having severe fatigue and nausea.  --> 1 dose neoadjuvant cis/gem Dec 2016  -4 cycles were planned to be administered initially", "Patient_id"].values[0]
chemo.loc[chemo["Patient_id"] == pt, "Date NAC start"] = chemo.loc[chemo["Patient_id"] == pt, "Date of path"] - pd.DateOffset(weeks=16) # 4 weeks of recovery, and 3 cycles of chemo
chemo.loc[chemo["Patient_id"] == pt, "Date NAC discontinuation"] = chemo.loc[chemo["Patient_id"] == pt, "Date NAC start"] + pd.DateOffset(days=1) # pt only received dose 1 of chemo

# Now we do Adjuvant Therapy
chemo["Date discontinuation adjuvant"] = pd.to_datetime(chemo["Date last adjuvant cycle"]) + pd.DateOffset(weeks=3)

# if a pt didn't receive neither NAC nor adjuvant, drop them from the df.
mask = (chemo["NAC"].isna()) & (chemo["Date last adjuvant cycle"].isna())
chemo = chemo[~mask]

# A little bit of renaming so that it matches the previous treatment dataset
chemo.replace('', np.nan, inplace=True)
chemo.dropna(axis=1, how='all', inplace=True)

# from wide to long format
# NAC df
chemo_nac = chemo[['Patient_id', 'NAC', 'Date NAC start', 'Date NAC discontinuation', 'Comments']]
chemo_nac = chemo_nac[~pd.isnull(chemo_nac["NAC"])].reset_index(drop = True)
chemo_nac["Treatment"] = "NAC"
chemo_nac = chemo_nac.rename(columns = {"NAC": "Drug", "Date NAC start": "Date start", "Date NAC discontinuation": "Date discontinuation"})
chemo_nac["Treatment line"] = "NAC"
chemo_nac["Treatment"] = "Chemotherapy"
chemo_nac["Reason discontinuation"] = np.nan

# ADJUVANT DF
chemo_adj = chemo[['Patient_id', 'Adjuvant', 'Date start adjuvant', 'Date discontinuation adjuvant', 'Comments']]
chemo_adj = chemo_adj[~pd.isnull(chemo_adj["Adjuvant"])].reset_index(drop = True)
chemo_adj["Treatment"] = "Adjuvant"
chemo_adj = chemo_adj.rename(columns = {"Adjuvant": "Drug", "Date start adjuvant": "Date start", "Date discontinuation adjuvant": "Date discontinuation"})
chemo_adj["Treatment line"] = "Adjuvant"
chemo_adj["Treatment"] = "Chemotherapy"
chemo_adj["Reason discontinuation"] = np.nan

chemo_renaming_dict = {
    "Cisplatin, gemcitabine": 'Cisplatin + Gemcitabine',
    'Cisplatin, gemcitabine (GC)': 'Cisplatin + Gemcitabine',
    'Methotrexate, vinblastine, doxorubicin, cisplatin (MVAC)': 'Methotrexate + vinblastine + doxorubicin + cisplatin (MVAC)',
    "Cisplatin, gemcitabine": 'Cisplatin + Gemcitabine',
    'Carboplatin, gemcitabine': "Carboplatin + Gemcitabine"
}

chemo_combined = pd.concat([chemo_nac, chemo_adj]).reset_index(drop = True)
chemo_combined["Drug"] = chemo_combined["Drug"].map(chemo_renaming_dict)
treatment_df = pd.concat([mibc_treatments, maintenance_df, chemo_combined])
treatment_df["Date start"] = pd.to_datetime(treatment_df["Date start"])
treatment_df["Date discontinuation"] = pd.to_datetime(treatment_df["Date discontinuation"])
treatment_df = treatment_df[["Patient_id", "Treatment line", "Treatment", "Drug", "Date start", "Date discontinuation", "Reason discontinuation", "Comments", ]]

# Some manual additions that we confirmed with the data collectors
dict1 = {"Patient_id": "20-099", "Treatment line": "NAC", "Treatment": "Chemotherapy", "Drug": "Cisplatin + Gemcitabine", "Reason discontinuation": "Completed intended treatment", "Date discontinuation": "2021-Feb-12", "Comments": "np.nan"}
dict2 = {"Patient_id": "21-380", "Treatment line": "NAC", "Treatment": "Chemotherapy", "Drug": "Cisplatin + Gemcitabine", "Reason discontinuation": "Completed intended treatment", "Date discontinuation": "2021-Mar-31", "Comments": "np.nan"}
dict3 = {"Patient_id": "22-177", "Treatment line": "NAC", "Treatment": "Chemotherapy", "Drug": "Cisplatin + Gemcitabine", "Reason discontinuation": "Completed intended treatment", "Date discontinuation": "2021-Apr-27", "Comments": "np.nan"}
dict4 = {"Patient_id": "23-088", "Treatment line": "NAC", "Treatment": "Chemotherapy", "Drug": "Cisplatin + Gemcitabine", "Reason discontinuation": "Completed intended treatment", "Date discontinuation": "2023-Jan-06", "Comments": "np.nan"}
df = pd.DataFrame([dict1, dict2, dict3, dict4])
df["Date discontinuation"] = pd.to_datetime(df["Date discontinuation"])
df["Date start"] = df["Date discontinuation"] - pd.DateOffset(weeks=12)
df["Date start"] = pd.to_datetime(df["Date start"])

treatment_df = pd.concat([treatment_df, df]).reset_index(drop = True)
treatment_df["Date start"] = treatment_df["Date start"].dt.strftime('%Y-%m-%d')
treatment_df["Date discontinuation"] = treatment_df["Date discontinuation"].dt.strftime('%Y-%m-%d')
treatment_df = treatment_df[treatment_df["Treatment"] != "Palliative / best supportive care"].reset_index(drop = True)
del treatment_df["Comments"]
treatment_df.to_csv(os.path.join(DIR_output, "treatment.csv"), index = False)