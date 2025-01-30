"""
Makes the table 1 for bladder
"""

import pandas as pd
import numpy as np
import os

enrollment = pd.read_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/clinical_data.csv")
enrollment = enrollment[enrollment["First sample?"] == True]

# Sex
n_total                             = enrollment.shape[0]
n_male                              = enrollment[enrollment["Sex"] == "Male"].shape[0]
n_female                            = enrollment[enrollment["Sex"] == "Female"].shape[0]
perc_male                           = round((n_male/n_total)*100, 1)
perc_female                         = round((n_female/n_total)*100, 1)

# Stage at initial diagnosis
stage                               = enrollment["Disease stage at initial diagnosis"].value_counts(dropna = False).reset_index()
n_nmibc                             = stage[stage["index"] == "Localized NMIBC"]["Disease stage at initial diagnosis"].iloc[0]
n_mibc                              = stage[stage["index"] == "Localized MIBC"]["Disease stage at initial diagnosis"].iloc[0]
n_muc                               = stage[stage["index"] == "Metastatic"]["Disease stage at initial diagnosis"].iloc[0]
perc_nmibc                          = round(n_nmibc/n_total*100, 1)
perc_mibc                           = round(n_mibc/n_total*100, 1)
perc_muc                            = round(n_muc/n_total*100, 1)

# Variant histology at dx
hist                                = enrollment["Histology of pathology specimen"].value_counts(dropna = False).reset_index()
n_pure_urothelial                   = hist[hist["index"] == "Pure urothelial"]["Histology of pathology specimen"].iloc[0]
n_mixed                             = hist[hist["index"] == "Mixed urothelial and variant"]["Histology of pathology specimen"].iloc[0]
n_pure_variant                      = hist[hist["index"] == "Pure variant"]["Histology of pathology specimen"].iloc[0]
n_unknown                           = hist[hist["index"].isna()]["Histology of pathology specimen"].iloc[0]

perc_pure_urothelial                = round(n_pure_urothelial/n_total*100, 1)
perc_mixed                          = round(n_mixed/n_total*100, 1)
perc_pure_variant                   = round(n_pure_variant/n_total*100, 1)
perc_unknown                        = round(n_unknown/n_total*100, 1)

# Age at initial diagnosis
ages                                = enrollment[["Patient_id", "Date of Birth", "Date of initial diagnosis", "Date of diagnosis of mUC"]]
ages["Date of initial diagnosis"]   = pd.to_datetime(ages['Date of initial diagnosis'])
ages["Date of Birth"]               = pd.to_datetime(ages['Date of Birth'])
ages["Date of diagnosis of mUC"]    = pd.to_datetime(ages['Date of diagnosis of mUC'])

ages["Age at Dx"]                   = round((ages["Date of initial diagnosis"] - ages["Date of Birth"]).dt.days/365.25)
median_age_at_dx                    = np.median(ages["Age at Dx"])
min_age_dx                          = int(ages["Age at Dx"].min())
max_age_dx                          = int(ages["Age at Dx"].max())
range_age_at_dx                     = f"{min_age_dx}-{max_age_dx}"

# Age at metastatic diagnosis
ages["Age at mUC"]                  = round((ages["Date of diagnosis of mUC"] - ages["Date of Birth"]).dt.days/365.25)
median_age_at_muc                   = np.median(ages["Age at mUC"])
min_age_muc                         = int(ages["Age at mUC"].min())
max_age_muc                         = int(ages["Age at mUC"].max())
range_age_at_muc                    = f"{min_age_muc}-{max_age_muc}"

# Smoking history
smoking                             = enrollment[["Patient_id", "Previous smoking history"]]
n_never                             = smoking[smoking["Previous smoking history"] == "Never smoker"].shape[0]
n_previous                          = smoking[smoking["Previous smoking history"] == "Previous smoker"].shape[0]
n_current                           = smoking[smoking["Previous smoking history"] == "Current smoker"].shape[0]
perc_never                          = round(n_never/n_total*100, 1)
perc_previous                       = round(n_previous/n_total*100, 1)
perc_current                        = round(n_current/n_total*100, 1)

# Treatment line details
n_L1 = enrollment[enrollment["Event"] == "1L"]["Patient_id"].shape[0]
n_L2 = enrollment[enrollment["Event"] == "2L"]["Patient_id"].shape[0]
n_L3 = enrollment[enrollment["Event"] == "3L"]["Patient_id"].shape[0]
n_maintenance = enrollment[enrollment["Event"] == "Maintenance"]["Patient_id"].shape[0]
n_palliative = enrollment[enrollment["Event"] == "Palliative / best supportive care"]["Patient_id"].shape[0]
n_other = enrollment[
    (enrollment["Event"] != "Palliative / best supportive care") &
    (enrollment["Event"] != "1L") &
    (enrollment["Event"] != "2L") &
    (enrollment["Event"] != "3L") & 
    (enrollment["Event"] != "Maintenance")]["Patient_id"].shape[0]

perc_L1 = round(n_L1/115*100, 1)
perc_L2 = round(n_L2/115*100, 1)
perc_L3 = round(n_L3/115*100, 1)
perc_maintenance = round(n_maintenance/115*100, 1)
perc_palliative = round(n_palliative/115*100, 1)
perc_other = round(n_other/115*100, 1)

# Treatment history prior to baseline blood collection, n (%)
mydict = {
    "Docetaxel": "Non-platinum chemotherapy", 
    "Enfortumab Vedotin": "Antibody-drug conjugate", 
    "Erdafitinib": "Targeted therapy"
    }

tx = enrollment[["Radical surgery performed", "Carboplatin", "Cisplatin", "Docetaxel", "Pembrolizumab", "Durvalumab", "Avelumab", "Enfortumab Vedotin", "Erdafitinib"]]
tx['Immunotherapy'] = tx['Pembrolizumab'] | tx['Durvalumab'] | tx['Avelumab']
tx['Platinum chemotherapy'] = tx['Cisplatin'] | tx['Carboplatin']
tx = tx.rename(columns = mydict)

tx = tx.drop(['Carboplatin', 'Cisplatin', 'Durvalumab', 'Pembrolizumab', 'Avelumab'], axis=1)
tx = tx.sum().reset_index().rename(columns = {"index": "Treatment", 0: "Counts"})

n_radical_surgery         = tx[tx["Treatment"] == "Radical surgery performed"]["Counts"].iloc[0]
n_platinum_chemo          = tx[tx["Treatment"] == "Platinum chemotherapy"]["Counts"].iloc[0]
n_immuno                  = tx[tx["Treatment"] == "Immunotherapy"]["Counts"].iloc[0]
n_adc                     = tx[tx["Treatment"] == "Antibody-drug conjugate"]["Counts"].iloc[0]
n_nonplatinum_chemo       = tx[tx["Treatment"] == "Non-platinum chemotherapy"]["Counts"].iloc[0]
n_targeted                = tx[tx["Treatment"] == "Targeted therapy"]["Counts"].iloc[0]
n_palliative              = n_total - np.sum([n_platinum_chemo, n_immuno, n_adc, n_nonplatinum_chemo, n_targeted])

perc_radical_surgery      = round((n_radical_surgery/n_total)*100, 2)
perc_platinum_chemo       = round((n_platinum_chemo/n_total)*100, 2)
perc_immuno               = round((n_immuno/n_total)*100, 2)
perc_adc                  = round((n_adc/n_total)*100, 2)
perc_nonplatinum_chemo    = round((n_nonplatinum_chemo/n_total)*100, 2)
perc_targeted             = round((n_targeted/n_total)*100, 2)
perc_palliative           = round((n_palliative/n_total)*100, 2)

# Information about follow up samples
path_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
sample_info = pd.read_csv(path_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
sample_info = sample_info[(sample_info["Diagnosis"] == "Bladder") & (sample_info["Timepoint"] == "During treatment")]
OT_sample_counts = sample_info["Patient_id"].value_counts().reset_index(drop = True).value_counts().reset_index().rename(columns = {"index": "Number of OT samples", "Patient_id": "Number of patients"})

n_1_OT_sample = OT_sample_counts[OT_sample_counts["Number of OT samples"] == 1]["Number of patients"].iloc[0]
n_2_OT_sample = OT_sample_counts[OT_sample_counts["Number of OT samples"] == 2]["Number of patients"].iloc[0]
n_3_OT_sample = OT_sample_counts[OT_sample_counts["Number of OT samples"] == 3]["Number of patients"].iloc[0]

# Mets information
# mets = enrollment[["Patient_id", "Bone mets", "Lymph node mets", "Lung mets", "Liver mets"]]
# n_bone_mets = mets[mets["Bone mets"] == "Yes"].shape[0]
# n_LN_mets = mets[mets["Lymph node mets"] == "Yes"].shape[0]
# n_lung_mets = mets[mets["Lung mets"] == "Yes"].shape[0]
# n_liver_mets = mets[mets["Liver mets"] == "Yes"].shape[0]
# n_unknown = mets[pd.isna(mets["Lung mets"])].shape[0]

# perc_bone_mets = round(n_bone_mets/n_total*100, 1)
# perc_LN_mets = round(n_LN_mets/n_total*100, 1)
# perc_lung_mets = round(n_lung_mets/n_total*100, 1)
# perc_liver_mets = round(n_liver_mets/n_total*100, 1)
# perc_unknown = round(n_unknown/n_total*100, 1)

# Make the table
lines = [
    f"Clinical characteristics of {n_total} patients with mUC", 
    "Biological sex, n (%)", 
    f"Male\t{n_male} ({perc_male})",
    f"Female\t{n_female} ({perc_female})",
    "Stage at initial diagnosis, n (%)",
    f"Localized NMIBC\t{n_nmibc} ({perc_nmibc})",
    f"Localized MIBC\t{n_mibc} ({perc_mibc})",
    f"Metastatic\t{n_muc} ({perc_muc})",
    f"Histology of the pathology specimen, n (%)",
    f"Pure urothelial\t{n_pure_urothelial} ({perc_pure_urothelial})",
    f"Mixed urothelial and variant\t{n_mixed} ({perc_mixed})",
    f"Pure variant\t{n_pure_variant} ({perc_pure_variant})",
    f"Unknown\t{n_unknown} ({perc_unknown})",
    f"Age at initial diagnosis, median (range)\t{median_age_at_dx} ({range_age_at_dx})",
    f"Age at metastatic diagnosis, median (range)\t{median_age_at_muc} ({range_age_at_muc})",
    "Smoking history, n (%)",
    f"Never smoker\t{n_never} ({perc_never})",
    f"Previous smoker\t{n_previous} ({perc_previous})",
    f"Current smoker\t{n_current} ({perc_current})",
    "Treatment line at time of baseline blood collection, n (%)"
    f"First-line\t{n_L1} ({perc_L1})",
    f"Second-line\t{n_L2} ({perc_L2})",
    f"Third-line\t{n_L3} ({perc_L3})",
    f"Maintenance\t{n_maintenance} ({perc_maintenance})",
    f"Palliative\t{n_palliative} ({perc_palliative})",
    f"Other\t{n_other} ({perc_other})",
    "Treatment history prior to baseline blood collection, n (%)",
    f"Cystectomy\t{n_radical_surgery} ({perc_radical_surgery})", 
    f"Platinum chemotherapy\t{n_platinum_chemo} ({perc_platinum_chemo})", 
    f"Immunotherapy\t{n_immuno} ({perc_immuno})", 
    f"Non-platinum chemotherapy\t{n_nonplatinum_chemo} ({perc_nonplatinum_chemo})",
    f"Targeted therapy\t{n_targeted} ({perc_targeted})",
    f"Antibody-drug conjugate\t{n_adc} ({perc_adc})",
    f"None/palliative\t{n_palliative} ({perc_palliative})",
    "Follow up samples (n)",
    f"1\t{n_1_OT_sample}", 
    f"2\t{n_2_OT_sample}",
    f"3\t{n_3_OT_sample}"
]

os.remove('/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/bladder_clinical_table.tsv')
# Open a file in write mode
with open('/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/bladder_clinical_table.tsv', 'w') as file:
    # Write each line to the file
    for line in lines:
        file.write(line + '\n')


