from scipy.stats import mannwhitneyu



"""
This script decides on which blood counts to use for the patients, and correlates them with CH status (pos and neg)
"""

# Step 1. Decide which GUBB draws can be used for blood counts analysis. They need to be never exposed to any treatment.

PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_treatment = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/treatment.csv"
PATH_date_surgery = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/operations.csv"
PATH_CHIP = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip_SSCS2_curated_complete.csv"

chemo_classification_dict = {
    "Carboplatin + Gemcitabine": "Platinum",
    "Cisplatin + Gemcitabine": "Platinum",
    "Paclitaxel chemotherapy": "Not platinum",
    "Carboplatin": "Platinum",
    "Docetaxel (taxane)": "Not platinum",
    "Platin-based + etoposide": "Platinum",
    "Etoposide": "Not platinum",
    "Methotrexate + vinblastine + doxorubicin + cisplatin (MVAC)": "Platinum"
    }

chip = pd.read_csv(PATH_CHIP)
chip_bladder = chip[(chip["Dependent"] == False) & (chip["Diagnosis"] == "Bladder")]

treatment_main = pd.read_csv(PATH_treatment)
treatment_main["Chemo_classification"] = treatment["Drug"].map(chemo_classification_dict)
treatment_main["Date start"] = pd.to_datetime(treatment_main["Date start"])
treatment_main["Date discontinuation"] = pd.to_datetime(treatment_main["Date discontinuation"])

sample_info_main = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
sample_info_main = sample_info_main[sample_info_main["Diagnosis"] == "Bladder"].reset_index(drop = True)
sample_info_main["Date_collected"] = pd.to_datetime(sample_info_main["Date_collected"], format='%Y%b%d')

merged_df = sample_info_main.merge(treatment_main, how = "left", on = "Patient_id")
merged_df["GUBB date minus date start"] = merged_df["Date_collected"] - merged_df["Date start"]
merged_df["GUBB date minus date discontinuation"] = merged_df["Date_collected"] - merged_df["Date discontinuation"]

# Initialize new columns to indicate the status of each blood sample
merged_df['Between_treatments'] = False
merged_df['Use for blood counts'] = False

results_list = []
for _, group in merged_df.groupby("Patient_id"):
    date_start_first_treatment = group["Date start"].dropna().min()  # Date the very first treatment is initiated
    group_filtered = group[group["Date_collected"] < date_start_first_treatment]  # Samples pre first therapy
    if not group_filtered.empty:
        results_list.append(group_filtered[["Patient_id", "Date_collected", "Timepoint"]].drop_duplicates())

# Concatenate the DataFrames in results_list
samples_to_use = pd.concat(results_list, ignore_index=True)

# Step 2. Explore blood counts for these patients. Things to consider:
# The blood draw should be as close as possible to the date_collected.
# The blood draw where the blood counts is coming from should be as close to the gubb blood draw as possible. 
# There needs to be at least 30 days between the post-surgery blood counts and the radical surgery  

# EXTRACTING BLOOD COUNTS FROM THE REDCAP OUTPUT
PATH_raw_table = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/UrothelialCancerData_DATA_LABELS_2023-12-12_1533.csv"

df_main = pd.read_csv(PATH_raw_table)
df_main = df_main[~((df_main["Repeat Instrument"] == "Pathology") | (df_main["Repeat Instrument"] == "Form 1 Urothelial Cancer And Treatment (old, read-only)"))] # drop old forms no longer used. This one 

blood_cols = ["STUDY ID for Urothelial Cancer ", 
              "Event Name",
              "What was the hemoglobin (Hb) value prior to [event-label] treatment initiation?",
              "Date of baseline hemoglobin (Hb) value",
              "What was the absolute neutrophil count prior to [event-label] treatment initiation?",
              "Date of baseline neutrophil value",
              "What was the absolute lymphocyte count prior to [event-label] treatment initiation?",
              "Date of baseline lymphocyte value",
              "What was the serum lactate dehydrogenase (LDH) value prior to [event-label] treatment initiation?",
              "What was the serum lactate dehydrogenase (LDH) upper limit of normal (ULN) at the laboratory where the LDH was recorded prior to [event-label] treatment initiation?",
              "Date of baseline serum lactate dehydrogenase (LDH) value",
              "What was the serum albumin value prior to [event-label] treatment initiation?",
              "What was the serum albumin lower limit of normal (LLN) at the laboratory where the albumin was recorded prior to [event-label] treatment initiation?",
              "Date of baseline serum albumin value"]

dict_rename = {"STUDY ID for Urothelial Cancer ": "Patient_id",              
               "What was the hemoglobin (Hb) value prior to [event-label] treatment initiation?": "Hb count",
               "Date of baseline hemoglobin (Hb) value": "Hb date",
               "What was the absolute neutrophil count prior to [event-label] treatment initiation?": "neutrophil count",
               "Date of baseline neutrophil value": "neutrophil date",
               "What was the absolute lymphocyte count prior to [event-label] treatment initiation?": "lymphocyte count",
               "Date of baseline lymphocyte value": "lymphocyte date",
               "What was the serum lactate dehydrogenase (LDH) value prior to [event-label] treatment initiation?": "LDH",
               "What was the serum lactate dehydrogenase (LDH) upper limit of normal (ULN) at the laboratory where the LDH was recorded prior to [event-label] treatment initiation?": "LDH ULN",
               "Date of baseline serum lactate dehydrogenase (LDH) value": "LDH date",
               "What was the serum albumin value prior to [event-label] treatment initiation?": "albumin",
               "What was the serum albumin lower limit of normal (LLN) at the laboratory where the albumin was recorded prior to [event-label] treatment initiation?": "albumin LLN",
               "Date of baseline serum albumin value": "albumin date"}

blood_df = df_main[blood_cols].rename(columns = dict_rename)
blood_df = blood_df.merge(samples_to_use, on = "Patient_id", how = "inner") # subset to pts in the cohort
# blood_df_dates.to_csv("/groups/wyattgrp/users/amunzur/COMPOST_BIN/test.csv")

blood_df_dates = blood_df[['Patient_id', 'Event Name', 'Timepoint', 'Hb date', 'neutrophil date', 'lymphocyte date', 'LDH date', 'albumin date', 'Date_collected']]
blood_df_dates = blood_df_dates.dropna(subset=['Hb date', 'neutrophil date', 'lymphocyte date', 'LDH date', 'albumin date']).reset_index(drop = True) # Drop rows with NaN values in specified columns

# Convert to date time format
date_cols = ['Hb date', 'neutrophil date', 'lymphocyte date', 'LDH date', 'albumin date']
blood_df_dates[date_cols] = blood_df_dates[date_cols].apply(pd.to_datetime)
blood_df_dates['Date_collected'] = pd.to_datetime(blood_df_dates['Date_collected'], format='%Y%b%d')

# Find the time difference between the date of blood draw for GUBB and the blood counts dates
blood_df_dates["Hb date difference"] = (blood_df_dates["Date_collected"] - blood_df_dates["Hb date"]).dt.days
blood_df_dates["neutrophil date difference"] = (blood_df_dates["Date_collected"] - blood_df_dates["neutrophil date"]).dt.days
blood_df_dates["lymphocyte date difference"] = (blood_df_dates["Date_collected"] - blood_df_dates["lymphocyte date"]).dt.days
blood_df_dates["LDH date difference"] = (blood_df_dates["Date_collected"] - blood_df_dates["LDH date"]).dt.days
blood_df_dates["albumin date difference"] = (blood_df_dates["Date_collected"] - blood_df_dates["albumin date"]).dt.days

# Also get the date of treatment for these patients, because some of them have blood counts taken after the GUBB. Need to make sure they weren't on treatment when the blood was drawn.
treatment = pd.read_csv(PATH_treatment)[["Patient_id", "Date start", "Drug"]].drop_duplicates().rename(columns = {"Date start": "Date start treatment"})
blood_df_dates = blood_df_dates.merge(treatment, how = "left", on = "Patient_id") # add treatment start dates to the blood counts df

# Now add the date of the surgeries
surgery_df = pd.read_csv(PATH_date_surgery)[["Patient_id", "Date of path"]].drop_duplicates()
blood_df_dates = blood_df_dates.merge(surgery_df, how = "left", on = "Patient_id") # add treatment start dates to the blood counts df

# Decide which blood counts to use, make sure to minimize the time between the blood counts and the GUBB draw.
my_dict_list = []
for patient_id, group in blood_df_dates.groupby("Patient_id"):
    patient_dict = {"Patient_id": patient_id}
    patient_dict["surgery date"] = date_surgery
    date_treatment_start = pd.to_datetime(group["Date start treatment"]).min()
    date_surgery = pd.to_datetime(group["Date of path"]).min()
    
    for count in ["Hb", "neutrophil", "lymphocyte", "LDH", "albumin"]:
        date_count = group[f"{count} date"].drop_duplicates()
        mask = (
            ((pd.notna(date_surgery)) & (date_count < date_treatment_start) & (date_count >= (date_surgery + timedelta(days=30))))  # Compare only when date_surgery is not null
            if pd.notna(date_surgery) else  # Only apply the condition if date_surgery is not null
            (date_count < date_treatment_start)  # If date_surgery is null, don't compare and use all date counts before treatment start
        )
        date_count_chosen = date_count[mask].min() if mask.any() else None
        patient_dict[f"{count} date"] = date_count_chosen
        
    
    my_dict_list.append(patient_dict)

dates_to_use = pd.DataFrame(my_dict_list)
sample_info_baselines = sample_info_main[sample_info_main["Timepoint"] == "Baseline"]
dates_to_use = dates_to_use.merge(sample_info_baselines, on = "Patient_id", how = "left")
dates_to_use = dates_to_use[["Patient_id", "Diagnosis", "Date_collected", "Timepoint", "surgery date", "Hb date", "neutrophil date", "lymphocyte date", "LDH date", "albumin date"]]

for count in ["Hb", "neutrophil", "lymphocyte", "LDH", "albumin"]:
    dates_to_use[f"GUBB date minus {count} date"] = dates_to_use["Date_collected"] - dates_to_use[f"{count} date"]


dates_to_use.dropna(subset = ['Hb date', 'neutrophil date', 'lymphocyte date', 'LDH date', 'albumin date'], inplace = True)
dates_to_use.to_csv("/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/blood_counts_dates_to_use.csv", index = False)

# START COLLECTING BLOOD COUNTS FROM THE DATES WE DETERMINED
dates_to_use = dates_to_use[['Patient_id', 'Hb date', 'neutrophil date', 'lymphocyte date', 'LDH date', 'albumin date']]
blood_df = blood_df[['Patient_id', 'Hb count', 'Hb date', 'neutrophil count','neutrophil date', 'lymphocyte count', 'lymphocyte date', 'LDH','LDH ULN', 'LDH date', 'albumin', 'albumin LLN', 'albumin date']]

for col in ['Hb date', 'neutrophil date', 'lymphocyte date', 'LDH date', 'albumin date']: 
    dates_to_use[col] = pd.to_datetime(dates_to_use[col])
    blood_df[col] = pd.to_datetime(blood_df[col])

blood_counts = dates_to_use.merge(blood_df, how = "left", on = ['Patient_id', 'Hb date', 'neutrophil date', 'lymphocyte date', 'LDH date', 'albumin date'])
blood_counts = blood_counts[['Patient_id', 'Hb count', 'neutrophil count', 'lymphocyte count', 'LDH', 'LDH ULN', 'albumin', 'albumin LLN']]
blood_counts["LDH ULN ratio"] = blood_counts["LDH"]/blood_counts["LDH ULN"]
blood_counts["albumin ULN ratio"] = blood_counts["albumin"]/blood_counts["albumin LLN"]

chip_status = annotate_mutation_status(chip_bladder, "Bladder", PATH_sample_information, annotate_what = "CHIP")
blood_counts = blood_counts.merge(chip_status)

# Separate the data for Positive and Negative patients
for count in ['Hb count', 'neutrophil count', 'lymphocyte count', 'LDH ULN ratio', 'albumin ULN ratio']:
    pos = blood_counts[blood_counts["CHIP status"] == "Positive"][count]
    neg = blood_counts[blood_counts["CHIP status"] == "Negative"][count]
    statistic, p_value = mannwhitneyu(pos, neg)
    age_vs_CH_presence(chip_bladder, count, PATH_sample_information, blood_counts, f"CHIP status and {count}", test_to_use = "MWU", figure_dir = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/bladder/blood_counts")

