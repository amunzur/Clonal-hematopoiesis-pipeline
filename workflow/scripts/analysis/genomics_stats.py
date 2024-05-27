"""
This script correlates mutation status (CH+ and CH-, PPM1D CH+ and CH- etc.) with various clinical data.
The main clinical data we use is the patient specific treatment landscape.

1. Prevalence of certain CH mutations and platinum chemo received 
2. Test for each gene in the panel. Any relationship with sex?


"""
import pandas as pd
import os
from statsmodels.stats.multitest import multipletests


PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_treatment = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/treatment.csv"
PATH_date_surgery = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/bladder/operations.csv"
PATH_CHIP = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip_SSCS2_curated_complete.csv"
PATH_panel_genes = "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/chip_panel_gene_categories.tsv"

PATH_kidney_clean = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/kidney_clin_clean.csv"
PATH_clinical_bladder = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Bladder_enrollment.csv"
########################################################
# 1. CH prevalence and platinum chemo
########################################################

gene = "TP53"
source_functions = "/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/visualization/UTILITIES_make_chip_plots.py"

with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

platinum_agents = [
    "Carboplatin + Gemcitabine",
    "Cisplatin + Gemcitabine",
    "Carboplatin",
    "Platin-based + etoposide",
    "Methotrexate + vinblastine + doxorubicin + cisplatin (MVAC)"
]

# Function to determine if a drug is platinum-based
def is_platinum(row):
    return any(row[agent] for agent in platinum_agents)

def convert_to_boolean(value):
    return value > 0


treatment_main = pd.read_csv(PATH_treatment)
treatment_main['Date start'] = pd.to_datetime(treatment_main['Date start'])
treatment_main['Date discontinuation'] = pd.to_datetime(treatment_main['Date discontinuation'])


sample_info_main = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
sample_info_main['Date_collected'] = pd.to_datetime(sample_info_main['Date_collected'], format='%Y%b%d')
sample_info_bladder = sample_info_main[sample_info_main["Diagnosis"] == "Bladder"].reset_index(drop = True)

# Annotate gene status for all samples
all_vars_chip = pd.read_csv(PATH_CHIP)
ch_bladder = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Dependent"] == False) & (all_vars_chip["Diagnosis"] == "Bladder")]
# chip_bladder = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Dependent"] == False) & (all_vars_chip["Diagnosis"] == "Bladder") & (all_vars_chip["VAF_n"] > 2)]

status_ch = annotate_mutation_status(ch_bladder, "Bladder", PATH_sample_information, annotate_what = "CHIP", annotate_gene = gene).merge(sample_info_main, how = "left")
status_ch = status_ch[["Patient_id", "Diagnosis", "Date_collected", "Timepoint", "CHIP status"]]
status_ch = status_ch[status_ch["Timepoint"] == "Baseline"]

# Annotate expose to certain drugs for all samples. To be true for exposure, the pt must be exposed to the drug for at least 3 weeks (one cycle)
results = []
for (patient_id, timepoint), group in sample_info_bladder.groupby(["Patient_id", "Timepoint"]):
    if timepoint == "Baseline":
        merged = group.merge(treatment_main, on = "Patient_id", how = "left")
        merged["Exposure to treatment"] = (merged["Date_collected"] - merged["Date start"]).dt.days > 14
        results.append(merged)

df = pd.DataFrame(pd.concat(results))[['Patient_id', 'Date_collected', 'Diagnosis', 'Timepoint', 'Drug','Exposure to treatment']].reset_index(drop = True)
df_wide = df.pivot_table(index=['Patient_id', 'Date_collected', 'Diagnosis', 'Timepoint'], columns='Drug', values='Exposure to treatment').reset_index()
df_final = df_wide.merge(status_ch, how = "left", on = ["Patient_id", "Diagnosis", "Date_collected", "Timepoint"]) # add gene CH status
# df_final = df_final.merge(ppm1d_status_chip, how = "left", on = ["Patient_id", "Timepoint"]) # add gene CHIP status

df_final_bool = df_final.applymap(lambda x: convert_to_boolean(x) if isinstance(x, float) else x)
df_final_bool = df_final_bool.rename(columns = {"CHIP status": f"{gene} status"})

# Prevalence of gene CH events in patients who received platinum chemo vs not. Not can mean nonplatinum chemo like docetaxel or no chemo at all
# Apply the function to each row to create a new column indicating platinum chemotherapy
df_final_bool['Received_platinum_chemo'] = df_final_bool.apply(is_platinum, axis=1)

contingency_table = pd.crosstab(df_final_bool[f'{gene} status'], df_final_bool['Received_platinum_chemo'])
odds_ratio, p_value = fisher_exact(contingency_table)
print(f'Odds Ratio: {odds_ratio}')
print(f'P-value: {p_value}')

# df_final_bool.to_csv("/groups/wyattgrp/users/amunzur/COMPOST_BIN/test.csv")

########################################################
# 2. Sex and CH prevalence
########################################################

# LOAD CHIP DATASETS
all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip_SSCS2_curated_complete.csv"))
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic_SSCS2_curated_complete.csv"))
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

# pull sex data
c1 = pd.read_csv(PATH_clinical_bladder)[["Patient_id", 'Sex']]
c2 = pd.read_csv(PATH_kidney_clean)[["Patient_id", 'Sex']]
sex_df = pd.concat([c1, c2]).reset_index(drop = True)

p_values = []
results = []

for gene in panel_genes:
    if all_vars_chip["Gene"].str.contains(gene).any():
        print(gene)
        p_value = check_mutation_prevalence_and_sex(all_vars_chip, sex_df, diagnosis, PATH_sample_information, annotate_what, annotate_gene=gene)
        if p_value is not None:
            p_values.append(p_value)
            results.append((gene, p_value))

# Apply Benjamini-Hochberg correction
reject, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

for i, (gene, p_value) in enumerate(results):
    if reject[i]:
        print(f"{gene}: Significant after correction")
        print(f"Corrected P-value: {pvals_corrected[i]}")
    else:
        print(f"{gene}: Not significant after correction")

