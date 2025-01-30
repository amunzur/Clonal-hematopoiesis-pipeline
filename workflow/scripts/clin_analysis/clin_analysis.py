import numpy as np 
import pandas as pd
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import statsmodels.api as sm
import matplotlib.patches as mpatches
from scipy.stats import kruskal
import scipy.stats as stats
from scipy.stats import linregress
from scipy.stats import fisher_exact
from sksurv.nonparametric import kaplan_meier_estimator

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix


with open(source_functions, 'r') as file:
    script_code = file.read()

exec(script_code)

run_fishers_exact(data, "CH status", col)

def run_fishers_exact(df, chip_variable, clinical_variable1): 
    """
    Given two variables in a df, make a contingency table and run fishers exact on it. 
    """
    df_per_pt = data[["Patient_id", chip_variable]]
    contingency_table = pd.crosstab(df_per_pt[chip_variable], df[clinical_variable1])
        
    if (contingency_table <= 2).any().any():
        print("Fisher's exact test is not performed due to cells with counts less than or equal to 2.")
        return None
    
    if contingency_table.shape[0] < 2 or contingency_table.shape[1] < 2:
        print("Fisher's exact test is not performed due to the contingency table being less than 2x2.")
        return None
    
    odds_ratio, p_value = fisher_exact(contingency_table)
    return[contingency_table, odds_ratio, p_value]

def run_mwu(df, categorical_variable, numerical_variable):
    g1_value = df[categorical_variable].unique()[0]
    g2_value = df[categorical_variable].unique()[1]
    g1 = df[df[categorical_variable] == g1_value][numerical_variable]
    g2 = df[df[categorical_variable] == g2_value][numerical_variable]
    stat, p_value = mannwhitneyu(g1, g2)
    return(p_value)

PATH_mutations = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip.csv"
PATH_clin_data = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/Supplementary tables - Clinical data - mRCC.tsv"
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"

sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Diagnosis", "Date_collected", "Timepoint"])
dict = pd.read_csv(PATH_clin_dict, "\t")
muts_main = pd.read_csv(PATH_mutations)
clin_main = pd.read_csv(PATH_clin_data, sep = "\t").rename(columns = {"GUBB ID": "Patient_id"})

clin = clin_main[["Patient_id", "irAE", "Age at GUBB draw"]]
clin = clin[~clin["irAE"].isna()].reset_index(drop = True)

ch = muts_main[
    (muts_main["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])) & 
    (muts_main["Timepoint"] == "Baseline") & 
    (muts_main['Diagnosis'] == 'Kidney') & 
    (muts_main["Dependent"] == False)].reset_index(drop = True)

# ch = ch.merge(sample_info["Patient_id"], how = "inner")

median_values = ch.groupby('Patient_id')['VAF_n'].median().reset_index()
n_mutations = ch.groupby('Patient_id')['VAF_n'].count().reset_index().rename(columns = {"VAF_n": "n_mutations"})
ch_status = annotate_mutation_status(ch, "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = False, drop_dependent = True).rename(columns = {"CHIP status": "CH status"})
chip_status = annotate_mutation_status(ch[ch["VAF_n"] >= 2].reset_index(drop = True), "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = False, drop_dependent = True)
ch10_status = annotate_mutation_status(ch[ch["VAF_n"] >= 10].reset_index(drop = True), "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = False, drop_dependent = True).rename(columns = {"CHIP status": "CH10 status"})
ch20_status = annotate_mutation_status(ch[ch["VAF_n"] >= 20].reset_index(drop = True), "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = False, drop_dependent = True).rename(columns = {"CHIP status": "CH20 status"})

ch_status = ch_status[ch_status["Timepoint"] == "Baseline"][["Patient_id", "CH status"]]
chip_status = chip_status[chip_status["Timepoint"] == "Baseline"][["Patient_id", "CHIP status"]]
ch10_status = ch10_status[ch10_status["Timepoint"] == "Baseline"][["Patient_id", "CH10 status"]]
ch20_status = ch20_status[ch20_status["Timepoint"] == "Baseline"][["Patient_id", "CH20 status"]]

data = clin.merge(ch_status, how = "left").merge(chip_status, how = "left").merge(ch10_status, how = "left").merge(ch20_status, how = "left").merge(median_values, how = "left").merge(n_mutations, how = "left")
data["VAF_n"].fillna(0, inplace = True)
data["n_mutations"].fillna(0, inplace = True)

for colname in ["CH status", "CHIP status", "CH10 status", "CH20 status"]:
    data[colname] = data[colname].map({"Negative": 0, "Positive": 1})

# data['CH_status_Charlson_score_interaction'] = data['CH status'] * data['Charlson_score']

data["irAe_binary"] = data["irAE"].map({True: 1, False: 0})
data["CH TET2 status binary"] = data["CH TET2 status"].map({"Negative": 0, "Positive": 1})

X = data[["CH TET2 status binary"]]  # Predictor variables
X_incl_const = sm.add_constant(X)  # Add constant for intercept
y = data['irAe_binary']  # Outcome variable

# Fit logistic regression model with interaction term
logit_model = sm.Logit(y, X_incl_const)
logit_result = logit_model.fit()
print(logit_result.summary())

# Plot log regression
ax.scatter(X["CH status"], logit_result.predict(sm.add_constant(X)))
fig.savefig("/groups/wyattgrp/users/amunzur/COMPOST_BIN/test2.png")



# TRYING RANDOM FOREST
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)

rf_classifier.fit(X_train, y_train) # Train the Random Forest classifier on the training data
y_pred = rf_classifier.predict(X_test) # Make predictions on the testing data

accuracy = accuracy_score(y_test, y_pred) # Evaluate the model performance
print("Accuracy:", accuracy)

print("\nClassification Report:") # Print classification report and confusion matrix
print(classification_report(y_test, y_pred))

print("\nConfusion Matrix:")
print(confusion_matrix(y_test, y_pred))

feature_importances = rf_classifier.feature_importances_
for i, importance in enumerate(feature_importances):
    print(f"Feature {i+1}: {importance}")

# muts = muts[['Patient_id', 'Sample_name_t', 'Date_collected', 'Timepoint', 'Diagnosis', 'Chrom', 'Gene', 'Consequence', 'Protein_annotation', 'VAF_n', 'Dependent']]
# ch_status = annotate_mutation_status(muts, "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = False, drop_dependent = True)
# data = clin.merge(ch_status, how="left", on = ["Patient_id"])
# data = data[(data["Timepoint"] == "Baseline") & ~pd.isnull(data["irAE"])][["Patient_id", "irAE", "CHIP status"]]
data["irAE bool"] = data["irAE"] > 0


# data["CHIP"] = data["VAF_n"] >= 2
# data["CH"] = ~pd.isnull(data["VAF_n"])
# data["CHIP_10"] = data["VAF_n"] >= 10
# data["CHIP_20"] = data["VAF_n"] >= 20


#=============================================
# IMMUNE ADVERSE EVENTS
# Does the prevalence of CHIP differ in patients with and without immune adverse events?
col = "irAE"
clin = clin_main[["Patient_id", "irAE", "Age at GUBB draw", "Charlson_score"]]
clin = clin[~pd.isnull(clin[col])]
clin[col] = clin[col].astype(int) > 0

# Clean up muts of interest
ch = muts_main[(muts_main["Timepoint"] == "Baseline") & (muts_main['Diagnosis'] == 'Kidney') & (muts_main['n_callers'] != 1) & (~muts_main['Patient_id'].str.startswith('VIP_'))].reset_index(drop = True)
ch = ch.merge(sample_info["Patient_id"], how = "inner")
median_values = ch.groupby('Patient_id')['VAF_n'].median().reset_index()
n_mutations = ch.groupby('Patient_id')['VAF_n'].count().reset_index().rename(columns = {"VAF_n": "n_mutations"})

median_values_dta = ch[ch["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])].groupby('Patient_id')['VAF_n'].median().reset_index().rename(columns = {"VAF_n": "DTA_VAF_n"})
n_mutations_dta = ch[ch["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])].groupby('Patient_id')['VAF_n'].count().reset_index().rename(columns = {"VAF_n": "DTA_n_mutations"})

median_values_tet2 = ch[ch["Gene"] == "TET2"].groupby('Patient_id')['VAF_n'].median().reset_index().rename(columns = {"VAF_n": "TET2_VAF_n"})
n_mutations_tet2 = ch[ch["Gene"] == "TET2"].groupby('Patient_id')['VAF_n'].count().reset_index().rename(columns = {"VAF_n": "TET2_n_mutations"})

# Now annotate CH status based on various criteria
ch_status = annotate_mutation_status(ch, "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = False, drop_dependent = True).rename(columns = {"CHIP status": "CH status"})
chip_status = annotate_mutation_status(ch[ch["VAF_n"] >= 2].reset_index(drop = True), "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = False, drop_dependent = True)
ch_dta_status = annotate_mutation_status(ch, "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = ["DNMT3A", "TET2", "ASXL1"], drop_dependent = True).rename(columns = {"CHIP status": "CH DTA status"})
chip_dta_status = annotate_mutation_status(ch[ch["VAF_n"] >= 2].reset_index(drop = True), "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = ["DNMT3A", "TET2", "ASXL1"], drop_dependent = True).rename(columns = {"CHIP status": "CHIP DTA status"})
ch_dnmt3a_status = annotate_mutation_status(ch, "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = "DNMT3A", drop_dependent = True).rename(columns = {"CHIP status": "CH DNMT3A status"})
ch_tet2_status = annotate_mutation_status(ch, "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = "TET2", drop_dependent = True).rename(columns = {"CHIP status": "CH TET2 status"})
ch_asxl1_status = annotate_mutation_status(ch, "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = "ASXL1", drop_dependent = True).rename(columns = {"CHIP status": "CH ASXL1 status"})
ch_tp53_status = annotate_mutation_status(ch, "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = "TP53", drop_dependent = True).rename(columns = {"CHIP status": "CH TP53 status"})

ch_status = ch_status[ch_status["Timepoint"] == "Baseline"][["Patient_id", "CH status"]]
chip_status = chip_status[chip_status["Timepoint"] == "Baseline"][["Patient_id", "CHIP status"]]
ch_dta_status = ch_dta_status[ch_dta_status["Timepoint"] == "Baseline"][["Patient_id", "CH DTA status"]]
chip_dta_status = chip_dta_status[chip_dta_status["Timepoint"] == "Baseline"][["Patient_id", "CHIP DTA status"]]
ch_dnmt3a_status = ch_dnmt3a_status[ch_dnmt3a_status["Timepoint"] == "Baseline"][["Patient_id", "CH DNMT3A status"]]
ch_tet2_status = ch_tet2_status[ch_tet2_status["Timepoint"] == "Baseline"][["Patient_id", "CH TET2 status"]]
ch_asxl1_status = ch_asxl1_status[ch_asxl1_status["Timepoint"] == "Baseline"][["Patient_id", "CH ASXL1 status"]]
ch_tp53_status = ch_tp53_status[ch_tp53_status["Timepoint"] == "Baseline"][["Patient_id", "CH TP53 status"]]

data = clin.merge(ch_status, how = "left").\
    merge(chip_status, how = "left").\
    merge(ch_dta_status, how = "left").\
    merge(chip_dta_status, how = "left").\
    merge(ch_dnmt3a_status, how = "left").\
    merge(ch_tet2_status, how = "left").\
    merge(ch_asxl1_status, how = "left").\
    merge(ch_tp53_status, how = "left").\
    merge(median_values, how = "left").\
    merge(n_mutations, how = "left").\
    merge(median_values_dta, how = "left").\
    merge(n_mutations_dta, how = "left").\
    merge(median_values_tet2, how = "left").\
    merge(n_mutations_tet2, how = "left")

data["VAF_n"].fillna(0, inplace = True)
data["n_mutations"].fillna(0, inplace = True)
data["DTA_VAF_n"].fillna(0, inplace = True)
data["DTA_n_mutations"].fillna(0, inplace = True)
data["TET2_VAF_n"].fillna(0, inplace = True)
data["TET2_n_mutations"].fillna(0, inplace = True)

# Any CH in any gene, vs CH in DTA genes
run_fishers_exact(data, "CH status", col)
run_fishers_exact(data, "CHIP status", col)
run_fishers_exact(data, "CH DTA status", col)
run_fishers_exact(data, "CHIP DTA status", col)
run_fishers_exact(data, "CH DNMT3A status", col) 
run_fishers_exact(data, "CH TET2 status", col)
run_fishers_exact(data, "CH ASXL1 status", col)
run_fishers_exact(data, "CH TP53 status", col)

# None significant
run_mwu(data, "irAE", "VAF_n")
run_mwu(data, "irAE", "n_mutations")
run_mwu(data, "irAE", "DTA_VAF_n")
run_mwu(data, "irAE", "DTA_n_mutations")
run_mwu(data, "irAE", "TET2_VAF_n")
run_mwu(data, "irAE", "TET2_n_mutations")
run_mwu(data, "irAE", "Age at GUBB draw") # confounding
run_mwu(data, "irAE", "Charlson_score") # confounding
#-------------------------------------------



#=============================================
# TREATMENT RECEIVED
# Does the prevalence of CHIP differ in patients who received immunotherapy vs another type of treatment?
df = data.merge(clin_main[["Patient_id", "Type of treatment"]])
df["immunotherapy_received"] = df["Type of treatment"].apply(lambda x: x in ["1", "2", "3", "4", "5", "8"])
col = "immunotherapy_received"
run_fishers_exact(df, "CH status", col)
run_fishers_exact(df, "CHIP status", col)
run_fishers_exact(df, "CH DTA status", col)
run_fishers_exact(df, "CHIP DTA status", col)
run_fishers_exact(df, "CH DNMT3A status", col) 
run_fishers_exact(df, "CH TET2 status", col)
run_fishers_exact(df, "CH ASXL1 status", col)

# Does the prevalence of CHIP differ in patients who received a certain type of treatment or not? 
df = data.merge(clin_main[["Patient_id", "Type of treatment"]])
treatment_dict = {"Ipilimumab-nivolumab": 1, "Pembrolizumab-axitinib": 2, "Avelumab-axtinib": 3, "Nivolumab-cabozantinib": 4, "Pembrolizumab-lenvatinib": 5, "sunitnib": 6, "pazopanib": 7, "quavonlimab-pembrolizumab-lenvatinib": 8}
for key, value in treatment_dict.items():
    print("Testing the relationship between CH, CHIP, CHIP > 10 and " + key + " usage")
    df[key] = df["Type of treatment"].apply(lambda x: x in [value])
    run_fishers_exact(df, "CH status", key)
    run_fishers_exact(df, "CHIP status", key)
    run_fishers_exact(df, "CH DTA status", key)
    run_fishers_exact(df, "CHIP DTA status", key)
    run_fishers_exact(df, "CH DNMT3A status", key) 
    run_fishers_exact(df, "CH TET2 status", key)
    run_fishers_exact(df, "CH ASXL1 status", key)
    print(" ")
#=============================================


#=============================================
# COMORBIDITIES
# Does the prevalence of CHIP differ in patients with a given comorbidity or not? 
df = data.merge(clin_main[["Patient_id", 'CCI_MI', 'CCI_CHF', 'CCI_PVD', 'CCI_stroke', 'CCI_Dementia','CCI_Lung', 'CCI_Rheum', 'CCI_Ulcer', 'CCI_Liver', 'CCI_Diabetes','CCI_Hemiplegia', 'CCI_Renal', 'CCI_Malignancy', 'CCI_Lymphoma','CCI_Leukemia', 'CCI _AIDS']])
cci_columns = [col for col in df.columns if col.startswith("CCI")]
df[cci_columns] = df[cci_columns].applymap(lambda x: True if int(x) > 0 else False)
for column in cci_columns:
    print("Testing the relationship between CH, CHIP, CHIP > 10 and " + column)
    run_fishers_exact(df, "CH status", column)
    run_fishers_exact(df, "CHIP status", column)
    run_fishers_exact(df, "CH DTA status", column)
    run_fishers_exact(df, "CHIP DTA status", column)
    run_fishers_exact(df, "CH DNMT3A status", column) 
    run_fishers_exact(df, "CH TET2 status", column)
    run_fishers_exact(df, "CH ASXL1 status", column)
    print(" ")
    
# Does the prevalence of CHIP differ in patients with Charlson_score above and below 4? Above 4: 0.983^(e^(0.9*4)) the 10 year survival probability halves.
df = data.copy()
df["Charlson_score_bool"] = df["Charlson_score"] >= 4
run_fishers_exact(df, "CH status", "Charlson_score_bool")
run_fishers_exact(df, "CHIP status", "Charlson_score_bool")
run_fishers_exact(df, "CH DTA status", "Charlson_score_bool")
run_fishers_exact(df, "CHIP DTA status", "Charlson_score_bool")
run_fishers_exact(df, "CH DNMT3A status", "Charlson_score_bool") 
run_fishers_exact(df, "CH TET2 status", "Charlson_score_bool") 
run_fishers_exact(df, "CH ASXL1 status", "Charlson_score_bool")
#=============================================


#=============================================
# METASTASIS
# Does the prevalence of CHIP differ in patients with and without plurimetastasis? 
df = data.merge(clin_main[["Patient_id", "Mets_at_tx_start"]])
df = df[~pd.isnull(df["Mets_at_tx_start"])] # dropping patients with missing records for mets
df["Plurimetastasis_bool"] = df["Mets_at_tx_start"] >= 4
run_fishers_exact(df, "CH status", "Plurimetastasis_bool")
run_fishers_exact(df, "CHIP status", "Plurimetastasis_bool")
run_fishers_exact(df, "CH DTA status", "Plurimetastasis_bool")
run_fishers_exact(df, "CHIP DTA status", "Plurimetastasis_bool")
run_fishers_exact(df, "CH DNMT3A status", "Plurimetastasis_bool") 
run_fishers_exact(df, "CH TET2 status", "Plurimetastasis_bool") 
run_fishers_exact(df, "CH ASXL1 status", "Plurimetastasis_bool")

# Does the prevalence of CHIP differ in patients with and metastasis in given regions? 
df = data.merge(clin_main[["Patient_id", "Adrenal_Met", "Lung_Met", "LN_Met", "Bone_Met", "Liver_Met", "Skin_Met", "Soft_tissue_Met", "Brain_Met", "Spinal_Met", "CNS_Met", "Thyroid_Met", "Pancreas_Met", "Nephrectomy bed_Met"]])
met_columns = [col for col in df.columns if "Met" in col]
df[met_columns] = df[met_columns].applymap(lambda x: True if x == 1 else False if x == 0 else pd.NA)
for column in met_columns:
    print("Testing the relationship between CH and" + column)
    df = df[~pd.isnull(df[column])]
    run_fishers_exact(df, "CH status", column)
    run_fishers_exact(df, "CHIP status", column)
    run_fishers_exact(df, "CH DTA status", column)
    run_fishers_exact(df, "CHIP DTA status", column)
    run_fishers_exact(df, "CH DNMT3A status", column) 
    run_fishers_exact(df, "CH TET2 status", column) 
    run_fishers_exact(df, "CH ASXL1 status", column)
    print(" ")
#=============================================