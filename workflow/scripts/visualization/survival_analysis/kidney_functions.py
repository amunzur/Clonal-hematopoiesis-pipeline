
def prepare_survival_data(clin_df):
    """Prepare survival data merged with patient data."""
    surv_df = clin_df[["Patient_id", "Death", "OS from cfDNA collection (mo)", "Progression", "PFS (mo)"]]
    return surv_df

def prepare_clinical_data(clin_df):
    """Prepare clinical data including sex, age, CCI, mets, IMDC risk, ctDNA status, and tumor subtype."""
    sex_and_age_df = clin_df[["Patient_id", "Sex", "Age at GUBB draw"]]
    cci_df = prepare_cci_data(clin_df)
    mets_df = prepare_mets_data(clin_df)
    imdc_df = prepare_imdc_data(clin_df)
    subtype_df = prepare_subtype_data(clin_df)
    surv_df = prepare_survival_data(clin_df)
    return sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df, surv_df

def prepare_cci_data(clin_df):
    """Prepare Charlson Comorbidity Index data."""
    cci_df = clin_df[["Patient_id", "Charlson_score"]]
    dict_to_replace = {
        "4*see comment, many comorbidities not meeting CCI": "4",
        "1** (see comment)": "1", 
        "3**see comment": "3", 
        "5* see comment": "5", 
        "1*see comment": "1", 
        "6*see comment": "6"
    }
    cci_df["Charlson_score"] = cci_df["Charlson_score"].replace(dict_to_replace).astype(int)
    median_value = cci_df["Charlson_score"].median()
    cci_df["CCI above median"] = cci_df["Charlson_score"] > median_value
    return cci_df[["Patient_id", "CCI above median"]]

def prepare_mets_data(clin_df):
    """Prepare visceral metastasis data."""
    col_names = []
    mets_df = clin_df[["Patient_id", "Adrenal_Met", "Lung_Met", "Liver_Met"]]
    mets_df[["Adrenal_Met", "Lung_Met", "Liver_Met"]] = mets_df[["Adrenal_Met", "Lung_Met", "Liver_Met"]].fillna(0).astype(int)
    mets_df["summed"] = mets_df[["Adrenal_Met", "Lung_Met", "Liver_Met"]].sum(axis=1)
    mets_df["Visceral mets"] = mets_df["summed"] > 0
    mets_df = mets_df[["Patient_id", "Visceral mets"]]
    return(mets_df)

def prepare_imdc_data(clin_df):
    """Prepare IMDC risk data."""
    imdc_df = clin_df[["Patient_id", "IMDC coded"]]
    imdc_df["IMDC_poor_to_intermediate_risk"] = imdc_df["IMDC coded"].isin(["2", "2a", "2b", "3"])
    return imdc_df[["Patient_id", "IMDC_poor_to_intermediate_risk"]]

def prepare_ctDNA_status(muts_df, diagnosis, PATH_sample_information, annotate_gene = False):
    """Prepare ctDNA status data."""
    ctDNA_status = annotate_mutation_status(muts_df, diagnosis, PATH_sample_information, annotate_what="ctDNA", annotate_gene = annotate_gene)
    ctDNA_status = ctDNA_status[ctDNA_status["Timepoint"] == "Baseline"]
    ctDNA_status["ctDNA positive"] = ctDNA_status["ctDNA status"] == "Positive"
    return ctDNA_status[["Patient_id", "ctDNA positive"]]

def prepare_chip_status(muts_df, diagnosis, PATH_sample_information, min_vaf_threshold = 0, annotate_gene = False):
    """Prepare CHIP status data."""
    muts_df = muts_df[muts_df["VAF_n"] >=min_vaf_threshold]
    chip_status = annotate_mutation_status(muts_df, diagnosis, PATH_sample_information, annotate_what = "CHIP", annotate_gene = annotate_gene)
    chip_status = chip_status[chip_status["Timepoint"] == "Baseline"][["Patient_id", "CHIP status"]]
    chip_status["CHIP positive"] = chip_status["CHIP status"] == "Positive"
    chip_status = chip_status[["Patient_id", "CHIP positive"]]
    return(chip_status)

def prepare_subtype_data(clin_df):
    """Prepare tumor subtype data."""
    subtype_df = clin_df[["Patient_id", "subtype"]]
    subtype_df["Clear cell"] = subtype_df["subtype"].isin(["0", "1"])
    return subtype_df[["Patient_id", "Clear cell"]]

def prepare_multivars_df(sex_and_age_df, cci_df, mets_df, ctDNA_status, subtype_df, chip_status, imdc_df):
    """Combine all variables into a single DataFrame and binarize categorical variables."""
    multivars_df = sex_and_age_df.merge(cci_df).merge(mets_df).merge(ctDNA_status).merge(subtype_df).merge(chip_status).merge(imdc_df)
    multivars_df_binarized = pd.get_dummies(multivars_df, columns=["CCI above median", "Sex", "Visceral mets", "ctDNA positive", "Clear cell", "CHIP positive", "IMDC_poor_to_intermediate_risk"], drop_first=True)
    multivars_df_binarized.rename(columns={
        "Age at GUBB draw": "Age at draw", 
        "CCI above median_True": "CCI>median", 
        "Sex_Male": "Male", 
        "Visceral mets_True": "Visceral mets", 
        "ctDNA positive_True": "ctDNA+", 
        "Clear cell_True": "Clear cell", 
        "CHIP positive_True": "CH+", 
        "IMDC_poor_to_intermediate_risk_True": "IMDC poor/int. risk"
    }, inplace=True)
    return multivars_df_binarized

def prepare_multivars_df_bladder(sex_and_age_df, smoking_status, ctDNA_status, chip_status):
    """Combine all variables into a single DataFrame and binarize categorical variables."""
    smoking_status.loc[smoking_status["Previous smoking history"].isin(["Current smoker", "Previous smoker"]), "Previous smoking history"] = "curr/prev smoker"
    multivars_df = sex_and_age_df_bladder.merge(chip_status).merge(smoking_status).merge(ctDNA_status)
    multivars_df_binarized = pd.get_dummies(multivars_df, columns=["Sex", "ctDNA positive", "CHIP positive", "Previous smoking history"], drop_first=True)
    multivars_df_binarized.rename(columns={
        "Age at baseline blood draw": "Age at draw", 
        "Sex_Male": "Male", 
        "ctDNA positive_True": "ctDNA+", 
        "CHIP positive_True": "CH+",
        "Previous smoking history_curr/prev smoker": "curr/prev smoker"
    }, inplace=True)
    return multivars_df_binarized


def prepare_multivars_df_with_3_genomic_vars(sex_and_age_df, cci_df, mets_df, subtype_df, imdc_df, genomics_df):
    """Combine all variables into a single DataFrame and binarize categorical variables."""
    multivars_df = sex_and_age_df.merge(cci_df).merge(mets_df).merge(subtype_df).merge(imdc_df).merge(genomics_df)
    multivars_df_binarized = pd.get_dummies(multivars_df, columns=["CCI above median", "Sex", "Visceral mets", "Clear cell", "IMDC_poor_to_intermediate_risk", "CH only", "ctDNA only", "CH and ctDNA"], drop_first=True)
    multivars_df_binarized.rename(columns={
        "Age at GUBB draw": "Age at draw", 
        "CCI above median_True": "CCI>median", 
        "Sex_Male": "Male", 
        "Visceral mets_True": "Visceral mets", 
        "Clear cell_True": "Clear cell", 
        "IMDC_poor_to_intermediate_risk_True": "IMDC poor/int. risk", 
        "CH only_True": "CH+", 
        "ctDNA only_True": "ctDNA+", 
        "CH and ctDNA_True": "CH+ and ctDNA+"
    }, inplace=True)
    return multivars_df_binarized

def make_survival_df(Path_clinical_data, PATH_sample_information, diagnosis):
    surv_df = pd.read_csv(Path_clinical_data)[["Patient_id", 'Date of last follow-up or death', 'Death at last follow up', 'Age at last follow-up or death']]
    pts = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])  
    
    pts = pts[(pts["Diagnosis"] == diagnosis) & (pts["Timepoint"] == "Baseline")].reset_index(drop = True)  
    
    surv_df = surv_df.merge(pts, how = "inner")
    surv_df["Date of last follow-up or death"] = pd.to_datetime(surv_df["Date of last follow-up or death"])
    surv_df["Date_collected"] = pd.to_datetime(surv_df["Date_collected"], format='%Y%b%d', errors='coerce')
    surv_df["Overall survival"] = (surv_df["Date of last follow-up or death"] - surv_df["Date_collected"]).astype('timedelta64[M]')
    surv_df["Overall survival"] = surv_df["Overall survival"].round().astype(int)
    
    return(surv_df)

def prepare_ctDNA_status_high_low(baseline_somatic, PATH_sample_information):
    """
    Groups patients into low and high ctDNA fraction based on interquartiles.
    """
    # Get ctDNA neg samples
    diagnosis = baseline_somatic["Diagnosis"].unique()[0]
    ctDNA_status = annotate_mutation_status(baseline_somatic, diagnosis, PATH_sample_information, annotate_what = "ctDNA")
    ctDNA_status = ctDNA_status[(ctDNA_status["Timepoint"] == "Baseline") & (ctDNA_status["ctDNA status"] == "Negative")]
    ctDNA_status["ctDNA status"] = "low"
    ctDNA_status = ctDNA_status[["Patient_id", "ctDNA status"]].reset_index(drop = True)
    
    # calculate the 25th percentile. 
    q25 = np.percentile(baseline_somatic["VAF_t"], 25)
    baseline_somatic = baseline_somatic.groupby("Patient_id")["VAF_t"].max().reset_index()
    baseline_somatic["ctDNA status"] = ["low" if vaf < q25 else "high" for vaf in baseline_somatic["VAF_t"]]
    baseline_somatic = baseline_somatic[["Patient_id", "ctDNA status"]]
    baseline_somatic = pd.concat([baseline_somatic, ctDNA_status])
    
    baseline_somatic = baseline_somatic.rename(columns = {"ctDNA status": "ctDNA positive"})
    baseline_somatic["ctDNA positive"] = baseline_somatic["ctDNA positive"].replace({"high": True, "low": False})
    
    return(baseline_somatic)
# mpl.rcParams['font.size'] = 10
# mpl.rcParams['text.color'] = 'k'
# mpl.rcParams['legend.fontsize'] = 10
# mpl.rcParams['legend.handletextpad'] = '0.8'
# mpl.rcParams['legend.labelspacing'] = '0.4'
# mpl.rcParams['pdf.fonttype'] = 42
# mpl.rcParams['ps.fonttype'] = 42
# mpl.rcParams['axes.linewidth'] = 1
# mpl.rcParams['xtick.labelsize'] = 10
# mpl.rcParams['ytick.labelsize'] = 10
# mpl.rcParams['axes.labelsize'] = 10

# DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
# PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
# PATH_gene_categories = os.path.join(DIR_working, "resources/panel/chip_panel_gene_categories.tsv")
# PATH_kidney_clinical = os.path.join(DIR_working, "resources/clinical_data/Kidney/mRCC clinical Data_clean.csv")
# PATH_clinical_bladder = os.path.join(DIR_working, "resources/clinical_data/Bladder_enrollment.csv")
# PATH_treatment_landscape = os.path.join(DIR_working, "resources/clinical_data/bladder/treatment.csv")
# PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv"
# PATH_kidney_clean = "/groups/wyattgrp/users/amunzur/pipeline/resources/clinical_data/kidney_clin_clean.csv"

# figure_dir = os.path.join(DIR_working, "results/figures/amazing_figures")
# source_functions = os.path.join(DIR_working, "workflow/scripts/visualization/UTILITIES_make_chip_plots.py")

# color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}

# with open(source_functions, 'r') as file:
#     script_code = file.read()

# exec(script_code)

# # LOAD CHIP DATASETS
# all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
# sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
# all_vars_chip = sample_info.merge(all_vars_chip, how = "inner")

# base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
# prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
# base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
# prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

# # LOAD SOMATIC DATASETS
# all_vars_somatic = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic.csv"))
# all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
# base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
# prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
# base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
# prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

# # SURVIVAL ANALYSIS
# # Multivar cox analysis with the following in the RCC cohort:
# # charlson score > median
# # has mets at TX start
# # tumor subtype
# # age at blood draw
# # sex
# # has detectable ctDNA or not 

# clin_df = pd.read_csv(PATH_kidney_clinical)

# sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df = prepare_clinical_data(clin_df)
# ctDNA_status = prepare_ctDNA_status(base_kidney_somatic, diagnosis, PATH_sample_information)
# chip_status = prepare_chip_status(base_kidney_chip, diagnosis, PATH_sample_information)
# multivars_df_binarized = prepare_multivars_df(sex_and_age_df, cci_df, mets_df, ctDNA_status, subtype_df, chip_status, imdc_df)
# surv_df = make_survival_df(PATH_kidney_clean, PATH_sample_information)
# merged_df = surv_df.merge(multivars_df_binarized)
# merged_df = merged_df[['Overall survival', 'Death at last follow up',
#                        'Age at blood draw', 'CCI > median', 'Male', 'Visceral mets', 'ctDNA positive', 
#                        'Clear cell subtype', 'CHIP positive', 'IMDC poor/intermediate risk']]   

# fig = plt.figure(figsize=(8, 4))
# gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.7], hspace=0.4)
# ax_km = fig.add_subplot(gs[0])
# ax_cox = fig.add_subplot(gs[1]) 

# ax_km, ax_cox = make_survival_curve_with_forest(merged_df, "CHIP positive", ax_km, ax_cox, [0.52, 0.48], [0.35, 0.22])
# ax_km.set_title("Overall survival in mRCC cohort (Any CH)")

# gs.tight_layout(fig)
# fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/survival_figure3/kidney_ch.png")

# ##############################################3
# # CH > 2%
# clin_df = pd.read_csv(PATH_kidney_clinical)

# sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df = prepare_clinical_data(clin_df)
# ctDNA_status = prepare_ctDNA_status(base_kidney_somatic, diagnosis, PATH_sample_information)
# chip_status = prepare_chip_status(base_kidney_chip, diagnosis, PATH_sample_information, min_vaf_threshold = 2)
# multivars_df_binarized = prepare_multivars_df(sex_and_age_df, cci_df, mets_df, ctDNA_status, subtype_df, chip_status, imdc_df)
# surv_df = make_survival_df(PATH_kidney_clean, PATH_sample_information)
# merged_df = surv_df.merge(multivars_df_binarized)
# merged_df = merged_df[['Overall survival', 'Death at last follow up',
#                        'Age at blood draw', 'CCI > median', 'Male', 'Visceral mets', 'ctDNA positive', 
#                        'Clear cell subtype', 'CHIP positive', 'IMDC poor/intermediate risk']]   

# fig = plt.figure(figsize=(8, 4))
# gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.7], hspace=0.4)
# ax_km = fig.add_subplot(gs[0])
# ax_cox = fig.add_subplot(gs[1]) 

# ax_km, ax_cox = make_survival_curve_with_forest(merged_df, "CHIP positive", ax_km, ax_cox, [0.52, 0.48], [0.35, 0.22])
# ax_km.set_title("Overall survival in mRCC cohort (CH>2)")

# gs.tight_layout(fig)
# fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/survival_figure3/kidney_chip.png")

# ##############################################3
# # Now stratify it into 3 groups: CH only, ctDNA only, CH and ctDNA
# clin_df = pd.read_csv(PATH_kidney_clinical)

# sex_and_age_df, cci_df, mets_df, imdc_df, subtype_df = prepare_clinical_data(clin_df)
# ctDNA_status = prepare_ctDNA_status(base_kidney_somatic, diagnosis, PATH_sample_information)
# chip_status = prepare_chip_status(base_kidney_chip, diagnosis, PATH_sample_information)

# genomics_df = chip_status.merge(ctDNA_status)
# genomics_df["CH only"] = (genomics_df["CHIP positive"] == True) & (genomics_df["ctDNA positive"] == False)
# genomics_df["ctDNA only"] = (genomics_df["CHIP positive"] == False) & (genomics_df["ctDNA positive"] == True)
# genomics_df["CH and ctDNA"] = (genomics_df["CHIP positive"] == True) & (genomics_df["ctDNA positive"] == True)
# genomics_df = genomics_df[["Patient_id", "CH only", "ctDNA only", "CH and ctDNA"]]

# multivars_df_binarized = prepare_multivars_df_with_3_genomic_vars(sex_and_age_df, cci_df, mets_df, subtype_df, imdc_df, genomics_df)
# surv_df = make_survival_df(PATH_kidney_clean, PATH_sample_information)
# merged_df = surv_df.merge(multivars_df_binarized)
# merged_df = merged_df[['Overall survival', 'Death at last follow up',
#                        'Age at blood draw', 'CCI > median', 'Male', 'Visceral mets',
#                        "CH only", "ctDNA only", "CH and ctDNA",
#                        'Clear cell subtype', 'IMDC poor/intermediate risk']]   


# fig = plt.figure(figsize=(8, 4))
# gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.7], hspace=0.4)
# ax_km = fig.add_subplot(gs[0])
# ax_cox = fig.add_subplot(gs[1]) 

# ax_km, ax_cox = make_survival_curve_with_forest_3_groups(merged_df, ax_km, ax_cox, add_legend = False)

# ax_km.set_title("Overall survival in mRCC cohort (CH>2)")

# gs.tight_layout(fig)
# fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/survival_figure3/kidney_3_categories.png")
