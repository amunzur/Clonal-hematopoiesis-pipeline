import pandas as pd
import numpy as np
import os 
import re
from scipy.stats import ttest_ind
from lifelines import KaplanMeierFitter
from datetime import datetime
from lifelines.statistics import logrank_test
from scipy.stats import fisher_exact
from lifelines import KaplanMeierFitter
from scipy.stats import mannwhitneyu
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr
from scipy.stats import chi2_contingency

# LOAD SOMATIC DATASETS
all_vars_somatic = pd.read_csv(os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic.csv"))
all_vars_somatic = all_vars_somatic[~all_vars_somatic["Patient_id"].isin(["20-313", "21-184", "21-430"])] # exclude some samples due to oxidative damage
all_vars_somatic = all_vars_somatic[all_vars_somatic["Dependent"] == False]
base_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "Baseline") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_somatic = all_vars_somatic[(all_vars_somatic["Timepoint"] == "During treatment") & (all_vars_somatic["Diagnosis"] == "Bladder")].reset_index(drop = True)

all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip.csv"))
all_vars_chip = all_vars_chip[all_vars_chip["Dependent"] == False]
base_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
prog_kidney_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Kidney")].reset_index(drop = True)
base_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "Baseline") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)
prog_bladder_chip = all_vars_chip[(all_vars_chip["Timepoint"] == "During treatment") & (all_vars_chip["Diagnosis"] == "Bladder")].reset_index(drop = True)

np.median(base_kidney_somatic["VAF_t"])
np.median(base_kidney_chip["VAF_t"])
mannwhitneyu(base_kidney_somatic["VAF_t"], base_kidney_chip["VAF_t"])

np.median(base_bladder_somatic["VAF_t"])
np.median(base_bladder_chip["VAF_t"])
mannwhitneyu(base_bladder_somatic["VAF_t"], base_bladder_chip["VAF_t"])