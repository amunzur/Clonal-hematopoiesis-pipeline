import pandas as pd
import numpy as np
import os 
import re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from itertools import cycle
from scipy.stats import ttest_ind
import seaborn as sns
from lifelines import KaplanMeierFitter
from adjustText import adjust_text
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.cm import get_cmap
from datetime import datetime
import matplotlib.cm as cm
from lifelines.statistics import logrank_test
from scipy.stats import fisher_exact
from lifelines import KaplanMeierFitter
import matplotlib.ticker as ticker

"""
This script makes gene violins based on the cleaned cnvkit results computed by the CN_analysis.py script.
"""

DIR_figures_cnvkit = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/cnvkit"
DIR_cnvkit_results = "/groups/wyattgrp/users/amunzur/pipeline/results/cnvkit/final_result" # cleaned cnvkit results
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
samples = pd.read_csv(PATH_sample_information, names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"], sep = "\t")
samples = samples[samples["Diagnosis"].str.contains("Bladder|Kidney", regex = True)]

for i, sample in samples.iterrows(): 
    patient_id = sample["Patient_id"]
    diagnosis = sample["Diagnosis"]
    date_collected = sample["Date_collected"]
    timepoint = sample["Timepoint"]
    figure_path = os.path.join(DIR_figures_cnvkit, diagnosis, patient_id + "_" + date_collected + ".png")
    print(f"Working on {patient_id}.")
    
    # identify file names that we will work with
    file_cfDNA_list = [os.path.join(DIR_cnvkit_results, file) for file in os.listdir(DIR_cnvkit_results) if "cfDNA" in file and patient_id in file and date_collected in file]
    file_WBC_list = [os.path.join(DIR_cnvkit_results, file) for file in os.listdir(DIR_cnvkit_results) if "WBC" in file and patient_id in file and date_collected in file]
    
    # Check to make sure there is only one file each
    if len(file_cfDNA_list) > 1 or len(file_WBC_list) > 1:
        raise ValueError("More than one file matched for cfDNA or WBC.")
    elif len(file_cfDNA_list) == 0 or len(file_WBC_list) == 0:
        print(f"No files found for {patient_id}. Continuing with other samples.")
    else:
        wbc_coverage_df = pd.read_csv(file_WBC_list[0], sep = "\t")
        cfDNA_coverage_df = pd.read_csv(file_cfDNA_list[0], sep = "\t")
        
        # Initialize the figure
        fig = plt.figure(figsize=(10, 8))
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace = 0.5)
        ax_t = plt.subplot(gs[0])
        ax_n = plt.subplot(gs[1])
        ax_t = make_violin_gene_cn(cfDNA_coverage_df, ax = ax_t, ax_title = "cfDNA")
        ax_n = make_violin_gene_cn(wbc_coverage_df, ax = ax_n, ax_title = "WBC")
        fig.suptitle(f"{diagnosis}, {patient_id}, {date_collected}, {timepoint}")
        fig.savefig(figure_path)





