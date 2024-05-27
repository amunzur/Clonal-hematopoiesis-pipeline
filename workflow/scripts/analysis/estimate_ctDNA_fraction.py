"""
Given a list of somatic mutations calculate the ctDNA fraction of the samples.
"""
# Estimate cfDNA using mutations
import pandas as pd
import numpy as np
import os
from scipy.stats import binom

def get_adj_vaf(vaf, depth):
    """
    Adjust the mutation VAF using binomial distribution.
    """
    alt = vaf*depth
    for p in np.arange(0,1,0.005):
        dist = binom(depth, p)
        if dist.cdf(alt) < 0.95:
            vaf = p; break;
    return(vaf)

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling"
PATH_mutations_curated = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic_SSCS2_curated_complete.csv"
PATH_sample_dict = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_all_samples = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_list_bladder_and_kidney.tsv"

PATH_mutation_ctfractions = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/ctDNA_fractions.csv" # output path
PATH_bladder_ctDNA = "/groups/wyattgrp/users/amunzur/pipeline/resources/ct_fractions/ERBB2 Supplemental Tables - ctDNA cohort.csv"


DIR_cnvkit = "/groups/wyattgrp/users/amunzur/pipeline/results/cnvkit/final_result"

sample_dict = pd.read_csv(PATH_sample_dict, "\t", names =["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
muts = pd.read_csv(PATH_mutations_curated)

# Gnenerate the median log2 ratio per gene per sample
mylist = []
for file in os.listdir(DIR_cnvkit):
    sample_name = file.replace(".tsv", "")
    print(sample_name)
    median_log2 = pd.read_csv(os.path.join(DIR_cnvkit, file), sep = "\t")[["gene", "log2"]]
    median_log2["median"] = median_log2.groupby("gene").transform("median")
    median_log2 = median_log2[["gene", "median"]].drop_duplicates().reset_index(drop = True)
    median_log2["Sample_name_t"] = sample_name
    mylist.append(median_log2)

gene_medians = pd.concat(mylist).reset_index(drop = True).rename(columns = {"gene": "Gene", "median": "Median log2"})
muts = muts.merge(gene_medians, how = "left")
muts["VAF_t"] = muts["VAF_t"]/100

muts = muts[(muts['Median log2'] > -0.3) & (muts['Median log2'] < 0.3) &
            (muts['Depth_t'] > 40) &
            ((muts['Diagnosis'] == "Kidney") | (muts['Diagnosis'] == "Bladder")) &
            (~muts['Chrom'].str.contains("X")) &
            (~muts['Chrom'].str.contains("Y")) & 
            (~muts["Patient_id"].isin(['20-313', '21-430', '21-184']))]

muts['Adj_allele_frequency'] = muts.apply(lambda row: get_adj_vaf(row['VAF_t'], row['Depth_t']), axis = 1)
muts["Max_Sample_adj_VAF"] = muts.groupby('Sample_name_t')["Adj_allele_frequency"].transform('max')
muts = muts[muts["Adj_allele_frequency"] == muts["Max_Sample_adj_VAF"]]
muts['row_max'] = muts[['Depth_t']].max(axis=1) #To break any ties for max VAF!
muts = muts.loc[muts.groupby('Sample_name_t')['row_max'].idxmax()]

muts['Mutation_ctDNA_fraction'] = 2/(1+1/(muts['Max_Sample_adj_VAF']))
muts = muts.merge(sample_dict, how = "outer") # some pts have no mutations available to call CT fractions, this retains them and allows a helpful comparison between CN based and mut based estimates across the whole cohort.
muts = muts[["Patient_id", "Diagnosis", "Sample_name_t", "Date_collected", "VAF_t", "Depth_t", "Gene", "Consequence", "Mutation_ctDNA_fraction"]]
muts = muts[~muts["Patient_id"].str.contains("Normal")] # remove VIP normals
muts["Mutation_ctDNA_fraction"] = muts["Mutation_ctDNA_fraction"].fillna(0)
muts["ctDNA_status"] = muts["Mutation_ctDNA_fraction"].apply(lambda x: "ctDNA negative" if x == 0 else "ctDNA positive")
muts['Date_collected'] = pd.to_datetime(muts['Date_collected'], format='%Y%b%d')

# add the sample name to the ctneg samples
all_samps = pd.read_csv(PATH_all_samples, sep = "\t").rename(columns = {"sample_names": "Sample_name_t"})
all_samps["Patient_id"] = all_samps["Sample_name_t"].str.split("_").str[0]
all_samps["Patient_id"] = all_samps["Patient_id"].str.replace("GU-", "")
all_samps["Date_collected"] = all_samps["Sample_name_t"].str.split("[-_]").str[-1]
cfDNA_samps = all_samps[all_samps["Sample_name_t"].str.contains("cfDNA")]
cfDNA_samps["Date_collected"] = pd.to_datetime(cfDNA_samps['Date_collected'], format='%Y%b%d')
del muts["Sample_name_t"]
muts = muts.merge(cfDNA_samps, how = "left")

# Add the ctDNA fraction estimations from the bladder project to compare
bladder_ctDNA = (
    pd.read_csv(PATH_bladder_ctDNA)[["Patient ID", "Sample ID", "Collection date", "ctDNA fraction (%)"]]
    .rename(columns={"Patient ID": "Patient_id", "Sample ID": "Sample_name_t", "Collection date": "Date_collected", "ctDNA fraction (%)": "ctDNA_fraction_bladder_project"})
    .assign(ctDNA_fraction_bladder_project=lambda x: x["ctDNA_fraction_bladder_project"] / 100)
    .assign(Date_collected=lambda x: pd.to_datetime(x['Date_collected'], format='%d-%b-%Y'))
    .assign(Date_collected=lambda x: x['Date_collected'].dt.strftime('%Y%b%d'))
    .query('Patient_id.str.startswith("GU-")', engine='python')  # Filter rows where Patient_id starts with 'GU-'
    .assign(Patient_id=lambda x: x['Patient_id'].str.replace('GU-', ''))
)

# combined = muts.merge(bladder_ctDNA, how = "outer", on = ["Patient_id", "Date_collected"])
muts.to_csv(PATH_mutation_ctfractions, index = False)