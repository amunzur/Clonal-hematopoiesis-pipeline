"""
This script generates some stats that we mentioned in the paper. Goal is to be able to easily reproduce those stats if needed.
Things like % of the DNMT3A mutations in the whole cohort etc.
"""

import pandas as pd
import numpy as np
import os

PATH_CHIP = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/chip_SSCS2_curated_complete.csv"
PATH_SOMATIC = "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/somatic_SSCS2_curated_complete.csv"

chip = pd.read_csv(PATH_CHIP)
somatic = pd.read_csv(PATH_SOMATIC)

baseline_chip = chip[chip["Timepoint"] == "Baseline"].reset_index(drop = True)
baseline_somatic = somatic[somatic["Timepoint"] == "Baseline"].reset_index(drop = True)

import pandas as pd

# TOP 5 MUTATED GENES (percentages given) AT BASELINE - CHIP 
gene_counts_CH_baseline = baseline_chip.groupby("Diagnosis")["Gene"].value_counts().reset_index(name='Count') # Calculate the count of each gene for each diagnosis
total_counts_CH_baseline = gene_counts_CH_baseline.groupby("Diagnosis")["Count"].sum().reset_index(name='Total') # Calculate the total count of genes for each diagnosis
merged_CH_baseline = pd.merge(gene_counts_CH_baseline, total_counts_CH_baseline, on="Diagnosis") # Merge the gene counts and total counts DataFrames
merged_CH_baseline['Percentage'] = (merged_CH_baseline['Count'] / merged_CH_baseline['Total']) * 100 # Calculate the percentage of each gene count
sorted_genes_CH_baseline = merged_CH_baseline.sort_values(by=['Diagnosis', 'Percentage'], ascending=[True, False]) # Sort the DataFrame by percentage in descending order
top_genes_CH_baseline = sorted_genes_CH_baseline.groupby('Diagnosis').head(5) # Get the top 5 genes for each diagnosis

# CHIP STATS WITHOUT GROUPING BY DIAGNOSIS
gene_counts_CH_no_group = baseline_chip["Gene"].value_counts().reset_index(name='Count') # Calculate the count of each gene
total_counts_CH_no_group = gene_counts_CH_no_group["Count"].sum() # Calculate the total count of genes
gene_counts_CH_no_group['Percentage'] = (gene_counts_CH_no_group['Count'] / total_counts_CH_no_group) * 100 # Calculate the percentage of each gene count
sorted_genes_CH_no_group = gene_counts_CH_no_group.sort_values(by='Percentage', ascending=False) # Sort the DataFrame by percentage in descending order
top_genes_CH_no_group = sorted_genes_CH_no_group.head(5) # Get the top 5 genes

# TOP 5 MUTATED GENES (percentages given) AT BASELINE - SOMATIC
gene_counts_somatic = baseline_somatic.groupby("Diagnosis")["Gene"].value_counts().reset_index(name='Count') # Calculate the count of each gene for each diagnosis
total_counts_somatic = gene_counts_somatic.groupby("Diagnosis")["Count"].sum().reset_index(name='Total') # Calculate the total count of genes for each diagnosis
merged_somatic = pd.merge(gene_counts_somatic, total_counts_somatic, on="Diagnosis") # Merge the gene counts and total counts DataFrames
merged_somatic['Percentage'] = (merged_somatic['Count'] / merged_somatic['Total']) * 100 # Calculate the percentage of each gene count
sorted_genes_somatic = merged_somatic.sort_values(by=['Diagnosis', 'Percentage'], ascending=[True, False]) # Sort the DataFrame by percentage in descending order
top_genes_somatic = sorted_genes_somatic.groupby('Diagnosis').head(5) # Get the top 5 genes for each diagnosis

# NUMBER OF PATIENTS AND THE NUMBER OF MUTATIONS
patient_counts = baseline_chip["Patient_id"].value_counts() # Count the occurrences of each unique patient ID
mutation_counts = patient_counts.value_counts().reset_index() # Count the occurrences of each count of mutations
mutation_counts.columns = ['Number of Mutations', 'Number of Patients']
total_patients = len(patient_counts) # Calculate the percentage of patients for each count of mutations
mutation_counts['Percentage'] = (mutation_counts['Number of Patients'] / total_patients) * 100

"""
Out of the total 723 CH mutations called at baseline samples, 331 (45.7%) originated from bladder samples and 
392 (54.2%) originated from kidney samples. Across the whole cohort the top 5 mutated genes were DNMT3A (35.4%), TET2 (12.2%), PPM1D (7.05%), 
ASXL1 (5.94%) and TP53 (4.14%). The mutation distribution revealed that 31% of the cohort had one CH mutation, 26.3% had 2 and 18.2% had 3 mutations. The remaining XX percent of the cohort
had 4 or more mutations. Notably we observed one hypermutated case in a sample from a patient with bladder cancer, showcasing 17 CH mutations at baseline. The median 
WBC VAF of CH mutations was XX in kidney and XX in bladder samples. 


"""