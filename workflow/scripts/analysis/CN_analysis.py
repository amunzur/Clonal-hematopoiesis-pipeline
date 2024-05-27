"""
This script uses cnvkit outputs and generates a table that indicates the median log value across all probes in each gene. 
Separate values are provided for cfDNA and WBC samples. 
conda activate pybedtools
"""
import pandas as pd
import pybedtools
import os

def return_probe_coverage(PATH_targets, PATH_coverage, DIR_output): 
    """
    Given a gene name and a path to a file that has the fixed coverage information, return the coverage of all probes
    that overlap the gene region.
    path_targets: bed file from the panel that also has the gene annotations
    path_coverage: coverage output for a single sample from cnvkit
    Output is a df that has the log2 ratios of all regions in the panel in a sample of interest.
    """
    # Prepare files for bedtools intersect
    coverage = pd.read_csv(PATH_coverage, sep = "\t").rename(columns = {"start": "probe_start", "end": "probe_end"})
    coverage = coverage[coverage["chromosome"].str.startswith("chr")].reset_index(drop = True)
    
    targets = pd.read_csv(PATH_targets, sep = "\t", names = ["chromosome", "target_start", "target_end", "gene"])
    
    # convert to bed
    coverage_bed = pybedtools.BedTool.from_dataframe(coverage[["chromosome", "probe_start", "probe_end"]])
    targets_bed = pybedtools.BedTool.from_dataframe(targets[["chromosome", "target_start", "target_end"]])
    
    # Intersect and add some columns
    sample_df = targets_bed.intersect(coverage_bed, wb = True, wa = True).to_dataframe()
    sample_df.columns = ["chromosome", "target_start", "target_end", "chromosome_del", "probe_start", "probe_end"]
    sample_df = sample_df.merge(targets, how = "left") # add gene annotations
    sample_df = sample_df.merge(coverage.drop("gene", axis = 1), how="left") # add log ratio
    del sample_df["chromosome_del"]
    
    path_output = os.path.join(DIR_output, os.path.basename(PATH_coverage).replace(".cnr", ".tsv"))
    sample_df.to_csv(path_output, index = False, sep = "\t")


# Bed file that has both genomic locations and the gene names corresponding to those positions
DIR_output = "/groups/wyattgrp/users/amunzur/pipeline/results/cnvkit/final_result"
PATH_targets = "/groups/wyattgrp/users/amunzur/pipeline/resources/panel/targets_and_genes.bed"

DIR_probes_cfDNA = "/groups/wyattgrp/users/amunzur/pipeline/results/cnvkit/coverage_fixed/cfDNA"
DIR_probes_WBC = "/groups/wyattgrp/users/amunzur/pipeline/results/cnvkit/coverage_fixed/WBC"

files_cfDNA = [os.path.join(DIR_probes_cfDNA, file) for file in os.listdir(DIR_probes_cfDNA)]
files_WBC = [os.path.join(DIR_probes_WBC, file) for file in os.listdir(DIR_probes_WBC)]

for coverage_file in files_cfDNA:
    print(coverage_file)
    return_probe_coverage(PATH_targets, coverage_file, DIR_output)

for coverage_file in files_WBC:
    print(coverage_file)
    return_probe_coverage(PATH_targets, coverage_file, DIR_output)
