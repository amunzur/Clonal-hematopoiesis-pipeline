# Downsamples each wbc to 1500, 1000, 750, 500, 250

import subprocess
import pandas as pd
import os 
import numpy as np
import time

DIR_bams = "/groups/wyattgrp/users/amunzur/pipeline/results/data/bam/SSCS2_final"
PATH_sample_information = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_information.tsv"
PATH_median_depth = "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/SSCS2/averaged_depths.txt"

sample_df = pd.read_csv(PATH_sample_information, sep="\t", header=None, names = ["Patient_id", "Date", "Diagnosis", "Timepoint"])
sample_df = sample_df[sample_df["Timepoint"] == "Baseline"].reset_index(drop = True)

# pts_remaining = [
# "GU-18-526",
# "GU-19-201",
# "GU-20-006",
# "GU-20-011",
# "GU-20-034",
# "GU-20-037",
# "GU-20-204",
# "GU-20-326",
# "GU-20-383",
# "GU-21-024",
# "GU-21-128",
# "GU-21-163",
# "GU-21-193",
# "GU-21-234",
# "GU-21-430",
# "GU-21-520",
# "GU-21-532",
# "GU-21-545",
# "GU-22-013",
# "GU-22-021",
# "GU-22-022",
# "GU-22-054",
# "GU-22-082",
# "GU-22-088",
# "GU-22-154",
# "GU-22-161",
# "GU-22-205",
# "GU-22-211",
# "GU-22-217",
# "GU-22-273"]

DIR_outputs = "/groups/wyattgrp/users/amunzur/pipeline/results/wbc_downsampling"

def return_subsample_fraction(PATH_median_depth, sample_name, target_depth):
    """
    Based on how much depth the sample got and the target depth we would like to subsample to, returns the subsampling factor.
    """
    median_depth_df = pd.read_csv(PATH_median_depth, header = None, sep = " ")
    median_depth_df.columns = ["sample_name", "av", "median"]
    acquired_depth = median_depth_df[median_depth_df["sample_name"] == sample_name]["median"].iloc[0]
    subsample_factor = target_depth/acquired_depth
    return(subsample_factor)

# Function to run commands in batches
# Function to manage active processes and run commands
def manage_active_processes(commands_list, batch_size):
    active_processes = []  
    idx = 0
    while idx < len(commands_list) or active_processes:
        # Start new processes as long as there are less than batch_size active processes
        while len(active_processes) < batch_size and idx < len(commands_list):
            command = commands_list[idx]
            print(f"Starting command in background: {command}")
            process = subprocess.Popen(command, shell=True)
            active_processes.append(process)
            idx += 1
        
        # Check if any processes have completed
        for process in list(active_processes):
            if process.poll() is not None:  # Process has completed
                active_processes.remove(process)
        
        time.sleep(1)  # Adjust sleep time as needed
    print("All commands have been started and are running.")

# Subsample the bam files
commands_list = []
for target_depth in [500]:
    for i, row in sample_df.iterrows():
        
        pt = row["Patient_id"]
        
        # Return sample name and the path to the bam file
        sample_name = [p for p in os.listdir(DIR_bams) if pt in p and "WBC" in p and p.endswith(".bam") and "aseline" in p][0].replace(".bam", "")
        path_sample = os.path.join(DIR_bams, sample_name + ".bam")
        subsample_factor = return_subsample_fraction(PATH_median_depth, sample_name, target_depth)
        
        # set the output path
        path_output = os.path.join(DIR_outputs, "bam", f"depth_{target_depth}/SSCS2_final", sample_name+".bam")
        os.makedirs(os.path.dirname(path_output), exist_ok=True)
        
        if subsample_factor >= 1:
            command = f"cp {path_sample} {path_output} && samtools index {path_output}"
            commands_list.append(command)
        else:
            command = f"samtools view -s {subsample_factor} -b {path_sample} > {path_output} && samtools index {path_output}"
            commands_list.append(command)

batch_size = 60
manage_active_processes(commands_list, batch_size)