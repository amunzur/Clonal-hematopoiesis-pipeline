import os
import pandas as pd
import numpy as np
import subprocess
import argparse

# === Utility Functions ===
def run_grep_in_groups(commands, group_size=100):
    """
    Runs grep commands in specified batches.
    """
    total_commands = len(commands)
    for i in range(0, total_commands, group_size):
        group = commands[i:i + group_size]
        print(f"Running group {i // group_size + 1} of {total_commands // group_size + 1}")
        processes = [subprocess.Popen(cmd, shell=True) for cmd in group]
        for p in processes:
            p.wait()
        print(f"Group {i // group_size + 1} completed.")

def get_read_support_from_grep(path_grepped_output, type_of_alteration, alt_base=None):
    """
    Returns read support for a given mutation from the grepped pileup.
    """
    try:
        df = pd.read_csv(path_grepped_output, sep="\t", names=["chrom", "position", "ref", "depth_red", "seq", "quality_scores"])
        if df.empty:
            return 0
        seq = df["seq"][0]
        if type_of_alteration.upper() == "DELETION":
            return seq.count("*")
        elif type_of_alteration.upper() == "INSERTION":
            return seq.count("+")
        elif type_of_alteration.upper() == "SNV" and alt_base is not None:
            return seq.count(alt_base.upper())
    except Exception:
        return 0
    return 0

def get_depth_from_pileup(path_grepped_output):
    """
    Calculates the total depth at a given position using the grepped pileup.
    """
    try:
        df = pd.read_csv(path_grepped_output, sep="\t", names=["chrom", "position", "ref", "depth_red", "seq", "quality_scores"])
        return df["depth_red"][0] if not df.empty else 0
    except Exception:
        return 0

# === Load data ===
def main(DIR_mpileups, DIR_filtered_mpileups, PATH_sample_info, PATH_mutations, PATH_output, force_run):
    
    os.makedirs(DIR_filtered_mpileups, exist_ok=True)
    
    sample_info = pd.read_csv(PATH_sample_info, sep="\t")
    mutations = pd.read_csv(PATH_mutations)
    mutations["Position"] = mutations["Position"].astype(str)   
    
    # Build a mapping: patient -> mpileup files
    patient_to_mpileups = {}
    for f in os.listdir(DIR_mpileups):
        for patient in sample_info["Patient_id"].unique():
            if patient in f:
                patient_to_mpileups.setdefault(patient, []).append(os.path.join(DIR_mpileups, f))
    
    # Step 1: Run greps only within patient samples ===
    grep_cmds = []
    for patient_id, mut_df in mutations.groupby("Patient_id"):
        mpileup_files = patient_to_mpileups.get(patient_id, []) # Locate patients' all pileup files
        
        for _, row in mut_df.iterrows():
            chrom = row["Chrom"]
            pos = row["Position"]
            
            for mpileup_path in mpileup_files:
                sample_name = os.path.basename(mpileup_path).replace(".mpileup", "")
                out_path = os.path.join(DIR_filtered_mpileups, f"{chrom}_{pos}_{sample_name}.mpileup")
                
                # Append the command if force_run is set to true, or if the file doesn't exist. 
                if force_run:
                    grep_cmds.append(f"grep '{chrom}.*{pos}' {mpileup_path} > {out_path}")
                else:
                    if not os.path.exists(out_path):
                        grep_cmds.append(f"grep '{chrom}.*{pos}' {mpileup_path} > {out_path}")  
    
    run_grep_in_groups(grep_cmds, group_size=100)   
    
    # === Step 2: Build final matrix ===
    results = []    
    
    for patient_id, mut_df in mutations.groupby("Patient_id"):
        mpileup_files = patient_to_mpileups.get(patient_id, [])
        for _, row in mut_df.iterrows():
            chrom = row["Chrom"]
            pos = row["Position"]
            alt = row["Alt"]
            ref = row["Ref"]
            variant_type = row["Type"]
            mutation_id = f"{chrom}_{pos}_{ref}_{alt}_{variant_type}"
            
            for mpileup_path in mpileup_files:
                sample = os.path.basename(mpileup_path).replace(".mpileup", "")
                grep_output = os.path.join(DIR_filtered_mpileups, f"{chrom}_{pos}_{sample}.mpileup")
                
                alt_count = get_read_support_from_grep(grep_output, variant_type, alt_base=alt)
                depth = get_depth_from_pileup(grep_output)
                vaf = (alt_count / depth) * 100 if depth > 0 else 0
                
                results.append({
                    "Mutation_ID": mutation_id,
                    "Patient_id": patient_id,
                    "Gene": row["Gene"],
                    "Chrom": chrom,
                    "Position": pos,
                    "Ref": ref,
                    "Alt": alt,
                    "Type": variant_type,
                    "Sample": sample,
                    "Alt_count": alt_count,
                    "Depth": depth,
                    "VAF(%)": vaf
                })  
    
    # Save result matrix ===
    df_matrix = pd.DataFrame(results)
    df_matrix.to_csv(PATH_output, index=False)
    print(f"Saved matrix to {PATH_output}")


def run_main():
    parser = argparse.ArgumentParser(description="Perform dependent mutation calling.")
    parser.add_argument("--DIR_mpileups", required=True, help="Directory where all raw mpileup files are located.")
    parser.add_argument("--DIR_filtered_mpileups", required=True, help="Filtered mpileups after running grep will be saved here.")
    parser.add_argument("--PATH_sample_info", required=True, help="Path to four column sample information.")
    parser.add_argument("--PATH_mutations", required=True, help="Mutations list.")
    parser.add_argument("--PATH_output", required=True, help="Only output file.")
    parser.add_argument("--force_run", type=lambda x: x.lower() == 'true', required=True, help="Force run grep? Provide true or false")
    
    args = parser.parse_args()
    
    DIR_mpileups = args.DIR_mpileups
    DIR_filtered_mpileups=args.DIR_filtered_mpileups
    PATH_sample_info = args.PATH_sample_info
    PATH_mutations = args.PATH_mutations
    PATH_output = args.PATH_output
    force_run = args.force_run
    
    main(
        DIR_mpileups=DIR_mpileups,
        DIR_filtered_mpileups=DIR_filtered_mpileups,
        PATH_sample_info=PATH_sample_info,
        PATH_mutations=PATH_mutations,
        PATH_output=PATH_output,
        force_run=force_run
    )

if __name__ == "__main__":
    run_main()

"""
Run example:

python get_read_support_matrix.py \
  --DIR_mpileups path/to/mpileups \
  --DIR_filtered_mpileups path/to/filtered_mpileups \
  --PATH_sample_info path/to/sample_info.tsv \
  --PATH_mutations path/to/mutations.csv \
  --PATH_output path/to/output/per_sample_read_support.csv \
  --force_run true
"""
