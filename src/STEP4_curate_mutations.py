import numpy as np
import pandas as pd 
import os
import re
import argparse

"""
After curating mutations by going through IGV snapshots, this script detects which mutations have been removed, and then adds them to the blacklist.
"""

def do_misc_filtering_function(df):
    """
    Performs additional filtering by removing variants based on ratio between cfDNA and WBC VAF,
    and removed mutations in LPAR6.
    """
    df = df[df["tumor_wbc_vaf_ratio"] < 20]
    df = df[df["Gene"] != "LPAR6"]
    
    df = df.reset_index(drop = True)
    
    return(df)

def curate(PATH_muts, DIR_curated_screenshots, PATH_retained_variants, PATH_excluded_variants, PATH_sample_information, do_misc_filtering): 
    
    all_muts = pd.read_csv(PATH_muts)
    
    # Subset to patients if sample information is provided
    if PATH_sample_information and os.path.exists(PATH_sample_information):
        
        # Load sample info df
        sample_info = pd.read_csv(
            PATH_sample_information, sep="\t",
            names=["Patient_id", "Date_collected", "Diagnosis", "Timepoint"]
        )[
            ["Patient_id", "Diagnosis"]
        ].drop_duplicates()

        pts = sample_info["Patient_id"].unique()
        all_muts = all_muts[all_muts["Patient_id"].isin(pts)]
        
    if "Sample_name_t" in all_muts.columns:
        all_muts["IGV_screenshot_name"] = all_muts.apply(lambda row: "_".join(map(str, row[['Gene', 'Protein_annotation', 'Chrom', 'Position', 'Sample_name_t']])), axis=1) + ".png" # all called muts
    else:
        all_muts["IGV_screenshot_name"] = all_muts.apply(lambda row: "_".join(map(str, row[['Gene', 'Protein_annotation', 'Chrom', 'Position', 'Sample_name']])), axis=1) + ".png" # all called muts
    
    curated_screenshots = [f for f in os.listdir(DIR_curated_screenshots) if f.endswith(".png")]
    curated_muts = pd.DataFrame({"IGV_status": "keep", "IGV_screenshot_name": curated_screenshots}) # muts we are keeping
    merged = all_muts.merge(curated_muts, how = "left")
    
    muts_to_keep = merged[merged["IGV_status"] == "keep"]
    muts_to_exclude = merged[pd.isnull(merged["IGV_status"])]
    
    muts_to_keep.loc[muts_to_keep["Gene"].str.contains("U2AF"), "Gene"] = "U2AF1"
    
    del muts_to_keep["IGV_status"]
    del muts_to_exclude["IGV_status"]
        
    if do_misc_filtering: 
        muts_to_keep = do_misc_filtering_function(muts_to_keep)
    
    muts_to_exclude.to_csv(PATH_excluded_variants, index = False)
    muts_to_keep.to_csv(PATH_retained_variants, index = False)
    
    print(f"Curated variants saved to {PATH_retained_variants}")
    print(f"Excluded {len(muts_to_exclude)} variants (not present in curated screenshots).")

def main():
    parser = argparse.ArgumentParser(description="Curate mutations after IGV review.")
    parser.add_argument("--do_misc_filtering", action="store_true", help="Perform additional filtering.")
    parser.add_argument("--DIR_curated_screenshots", required=True, help="")
    parser.add_argument("--PATH_muts", required=True, help="Path to variants df to curate.")
    parser.add_argument("--PATH_excluded_variants", required=True, help="Excluded variants will be written here.")
    parser.add_argument("--PATH_retained_variants", required=True, help="Retained variants after IGV curation will be saved here.")
    parser.add_argument("--PATH_sample_information", required=True, help="If provided only patients present in this file will be included in the final curated variants list.")
    
    args = parser.parse_args()
    
    do_misc_filtering = args.do_misc_filtering
    DIR_curated_screenshots=args.DIR_curated_screenshots
    PATH_muts = args.PATH_muts
    PATH_excluded_variants = args.PATH_excluded_variants
    PATH_retained_variants = args.PATH_retained_variants
    PATH_sample_information = args.PATH_sample_information
        
    curate(
        PATH_muts=PATH_muts,
        DIR_curated_screenshots=DIR_curated_screenshots,
        PATH_retained_variants=PATH_retained_variants,
        PATH_excluded_variants=PATH_excluded_variants,
        PATH_sample_information=PATH_sample_information,
        do_misc_filtering=do_misc_filtering,
    )

if __name__ == "__main__":
    main()

"""
Run example:

python STEP4_curate_mutations.py \
    --do_misc_filtering \
    --DIR_curated_screenshots /path/to/curated_chip_variants/ \
    --PATH_muts /path/to/min_alt_reads_5_CHIP_after_filtering.csv \
    --PATH_excluded_variants /path/to/excluded_variants.csv \
    --PATH_retained_variants /path/to/retained_variants.csv \
    --PATH_sample_information /path/to/sample_information.tsv
"""
