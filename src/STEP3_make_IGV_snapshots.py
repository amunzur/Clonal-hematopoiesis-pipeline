#!/usr/bin/env python3

"""
Generates IGV screenshots for manual variant curation.
"""

import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser("Generate IGV screenshots.")
parser.add_argument("PATH_variants_df", required=True, help="Path to variant calls to generate IGV screenshots.", type=str)
parser.add_argument("PATH_batch", required=True, help="IGV batch script will be saved here.", type=str)
parser.add_argument("DIR_snapshots", required=True, help="DIR to save IGV snapshots.", type=str)
parser.add_argument("given_range", required=False, default=100, help="IGV window range in basepairs.", type=int)
parser.add_argument("WBC_only", required=True, default=False, help="Set to False if calls are cfDNA and WBC paired. ", type=boolean)

args = parser.parse_args()

def make_igv_batch_script(PATH_variants_df, PATH_batch, DIR_snapshots, given_range, WBC_only):
	"""
	Generates an IGV batch script which can be used to automate IGV screenshot generation.
	PATH_variants_df should have these cols at minimum:
	"Patient_id", "Protein_annotation", "Gene", "Chrom", "Position", "Sample_name_t", "Path_bam_t", "Path_bam_n"
	when WBC_only is set to True Path_bam_t is not necessary.
	"""
	
	df= pd.read_csv(PATH_variants_df)
	
	for index, row in df.iterrows(): # iterate through each snv in the maf
		if WBC_only: 
			BAM_wbc = str(row["Path_bam_n"])
		else:
			BAM_wbc = str(row["Path_bam_n"])
			BAM_tumor = str(row["Path_bam_t"])
			
		position = int(row["Position"])
		
		if given_range: # if the user wants to see a given range, consider the end position of the range as well
			start_position = str(int(position - given_range/2))
			end_position = str(int(position + given_range/2))	
		
		os.makedirs(DIR_snapshots, exist_ok = True) # make the snapshost dir if it doesnt exist already
		if WBC_only:
			output_file_name = row["Gene"] + "_" + str(row["Protein_annotation"]) + "_" + row["Chrom"] + "_" + str(row["Position"]) + "_" + row["Sample_name"] + ".png" # one snapshot for each variant		
		else:
			output_file_name = row["Gene"] + "_" + str(row["Protein_annotation"]) + "_" + row["Chrom"] + "_" + str(row["Position"]) + "_" + row["Sample_name_t"] + ".png" # one snapshot for each variant		
		
		# Begin compiling the batch file
		with open(PATH_batch, 'a') as the_file:
			the_file.write('new\n')
			if index == 0: the_file.write('genome hg38\n') # only load the genome if it is the first time we are starting IGV
			if WBC_only:
				the_file.write(str("load " + BAM_wbc + '\n'))
			else:
				the_file.write(str("load " + BAM_tumor + '\n'))
				the_file.write(str("load " + BAM_wbc + '\n'))
			the_file.write(str("snapshotDirectory " + DIR_snapshots + '\n'))
			chrom = str(row["Chrom"])
			if given_range: the_file.write(str('goto ' + chrom + ":" + start_position + "-" + end_position + "\n"))
			else: the_file.write(str('goto ' + chrom + ":" + str(position) + "\n"))
			the_file.write('sort start location\n')
			the_file.write('sort base\n')
			the_file.write('maxPanelHeight -1\n')
			the_file.write(str('snapshot ' + output_file_name + '\n'))
			the_file.write("\n")
	
	with open(PATH_batch, 'a') as the_file: 
		the_file.write('exit') # append an exist statement at the end so that IGV closes on its own
	
	print("~/IGV_Linux_2.11.3/igv.sh --batch ", PATH_batch) # print the command needed to turn IGV to terminal

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate IGV screenshots.")
    parser.add_argument("PATH_variants_df", type=str, help="Path to variant calls to generate IGV screenshots.")
    parser.add_argument("PATH_batch", type=str, help="IGV batch script will be saved here.")
    parser.add_argument("DIR_snapshots", type=str, help="Directory to save IGV snapshots.")
    parser.add_argument("--given_range", type=int, default=100, help="IGV window range in basepairs.")
    parser.add_argument("--WBC_only", type=lambda x: x.lower() == 'true', default=False, help="Set to true if calls are WBC-only. Defaults to false.")

    args = parser.parse_args()

    if os.path.exists(args.PATH_batch):
        os.remove(args.PATH_batch)

    make_igv_batch_script(
        args.PATH_variants_df,
        args.PATH_batch,
        args.DIR_snapshots,
        args.given_range,
        args.WBC_only
    )

# Example use:
# python3 generate_igv_snapshots.py \
#     path/to/variants.csv \
#     path/to/igv_batch.txt \
#     path/to/output_screenshots \
#     --given_range 200 \
#     --WBC_only true
