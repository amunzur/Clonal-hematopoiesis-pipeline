#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 13:45:38 2021
@author: amunzur
"""
import os
import pandas as pd
import shutil
from shutil import copyfile

def make_igv_batch_script(variants_df, PATH_batch, DIR_snapshots, DIR_vancouver_bams, DIR_finland_bams, add_prefix, add_suffix, given_range):

	for index, row in variants_df.iterrows(): # iterate through each snv in the maf

		BAM_vancouver = os.path.join(DIR_vancouver_bams, str(add_prefix + row["Sample_name"] + add_suffix))
		BAM_finland = os.path.join(DIR_finland_bams, str(add_prefix + row["Sample_name_finland"] + add_suffix))
				
		position = int(row["Position"])
		if given_range: # if the user wants to see a given range, consider the end position of the range as well
			start_position = str(int(position - given_range/2))
			end_position = str(int(position + given_range/2))
			
		# os.remove(IGV_script)
		os.makedirs(DIR_snapshots, exist_ok=True) # make the snapshost dir if it doesnt exist already
		output_file_name = row["Sample_name"] + "_" + str(row["Chrom"]) + "_" + str(position) + ".png" # one snapshot for each snv
		
		with open(PATH_batch, 'a') as the_file:
			the_file.write('new\n')
			if index == 0: the_file.write('genome hg38\n') # only load the genome if it is the first time we are starting IGV
			the_file.write(str("load " + BAM_finland + '\n'))
			the_file.write(str("load " + BAM_vancouver + '\n'))
			the_file.write(str("snapshotDirectory " + DIR_snapshots + '\n'))
				
			chrom = str(row["Chrom"])
			if given_range: the_file.write(str('goto ' + chrom + ":" + start_position + "-" + end_position + "\n"))
			else: the_file.write(str('goto ' + chrom + ":" + str(position) + "\n"))
			the_file.write('sort start location\n')
			the_file.write('sort base\n')
			the_file.write(str('snapshot ' + output_file_name + '\n'))
			the_file.write("\n")

	with open(PATH_batch, 'a') as the_file: the_file.write('exit') # append an exist statement at the end so that IGV closes on its own

			
	print(PATH_batch) # print to terminal

given_range = 200
add_prefix = ""
add_suffix = ".bam"
cohort_name = "new_chip_panel"

PATH_varscan = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined", cohort_name, "varscan.csv")
PATH_vardict = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined", cohort_name, "vardict.csv")

DIR_vancouver_bams = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/data/bam", cohort_name, "SC_penalty")
DIR_finland_bams = "/groups/wyattgrp/data/bam/kidney"

varscan = pd.read_csv(PATH_varscan)
vardict = pd.read_csv(PATH_vardict)

combined = pd.concat([varscan, vardict], ignore_index = True)
combined = combined.drop_duplicates(subset = ["Sample_name", "Patient_ID", "Chrom", "Position", "Ref", "Alt"], ignore_index = True)
combined.loc[(combined["Variant"] == "deletion") | (combined["Variant"] == "insertion"), "Position"] = combined.Position + 1 # must add 1 so that IGV sorts correctly

# Alert samples
combined_ALERT = combined.loc[combined['Status'] == "ALERT"] 
IGV_script_ALERT = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/IGV_batch_scripts", cohort_name + "_ALERT.txt") # output path, the batch script to be used to run in IGV
DIR_snap_ALERT = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/figures/IGV_snapshots", cohort_name + "ALERT")

# OK samples 
combined_OK = combined.loc[combined['Status'] == "OK"] 
IGV_script_OK = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/IGV_batch_scripts", cohort_name + "_OK.txt") # output path, the batch script to be used to run in IGV
DIR_snap_OK = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/figures/IGV_snapshots", cohort_name + "OK")

# Great samples
combined_Great = combined.loc[combined['Status'] == "Great"] 
IGV_script_Great = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/IGV_batch_scripts", cohort_name + "_Great.txt") # output path, the batch script to be used to run in IGV
DIR_snap_Great = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/figures/IGV_snapshots", cohort_name + "Great")

make_igv_batch_script(combined_ALERT, IGV_script_ALERT, DIR_snap_ALERT, DIR_vancouver_bams, DIR_finland_bams, add_prefix, add_suffix, given_range)
make_igv_batch_script(combined_OK, IGV_script_OK, DIR_snap_OK, DIR_vancouver_bams, DIR_finland_bams, add_prefix, add_suffix, given_range)
make_igv_batch_script(combined_Great, IGV_script_Great, DIR_snap_Great, DIR_vancouver_bams, DIR_finland_bams, add_prefix, add_suffix, given_range)

# to run this script execute on the terminal
# /home/amunzur/IGV_Linux_2.11.3/igv.sh --batch /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/IGV_batch_scripts/new_chip_panel_ALERT.txt
# /home/amunzur/IGV_Linux_2.11.3/igv.sh --batch /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/IGV_batch_scripts/new_chip_panel_OK.txt
# /home/amunzur/IGV_Linux_2.11.3/igv.sh --batch /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/IGV_batch_scripts/new_chip_panel_Great.txt

# Now this next bit is for making snapshots of vardict alert samples, they look a bit suspicious. 
# Only retain ALERT vars called by vardict only.
variants_df = vardict[(vardict.Status == "ALERT") & (vardict.variant_caller == "Vardict")]
PATH_batch = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/IGV_batch_scripts", cohort_name + "_vardict_only_ALERT.txt")
DIR_snapshots = os.path.join("/groups/wyattgrp/users/amunzur/pipeline/results/figures/IGV_snapshots", cohort_name + "_vardict_only_ALERT")

make_igv_batch_script(variants_df, PATH_batch, DIR_snapshots, DIR_vancouver_bams, DIR_finland_bams, add_prefix, add_suffix, given_range)
# /home/amunzur/IGV_Linux_2.11.3/igv.sh --batch /groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/IGV_batch_scripts/new_chip_panel_vardict_only_ALERT.txt