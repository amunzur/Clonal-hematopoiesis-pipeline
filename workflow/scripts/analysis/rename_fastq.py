import os 
import re
import subprocess
import pandas as pd

'''
Given a directory that contains the newly sequenced samples, rename them by matching the barcodes to the sequencing IDs. 
'''

def find_fastq_files(DIR_data): 
    '''
    Given a directory, find all files ending with either fq.gz or fastq.gz. Save their absolute paths to a list.
    '''
    file_names = []
    for root, dirs, files in os.walk(DIR_data):
        for file in files: 
            if file.endswith(("fastq.gz", "fq.gz")):
                file_names.append(os.path.join(root, file))
    return file_names

def find_barcode(str_list):
    '''
    Given a file path to a fastq file with the barcode somewhere in its path, return the barcode.
    '''
    letters = ["A", "C", "G", "T", "-"]
    barcode = [x for x in str_list if all(letter in letters for letter in x)][0]
    return(barcode)

def get_reverse_compliment(some_fastq):
    '''
    Sometimes the barcode includes a portion that is the reverse compliment. Check sample sheet and determine if this is the case. 
    This function identifies the reverse compliment of the second part of the barcode and returns a complete new file name as string.
    '''
    seq_id = some_fastq.split("_")[3].split("-")[1] # to take reverse compliment of 
    id_reverse = seq_id[::-1] # reverse
    mydict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'} # dict to take the reverse 
    id_reverse_compliment = []
    
    for char in  id_reverse:
        for pattern, repl in mydict.items():
            if char == pattern: # only continue if the pattern and the character we are at match
                s = re.sub(pattern, repl, char)
                id_reverse_compliment.append(s)
    
    id_reverse_compliment = "".join(id_reverse_compliment)
    
    # replace the id in the fastq name with the reverse compliment 
    some_fastq = some_fastq.replace(seq_id, id_reverse_compliment)	
    return(some_fastq)

def get_sequencing_id(some_fastq, PATH_barcode_sheet, reverse_compliment = False):
    '''
    Given a file path to a fastq file with barcodes, scan the barcodes sheet to find its sequencing ID and return it as a string.
    '''
    df_barcodes = pd.read_csv(PATH_barcode_sheet)    
    barcode = find_barcode(some_fastq.split("_")) # barcodes as they appear in the identifier sheets, we'll filter for those
    pool = os.path.basename(some_fastq).split("_")[0]
    if reverse_compliment: 
        barcode = get_reverse_compliment(some_fastq)
    seq_id = df_barcodes[(df_barcodes["Barcode"] == barcode) & (df_barcodes["GSC Pool ID"] == pool)]["Sequencing_ID"].values[0]
    
    if len([seq_id]) > 1: 
        raise ValueError(f"The barcode {barcode} from {os.path.basename(some_fastq)} matched to more than one sample. Specify pool.")
    elif len([seq_id]) == 0:
        raise ValueError(f"The barcode {barcode} from {os.path.basename(some_fastq)} didn't match any sample.")
    
    # Now consider if this is the 1st or the 2nd pair    
    suffix = "_1.fq.gz" if "_1_" in os.path.basename(some_fastq) else "_2.fq.gz"
    seq_id += suffix
    
    return {'FastQ_path': some_fastq, 'barcode': barcode, 'seq_id': seq_id}

def write_script(file_names, DIR_output, PATH_script):
    '''
    Writes a bash script to rename move the fastq files into a new location. You would then simply run that script to rename the fastqs.
    file_names: the list of dictionaries, outputted by the get_sequencing_id function.
    DIR_output: the dir to save the renamed fastq files.
    PATH_script: The bash script that will be generated to rename and move files.
    '''
    if os.path.exists(PATH_script):
        try:
            os.remove(PATH_script)
        except OSError as e:
            print(f"Error: {e.filename} - {e.strerror}")            
    for my_dict in file_names:
        source=my_dict["FastQ_path"]
        destination=os.path.join(DIR_output, my_dict["seq_id"])
        command="mv " + source + " " + destination
        with open(PATH_script, 'a') as file:
            file.write(command)
            file.write("\n")    


# FEB 26TH 2024 WGS SAMPLES
for sequencing_date in ["wgs_26feb2024", "wgs_29feb2024", "wgs_06march2024"]:
    DIR_data=f"/groups/wyattgrp/data/fq/{sequencing_date}"
    DIR_output="/groups/wyattgrp/data/fq/deep_wgs/fq"
    PATH_barcode_sheet="/groups/wyattgrp/users/amunzur/wgs/resources/sample_maps.tsv"
    PATH_script=f"/groups/wyattgrp/users/amunzur/wgs/scripts/batch_scripts/{sequencing_date}.bash"
    file_names = find_fastq_files(DIR_data)   
    file_names = [get_sequencing_id(some_fastq, PATH_barcode_sheet) for some_fastq in file_names]
    write_script(file_names, DIR_output, PATH_script)

# Doing this one separately
DIR_data=f"/groups/wyattgrp/data/fq/wgs_08march2024"
DIR_output="/groups/wyattgrp/data/fq/deep_wgs/fq_08march2024"
PATH_barcode_sheet="/groups/wyattgrp/users/amunzur/wgs/resources/sample_maps.tsv"
PATH_script=f"/groups/wyattgrp/users/amunzur/wgs/scripts/batch_scripts/wgs_08march2024.bash"
file_names = find_fastq_files(DIR_data)   
file_names = [get_sequencing_id(some_fastq, PATH_barcode_sheet) for some_fastq in file_names]
write_script(file_names, DIR_output, PATH_script)





# FEB 29TH 2024 WGS SAMPLES
DIR_data="/groups/wyattgrp/data/fq/wgs_29feb2024"
DIR_output="/groups/wyattgrp/data/fq/deep_wgs/fq"
PATH_barcode_sheet="/groups/wyattgrp/users/amunzur/wgs/resources/sample_maps.tsv"
PATH_script="/groups/wyattgrp/users/amunzur/wgs/scripts/batch_scripts/wgs_29feb2024.bash"
file_names = find_fastq_files(DIR_data)   
file_names = [get_sequencing_id(some_fastq, PATH_barcode_sheet) for some_fastq in file_names]
write_script(file_names, DIR_output, PATH_script)

# MATCH 6TH 2024 WGS SAMPLES
DIR_data="/groups/wyattgrp/data/fq/wgs_26feb2024"
DIR_output="/groups/wyattgrp/data/fq/deep_wgs/fq"
PATH_barcode_sheet="/groups/wyattgrp/users/amunzur/wgs/resources/sample_maps.tsv"
PATH_script="/groups/wyattgrp/users/amunzur/wgs/scripts/batch_scripts/wgs_26feb2024.bash"
file_names = find_fastq_files(DIR_data)   
file_names = [get_sequencing_id(some_fastq, PATH_barcode_sheet) for some_fastq in file_names]
write_script(file_names, DIR_output, PATH_script)

# MARCH 8TH 2024 WGS SAMPLES
DIR_data="/groups/wyattgrp/data/fq/wgs_26feb2024"
DIR_output="/groups/wyattgrp/data/fq/deep_wgs/fq"
PATH_barcode_sheet="/groups/wyattgrp/users/amunzur/wgs/resources/sample_maps.tsv"
PATH_script="/groups/wyattgrp/users/amunzur/wgs/scripts/batch_scripts/wgs_26feb2024.bash"
file_names = find_fastq_files(DIR_data)   
file_names = [get_sequencing_id(some_fastq, PATH_barcode_sheet) for some_fastq in file_names]
write_script(file_names, DIR_output, PATH_script)


# POOL 6
DIR_data="/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/merged/GU-BCa_CHIP_28Apr2023_Pool6" # directory where the fastq files are saved, right out of the sequencer
DIR_output="/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/merged"
PATH_barcode_sheet="/groups/wyattgrp/users/amunzur/pipeline/resources/sample_maps/GU-BCa_CHIP_28Apr2023_Pool6_barcodes.tsv"
PATH_script="/groups/wyattgrp/users/amunzur/pipeline/resources/bash_scripts/GU-BCa_CHIP_28Apr2023_Pool6.bash"
file_names = find_fastq_files(DIR_data)   
file_names = [get_sequencing_id(some_fastq, PATH_barcode_sheet) for some_fastq in file_names]
write_script(file_names,  DIR_output, PATH_script)

# POOL 5
DIR_data="/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/merged/GU-BCa_CHIP_Pool5"
DIR_output="/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/merged"
PATH_barcode_sheet="/groups/wyattgrp/users/amunzur/pipeline/resources/sample_maps/GU-BCa_CHIP_Pool5_barcodes.tsv"
PATH_script="/groups/wyattgrp/users/amunzur/pipeline/resources/bash_scripts/GU-BCa_CHIP_Pool5.bash"
file_names = find_fastq_files(DIR_data)   
file_names = [get_sequencing_id(some_fastq, PATH_barcode_sheet) for some_fastq in file_names]
write_script(file_names,  DIR_output, PATH_script)

# POOLS 1,2,3,4
DIR_data="/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/merged/GU-BCa_CHIP_26Apr2023_Pool1_2_3_4"
DIR_output="/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/merged"
PATH_barcode_sheet="/groups/wyattgrp/users/amunzur/pipeline/resources/sample_maps/GU-BCa_CHIP_28Apr2023_Pool1_2_3_4_barcodes.tsv"
PATH_script="/groups/wyattgrp/users/amunzur/pipeline/resources/bash_scripts/GU-BCa_CHIP_Pool1_2_3_4.bash"
file_names = find_fastq_files(DIR_data)   
file_names = [get_sequencing_id(some_fastq, PATH_barcode_sheet) for some_fastq in file_names]
write_script(file_names,  DIR_output, PATH_script)

# Kidney OT samples
DIR_data="/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/kidney_OT"
DIR_output="/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/merged/kidney"
PATH_barcode_sheet="/groups/wyattgrp/users/amunzur/pipeline/resources/sample_maps/kidney_OT_Jan31_2024.tsv"
PATH_script="/groups/wyattgrp/users/amunzur/pipeline/resources/bash_scripts/kidneyOT.bash"
file_names = find_fastq_files(DIR_data)   
file_names = [get_sequencing_id(some_fastq, PATH_barcode_sheet) for some_fastq in file_names]
write_script(file_names,  DIR_output, PATH_script)