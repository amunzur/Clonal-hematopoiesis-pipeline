# Adds read suffixes to read names in FastQ files.
import os 
import sys

PATH_input = sys.argv[1]
PATH_output = os.path.join(os.path.dirname(PATH_input) + "_readsuffix", os.path.basename(PATH_input))

# PATH_input = "/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/trimmed/VIP-Normal-020_gDNA__IDT_2021Apr30_R2_masked_val_2.fq"
# PATH_output = "/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/trimmed/test_R2.fastq"

if not os.path.exists(os.path.dirname(PATH_output)):
    os.makedirs(os.path.dirname(PATH_output))

fastq = open(PATH_input, "r").readlines()
reads_names = fastq[::4]

if "R1" in PATH_input:
    reads_names_new = [x.replace(" ", "/1 ") for x in reads_names]
elif "R2" in PATH_input:
    reads_names_new = [x.replace(" ", "/2 ") for x in reads_names]

# To remove the suffix:
# if "R1" in PATH_input:
#     reads_names_new = [x.replace("/1 ", " ") for x in reads_names]
# elif "R2" in PATH_input:
#     reads_names_new = [x.replace("/2 ", " ") for x in reads_names]


fastq[::4] = reads_names_new

if os.path.exists(PATH_output): 
    os.remove(PATH_output)

with open(PATH_output, "a") as output:
    output.writelines(fastq)