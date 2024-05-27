# Goal of this script is to concat fastq files. Some samples were resequenced, the gDNA smaples, concatting them increases the depth.

import os
import glob
import pandas as pd
import fnmatch

sample_sheet = pd.read_csv("/groups/wyattgrp/users/amunzur/data/sample_sheet.csv")
DIR_samples = "/groups/wyattgrp/users/amunzur/pipeline/results/data/fastq/seq_runs"

DIR_seq1 = os.path.join(DIR_samples, "seqRun1")
DIR_seq2 = os.path.join(DIR_samples, "seqRun2")

os.listdir(DIR_seq1)
os.listdir(DIR_seq2)

matchers = sample_sheet[(sample_sheet["Sequencing batch"] == 2)]["Sequencing ID"].tolist() # list of patients that got resequenced
matchers = ["_".join(x.split("_")[0:3]) for x in samples] # strip off the date in case there is a mistake in the dates, being extra careful here

# the matching samples in seqDir1
matching = [x for x in A if any(b in x for b in matchers)]

matching = [x for x in os.listdir(DIR_seq1) if any(b in x for b in os.listdir(DIR_seq2))]



A = ["_".join(x.split("_")[0:3]) for x in os.listdir(DIR_seq1)] # strip off the date in case there is a mistake in the dates, being extra careful here
A= list(set(A))
