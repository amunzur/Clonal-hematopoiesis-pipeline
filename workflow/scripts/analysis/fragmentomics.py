"""
This script compares fragmentomics profiles of samples with and without CH.
"""

import pandas as pd
import numpy as np
import os
import pysam

# conda activate pysam_environment

def get_fragments_with_deletion(bam_file, chromosome, position):
    """
    Returns the fragment lengths of the reads that have a deletion event at a specific genomic position.
    """
    read_info_with_deletion = []  # List to store read information with deletion
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate over reads in the BAM file
        for read in bam.fetch(chromosome, position - 1, position):
            ref_pos = read.reference_start  # Initialize reference position
            fragment_id = read.query_name.split('/')[0]  # Extract fragment ID from read name
            fragment_length = None  # Placeholder for fragment length
            # Iterate over the CIGAR tuples to find deletions
            for op, length in read.cigartuples:
                if op == 2:  # Deletion
                    # Check if the specified position falls within this deletion
                    if ref_pos <= position - 1 < ref_pos + length:
                        # Get the fragment length if not already obtained
                        if fragment_length is None:
                            if read.is_paired:
                                fragment_length = abs(read.template_length)
                        # Add read information to the list
                        read_info_with_deletion.append({
                            'read_id': read.query_name,
                            'fragment_id': fragment_id,
                            'fragment_length': fragment_length
                        })
                        break  # Stop checking further CIGAR tuples for this read
                # Update the reference position
                if op in {0, 2, 3, 7, 8}:  # M, D, N, =, X
                    ref_pos += length
                elif op == 1 or op == 4:  # I or S
                    continue
    fragment_lengths = pd.DataFrame(read_info_with_deletion).drop_duplicates("fragment_id")["fragment_length"].tolist()
    return(fragment_lengths)

def get_fragments_with_snv(bam_file, chromosome, position, ref_base, alt_base):
    """
    Returns the fragment lengths of the reads that have an SNV event (either ref_base or alt_base)
    at a specific genomic position.
    """
    read_info_with_snv = []  # List to store read information with SNV
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate over reads in the BAM file
        for read in bam.fetch(chromosome, position - 1, position):
            ref_pos = read.reference_start  # Initialize reference position
            read_pos = 0  # Initialize read position
            fragment_id = read.query_name.split('/')[0]  # Extract fragment ID from read name
            fragment_length = None  # Placeholder for fragment length
            # Iterate over the CIGAR tuples to find matches and mismatches
            for op, length in read.cigartuples:
                if op == 0 or op == 7 or op == 8:  # Match or Mismatch
                    for i in range(length):
                        if ref_pos == position - 1:
                            read_base = read.query_sequence[read_pos]
                            if read_base == alt_base:
                                # Get the fragment length if not already obtained
                                if fragment_length is None:
                                    if read.is_paired:
                                        fragment_length = abs(read.template_length)
                                # Add read information to the list
                                read_info_with_snv.append({
                                    'read_id': read.query_name,
                                    'fragment_id': fragment_id,
                                    'fragment_length': fragment_length
                                })
                                break  # Stop checking further CIGAR tuples for this read
                        ref_pos += 1
                        read_pos += 1
                    if ref_pos > position - 1:
                        break  # If we've passed the position, no need to check further
                elif op == 2 or op == 3:  # Deletion or N (skip region in reference)
                    ref_pos += length
                elif op == 1:  # Insertion
                    read_pos += length
                elif op == 4 or op == 5:  # Soft clipping or hard clipping
                    read_pos += length
    fragment_lengths = pd.DataFrame(read_info_with_snv).drop_duplicates("fragment_id")["fragment_length"].tolist()
    return fragment_lengths

def get_fragments_with_insertion(bam_file, chromosome, position):
    """
    Returns the fragment lengths of the reads that have an insertion event at a specific genomic position.
    """
    read_info_with_insertion = []  # List to store read information with insertion
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate over reads in the BAM file
        for read in bam.fetch(chromosome, position - 1, position):
            ref_pos = read.reference_start  # Initialize reference position
            read_pos = 0  # Initialize read position
            fragment_id = read.query_name.split('/')[0]  # Extract fragment ID from read name
            fragment_length = None  # Placeholder for fragment length
            # Iterate over the CIGAR tuples to find insertions
            for op, length in read.cigartuples:
                if op == 1:  # Insertion
                    # Get the fragment length if not already obtained
                    if fragment_length is None:
                        if read.is_paired:
                            fragment_length = abs(read.template_length)
                    # Add read information to the list
                    read_info_with_insertion.append({
                        'read_id': read.query_name,
                        'fragment_id': fragment_id,
                        'fragment_length': fragment_length
                    })
                    break  # Stop checking further CIGAR tuples for this read
                # Update the reference position and read position
                if op in {0, 7, 8}:  # M, =, X
                    ref_pos += length
                    read_pos += length
                elif op == 2 or op == 3:  # D, N
                    ref_pos += length
                elif op == 1:  # I
                    read_pos += length
                elif op == 4 or op == 5:  # S, H
                    read_pos += length
    fragment_lengths = pd.DataFrame(read_info_with_insertion).drop_duplicates("fragment_id")["fragment_length"].tolist()
    return fragment_lengths

def get_fragments_without_variant(bam_file, chromosome, position, variant_type, ref_base=None, alt_base=None):
    """
    Returns the fragment lengths of the reads that do not have a variant (either SNV, deletion, or insertion)
    at a specific genomic position.
    """
    read_info_without_variant = []  # List to store read information without variant
    wild_type_count = 0  # Counter
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate over reads in the BAM file
        for read in bam.fetch(chromosome, position - 1, position):
            ref_pos = read.reference_start  # Initialize reference position
            read_pos = 0  # Initialize read position
            fragment_id = read.query_name.split('/')[0]  # Extract fragment ID from read name
            fragment_length = None  # Placeholder for fragment length
            has_variant = False  # Flag to check if read has variant
            # Iterate over the CIGAR tuples to find matches and mismatches
            for op, length in read.cigartuples:
                if op == 0 or op == 7 or op == 8:  # Match or Mismatch
                    for i in range(length):
                        if ref_pos == position - 1:
                            read_base = read.query_sequence[read_pos]
                            if variant_type.lower() == "snv":
                                if read_base == alt_base or read_base != ref_base:
                                    has_variant = True
                                    break
                                break
                        ref_pos += 1
                        read_pos += 1
                    if ref_pos > position - 1:
                        break  # If we've passed the position, no need to check further
                elif op == 2 or op == 3:  # Deletion or N (skip region in reference)
                    if ref_pos <= position - 1 < ref_pos + length:
                        if variant_type.lower() == "deletion":
                            has_variant = True
                            break
                    ref_pos += length
                elif op == 1:  # Insertion
                    # Check if the insertion spans the specified position
                    if ref_pos <= position - 1 < ref_pos + length:
                        if variant_type.lower() == "insertion":
                            has_variant = True
                            break
                    read_pos += length  # Move to the next position in the read
                elif op == 4 or op == 5:  # Soft clipping or hard clipping
                    read_pos += length
            if not has_variant:
                wild_type_count += 1  
                # Get the fragment length if not already obtained
                if fragment_length is None:
                    if read.is_paired:
                        fragment_length = abs(read.template_length)
                # Add read information to the list
                read_info_without_variant.append({
                    'read_id': read.query_name,
                    'fragment_id': fragment_id,
                    'fragment_length': fragment_length
                })
    fragment_lengths = pd.DataFrame(read_info_without_variant).drop_duplicates("fragment_id")["fragment_length"].tolist()
    return fragment_lengths

DIR_working = "/groups/wyattgrp/users/amunzur/pipeline"
PATH_sample_information = os.path.join(DIR_working, "resources/sample_lists/sample_information.tsv")
DIR_output = "/groups/wyattgrp/users/amunzur/pipeline/results/fragmentomics"

sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
sample_info = sample_info[sample_info["Diagnosis"] != "Healthy"]

all_vars_chip = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/chip_SSCS2_curated_complete.csv"))
all_vars_chip = sample_info.merge(all_vars_chip, how = "inner")

all_vars_ctDNA = pd.read_csv(os.path.join(DIR_working, "results/variant_calling/somatic_SSCS2_curated_complete.csv"))
all_vars_ctDNA = sample_info.merge(all_vars_ctDNA, how = "inner")

#######################################################################
# CH mutations
#######################################################################
# Get the fragment lengths from CH mutations
# these lists keep info from the entire cohort. In the loop we also collect info on a per sample and per position basis.
deletion_fragments = []
insertion_fragments = []
snv_fragments = []
row_list = []
for i, row in all_vars_chip.iterrows():
    if not row["Dependent"]:
        print(i)
        chromosome = row["Chrom"]
        position = row["Position"]
        event_type = row["Type"]
        sample_name = row["Sample_name_t"]
        bam_file = os.path.join(DIR_working, "results/data/bam/SSCS2_final", sample_name + ".bam")
        if event_type == "SNV":
            ref_base = row["Ref"]
            alt_base = row["Alt"]
            ch_fragments = get_fragments_with_snv(bam_file, chromosome, position, ref_base, alt_base)
            snv_fragments.extend(ch_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position, variant_type=event_type, ref_base=ref_base, alt_base=alt_base)
        elif event_type == "Deletion":
            ch_fragments = get_fragments_with_deletion(bam_file, chromosome, position)
            deletion_fragments.extend(ch_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position, variant_type=event_type)
        else:
            ch_fragments = get_fragments_with_insertion(bam_file, chromosome, position)
            insertion_fragments.extend(ch_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position+1, variant_type=event_type)
        # Add the fragment lengths to the row
        row["CH fragments"] = ch_fragments
        row["WT fragments"] = wt_fragments
        row_list.append(row)      

df = pd.DataFrame(row_list)
df.to_csv(os.path.join(DIR_output, "CH.csv"))

#######################################################################
# ctDNA mutations
#######################################################################
# Get the fragment lengths from ctDNA mutations
ctDNA_deletion_fragments = []
ctDNA_insertion_fragments = []
ctDNA_snv_fragments = []
row_list = []
for i, row in all_vars_ctDNA.iterrows():
    if not row["Dependent"]:
        print(i)
        chromosome = row["Chrom"]
        position = row["Position"]
        event_type = row["Type"]
        sample_name = row["Sample_name_t"]
        bam_file = os.path.join(DIR_working, "results/data/bam/SSCS2_final", sample_name + ".bam")
        if row["Sample_name_t"] == "GU-23-022_cfDNA-Baseline-IDT-2023Jan11" and position == 26696588:
            position += 1
        if event_type == "SNV":
            ref_base = row["Ref"]
            alt_base = row["Alt"]
            ctDNA_fragments = get_fragments_with_snv(bam_file, chromosome, position, ref_base, alt_base)
            ctDNA_snv_fragments.extend(ctDNA_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position, variant_type=event_type, ref_base=ref_base, alt_base=alt_base)
        elif event_type == "Deletion":
            ctDNA_fragments = get_fragments_with_deletion(bam_file, chromosome, position)
            ctDNA_deletion_fragments.extend(ctDNA_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position, variant_type=event_type)
        else:
            ctDNA_fragments = get_fragments_with_insertion(bam_file, chromosome, position)
            ctDNA_insertion_fragments.extend(ctDNA_fragments)
            wt_fragments = get_fragments_without_variant(bam_file, chromosome, position+1, variant_type=event_type)
        # Add the fragment lengths to the row
        row["ctDNA fragments"] = ctDNA_fragments
        row["WT fragments"] = wt_fragments
        row_list.append(row)      

df = pd.DataFrame(row_list)
df.to_csv(os.path.join(DIR_output, "ctDNA.csv"))

