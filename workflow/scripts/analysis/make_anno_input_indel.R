# This script modifies the VarScan output to make an annovar input file

library(dplyr)
library(argparse)

parser <- ArgumentParser(description = "Make the input files for running ANNOVAR on indel files from VarScan2.")

# passed from Snakemake into the script
parser$add_argument('--PATH_VarScan_indel', metavar='FILE', type='character', help="The output directory from VarScan2 that contains the indels. (input of this script)")
parser$add_argument('--ANNOVAR_indel_input', metavar='FILE', type='character', help="The input file we make to run ANNOVAR on the indel file outputted by VarScan2. (output of this script)")

args <- parser$parse_args()

# DIR_varscan_indel <- "/groups/wyattgrp/users/amunzur/chip_project/VarScan2_results/WBC_only/new_chip_panel_CP/indel"
# PATH_VarScan_indel <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/VarScan2/indel/batch5/GUBB-19-407_cfDNA_Baseline_IDT_2019Jun05.vcf"

# works on only one path. parallelize using lapply.
make_anno_input <- function(PATH_VarScan_indel, ANNOVAR_indel_input) {

# come back to this later to make this piece of code more readable
	df_main <- varscan_df %>%
		mutate(
			Var_type = case_when(
				grepl("-", Alt) ~ "DEL", 
				grepl("+", Alt) ~ "INS",
				TRUE ~ "FUCK"), 
			Alt = gsub("\\+|-", "", Alt), # remove - and + signs
			Alt = ifelse(Var_type == "DEL", "-", Alt), # if DEL, alt allele must be -. Otherwise (insertion) alt column doesnt change
			Ref = ifelse(Var_type == "INS", "-", Ref), # if INS the ref becomes -. Otherwise (DEL) it stays the same.
			Start = ifelse(Var_type == "DEL", as.numeric(Start) + 1, as.numeric(Start)), # if DEL, add 1 to pos
			End = ifelse(Var_type == "DEL", as.numeric(Start) + nchar(Alt) - 1, as.numeric(Start))) %>% # if DEL we add 1 to the start  
		rename(Chrom = Chr) %>%
		select(Chrom, Start, End, Ref, Alt)

	
	df_main <- as.data.frame(read.delim(PATH_VarScan_indel)) %>%
			mutate(
				Var_type = case_when(
					grepl("\\-", VarAllele) ~ "DEL", 
					grepl("\\+", VarAllele) ~ "INS",
					TRUE ~ "FUCK"), 
				VarAllele = gsub("\\+|-", "", VarAllele))

	# DELETIONS
	df_del <- df_main %>%
		filter(Var_type == "DEL") %>%
		mutate(End = Position) %>%
		select(Chrom, Position, End, VarAllele, Var_type, VarFreq) %>%
		mutate(
			Position = Position + 1, 
			End = as.numeric(Position) + nchar(VarAllele) - 1, 
			Var = "-", 
			VarFreq = readr::parse_number(as.character(VarFreq))) %>%
		arrange(Position, -VarFreq) %>%
  		filter(duplicated(Position) == FALSE) %>% # If position is duplicated, only keep the variant with the higher allele frequency
		rename(
			Start = Position,
			Ref = VarAllele) %>%
		select(Chrom, Start, End, Ref, Var)

	# INSERTIONS
	df_ins <- df_main %>%
		filter(Var_type == "INS") %>%
		mutate(
			End = Position, 
			Ref = "-", 
			VarFreq = readr::parse_number(as.character(VarFreq))) %>%
		select(Chrom, Position, End, Ref, VarAllele, VarFreq) %>%
		arrange(Position, -VarFreq) %>%
  		filter(duplicated(Position) == FALSE) %>% # If position is duplicated, only keep the variant with the higher allele frequency
		rename(
			Start = Position, 
			Var = VarAllele) %>%
		select(Chrom, Start, End, Ref, Var)

	# COMBINED
	df_main_combined <- rbind(df_ins, df_del) %>%
		arrange(Chrom, Start)

	write.table(df_main_combined, ANNOVAR_indel_input, sep="\t",  col.names=FALSE, row.names=FALSE, quote = FALSE)
}

# lapply(as.list(list.files(DIR_varscan_indel, full.names = TRUE, pattern = "\\.vcf$")), make_anno_input)

make_anno_input(PATH_VarScan_indel = args$PATH_VarScan_indel,
	ANNOVAR_indel_input = args$ANNOVAR_indel_input)