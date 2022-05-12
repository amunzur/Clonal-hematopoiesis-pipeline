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

	df_main <- as.data.frame(read.delim(PATH_VarScan_indel)) %>%
		mutate(
			Var_type = case_when(
				grepl("-", VarAllele) ~ "DEL", 
				grepl("+", VarAllele) ~ "INS",
				TRUE ~ "FUCK"), 
			VarAllele = gsub("\\+|-", "", VarAllele), # remove - and + signs
			VarAllele = ifelse(Var_type == "DEL", "-", VarAllele), # if DEL, alt allele must be -. Otherwise (insertion) alt column doesnt change
			Ref = ifelse(Var_type == "INS", "-", Ref), # if INS the ref becomes -. Otherwise (DEL) it stays the same.
			Start = Position, 
			Start = ifelse(Var_type == "DEL", as.numeric(Position) + 1, as.numeric(Position)), # if DEL, add 1 to pos
			End = ifelse(Var_type == "DEL", as.numeric(Start) + nchar(VarAllele) - 1, as.numeric(Start))) %>% # if DEL we add 1 to the start  
		rename(Alt = VarAllele) %>%
		select(Chrom, Start, End, Ref, Alt)

	write.table(df_main, ANNOVAR_indel_input, sep="\t", col.names=FALSE, row.names=FALSE, quote = FALSE)
}

make_anno_input(PATH_VarScan_indel = args$PATH_VarScan_indel,
	ANNOVAR_indel_input = args$ANNOVAR_indel_input)