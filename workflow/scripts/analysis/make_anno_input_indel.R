# This script modifies the VarScan output to make an annovar input file

library(dplyr)
library(argparse)

parser <- ArgumentParser(description = "Make the input files for running ANNOVAR on indel files from VarScan2.")

# passed from Snakemake into the script
parser$add_argument('--PATH_VarScan_indel', metavar='FILE', type='character', help="The output directory from VarScan2 that contains the indels. (input of this script)")
parser$add_argument('--ANNOVAR_indel_input', metavar='FILE', type='character', help="The input file we make to run ANNOVAR on the indel file outputted by VarScan2. (output of this script)")

args <- parser$parse_args()

# DIR_varscan_indel <- "/groups/wyattgrp/users/amunzur/chip_project/VarScan2_results/WBC_only/new_chip_panel_CP/indel"
# PATH_VarScan_indel <- "/groups/wyattgrp/users/amunzur/chip_project/VarScan2_results/WBC_only/new_chip_panel_CP/indel/GUBB-17-267-gDNA-Baseline-IDT-2017Oct04_S16.vcf"

# works on only one path. parallelize using lapply.
make_anno_input <- function(PATH_VarScan_indel, ANNOVAR_indel_input) {

	df <- as.data.frame(read.delim(PATH_VarScan_indel)) %>%
			mutate(Position = as.numeric(Position), 
					End = Position, 
					VarAllele = as.character(VarAllele), 
					Ref = as.character(Ref)) %>%
			select(Chrom, Position, End, Ref, VarAllele)

	names(df) <- c("Chrom", "Start", "End", "Ref", "Var")

	# deletions 
	idx <- grep("\\-", df$Var)
	df$Var[idx] <- "-"
	df$Var <- gsub("\\+", "", df$Var)

	# insertions
	idx <- grep("\\+", df$Var)
	df$Ref[idx] <- "-"

	write.table(df, ANNOVAR_indel_input, sep="\t",  col.names=FALSE, row.names=FALSE, quote = FALSE)
}

# lapply(as.list(list.files(DIR_varscan_indel, full.names = TRUE, pattern = "\\.vcf$")), make_anno_input)

make_anno_input(PATH_VarScan_indel = args$PATH_VarScan_indel,
	ANNOVAR_indel_input = args$ANNOVAR_indel_input)