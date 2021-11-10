library(tidyverse)
library(argparse)
library(stringr)

parser <- ArgumentParser(description = "Make the input files for running ANNOVAR on indels and snvs from Vardict.")

# passed from Snakemake into the script
parser$add_argument('--PATH_Vardict_output', metavar='FILE', type='character', help="The output file from Vardict that contains the variants. (input of this script)")
parser$add_argument('--PATH_ANNOVAR_input', metavar='FILE', type='character', help="The input file we make to run ANNOVAR on the variants file outputted by Vardict. (output of this script)")

args <- parser$parse_args()

# PATH_Vardict_output <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/new_chip_panel/GUBB-17-012-gDNA-Baseline-IDT-2017Jan17_S1.vcf"
# PATH_ANNOVAR_input <- "/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_outputs/Vardict/new_chip_panel/GUBB-17-012-gDNA-Baseline-IDT-2017Jan17_S1_anno.tsv"

make_anno_input <- function(PATH_Vardict_output, PATH_ANNOVAR_input) {

	tmp_vcf <- readLines(PATH_Vardict_output)
	tmp_vcf_data <- read.table(PATH_Vardict_output, stringsAsFactors = FALSE)

	# filter for the columns names
	tmp_vcf <- tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
	vcf_names <- unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
	names(tmp_vcf_data) <- vcf_names

	df <- tmp_vcf_data %>%
			mutate(End = POS, 
					variant_type = str_extract(INFO, "Deletion|Insertion|SNV|Complex")) %>%
			select("#CHROM", "POS", "End", "REF", "ALT", "variant_type")

	names(df) <- c("Chrom", "Start", "End", "Ref", "Var", "variant_type")

	# deletions 
	idx <- grep("Deletion", df$variant_type)
	df$Var[idx] <- "-"

	# fix the positions
	df$End <- df$Start + unlist(lapply(as.list(df$Ref), nchar))
	df$End <- df$End -1

	# insertions
	# idx <- grep("Insertion", df$variant_type)
	# df$Ref[idx] <- "-"

	df$variant_type <- NULL
	write.table(df, PATH_ANNOVAR_input, sep="\t",  col.names=FALSE, row.names=FALSE, quote = FALSE)
}

make_anno_input(PATH_Vardict_output = args$PATH_Vardict_output,
	PATH_ANNOVAR_input = args$PATH_ANNOVAR_input)