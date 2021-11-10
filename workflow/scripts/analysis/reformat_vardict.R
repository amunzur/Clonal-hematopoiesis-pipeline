library(tidyverse)
library(stringr)
library(argparse)

parser <- ArgumentParser(description = "Make the input files for running ANNOVAR on indel files from VarScan2.")

# passed from Snakemake into the script
parser$add_argument('--PATH_Vardict_output', metavar='FILE', type='character', help="The output directory from Vardict that contains the indels. (input of this script)")
parser$add_argument('--PATH_Vardict_reformatted', metavar='FILE', type='character', help="The dir that contains the reformatted Vardict outputs, output of this dir.")

args <- parser$parse_args()


# PATH_Vardict_output <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/new_chip_panel/GUBB-17-185-gDNA-Baseline-IDT-2017Aug30_S13.vcf"

reformat_vardict_output <- function(PATH_Vardict_output, PATH_Vardict_reformatted){

	tmp_vcf <- readLines(PATH_Vardict_output)
	tmp_vcf_data <- read.table(PATH_Vardict_output, stringsAsFactors = FALSE)

	# filter for the columns names
	tmp_vcf <- tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
	vcf_names <- unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
	names(tmp_vcf_data) <- vcf_names

	info_col_names <- c("SAMPLE", "TYPE", "DP", "VD", "AF", "BIAS", "REFBIAS", "VARBIAS", "PMEAN", "PSTD", "QUAL", "QSTD", "SBF", "ODDRATIO", "MQ", "SN", "HIAF", "ADJAF", "SHIFT3", "MSI", "MSILEN", "NM", "HICNT", "HICOV", "LSEQ", "RSEQ", "DUPRATE", "SPLITREAD", "SPANPAIR")

	# make a df with the info cols only
	info_df <- tmp_vcf_data %>%
			separate(col = INFO, sep = ";",	into = info_col_names, remove = FALSE) %>%
			select(info_col_names)

	# string replacement in the new df with the info information about variants 
	info_df <- as.data.frame(apply(info_df, 2, function(x) gsub(".*=", "", x)))

	# combine the two data frames
	var_df <- select(tmp_vcf_data, "#CHROM", "POS", "REF", "ALT", "FILTER")
	df <- as.data.frame(cbind(var_df, info_df)) %>%
			select("SAMPLE", "#CHROM", "POS", "REF", "ALT", "AF", "DP", "VD", "MQ", "TYPE", "FILTER")
	
	names(df) <- c("Sample_name", "Chr", "Start", "Ref", "Alt", "VarFreq", "Reads1", "Reads2", "Mapping_quality", "variant", "FILTER")
	df$variant <- tolower(df$variant)

	write_delim(df, PATH_Vardict_reformatted, delim = "\t")
}

reformat_vardict_output(PATH_Vardict_output = args$PATH_Vardict_output, 
						PATH_Vardict_reformatted = args$PATH_Vardict_reformatted)