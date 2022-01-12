library(tidyverse)
library(stringr)
library(argparse)

parser <- ArgumentParser(description = "Make the input files for running ANNOVAR on files from GATK.")

# passed from Snakemake into the script
parser$add_argument('--PATH_GenotypeVCFs_output', metavar='FILE', type='character', help="GATK output")
parser$add_argument('--PATH_GenotypeVCFs_reformatted', metavar='FILE', type='character', help="The dir that contains the reformatted Vardict outputs, output of this dir.")
parser$add_argument('--PATH_GenotypeVCFs_annovar_input', metavar='FILE', type='character', help="The dir that contains the reformatted Vardict outputs, output of this dir.")

args <- parser$parse_args()

PATH_GenotypeVCFs_output <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Haplotype_caller/GenotypeGVCFs_table/batch5.tsv"

reformat_gatk_output <- function(PATH_GenotypeVCFs_output, PATH_GenotypeVCFs_reformatted, PATH_GenotypeVCFs_annovar_input){

	df <- read_delim(PATH_GenotypeVCFs_output, delim = "\t")
	df_AD <- df[, grep("AD", names(df))] # slice of df with allele count info
	variant_info <- df[, c(1, 2, 3, 4, 5)]	# info about variants: chrom, pos etc
	df <- as.data.frame(cbind(variant_info, df_AD)) # put back together 
	# names(df) <- gsub("\\.GT|\\.AD|\\.DP|\\.GQ|\\.PL", "", names(df)) # remove the ".AD" at the end of the column names

	df <- df %>%
			pivot_longer(cols = ends_with(".AD"), # put into long format, easier to parse and understand.
						names_to = "sample_names",
						values_to = "read_depth") %>%
			mutate(read_depth = sub(",", "_", read_depth), 
					ALT = gsub("\\*", "DEL", ALT), 
					sample_names = gsub("\\.AD", "", sample_names)) %>%
			separate(col = read_depth, 
					into = c("Ref_reads", "Alt_reads"), 
					sep = "_") %>%
			separate_rows(ALT, Alt_reads) %>%
			filter(Alt_reads != 0, Ref_reads != 0) %>%
			mutate(Alt_reads = as.numeric(Alt_reads), 
				   Ref_reads = as.numeric(Ref_reads), 
				   VAF = Alt_reads/(Alt_reads + Ref_reads))

	# Now generate the annovar input
	anno_df <- df %>% 
			mutate(END_POS = POS + nchar(REF)) %>%
			select(CHROM, POS, END_POS, REF, ALT)
	
	write_csv(df, PATH_GenotypeVCFs_reformatted)
	write.table(anno_df, PATH_GenotypeVCFs_annovar_input, sep = "\t",  col.names = FALSE, row.names = FALSE)
			
}

reformat_gatk_output(PATH_GenotypeVCFs_output = args$PATH_GenotypeVCFs_output, 
					PATH_GenotypeVCFs_reformatted = args$PATH_GenotypeVCFs_reformatted, 
					PATH_GenotypeVCFs_annovar_input = args$PATH_GenotypeVCFs_annovar_input)