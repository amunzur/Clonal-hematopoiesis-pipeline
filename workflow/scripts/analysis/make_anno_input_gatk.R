library(tidyverse)
library(argparse)

parser <- ArgumentParser(description = "Make the input files for running ANNOVAR on files from GATK.")

# passed from Snakemake into the script
parser$add_argument('--PATH_vcf', metavar='FILE', type='character', help="PATH to the vcf file from Haplotype caller in zipped format")
parser$add_argument('--PATH_vcf_table', metavar='FILE', type='character', help="Haplotype caller output in table format")
parser$add_argument('--PATH_ANNOVAR_input', metavar='FILE', type='character', help="The path that contains the inputs for running ANNOVAR")

args <- parser$parse_args()

# reorganize such that in each row there is only one variant reported
parse_locus_info <- function(PATH_vcf){
	
	# reading the vcf file twice, once for col names and once for the data
	tmp_vcf <- readLines(PATH_vcf)
	tmp_vcf_data <- read.table(PATH_vcf, stringsAsFactors = FALSE)

	# filter for the columns names
	tmp_vcf <- tmp_vcf[-(grep("#CHROM", tmp_vcf) +1) :- (length(tmp_vcf))]
	vcf_names <- unlist(strsplit(tmp_vcf[length(tmp_vcf)], "\t"))
	names(tmp_vcf_data) <- vcf_names

	# filter for the cols we are interested in keeping
	vars <- tmp_vcf_data %>%
		rename(CHROM = "#CHROM") %>%
			select(CHROM, POS, REF, ALT)
			
	# find which rows have more than one variant, do separate_rows on them
	idx <- grep(",", vars$ALT)
	x <- vars[idx, ]

	y <- mapply(rep, times = str_count(x$ALT, ",") + 1, x = as.vector(x$REF)) # repeat the ref as many times as the alt is repeated in a given row
	x$REF <- unlist(lapply(y, paste, collapse = ",")) # collapse in the same format as the alt
	
	x <- x %>%
		mutate(ALT = gsub("\\*", 0, ALT)) %>%
		separate_rows(REF, ALT, convert = TRUE)

	# delete the the idx rows from the main df, instead append the x to the dataframe
	z <- rbind(vars[-idx, ], x)
	z$POS <- as.numeric(z$POS)
	return(z)
}

# reformat such that each variant is associated with a sample, vaf etc
parse_variant_info <- function(PATH_vcf_table){
	df <- read_delim(PATH_vcf_table, delim = "\t")
	df_AD <- df[, grep("AD", names(df))] # slice of df with allele count info
	variant_info <- df[, c(1, 2, 3, 4, 5)]	# info about variants: chrom, pos etc
	df <- as.data.frame(cbind(variant_info, df_AD)) # put back together 
	# names(df) <- gsub("\\.GT|\\.AD|\\.DP|\\.GQ|\\.PL", "", names(df)) # remove the ".AD" at the end of the column names

	df <- df %>%
			gather("sample_names", "read_depth", ends_with(".AD")) %>%
			mutate(read_depth = sub(",", "_", read_depth), ALT = gsub("\\*", 0, ALT), sample_names = gsub("\\.AD", "", sample_names)) %>%
			separate(col = read_depth, into = c("Ref_reads", "Alt_reads"), sep = "_") %>%
			separate_rows(ALT, Alt_reads) %>%
			mutate(Alt_reads = as.numeric(Alt_reads), Ref_reads = as.numeric(Ref_reads), VAF = Alt_reads/(Alt_reads + Ref_reads))

	return(df)
}

# Parse the locus info to make ANNOVAR inputs 
make_anno_input <- function(locus_info, ANNOVAR_input) {

		locus_info$POS <- as.numeric(locus_info$POS)
		df <- locus_info %>%
			mutate(
				nchar_REF = nchar(REF), 
				nchar_ALT = nchar(ALT), 
				nchar_dif = nchar_REF - nchar_ALT, 
				TYPE = case_when(
					nchar_dif == 0 ~ "SNP", 
					nchar_dif < 0 ~ "INS",
					nchar_dif > 0 ~ "DEL"),
				POS_new = case_when(
					nchar_REF == 1 ~ POS, # SNP 
					nchar_REF > 1 ~ POS + 1, # more than 1 base at REF, add 1 to pos
					TRUE ~ 9999), 
				END_new = case_when(
					nchar_REF == 1 ~ POS, # SNP
					nchar_REF > 1 ~ POS + nchar_REF - 1, # more than 
					TRUE ~ 9999), 
				REF = case_when(
					TYPE == "SNP" ~ REF, 
					TYPE != "SNP" ~ substring(REF, 2)),
				ALT = case_when(
					TYPE == "SNP" ~ ALT, 
					TYPE != "SNP" ~ substring(ALT, 2)), 
				REF = case_when(
					nchar(REF) == 0 ~ "-",
					nchar(REF) != 0 ~ REF),
				ALT = case_when(
					nchar(ALT) == 0 ~ "-",
					nchar(ALT) != 0 ~ ALT)) %>%
			select(CHROM, POS_new, END_new, REF, ALT)

	write.table(df, ANNOVAR_input, sep="\t",  col.names = FALSE, row.names = FALSE, quote = FALSE)
}

MAIN <- function(PATH_ANNOVAR_input, PATH_vcf, PATH_vcf_table){

	locus_info <- parse_locus_info(PATH_vcf) # each row is a variant 
	sample_info <- parse_variant_info(PATH_vcf_table) # each row is a sample

	locus_info <- inner_join(locus_info, sample_info[, c("CHROM", "POS", "REF", "ALT", "TYPE")]) # add the type col to locus info
	locus_info <- locus_info[!duplicated(locus_info), ] # remove duplicated rows

	make_anno_input(locus_info, PATH_ANNOVAR_input)
}

# INPUTS
# cohort_name <- "batch4"
# PATH_vcf <- file.path("/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Haplotype_caller/GenotypeGVCFs/test") # unzipped results
# PATH_vcf_table <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Haplotype_caller/GenotypeGVCFs_table/batch4_OLD.tsv" # HaplotypeCaller results in table format
# DIR_ANNOVAR_input <- "/groups/wyattgrp/users/amunzur/pipeline/results/data/annovar_inputs/haplotype_caller/"

MAIN(
	PATH_vcf = args$PATH_vcf, 
	PATH_vcf_table = args$PATH_vcf_table, 
	PATH_ANNOVAR_input = args$PATH_ANNOVAR_input){
