library(tidyverse)

##################################################
# CLEAN UP BETS
##################################################
process_bets <- function(path_to_raw, path_to_cleaned) {
	
	df <- as.data.frame(read.delim(path_to_raw))
	names(df)[9:ncol(df)] <- paste(names(df)[9:ncol(df)], "bam", sep=".") # add the bam extension so it is easier to find

	df <- pivot_longer(data = df, cols = ends_with(".bam"), # put into long format, easier to parse and understand.
							names_to = "sample_names", 
							values_to = "variants")

	df <- df[grep("[*]", df$variants), ] # keep the rows with the asterix only
	df[, 2] <- NULL

	df$sample_names <- gsub('\\.', '-', df$sample_names)
	df$sample_names <- gsub('-bam', '.bam', df$sample_names)

	message(c("There are ", dim(df)[1], " mutations."))
	message(c("Saving to ", path_to_cleaned))
	
	write_delim(df, path_to_cleaned, delim = "\t")
}



path_to_raw <- "/groups/wyattgrp/users/amunzur/pipeline/resources/betastasis/RAW_mutations_no_germline_filter.tsv"
path_to_cleaned <- "/groups/wyattgrp/users/amunzur/pipeline/resources/betastasis/CLEANED_mutations_no_germline_filter.tsv"

process_bets(path_to_raw, path_to_cleaned)