library(tidyverse)
library(gridExtra)
library(cowplot)
library(scales)

cool_theme <-

  theme(panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", size = 1),
      axis.ticks = element_line(colour = "black", size = 2),
      axis.text = element_text(size=10),
      axis.text.x = element_text(vjust=0.5, colour = "black", size=12),
      axis.text.y = element_text(vjust=0.5, colour = "black", size=12),
      axis.title = element_text(size=14, face="bold"),
      legend.title = element_text(color = "black", size = 12),
      legend.text = element_text(color = "black", size = 12),
      axis.ticks.length=unit(0.15, "cm"),
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 20)),
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

# Stacked bar plot showing gene function and number of different variants
MAKE_variant_type_plot <- function(df_main) {

	df <- df_main %>% select(Variant, Function)
	df <- df[!df$Function == ".", ]

	p <- ggplot(data = df, aes(x = Variant, fill = Function)) + 
		geom_bar(stat = "count") + 
		scale_y_continuous(breaks = pretty_breaks(10), name = "Count") + 
		scale_fill_brewer(palette = "Set1") +
		cool_theme + 
		ggtitle("Number of variants and functions") 

	return(p)

}

# counts for each gene, only considering exonic variants
MAKE_gene_counts_plot <- function(df_main) {

	df <- df_main %>% 
			select(Function, Gene, Variant) %>%
			filter(Function == "exonic")

	df <- as.data.frame(table(df)) # we need the counts information to sort accordingly
	df <- df[order(df$Freq), ]
	df$Gene <- factor(df$Gene, levels = unique(df$Gene))

	p <- ggplot(data = df, aes(x = Gene, y = Freq, fill = Variant)) + 
		geom_bar(stat = "identity") + 
		scale_y_continuous(breaks = pretty_breaks(10), name = "Count") + 
		cool_theme +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
		ggtitle("Distribution of exonic variants")

	y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y.major_source

	p <- p + geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "black")

	return(p)

}

# make vaf plot, patient ID on y, color variants by function
MAKE_vaf_plot_function <- function(df_main) {
	
	df <- select(df_main, Sample_name, Patient_ID, VAF, Function)
	df <- df[!df$Function == ".", ]

	# x <- str_split(df$Sample_name, "-") # split the sample name
	# df$patient_id <- paste(lapply(x, "[", 2), lapply(x, "[", 3), sep = "-") # paste 2nd and 3rd elements to create the sample name 
	df <- df[order(df$VAF), ]
	df$Patient_ID <- factor(df$Patient_ID, levels = unique(df$Patient_ID))

	p <- ggplot(data = df, aes(x = VAF, y = Patient_ID)) + 
			geom_point(aes(color = Function), size = 0.00001) +
			scale_x_continuous(breaks = pretty_breaks(15), name = "Variant Allele Frequency") + 
			scale_color_brewer(palette = "Set1") +
			ylab("Patient ID") + 
			cool_theme + 
			ggtitle("Showing variants per patient, colored according to function")


	y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y.major_source

	p <- p + geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "black") + 
			geom_point(aes(color = Function), size = 2)

	return(p)

}

# make vaf plot, patient ID on y, color variants by variant type
MAKE_vaf_plot_variant <- function(df_main) {
	
	df <- select(df_main, Sample_name, Patient_ID, VAF, Function, Variant)

	# x <- str_split(df$Sample_name, "-") # split the sample name
	# df$patient_id <- paste(lapply(x, "[", 2), lapply(x, "[", 3), sep = "-") # paste 2nd and 3rd elements to create the sample name 
	df <- df[order(df$VAF), ]
	df$Patient_ID <- factor(df$Patient_ID, levels = unique(df$Patient_ID))

	p <- ggplot(data = df, aes(x = VAF, y = Patient_ID)) + 
			geom_point(aes(color = Variant), size = 0.00001) +
			scale_x_continuous(breaks = pretty_breaks(15), name = "Variant Allele Frequency") + 
			ylab("Patient ID") + 
			cool_theme + 
			ggtitle("Showing variants per patient, colored according to variant type")

	y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y.major_source

	p <- p + geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "black") + 
			geom_point(aes(color = Variant), size = 2)

	return(p)

}

# shows the averaged depth from each sample in a bar plot, ordered from low to high. 
# vaf_plot_variant is needed because we present the patients in the same order
MAKE_average_depth_per_patient <- function(df_average_depth, vaf_plot_variant) {

	df_average_depth <- read.delim(PATH_to_averaged_depth, sep = " ", header = FALSE)
	names(df_average_depth) <- c("bam_file", "average_depth")

	y.labels <- ggplot_build(vaf_plot_variant)$layout$panel_params[[1]]$y.labels

	df <- df_average_depth[order(df_average_depth$average_depth), ]
	x <- str_split(df$bam_file, "-") # split the sample name
	df$Patient_ID <- paste(lapply(x, "[", 1), lapply(x, "[", 2), lapply(x, "[", 3), sep = "-") # paste 2nd and 3rd elements to create the sample name 
	df$Patient_ID <- factor(df$Patient_ID, levels = y.labels)
	df <- df[!is.na(df$Patient_ID), ] # drop NAs

	p <- ggplot(data = df, aes(x = Patient_ID, y = average_depth)) + 
			geom_bar(stat = "identity") + 
			cool_theme +
			scale_y_continuous(breaks = pretty_breaks(15), name = "Average depth") + 
			xlab("Patient ID") + 
			cool_theme + 
			ggtitle("Showing average depth across all positions per patient") + 
			coord_flip()

	x.ticks <- ggplot_build(p)$layout$panel_params[[1]]$x.major_source

	p <- p + geom_hline(yintercept = x.ticks, size = 0.5, linetype = "dotted", color = "black")
	return(p)

}

# patients that didnt have any variants will have a NA, we will make a separate plot for them. 
MAKE_average_depth_per_patient2 <- function(df_average_depth, vaf_plot_variant) {

	df_average_depth <- read.delim(PATH_to_averaged_depth, sep = " ", header = FALSE)
	names(df_average_depth) <- c("bam_file", "average_depth")

	y.labels <- ggplot_build(vaf_plot_variant)$layout$panel_params[[1]]$y.labels

	df <- df_average_depth[order(df_average_depth$average_depth), ]
	x <- str_split(df$bam_file, "-") # split the sample name
	df$Patient_ID <- paste(lapply(x, "[", 1), lapply(x, "[", 2), lapply(x, "[", 3), sep = "-") # paste 2nd and 3rd elements to create the sample name 
	df$Patient_ID <- factor(df$Patient_ID, levels = y.labels)

	# patients that didnt have any variants will have a NA in the Patient_ID column, we will make a separate plot for them. 
	df_no_var <- df[is.na(df$Patient_ID), ] # identify the NA patients
	x <- str_split(df_no_var$bam_file, "-") # split the sample name
	df_no_var$Patient_ID <- paste(lapply(x, "[", 1), lapply(x, "[", 2), lapply(x, "[", 3), sep = "-") # paste 2nd and 3rd elements to create the sample name 
	df_no_var <- df_no_var[order(df_no_var$average_depth), ]
	df_no_var$Patient_ID <- factor(df_no_var$Patient_ID, levels = unique(df_no_var$Patient_ID))

	p <- ggplot(data = df_no_var, aes(x = Patient_ID, y = average_depth)) + 
			geom_bar(stat = "identity") + 
			cool_theme +
			scale_y_continuous(breaks = pretty_breaks(15), name = "Average depth") + 
			xlab("Patient ID") + 
			cool_theme + 
			ggtitle("Showing average depth across all positions for patients without any detected variants") + 
			coord_flip()

	x.ticks <- ggplot_build(p)$layout$panel_params[[1]]$x.major_source

	p <- p + geom_hline(yintercept = x.ticks, size = 0.5, linetype = "dotted", color = "black")
	return(p)

}

# show the depth at each snv. for indels, the position on the left, just before the indel starts, is given.
# we use vaf_plot_variant to make sure the order of samples here is the same as the previous plot
MAKE_depth_at_variants <- function(df_main, vaf_plot_variant){

	df <- df_main

	# adding suffixes to patient ids with more than one mutation 
	repeated <- unname(table(df$Patient_ID))
	suffixes <- unlist(mapply(seq, rep(1, length(repeated)), repeated))

	# some patients only appear once, remove the suffix from them
	df$Mutations <- paste(df_main$Patient_ID, unlist(mapply(seq, rep(1, length(repeated)), repeated)), sep = "_")

	y.labels <- ggplot_build(vaf_plot_variant)$layout$panel_params[[1]]$y.labels
	df$Patient_ID <- factor(df$Patient_ID, levels = y.labels)

	p <- ggplot(data = df, aes(x = Mutations, y = depth)) + 
			geom_bar(stat = "identity") + 
			cool_theme +
			scale_y_continuous(breaks = pretty_breaks(15), name = "Depth at position") + 
			xlab("Patient ID") + 
			ggtitle("Showing depth at the position of detected variants") + 
			coord_flip()

	return(p)

}

SAVE_plot <- function(PATH_figure_main, NAME_p, height_cm, width_cm, p) {

	message(c("Saving to ", file.path(PATH_figure_main, NAME_p)))
	ggsave(file.path(PATH_figure_main, NAME_p), p, width = width_cm, height = height_cm, units = "cm")

}

MAIN <- function(PATH_figure_main, PATH_to_variants, PATH_to_averaged_depth) {

	df_main <- read_csv(PATH_to_variants)

	p1 <- MAKE_variant_type_plot(df_main)
	p2 <- MAKE_gene_counts_plot(df_main)
	p3 <- MAKE_vaf_plot_function(df_main) 
	p4 <- MAKE_vaf_plot_variant(df_main)  
	p5 <- MAKE_average_depth_per_patient(PATH_to_averaged_depth, p4) 
	p6 <- MAKE_average_depth_per_patient2(PATH_to_averaged_depth, p4)
	p7 <- MAKE_depth_at_variants(df_main, p4)

	dir.create(PATH_figure_main)

	SAVE_plot(PATH_figure_main, "variant_type.png", 15, 15, p1)
	SAVE_plot(PATH_figure_main, "gene_counts.png", 15, 20, p2)
	SAVE_plot(PATH_figure_main, "vaf_function.png", 20, 30, p3)
	SAVE_plot(PATH_figure_main, "vaf_variant.png", 20, 30, p4)
	SAVE_plot(PATH_figure_main, "average_depth1.png", 20, 30, p5)
	SAVE_plot(PATH_figure_main, "average_depth2.png", 8, 22, p6)
	SAVE_plot(PATH_figure_main, "depth_at_variants.png", 22, 30, p7)

}

PATH_variants <- "/groups/wyattgrp/users/amunzur/chip_project/variant_lists/chip_variants.csv" # Before subsetting to IGV\
PATH_figure_main <- "/groups/wyattgrp/users/amunzur/chip_project/figures/main_figures/new_chip_panel/before_manual_curation"

PATH_variants_igv <- "/groups/wyattgrp/users/amunzur/chip_project/variant_lists/chip_variants_igv_subsetted.csv" # Before subsetting to IGV\ # After subsetting to IGV
PATH_figure_main_igv <- "/groups/wyattgrp/users/amunzur/chip_project/figures/main_figures/new_chip_panel/after_manual_curation"

PATH_to_averaged_depth <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/new_chip_panel_CP.txt"
DIR_individual_depths <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/individual_samples/new_chip_panel_CP"


MAIN(PATH_figure_main, PATH_variants, PATH_to_averaged_depth)
MAIN(PATH_figure_main_igv, PATH_variants_igv, PATH_to_averaged_depth)

# combined_plots <- plot_grid(p1, p2, p3, p4, p5, p6, p7, labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G'), label_size = 12)