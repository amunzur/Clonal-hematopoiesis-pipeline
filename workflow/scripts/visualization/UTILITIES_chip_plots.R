# library(tidyverse)
# library(gridExtra)
# library(cowplot)
# library(scales)

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
MAKE_variant_type_plot <- function(df_main, colors_vector) {

	df <- df_main %>% select(Variant, Function)
	df <- df[!df$Function == ".", ]

	p <- ggplot(data = df, aes(x = Variant, fill = Function)) + 
		geom_bar(stat = "count") + 
		scale_y_continuous(breaks = pretty_breaks(10), name = "Count") + 
		scale_fill_manual(values = colors_vector) +
		cool_theme + 
		ggtitle("Number of variants and functions") 

	return(p)

}

# Simple bar plot showing the number of chips per patient
MAKE_chips_per_patient <- function(df_main, color_what, colors_vector, plot_title) {

	# Total number of chips per patient, to merge later on
	df_total_counts <- df_main %>%
		select(Patient_ID) %>%
		group_by(Patient_ID) %>%
		summarise(total_counts = n())

	df_variants <- df_main %>% 
		select(Patient_ID, Gene, all_of(color_what)) %>%
		group_by_at(c("Patient_ID", color_what)) %>%
		summarise(variant_counts = n())

	df <- left_join(df_variants, df_total_counts, by = "Patient_ID") %>%
			arrange(total_counts)
	
	df$Patient_ID <- factor(df$Patient_ID, levels = unique(df$Patient_ID))

	p <- ggplot(data = df, aes_string(x = "Patient_ID", y = "variant_counts", fill = color_what)) + 
		geom_bar(stat = "identity") + 
		scale_y_continuous(breaks = pretty_breaks(10), name = "Count") + 
		scale_fill_manual(values = colors_vector) + 
		cool_theme +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
		ggtitle(plot_title)

	y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y.major_source

	p <- p + geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "gray62")

	return(p)

}

# Make a plot to show the number of chips identified by the two variant callers.
make_variant_caller_counts <- function(df_main, colors_vector){

	df <- df_main %>%
		select(variant_caller, Variant) %>%
		group_by(variant_caller, Variant) %>%
		summarize(n = n()) %>%
		group_by(variant_caller)

	p <- ggplot(data = df, aes(x = variant_caller, y = n, fill = Variant)) + 
		geom_bar(position="stack", stat="identity") + 
		scale_y_continuous(breaks = pretty_breaks(10), name = "Count") + 
		scale_fill_manual(values = colors_vector) + 
		xlab("Variant caller") + 
		cool_theme +
		theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5)) +
		ggtitle("Distribution of CHIP variants across variant callers")

	return(p)

}

MAKE_variant_type_plot <- function(df_main, colors_vector) {

	df <- df_main %>% select(Variant, Function)
	df <- df[!df$Function == ".", ]

	p <- ggplot(data = df, aes(x = Variant, fill = Function)) + 
		geom_bar(stat = "count") + 
		scale_y_continuous(breaks = pretty_breaks(10), name = "Count") + 
		scale_fill_manual(values = colors_vector) +
		cool_theme + 
		ggtitle("Number of variants and functions") 

	return(p)

}

# counts for each gene, only considering exonic variants
MAKE_gene_counts_plot <- function(df_main, colors_vector) {

	# Counts how many total number of snvs and indels a gene has in total. Only contains one occurrence of each gene. 
	df_total_counts <- df_main %>% 
			select(Function, Gene, Variant) %>%
			filter(Function == "exonic") %>%
			group_by(Gene) %>%
			mutate(total_count = table(Gene)) %>% 
			distinct(Gene, .keep_all = TRUE) %>%
			select(Gene, total_count)

	df_variants <- df_main %>% 
			select(Function, Gene, Variant) %>%
			filter(Function == "exonic")

	df_variants <- as.data.frame(table(df_variants)) # we need the counts information to sort accordingly
	df <- left_join(df_variants, df_total_counts, by = c("Gene"))
	df <- df[order(df$total_count), ]
	df$Gene <- factor(df$Gene, levels = unique(df$Gene))

	p <- ggplot(data = df, aes(x = Gene, y = Freq, fill = Variant)) + 
		geom_bar(stat = "identity") + 
		scale_y_continuous(breaks = pretty_breaks(10), name = "Count") + 
		scale_fill_manual(values = colors_vector) +
		cool_theme +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
		ggtitle("Distribution of exonic variants")

	y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y.major_source

	p <- p + geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "gray62")

	return(p)

}

# make vaf plot, patient ID on y, color variants by function
MAKE_vaf_plot_function <- function(df_main, colors_vector) {
	
	df <- select(df_main, Sample_name, Patient_ID, VAF, Function, Depth)
	df <- df[!df$Function == ".", ]

	# x <- str_split(df$Sample_name, "-") # split the sample name
	# df$patient_id <- paste(lapply(x, "[", 2), lapply(x, "[", 3), sep = "-") # paste 2nd and 3rd elements to create the sample name 
	df <- df[order(df$VAF), ]
	df$Patient_ID <- factor(df$Patient_ID, levels = unique(df$Patient_ID))

	p <- ggplot(data = df, aes(x = VAF, y = Patient_ID, label = Depth)) + 
			geom_point(aes(color = Function), size = 0.00001) +
			scale_x_continuous(breaks = pretty_breaks(15), name = "Variant Allele Frequency") + 
			scale_color_manual(values = colors_vector) +
			ylab("Patient ID") + 
			cool_theme + 
			ggtitle("Showing variants per patient, colored according to function")


	y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y.major_source

	p <- p + geom_hline(yintercept = y.ticks, size = 0.2, linetype = "dotted", color = "gray") + 
			geom_point(aes(color = Function), size = 1)

	return(p)

}

# make vaf plot, patient ID on y, color variants by variant type
MAKE_vaf_plot_variant <- function(df_main, colors_vector) {
	
	df <- select(df_main, Sample_name, Patient_ID, VAF, Function, Variant, Depth)

	df <- df[order(df$VAF), ]
	df$Patient_ID <- factor(df$Patient_ID, levels = unique(df$Patient_ID))

	p <- ggplot(data = df, aes(x = VAF, y = Patient_ID, label = Depth)) + 
			geom_point(aes(color = Variant), size = 0.00001) +
			scale_x_continuous(breaks = pretty_breaks(15), name = "Variant Allele Frequency") + 
			scale_color_manual(values = colors_vector) +
			ylab("Patient ID") + 
			cool_theme + 
			ggtitle("Showing variants per patient, colored according to variant type")

	y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y.major_source

	p <- p + geom_hline(yintercept = y.ticks, size = 0.2, linetype = "dotted", color = "gray62") + 
			geom_point(aes(color = Variant), size = 1)

	return(p)

}

# shows the averaged depth from each sample in a bar plot, ordered from low to high. 
# vaf_plot_variant is needed because we present the patients in the same order
MAKE_depth_plots <- function(df_main, vaf_plot_variant) {

	df <- df_main %>%
		select(Patient_ID, Averaged_depth, Median_depth) %>%
		distinct(Patient_ID, Averaged_depth, Median_depth, .keep_all = TRUE) %>%
		gather(Depth, Values, "Averaged_depth", "Median_depth")


	y.labels <- ggplot_build(vaf_plot_variant)$layout$panel_params[[1]]$y.labels
	df$Patient_ID <- factor(df$Patient_ID, levels = y.labels)

	df <- df[!is.na(df$Patient_ID), ] # drop NAs

	p <- ggplot(data = df, aes(x = Patient_ID, y = Values, fill = Depth)) + 
			geom_bar(position="dodge", stat="identity") + 
			cool_theme +
			scale_y_continuous(breaks = pretty_breaks(15), name = "Depth") + 
			xlab("Patient ID") + 
			cool_theme + 
			ggtitle("Median and mean depth of the samples") + 
			coord_flip()

	x.ticks <- ggplot_build(p)$layout$panel_params[[1]]$x.major_source

	p <- p + geom_hline(yintercept = x.ticks, size = 0.5, linetype = "dotted", color = "gray62")

	return(p)

}


SAVE_plot <- function(PATH_figure_main, NAME_p, height_cm, width_cm, p) {

	message(c("Saving to ", file.path(PATH_figure_main, NAME_p)))
	ggsave(file.path(PATH_figure_main, NAME_p), p, width = width_cm, height = height_cm, units = "cm")

}

MAIN <- function(PATH_figure_main, PATH_to_variants) {

	# Helps make sure to maintain the same order in each plot
	df_main <- read_csv(PATH_to_variants)
	df_main$Function <- factor(df_main$Function, levels = sort(unique(df_main$Function), decreasing = TRUE))
	df_main$Variant <- factor(df_main$Variant, levels = sort(unique(df_main$Variant), decreasing = FALSE))

	p1 <- MAKE_variant_type_plot(df_main)
	p2 <- MAKE_gene_counts_plot(df_main)
	p3 <- MAKE_vaf_plot_function(df_main) 
	p4 <- MAKE_vaf_plot_variant(df_main)  
	# p5_positional_depth <- MAKE_average_depth_per_patient(PATH_to_averaged_depth, "Depth", "Depth at variants", "Depth at variant positions", p4) 
	p6_mean_depth <- MAKE_depth_plots(df_main, "Averaged_depth", "Averaged sample depth", "Averaged sample depth in patients with at least one variant", p4) 
	p7_median_depth <- MAKE_depth_plots(df_main, "Median_depth", "Median sample depth", "Median sample depth in patients with at least one variant", p4)
	p7 <- MAKE_chips_per_patient(df_main, "Variant", variant_colors, "Distribution of CHIP variants per patient, colored according to variant type.") # barplot showing chips per patient, colored according to snv, and indel
	p8 <- MAKE_chips_per_patient(df_main, "Function", function_colors, "Distribution of CHIP variants per patient, colored according to variant functions.") # barplot showing chips per patient, colored according to function (exonic, splicing etc)
	p9 <- make_variant_caller_counts(df_main)

	dir.create(PATH_figure_main)

	SAVE_plot(PATH_figure_main, "variant_type.png", 16, 16, p1)
	SAVE_plot(PATH_figure_main, "gene_counts.png", 16, 30, p2)
	SAVE_plot(PATH_figure_main, "vaf_function.png", 20, 30, p3)
	SAVE_plot(PATH_figure_main, "vaf_variant.png", 20, 30, p4)

	SAVE_plot(PATH_figure_main, "average_depths.png", 20, 30, p6_mean_depth)
	SAVE_plot(PATH_figure_main, "median_depth.png", 20, 30, p7_median_depth)
	SAVE_plot(PATH_figure_main, "chips_per_patient_Variant.png", 22, 30, p7)
	SAVE_plot(PATH_figure_main, "chips_per_patient_Function.png", 22, 30, p8)
	SAVE_plot(PATH_figure_main, "variant_caller_counts.png", 15, 15, p9)
}

PATH_to_variants <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined/combined.csv" # Before subsetting to IGV
PATH_figure_main <- "/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/new_chip_panel/before_manual_curation"

# PATH_variants_igv <- "/groups/wyattgrp/users/amunzur/chip_project/variant_lists/chip_variants_igv_subsetted.csv" # Before subsetting to IGV\ # After subsetting to IGV
# PATH_figure_main_igv <- "/groups/wyattgrp/users/amunzur/chip_project/figures/main_figures/new_chip_panel/after_manual_curation"

# PATH_to_averaged_depth <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/new_chip_panel_CP.txt"
# DIR_individual_depths <- "/groups/wyattgrp/users/amunzur/chip_project/metrics/coverage_information/individual_samples/new_chip_panel_CP"

# MAIN(PATH_figure_main, PATH_to_variants)
# MAIN(PATH_figure_main_igv, PATH_variants_igv)

# combined_plots <- plot_grid(p1, p2, p3, p4, p5, p6, p7, labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G'), label_size = 12)