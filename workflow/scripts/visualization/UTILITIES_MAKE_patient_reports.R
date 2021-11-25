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
MAKE_chips_per_patient <- function(df_main, colors_vector, plot_title) {

	# Total number of chips per patient, to merge later on
	df <- df_main %>%
				select(Chrom, VAF, Variant, Function, Depth) %>%
				mutate(VAF_percent = (as.numeric(VAF))*100, 
						y.labels = percent(as.numeric(VAF)))

	p <- ggplot(data = df, aes(x = Function, y = VAF_percent, color = Variant, label = Depth)) + 
		geom_point(size = 3) + 
		scale_y_continuous(breaks = pretty_breaks(10), name = "VAF", labels = function(x) paste0(x, "%")) + 
		scale_fill_manual(values = colors_vector) + 
		cool_theme +
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
		scale_y_continuous(breaks = pretty_breaks(5), name = "Count") + 
		scale_fill_manual(values = colors_vector) +
		cool_theme +
		ggtitle("Distribution of exonic variants")

	y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y.major_source

	p <- p + geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "gray62")

	return(p)

}

