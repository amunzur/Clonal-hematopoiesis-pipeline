library(tidyverse)
library(scales)
library(RColorBrewer)
library(gridExtra)

PATH_file <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined/curated.csv"
DIR_figure <- "/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/figures"

# Only keep unique variants and remove synonymous mutations
DF <- as.data.frame(read_csv(PATH_file)) %>%
    distinct(AAchange, .keep_all = TRUE) %>%
    filter(!grepl("Synonymous", Effects, ignore.case = TRUE))

# Find the longest splicing variant
protein_annots <- str_split(DF$Protein_annotation, ":")

idx_list <- list()
for (substr in protein_annots) {
    idx_list <- append(idx_list, which.max(as.numeric(str_match(substr, "\\d+")))) # the idx of the annotations we keep for each unique variant, the largest splice variant
}
DF$Protein_annotation <- mapply("[", protein_annots, idx_list) # choose the splice variant with the highest value

# Plotting utilities
cool_theme <-

  theme(panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", size = 0.5),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(vjust=0.5, colour = "black", size=8),
      axis.text.y = element_text(vjust=0.5, colour = "black", size=8),
      axis.title = element_text(size = 10, face="bold"),
      legend.title = element_text(color = "black", size = 12),
      legend.text = element_text(color = "black", size = 12),
      axis.ticks.length=unit(0.15, "cm"),
      axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 20)),
      axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

# Gene counts
gene_counts <- function(DF, DIR_figure, fig_height, fig_width) {

    df_total_counts <- DF %>% 
        select(
            Gene, 
            Variant, 
            AAchange) %>%
	    group_by(Gene) %>%
	    mutate(total_count = table(Gene)) %>% 
	    select(
            Gene, 
            total_count) %>%
        distinct(
            Gene, 
            .keep_all = TRUE)

    df_variants <- as.data.frame(table(select(DF, Gene, Variant)))
    df <- left_join(df_variants, df_total_counts, by = c("Gene"))

    df <- df[order(df$total_count), ] # order the genes based on the total number of indels and SNVs it has
    df$Gene <- factor(df$Gene, levels = unique(df$Gene))

    p <- ggplot(data = df, aes(fill = Variant, x = Gene, y = Freq)) + 
         geom_bar(position="stack", stat="identity") +
         scale_y_continuous(breaks = pretty_breaks(10), name = "Count") + 
         cool_theme +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
         ylab("Number of variants")

    y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y$minor_breaks
    p <- p + geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "gray62")

    ggsave(file.path(DIR_figure, "gene_counts.jpg"), p, width = fig_width, height = fig_height, units = "cm")

}

tumor_wbc_scatter <- function(DF, DIR_figure, variant_type, fig_width, fig_height){

    df_main <- DF %>%
        select(VAF_t, VAF_n, Variant) %>%
        mutate(VAF_t = VAF_t*100, 
               VAF_n = VAF_n*100) %>%
        rename(Tumor = VAF_t,
               WBC = VAF_n) %>%
        filter(Variant == variant_type)

    p <- ggplot(data = df_main, aes(x = WBC, y = Tumor)) + 
         geom_point(size = 0.5) +
         scale_x_continuous(breaks = pretty_breaks(10)) + 
         scale_y_continuous(breaks = pretty_breaks(10)) + 
         geom_smooth(method = "lm", size = 0.5) +
         xlab("WBC VAF (%)") + 
         ylab("Tumor VAF (%)") +
         cool_theme +
         ggtitle(str_to_title(variant_type))

    ggsave(file.path(DIR_figure, paste(str_to_title(variant_type), "VAF.jpg", sep = "_")), p, width = fig_width, height = fig_height, units = "cm")

    return(p)

}

###################################################
# Tumor - WBC dot plot. Separated into insertion, deletion and SNV
###################################################

tumor_wbc_dotplot <- function(DF, variant_type, tumor_threshold_lower, tumor_threshold_upper, DIR_figure, fig_height, fig_width){

    df_main <- DF %>%
        distinct(
            Protein_annotation, 
            .keep_all = TRUE) %>% # only unique variants
        select(
            Protein_annotation, 
            VAF_t, 
            VAF_n, 
            Variant, 
            Gene) %>%
        mutate(
            VAF_t = VAF_t*100, 
            VAF_n = VAF_n*100, 
            Annotation = paste(Gene, Protein_annotation)) %>%
        select(
            -Protein_annotation, 
            -Gene) %>%
        rename(
            Tumor = VAF_t,
            WBC = VAF_n) %>%
        filter(
            Variant == variant_type, 
            Tumor >= tumor_threshold_lower & Tumor < tumor_threshold_upper) %>%
        arrange(Tumor) %>%
        gather(
            Sample_type, 
            VAF, 
            1:2) %>%
        mutate(Annotation = factor(Annotation, levels = unique(Annotation)))

    p <- ggplot(data = df_main, aes(color = Sample_type, y = Annotation, x = VAF)) + 
             geom_point() +
             scale_x_continuous(breaks = pretty_breaks(15)) + 
             xlab("VAF (%)") +
             geom_smooth(method = "lm") +
             cool_theme + 
             geom_hline(yintercept = 1:(dim(df_main)[1]/2), size = 0.3, linetype = "dotted", color = "gray62") + 
             ggtitle(paste("CHIP mutations -", str_to_title(variant_type))) + 
             theme(legend.position = "none")
    
    plot_name <- paste0(paste(variant_type, tumor_threshold_lower, tumor_threshold_upper, sep = "_"), ".jpg")
    ggsave(file.path(DIR_figure, plot_name), p, width = fig_width, height = fig_height, units = "cm")

}

################################################
# LOLLIPOP PLOTS TO COMPARE TUMOR AND WBC VAF AT THE SAME TIME, SEPARATED BY VARIANT TYPE
################################################
tumor_wbc_dotplot(DF = DF, variant_type = "deletion", tumor_threshold_lower = 1.5, tumor_threshold_upper = 100, DIR_figure = DIR_figure, fig_height = 7, fig_width = 15)
tumor_wbc_dotplot(DF = DF, variant_type = "deletion", tumor_threshold_lower = 0, tumor_threshold_upper = 1.5, DIR_figure = DIR_figure, fig_height = 20, fig_width = 15)
tumor_wbc_dotplot(DF = DF, variant_type = "insertion", tumor_threshold_lower = 0, tumor_threshold_upper = 100, DIR_figure = DIR_figure, fig_height = 5, fig_width = 15)
tumor_wbc_dotplot(DF = DF, variant_type = "snv", tumor_threshold_lower = 0, tumor_threshold_upper = 2, DIR_figure = DIR_figure, fig_height = 20, fig_width = 15)
tumor_wbc_dotplot(DF = DF, variant_type = "snv", tumor_threshold_lower = 2, tumor_threshold_upper = 100, DIR_figure = DIR_figure, fig_height = 20, fig_width = 15)

################################################
# SCATTER PLOTS TO COMPARE TUMOR AND WBC VAFs, SEPARATED BY VARIANT TYPE
################################################
p1 <- tumor_wbc_scatter(DF = DF, DIR_figure = DIR_figure, variant_type = "deletion", fig_height = 15, fig_width = 12)
p2 <- tumor_wbc_scatter(DF = DF, DIR_figure = DIR_figure, variant_type = "insertion", fig_height = 15, fig_width = 12)
p3 <- tumor_wbc_scatter(DF = DF, DIR_figure = DIR_figure, variant_type = "snv", fig_height = 15, fig_width = 12)

combined_VAF <- grid.arrange(p1, p2, p3, nrow = 1)
ggsave(file.path(DIR_figure, "combined_VAF"), combined_VAF, width = 36, height = 15, units = "cm")


################################################
# GENE COUNTS
################################################
gene_counts(DF = DF, DIR_figure = DIR_figure, fig_height = 15, fig_width = 20)