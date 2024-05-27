library(tidyverse)
library(scales)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(stringi)


cool_theme <-

  theme(panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", size = 0.5),
      axis.ticks = element_line(colour = "black", size = 0.5),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(vjust=0.5, colour = "black", size = 10),
      axis.text.y = element_text(vjust=0.5, colour = "black", size = 10),
      axis.title = element_text(size = 15, face="bold"),
      legend.title = element_text(color = "black", size = 12),
      legend.text = element_text(color = "black", size = 12),
      axis.ticks.length=unit(0.15, "cm"),
      axis.title.x = element_text(size = 15, margin = margin(t = 20, r = 0, b = 0, l = 20)),
      axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 20, b = 0, l = 0)))

# PATH_variants <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined/combined_vars_ALL.csv"
PATH_variants <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/Vardict/finalized/SSCS1/chip_variants.csv"
DEPTH_info <- "/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/combined/combined_vars_ALL_DEPTH.csv"

DIR_figure <- "/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/consensus_figures"
depth_df <- read_csv(DEPTH_info) %>% select(Sample_name, Chrom, Protein_annotation, Ref_reads_n, Alt_reads_n, Ref_reads_t, Alt_reads_t)
DF <- as.data.frame(read_csv(PATH_variants)) %>% 
            filter(
                !grepl('VIP', Sample_name_t), 
                Effects != "synonymous") %>%
            group_by(Protein_annotation) %>%
            mutate(Duplicated_variant = n(), Patient_ID = gsub("_gDNA.*", "", Sample_name_n)) %>%
            ungroup() %>%
            mutate(
                Effects = str_replace_all(Effects, c("nonframeshift_deletion" = "Nonframeshift deletion", "nonframeshift_delins" = "Nonframeshift insertion", "stop_gain" = "Stop gain", "missense" = "Missense", "synonymous" = "Synonymous", "start_loss" = "Start loss", "NA" = "NA", "frameshift" = "Frameshift", "splicing" = "Splicing")), 
                Function = str_replace_all(Function, c("exonic" = "Exonic", "ncRNA_intronic" = "ncRNA intronic", "splicing" = "Splicing", "upstream" = "Upstream")), 
                Variant = str_replace_all(Variant, c("deletion" = "Deletion", "insertion" = "Insertion", "snv" = "SNV"))) %>%
            filter(!Protein_annotation %in% c("p.P459R", "p.P459delinsPDAPADPDSGAAR", "p.P471R", "p.D448A", "p.E472D", "p.R383Q", "p.V1956G", "p.G643fs", "p.G394fs", "p.752_753del", "p.186_187del", "p.186_187del", "p.P1326delinsPQ"))

DF <- left_join(DF, depth_df, by = c("Sample_name", "Chrom", "Protein_annotation", "Alt_reads_n", "Alt_reads_t"))

# Gene counts
gene_counts <- function(DF, plot_fill, DIR_figure, fig_height, fig_width, color_palette, y_title, file_name) {   
    
    df_gene_counts <- DF %>% 
        select(Patient_ID, Gene, !!plot_fill) %>%
        drop_na(!!plot_fill) %>%
        arrange(Gene) %>%
	    group_by(Gene) %>%
        select(-Patient_ID) %>%
        add_count(Gene) %>%
        arrange(n)

    df_gene_counts$Gene <- factor(df_gene_counts$Gene, levels = unique(df_gene_counts$Gene))

    p <- ggplot(data = df_gene_counts, aes(x = Gene, fill = get(plot_fill))) + 
         geom_bar(position="stack") +
         scale_y_continuous(breaks = pretty_breaks(8), name = y_title) + 
         cool_theme +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
         scale_fill_brewer(name = plot_fill, palette = color_palette, direction=1)

    y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y$breaks

    p <- p + geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "gray62")

    ggsave(file.path(DIR_figure, file_name), p, width = fig_width, height = fig_height, units = "cm")

}

lm_eqn <- function(df){
    m <- lm(WBC ~ Tumor, df)
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

tumor_wbc_scatter <- function(DF, DIR_figure, fig_width, fig_height, filename, title, xmax, ymax, xpos_equation, ypos_equation){
    
    df_main <- DF %>%
        filter(Function == "Exonic") %>%
        select(Gene, VAF_t, VAF_n, Variant, Protein_annotation) %>%
        mutate(VAF_t = as.numeric(VAF_t)*100, VAF_n = as.numeric(VAF_n)*100) %>%
        rename(Tumor = VAF_t, WBC = VAF_n) %>%
        unite(Annotation, c("Gene", "Protein_annotation"), sep = " ")

    df_main_filtered <- filter(df_main, Tumor <= xmax, WBC <= ymax)

    p <- ggplot(data = df_main, aes(x = WBC, y = Tumor, label = Annotation)) + 
         geom_point(size = 1, color = "red") + 
         geom_text_repel(max.overlaps = 20) +
         scale_x_continuous(breaks = pretty_breaks(10), limits = c(0, xmax)) + 
         scale_y_continuous(breaks = pretty_breaks(10), limits = c(0, ymax)) + 
         geom_smooth(method = 'lm', formula = y ~ x, color = "black") +
         geom_text(x = xpos_equation, y = ypos_equation, label = lm_eqn(df_main_filtered), parse = TRUE, color = "brown1") +
         xlab("WBC VAF (%)") + 
         ylab("cfDNA VAF (%)") +
         cool_theme + 
         theme(strip.background = element_blank(), 
               strip.text.x = element_text(size = 15)) + 
         labs(title = title, 
            subtitle = "Only showing non-synonymous exonic variants")

    y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y$breaks
    x.ticks <- ggplot_build(p)$layout$panel_params[[1]]$x$breaks

    p <- p + 
          geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "gray62") + 
          geom_vline(xintercept = x.ticks, size = 0.5, linetype = "dotted", color = "gray62") 

    ggsave(file.path(DIR_figure, filename), p, width = fig_width, height = fig_height, units = "cm")

}

tumor_wbc_dotplot <- function(DF, DIR_figure, fig_height, fig_width){

    df_main <- DF %>%
        filter(Function == "Exonic") %>%
        select(Patient_ID, Protein_annotation, VAF_t, VAF_n, Variant, Gene) %>%
        add_count(Gene) %>%
        arrange(n) %>%
        unite(Annotation, c("Patient_ID", "Gene", "Protein_annotation"), sep = " ") %>%
        mutate(Annotation = factor(Annotation, levels = unique(Annotation)), VAF_t = VAF_t*100, VAF_n = VAF_n*100) %>%
        rename(cfDNA = VAF_t, WBC = VAF_n) %>%
        arrange(cfDNA) %>%
        gather(Sample_type, VAF, cfDNA, WBC) %>%
        rename(Type = Sample_type)

    p <- ggplot(data = df_main, aes(color = Type, y = Annotation, x = VAF)) + 
             geom_point(size = 2) +
             scale_x_continuous(breaks = pretty_breaks(15)) + 
             xlab("VAF (%)") +
             cool_theme + 
             labs(title = "Protein annotations", 
                  subtitle = "Only showing non-synonymous exonic variants")

    
    y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y$breaks
    p <- p + geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "gray62")

    ggsave(file.path(DIR_figure, "variants_dot_plot.jpg"), p, width = fig_width, height = fig_height, units = "cm")

}

patient_prevalence <- function(DF, DIR_figure, fig_height, fig_width){

    df_main <- DF %>%
        select(Patient_ID, Gene) %>%
        distinct(Patient_ID, Gene, .keep_all = TRUE) %>%
        group_by(Gene) %>%
        summarise(cnt = n()) %>%
        mutate(Prevalence = as.numeric(cnt / sum(cnt))*100) %>%
        arrange(Prevalence) %>%
        mutate(Gene = factor(Gene, levels = unique(Gene)))

    p <- ggplot(data = df_main, aes(x = Gene, y = Prevalence)) + 
             geom_bar(stat = 'identity') +
             scale_y_continuous(breaks = pretty_breaks(10), name = "% of unique patients that have at least one mutation") + 
             cool_theme + 
             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
             ggtitle("CHIP prevalence") + 
             theme(legend.position = "none")             
             

    y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y$minor_breaks
    p <- p + geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "gray62")
    
    ggsave(file.path(DIR_figure, "prevalence.jpg"), p, width = fig_width, height = fig_height, units = "cm")

}

depth <- function(DF, DIR_figure, fig_height, fig_width, y_title, file_name) {
    
    df_main <- DF %>%
        filter(Function == "Exonic") %>%
        select(Gene, Ref_reads_n, Alt_reads_n, Ref_reads_t, Alt_reads_t, Variant, Protein_annotation) %>%
        mutate(depth_n = Ref_reads_n + Alt_reads_n, depth_t = Ref_reads_t + Alt_reads_t)

    p <- ggplot(data = df_main, aes(x = depth_n, y = depth_t, label = Gene)) + 
         geom_point(size = 1, color = "red") + 
         geom_text_repel(max.overlaps = 5) +
         scale_x_continuous(breaks = pretty_breaks(10), limits = c(0, 2500)) + 
         scale_y_continuous(breaks = pretty_breaks(10), limits = c(0, 2500)) + 
         geom_smooth(method = "lm", size = 0.5, col = "black") +
         facet_wrap(~Variant) +
         xlab("WBC depth") + 
         ylab("cfDNA depth") +
         cool_theme + 
         theme(strip.background = element_blank(), 
               strip.text.x = element_text(size = 15)) + 
         labs(title = "cfDNA - WBC depth comparison", 
            subtitle = "Only showing non-synonymous exonic variants")

    y.ticks <- ggplot_build(p)$layout$panel_params[[1]]$y$breaks
    x.ticks <- ggplot_build(p)$layout$panel_params[[1]]$x$breaks

    p <- p + 
          geom_hline(yintercept = y.ticks, size = 0.5, linetype = "dotted", color = "gray62") + 
          geom_vline(xintercept = x.ticks, size = 0.5, linetype = "dotted", color = "gray62") 

    ggsave(file.path(DIR_figure, file_name), p, width = fig_width, height = fig_height, units = "cm")

}

plot_consensus_depth <- function(DF, DIR_figure, fig_height, fig_width, y_title, file_name) {
    
    df_main <- DF %>%
        filter(Function == "Exonic") %>%
        select(Patient_ID, Gene, Protein_annotation) %>%
        arrange(Patient_ID) %>%
        group_by(Patient_ID) %>%
        mutate(n_muts = n()) %>%
        distinct(Patient_ID, .keep_all = TRUE) %>%
        ungroup() %>%
        group_by(n_muts) %>%
        mutate(mutations = n()) %>%
        arrange(n_muts) %>%
        distinct(mutations, .keep_all = TRUE)

    p <- ggplot(data = df_main, aes(x = n_muts, y = mutations)) + 
         geom_bar(stat="identity") + 
         geom_text(aes(label = mutations), position=position_dodge(width=0.9), vjust=-0.25) +
         xlab("Number of mutations") + 
         ylab("Number of patients") +
         cool_theme + 
         theme(strip.background = element_blank(), 
               strip.text.x = element_text(size = 15))

    ggsave(file.path(DIR_figure, file_name), p, width = fig_width, height = fig_height, units = "cm")

}



patient_mut_plot <- function(DF, DIR_figure, fig_height, fig_width, y_title, file_name) {
    
    df_main <- DF %>%
        filter(Function == "Exonic") %>%
        select(Patient_ID, Gene, Protein_annotation) %>%
        arrange(Patient_ID) %>%
        group_by(Patient_ID) %>%
        mutate(n_muts = n()) %>%
        distinct(Patient_ID, .keep_all = TRUE) %>%
        ungroup() %>%
        group_by(n_muts) %>%
        mutate(mutations = n()) %>%
        arrange(n_muts) %>%
        distinct(mutations, .keep_all = TRUE)

    p <- ggplot(data = df_main, aes(x = n_muts, y = mutations)) + 
         geom_bar(stat="identity") + 
         geom_text(aes(label = mutations), position=position_dodge(width=0.9), vjust=-0.25) +
         xlab("Number of mutations") + 
         scale_x_continuous(breaks = pretty_breaks(7), name = "Number of patients") + 
         cool_theme + 
         theme(strip.background = element_blank(), 
               strip.text.x = element_text(size = 15))

    ggsave(file.path(DIR_figure, file_name), p, width = fig_width, height = fig_height, units = "cm")

}

patient_mut_plot(DF = DF, DIR_figure = DIR_figure, fig_height = 15, fig_width = 15, y_title, "patient_mut_plot.jpg") 



################################################
# LOLLIPOP PLOTS TO COMPARE TUMOR AND WBC VAF AT THE SAME TIME, SEPARATED BY VARIANT TYPE
################################################
tumor_wbc_dotplot(DF = DF, DIR_figure = DIR_figure, fig_height = 25, fig_width = 30)

################################################
# SCATTER PLOTS TO COMPARE TUMOR AND WBC VAFs, SEPARATED BY VARIANT TYPE
################################################
tumor_wbc_scatter(DF = DF, DIR_figure = DIR_figure, fig_height = 15, fig_width = 25, filename = "tumor_wbc_scatter_2.jpg", title = "cfDNA - WBC VAF comparison", xmax = 40, ymax = 40, xpos_equation = 5, ypos_equation = 35)

################################################
# GENE COUNTS
################################################
gene_counts(DF = DF, plot_fill = "Function", DIR_figure = DIR_figure, fig_height = 15, fig_width = 25, color_palette = "BrBG", file_name = "function.jpg", y_title = "Number of mutations")
gene_counts(DF = DF, plot_fill = "Effects", DIR_figure = DIR_figure, fig_height = 15, fig_width = 25, color_palette = "Paired", file_name = "exonic_effects.jpg", y_title = "Number of exonic mutations")
gene_counts(DF = DF, plot_fill = "Variant", DIR_figure = DIR_figure, fig_height = 15, fig_width = 25, color_palette = "Set2", file_name = "variants.jpg", y_title = "Number of exonic mutations")

################################################
# GENE COUNTS
################################################
patient_prevalence(DF = DF, DIR_figure = DIR_figure, fig_height = 20, fig_width = 20)

depth(DF = DF, DIR_figure = DIR_figure, fig_height = 15, fig_width = 45, y_title = "WBC and cfDNA total depth comparison", file_name = "WBC_cfDNA_depth.jpg")


############################################
# Looking for DNMT3A R882H mutations
df_main <- as.data.frame(read_csv("/groups/wyattgrp/users/amunzur/COMPOST_BIN/vardict_before_filtering.csv"))

annot <- "p.V617F"
df <- df_main %>%
        filter(Protein_annotation == annot) %>%
        select(Sample_name, VarFreq, Reads1, Reads2) %>%
        mutate(Patient_ID = gsub("_gDNA.*", "", Sample_name), VarFreq = 100 * VarFreq)

p_vaf <- ggplot(data = df, aes(x = Patient_ID, y = VarFreq)) + 
         geom_bar(stat="identity") + 
         geom_text(aes(label = Reads2), position=position_dodge(width=0.9), vjust=-0.25) +
         xlab("Patients") + 
         ylab("VAF (%)") +
         cool_theme + 
         theme(strip.background = element_blank(), 
               strip.text.x = element_text(size = 15)) +
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(file.path(DIR_figure, paste0("DNMT3A_", annot, ".jpg")), p_vaf, width = 9, height = 15, units = "cm")

# DEPTH PLOTS TO COMPARE CFDNA AND WBC SAMPLES
path_singlex <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/SSCS3/averaged_depths.txt"
path_duplex <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/DCS/averaged_depths.txt"

singlex <- as.data.frame(read_delim(path_singlex, " ")) %>% 
            filter(str_detect(Sample_name, "gDNA")) %>% 
            mutate(Patient_ID = gsub("_gDNA.*", "", Sample_name),
                   Type = "WBC") %>%
            select(Patient_ID, median_depth, Type)

duplex <- as.data.frame(read_delim(path_duplex, " ")) %>% 
            mutate(Patient_ID = gsub("_cfDNA.*", "", Sample_name),
                   Type = "cfDNA") %>%
            select(Patient_ID, median_depth, Type)

df <- rbind(singlex, duplex) %>% 
        pivot_wider(names_from = Type, values_from = median_depth)

p_depth <- ggplot(df, aes(x = WBC, y = cfDNA)) +
                geom_point() + 
                xlim(0, 2500) + ylim(0, 1000) + 
                geom_smooth(method = "lm", se = TRUE, col = "black") + 
                cool_theme + 
                xlab("WBC singlex depth") + 
                ylab("cfDNA duplex depth")

ggsave(file.path(DIR_figure, "consensus_depth_zoomed.jpg"), p_depth, height = 10, width = 10, units = "cm")
