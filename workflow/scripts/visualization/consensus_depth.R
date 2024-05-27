library(tidyverse)
library(scales)
library(ggrepel)

PATH_depth_SSCS1 <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/SSCS1/averaged_depths.txt"
PATH_depth_DCS <- gsub("SSCS1", "DCS", PATH_depth_SSCS1)
PATH_read_counts <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/read_counts/merged/read_counts.txt"
PATH_sample_sheet <- "/groups/wyattgrp/users/amunzur/data/sample_sheet.csv"

DIR_figure <- "/groups/wyattgrp/users/amunzur/pipeline/results/figures/main_figures/consensus_figures"

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

# SEQUENCING BATCH 3 - POOL 1
DF_reads <- read_delim(PATH_read_counts, "\t") %>% mutate(Reads_count = round(Reads_count / 1000000, 3)) %>% rename(Reads_count_millions = Reads_count)
DF_depth_DCS <- read_delim(PATH_depth_DCS, " ")
DF_depth_SSCS1 <- read_delim(PATH_depth_SSCS1, " ") %>% mutate(Patient_ID = gsub("_cfDNA.*|_gDNA.*", "", Sample_name), Sample_type = str_extract(Sample_name, "cfDNA|gDNA"))
DF_sample_info <- read_csv(PATH_sample_sheet) %>% select("Sequencing ID", "Sequencing batch", Pool) %>% rename(Sample_name = "Sequencing ID", Sequencing_batch = "Sequencing batch")
DF <- inner_join(DF_depth_SSCS1, DF_sample_info, by = "Sample_name")

make_consensus_depth_plots <- function(DF, PATH_read_counts, consensus_type, plot_name, height, width, DIR_figure, title) {

    DF <- DF %>%
          filter(Sequencing_batch == 3) %>%
          mutate(median_depth_SSCS1 = as.numeric(median_depth_SSCS1), 
                 median_depth_DCS = as.numeric(median_depth_DCS), 
                 Reads_count_millions = Reads_count_millions*2)

    p <- ggplot(DF, aes_string(x = "Reads_count_millions", y = paste0("median_depth_", consensus_type))) + 
         geom_point(size = 2) + 
         geom_text_repel(aes(label = Patient_ID)) +
         cool_theme + 
         scale_x_continuous(breaks = pretty_breaks(8), name = "Read count (millions)") + 
         scale_y_continuous(breaks = pretty_breaks(8), name = "Depth") + 
         ggtitle(title)

    ggsave(file.path(DIR_figure, paste0(consensus_type, "_depth_vs_read_count.jpg")), p, height = 20, width = 20, units = "cm")

}

make_consensus_depth_plots(DF = DF, consensus_type = "SSCS1", plot_name = "SSCS1_depth.jpg", height = 20, width = 80, DIR_figure = DIR_figure, title = "SINGLEX CONSENSUS DEPTH")
make_consensus_depth_plots(DF = DF, consensus_type = "DCS", plot_name = "DCS_depth.jpg", height = 20, width = 80, DIR_figure = DIR_figure, title = "DUPLEX CONSENSUS DEPTH")

# COMPARING THE DEPTH OF ALL SSCS1 WBC AND CFDNA SAMPLES
PATH_depth <- "/groups/wyattgrp/users/amunzur/pipeline/results/metrics/averaged_depth/SSCS1/averaged_depths.txt"

DF <- read_delim(PATH_depth, " ") %>%
        as.data.frame() %>%
        mutate(Sample_type = str_extract(Sample_name, "cfDNA|gDNA"),
               Patient_ID = gsub("_cfDNA.*|_gDNA.*", "", Sample_name),
               Patient_ID = gsub("GUBB", "GU", Patient_ID)) %>%
        left_join(x = .,  y = DF_sample_info) %>%
        distinct(Sample_name, .keep_all = TRUE) %>%
        mutate(label_color = case_when(Sequencing_batch == 3 ~ "blue",
                                       TRUE ~ "black"))


make_depth_plot <- function(DF, sample_type, DIR_figure, title, height, width) {
    
    DF <- DF %>%
            filter(Sample_type == paste(sample_type), 
                   !grepl("VIP", Patient_ID)) %>%
            arrange(median_depth) %>%
            mutate(Patient_ID = factor(Patient_ID, levels = unique(Patient_ID)))

    p <- ggplot(DF, aes(x = Patient_ID, y = median_depth)) + 
         geom_bar(stat="identity") +
         cool_theme + 
         scale_y_continuous(breaks = pretty_breaks(8), name = "Depth") + 
         xlab("Patient ID") +
         geom_hline(yintercept = median(DF$median_depth), color = "red", linetype = 2, size = 1) +
         geom_text(aes(7, median(DF$median_depth), label = paste("median:", median(DF$median_depth)), vjust = -1), color = "red", size = 5, show.legend = FALSE) +
         ggtitle(title) + 
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = DF$label_color))

    ggsave(file.path(DIR_figure, paste0(sample_type, "_SSCS1_depth.jpg")), p, height = height, width = width, units = "cm")

}

make_depth_plot(DF, sample_type = "cfDNA", DIR_figure, title = "cfDNA SSCS1 DEPTH", height = 20, width = 50)
make_depth_plot(DF, sample_type = "gDNA", DIR_figure, title = "WBC SSCS1 DEPTH", height = 20, width = 50)

# All samples starting with 22, comparing cfDNA depth to gDNA depth
make_depth_plot2 <- function(DF, sample_type, DIR_figure, title, height, width) {
    
    DF <- DF %>%
            arrange(median_depth) %>%
            mutate(Patient_ID = factor(Patient_ID, levels = unique(Patient_ID)))

    p <- ggplot(DF, aes(x = Patient_ID, y = median_depth, fill = Sample_type)) + 
         geom_bar(position="dodge", stat="identity") +
         cool_theme + 
         scale_y_continuous(breaks = pretty_breaks(8), name = "Depth") + 
         xlab("Patient ID") +
         geom_hline(yintercept = median(DF$median_depth), color = "red", linetype = 2, size = 1) +
         geom_text(aes(7, median(DF$median_depth), label = paste("median:", median(DF$median_depth)), vjust = -1), color = "red", size = 5, show.legend = FALSE) +
         ggtitle("Depth - last batch of samples")
         
    ggsave(file.path(DIR_figure, "_SSCS1_depth_last_batch.jpg"), p, height = height, width = width, units = "cm")

}

DF_depth_SSCS1 <- DF_depth_SSCS1 %>%
                  filter(startsWith(Sample_name, "GUBB-22"))