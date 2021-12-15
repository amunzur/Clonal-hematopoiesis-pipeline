# Bioinfo_pipeline

This pipeline is developed with the purpose of scanning targeted sequencing data from cancer patients in order to identify CHIP mutations. 
It uses Snakemake as a framework with a combination of R, Python and Bash. 

It consists of the following steps:  
- Merging FastQ files from different lanes into R1 and R2 separately for each sample  
- Processing the merged FastQ files (trimming, masking the low quality reads with a base quality less than 30)  
- Generating FastQC reports before and afte processing  
- Aligning the FastQ files against the Hg38  
- Processing the BAM files (sorting, indexing, marking duplicates, adding read groups, 3' and 5' clipping, filtering)  
- Running metrics such as insert size distribution, pileup files and coverage  
- Variant calling using VarScan2 and VarDict  
- Visualization and generating patient reports  

`workflow/scripts/analysis` contains the following scripts to further process the called variants:   
- **compute_averaged_depth.sh**: Computes the average depth of the BAM files, only considers regions in the panel.  
- **filter_tnvstats.sh**: Filters the tnvstats to only retain the regions in the panel.  
- **filter_vardict.R**: Combines and filters the variants called by VarDict across all samples in cohort based on the user defined parameters and thresholds.  
- **filter_VarScan.R**: Combines and filters the variants called by VarScan across all samples in cohort based on the user defined parameters and thresholds  
- **make_anno_input_indel.R**: Modifies the indel outputs to put in a format acceptable by ANNOVAR.   
- **make_anno_input_vardict.R**: Modifies the outputs by VarDict to put in a format acceptable by ANNOVAR.  
- **pdf_to_png.sh**: Converts the PDF images to PNG, mainly used for visuzalizing the insert size plots.  
- **process_bets.R**: Cleans up and melts the betastasis tables to retain only the called variants.  
- **reformat_vardict.R**: Minor reformatting  
- **UTILITIES_filter_vardict.R**: Utilities functions for `filter_vardict.R`  
- **UTILITIES_filter_varscan.R**: Utilities functions for `filter_VarScan.R` 
- **UTILITIES.R**: Utilities functions for both `filter_vardict.R` and `filter_VarScan.R` 
