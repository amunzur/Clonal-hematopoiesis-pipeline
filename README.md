# Bioinfo_pipeline

This pipeline is developed with the purpose of scanning targeted sequencing data from cancer patients in order to identify CHIP mutations.
It uses Snakemake as a framework with a combination of R, Python and Bash scripts. Note that this project is still in development! 

It consists of the following steps:
- Merging FastQ files from different lanes into R1 and R2 separately for each sample
- Processing the merged FastQ files (trimming, masking the low quality reads with a base quality less than 30)
- Generating FastQC reports before and afte processing
- Aligning the FastQ files against the Hg38
- Processing the BAM files (sorting, indexing, marking duplicates, adding read groups, 3' and 5' clipping, filtering)
- Running metrics such as insert size distribution, pileup files and coverage
- Variant calling using VarScan2, VarDict and HaplotypeCaller from GATK
- Visualization and generating patient reports

`workflow/scripts/analysis` contains the following scripts to further process the called variants. They should be run in the given order after processing the individual samples through the pipeline:
- **compute_averaged_depth.sh**: Computes the average depth of the BAM files, only considers regions in the panel.
- **seq_quality_metrics.py**: Combined sequencing quality metrics such as number of reads taken from raw FastQ files, average and median depth of aligned and processed bams and percentage of duplicates. Should be run for every cohort outside of the framework of Snakemake.
- **filter_vardict.R**, **filter_VarScan.R** and **filter_gatk.R** combines and filters the variants called by 3 variant callers across all samples in cohort based on the user defined parameters and thresholds. As output, it generates batch specific files.
- **combine_variant_callers.R**: Combines the output from all 3 variant callers into one large file where one unique variant is found in each line. It also adds 3 new columns to indicate which variant callers identified the variant. It also counts how many callers identified the variant of interest.
- **compare_with_tumor.R**: In the cases where both WBC and tumor samples are available, outputs a table only with variants with read support in both sample types. It also combined all the batch files into one large file, so no more batch specific results!

Utilities functions used by some scripts above are found here:
- **UTILITIES_filter_vardict.R**: Utilities functions for `filter_vardict.R`
- **UTILITIES_filter_varscan.R**: Utilities functions for `filter_VarScan.R`
- **UTILITIES_filter_gatk.R**: Utilities functions for `filter_gatk.R`
- **UTILITIES.R**: Utilities functions for both `filter_vardict.R` and `filter_VarScan.R`

These scripts are used in various steps of the pipeline:
- **make_anno_input_indel.R**: Modifies the indel outputs to put in a format acceptable by ANNOVAR.
- **make_anno_input_vardict.R**: Modifies the outputs by VarDict to put in a format acceptable by ANNOVAR.
- **pdf_to_png.sh**: Converts the PDF images to PNG, mainly used for visuzalizing the insert size plots.
- **process_bets.R**: Cleans up and melts the betastasis tables to retain only the called variants.
- **reformat_vardict.R**: Minor reformatting
- **filter_tnvstats.sh**: Filters the tnvstats to only retain the regions in the panel.

### Before running the pipeline
The pipeline needs a couple files to be present before it can start processing the files. These files need to be present for each batch, and they are: 
1. A text outlining the sample names, each sample name should be in a new line. 
2. Dummy files for raw FastQs. Sometimes the files come already merged from the sequencer, but these empty dummy files need to be created so that Snakemake doesn't complain about missing input files. `workflow/scripts/analysis/create_dummy_files.py` uses the identifier excel sheets to create these files. 

### Renaming FastQ files
Before we can do any analysis on the data, the raw FastQ files from the sequencer need to be renamed, that is the molecular identifiers in the file names need to be
match with the sample names. We use the identifier sheets found in `workflow/results/identifiers`. Scripts that are not a part of the pipeline (yet!) process the 
molecular IDs and rename the files accordingly. This process consists of finding the reverse compliment of the second barcode and matching it with the sequencing id. 