# Bioinfo_pipeline

This pipeline is developed with the purpose of scanning targeted sequencing data from cancer patients in order to identify CHIP mutations. 
It consists of the following steps: 
- Merging FastQ files from different lanes into R1 and R2 separately for each sample 
- Processing the merged FastQ files (trimming, masking the low quality reads with a base quality less than 30)
- Generating FastQC reports before and afte processing
- Aligning the FastQ files against the Hg38
- Processing the BAM files (sorting, indexing, marking duplicates, adding read groups, 3' and 5' clipping, filtering)
- Running metrics such as insert size distribution, pileup files and coverage
- Variant calling using VarScan2 and VarDict
- Visualization and generating patient reports