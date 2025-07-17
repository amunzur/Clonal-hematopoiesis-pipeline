# Comprehensive Snakemake Workflow for Genomic Data Processing and Variant Calling

## Table of Contents

- [Overview](#overview)  
- [Directory Structure and Variables](#directory-structure-and-variables)  
- [Pipeline Workflow](#pipeline-workflow)  
  - [Fastq Processing](#fastq-processing)  
  - [Bam Processing and Consensus BAM Generation](#bam-processing-and-consensus-bam-generation)  
  - [Quality Control](#quality-control)  
  - [Depth, Coverage and Metrics](#depth-coverage-and-metrics)  
  - [Variant Calling](#variant-calling)  
    - [Somatic Variant Calling](#somatic-variant-calling)  
    - [CHIP Variant Calling](#chip-variant-calling)  
  - [VCF Post-processing](#vcf-post-processing)  
- [Tools and Environments](#tools-and-environments)  
- [Usage](#usage)  
- [Contact](#contact)  

---

## Overview

This Snakemake pipeline is designed for processing paired-end sequencing data to generate high-quality BAM files and perform comprehensive variant calling. It covers:

- Fastq quality control, trimming, and read count metrics  
- BAM file generation from fastq, including consensus read calling using UMIs  
- Alignment, indel realignment, base quality score recalibration  
- Depth and insert size metrics calculations  
- Variant calling using multiple callers: Mutect2, VarDict, FreeBayes  
- VCF normalization, sorting, indexing, and decomposition for downstream analysis
---

## Directory Structure and Variables

The workflow expects several directory variables to be set in the `config.yaml` file:

- `DIR_trimmed_fastq`: directory with trimmed fastq files  
- `DIR_fastq`: directory with raw or merged fastq files  
- `DIR_bams`: directory for all BAM outputs (including intermediate and final BAMs, subset BAM directories will be automatically generated)  
- `DIR_results`: main results directory for variant calling, metrics, and other outputs  
- `DIR_umi_metrics`: directory for UMI family size histograms  
- `DIR_readcounts_metrics`: directory for read count metrics  
- `DIR_depth_metrics`: directory for depth of coverage outputs  
- `DIR_insertsize_metrics` & `DIR_insertsize_figures`: insert size metrics and plots  
- `DIR_HS_metrics`: hybrid selection metrics output  
- `DIR_metrics`: directory for other picard metrics (alignment summary etc.)  
- `PATH_hg38`: path to the hg38 reference fasta file  
- `PATH_bed`: target regions BED file  
- `PATH_baits` & `PATH_bed_intervals`: bait and intervals BED files for hybrid capture  
- `PATH_known_indels`, `PATH_gold_std_indels`, `PATH_SNP_db`: known sites for base recalibration  

---

## Pipeline Workflow

### Fastq Processing

- `run_fastqc_merged`: Run FastQC on merged raw fastq files, producing `.zip` and `.html` reports.  
- `fastq_read_counts`: Calculate read counts from fastq files and output sample-wise counts as text.  
- `trim_fastq`: Trim adapters and low-quality bases using `fastp`, generating trimmed fastq with reports (HTML/JSON).  
- `run_fastqc_trimmed`: Run FastQC again on trimmed fastq files for quality confirmation.

You can then run MultiQC on these files to generate one QC report.

---

### BAM Processing and Consensus BAM Generation

- `FastqToBam`: Converts paired fastq files into unaligned BAMs (uBAM) using `fgbio FastqToBam` preserving UMIs.  
- `BamtoFastq`: Convert uBAMs back to fastq for remapping.  
- `mapBAM`: Map fastq reads to reference genome using `bwa mem`.  
- `MergeBamAlignment`: Merge aligned BAM with original unaligned BAM to retain original info using Picard.  
- `indel_realignment`: Realign indels using ABRA2.  
- `fixmate`: Fix mate information and sort BAMs with Picard.  
- `recalibrate_bases` and `apply_base_scores`: Perform GATK base quality score recalibration.  
- `GroupReadsByUmi`: Group reads by UMI to identify consensus families.  
- `CallMolecularConsensusReads`: Generate consensus reads (consensus BAMs) based on UMI families.

---

### Quality Control

- `index_bams`: Index BAM files using `samtools index`.  
- `run_depth`: Calculate depth of coverage across target regions with `samtools depth`.  
- `run_mpileup`: Generate pileup files with `samtools mpileup`.  
- `run_insert_size`: Collect insert size metrics and plots using Picard.  
- `run_read_counts`: Calculate read counts for merged fastq files.  
- `hs_metrics`: Collect hybrid selection (capture) metrics using Picard.  
- `alignment_summary_metrics`: Generate alignment summary metrics using Picard.

---

### Variant Calling

#### Somatic Variant Calling

- `run_VarDict_somatic`: Run VarDict paired variant caller on cfDNA and matched WBC BAMs.  
- `run_mutect2_somatic`: Run GATK Mutect2 somatic variant caller on tumor and matched normal BAMs.  
- `run_freebayes_somatic`: Run FreeBayes variant caller on paired BAMs.  
- Post-process VCF files by compressing, indexing, sorting, normalizing, and decomposing using `bgzip`, `tabix`, `bcftools`, and `vt`.

#### Clonal Hematopoiesis Variant Calling

- `run_mutect2`: Run Mutect2 on consensus BAMs.  
- `unzip_mutect`: Decompress Mutect2 VCFs.  
- `run_VarDict_chip`: Run VarDict on consensus BAMs with lower thresholds appropriate for CHIP detection.  
- `run_freebayes_chip`: Run FreeBayes on consensus BAMs.  
- Post-processing steps similar to somatic workflow for compressing, indexing, sorting, normalizing, and decomposing VCFs.

---

### VCF Post-processing

- Compression of VCFs with `bgzip`.  
- Indexing of compressed VCFs using `tabix`.  
- Sorting of VCFs with `bcftools sort`.  
- Normalizing variants (splitting multiallelics, left-aligning indels) with `bcftools norm`.  
- Decomposing complex block substitutions into simpler forms with `vt decompose_blocksub`.

---

## Tools and Environments

- **Alignment & BAM manipulation:** `bwa mem`, `picard`, `samtools`, `fgbio`, `abra2`  
- **Quality Control:** `fastp`, `FastQC`, Picard tools  
- **Variant Calling:** `GATK Mutect2`, `VarDictJava`, `FreeBayes`  
- **VCF processing:** `bcftools`, `tabix`, `bgzip`, `vt`  
- **Base Recalibration:** `GATK BaseRecalibrator`, `ApplyBQSR`  
- **UMI consensus:** `fgbio` GroupReadsByUmi and CallMolecularConsensusReads  

Conda environment YAML files are specified for reproducibility and are provided in the `envs` directory. Snakemake will automatically generate the required environments.

---

## Post-processing tools

Following variant calling, a series of custom R and Python scripts are used to filter, combine, curate, and visualize variants for downstream interpretation. Please refer to each individual script for example usage.

### STEP1_filter_variants.R  
Applies initial filtering criteria to raw variant call files to remove low-quality and likely false-positive calls.

### STEP2_combine_variant_callers.R  
Merges variant calls from multiple callers into a unified dataset, resolving overlaps and conflicts.

### STEP3_make_IGV_snapshots.py  
Automates generation of IGV screenshots around variants of interest to facilitate manual review. After running this script the user needs to go through IGV screenshots and delete those that fail manual inspection.

### STEP4_curate_mutations.py  
Performs additional mutation curation based on manually curated IGV screenshots.

### STEP5_perform_dependent_calling.py  
Executes dependent or secondary variant calling steps that leverage curated variant lists or additional input data.

### combine_variant_callers_UTILITIES.R & UTILITIES.R  
Contain helper functions and utilities used across the above scripts for data manipulation, annotation, and quality control.

---

## Usage

1. Configure your environment variables and directory paths (e.g., `DIR_bams`, `DIR_results`, `PATH_hg38`, `PATH_bed`, etc.) in the Snakemake config or directly in the workflow script.  
2. Prepare input data: paired-end FASTQ files.  
3. Generate the snakemake environment using the .yaml file in `envs/snakemake_env.yaml`
4. Run the pipeline with Snakemake, specifying number of cores, e.g.:  
   ```bash
   snakemake --cores 12

