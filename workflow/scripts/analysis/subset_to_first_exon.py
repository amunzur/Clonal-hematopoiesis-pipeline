#!/bin/bash
import os
import subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

# Directories and file paths
DIR_bams="/groups/wyattgrp/users/amunzur/pipeline/results/data/bam/SSCS2_final"
PATH_bed="/groups/wyattgrp/users/amunzur/pipeline/resources/panel/chip_first_exon_bed.tsv"

DIR_output_bams="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/tss_bam_dir"
DIR_insert_sizes="/groups/wyattgrp/users/amunzur/pipeline/results/metrics/tss_insert_sizes"
DIR_slurm_scripts="/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/batch_scripts"

bam_list = sorted([os.path.join(DIR_bams, bam) for bam in os.listdir(DIR_bams) if bam.endswith(".bam")])
bed = pd.read_csv(PATH_bed, sep = "\t", names = ["chrom", "tss_start", "tss_end", "gene"])

# Subset the bams
with ThreadPoolExecutor(max_workers=200) as executor:
    futures = []
    for bam in bam_list:
        if "cfDNA" in bam:
            print(bam)
            for i, row in bed.iterrows():
                chrom = row["chrom"]
                tss_start = row["tss_start"]
                tss_end = row["tss_end"]
                gene = row["gene"]
                tss_bam = os.path.join(DIR_output_bams, os.path.basename(bam).replace(".bam", f"_{gene}_TSS.bam"))
                samtools_command = f"samtools view -b {bam} {chrom}:{tss_start}-{tss_end} -o {tss_bam} && samtools index {tss_bam}"
                futures.append(executor.submit(subprocess.run, samtools_command, shell=True))

# Wait for all futures to complete
for future in futures:
    future.result()

# Run insert size in subsetted bams
with ThreadPoolExecutor(max_workers=150) as executor:
    futures = []
    for bam in bam_list:
        if "cfDNA" in bam:
            print(bam)
            for i, row in bed.iterrows():
                chrom = row["chrom"]
                tss_start = row["tss_start"]
                tss_end = row["tss_end"]
                gene = row["gene"]
                tss_bam = os.path.join(DIR_output_bams, os.path.basename(bam).replace(".bam", f"_{gene}_TSS.bam"))
                insert_size_metrics = os.path.join(DIR_insert_sizes, os.path.basename(bam).replace(".bam", f"_{gene}_TSS.txt"))
                insert_size_metrics_pdf = insert_size_metrics.replace(".txt", ".pdf")
                picard_command = f"picard CollectInsertSizeMetrics I={tss_bam} O={insert_size_metrics} H="{insert_size_metrics%.txt}.pdf" M=0.5"
                futures.append(executor.submit(subprocess.run, samtools_command, shell=True))





# Read the BED file line by line
while IFS=$'\t' read -r chrom start end gene; do
    # Loop through each BAM file in the directory
    for bam_file in "$DIR_bams"/*.bam; do
        # Extract the file name without the path and extension
        bam_filename=$(basename "$bam_file" .bam)
        
        # Define the output BAM file path
        output_bam="$DIR_output_bams/${bam_filename}_${gene}_first_exon.bam"
        
        # Define the insert size metrics output file path
        insert_size_metrics="$DIR_insert_sizes/${bam_filename}_${gene}_insert_size_metrics.txt"
        
        # Define the SLURM script path
        slurm_script="$DIR_slurm_scripts/${bam_filename}_${gene}_process_bam.slurm"
        
        # Write the SLURM script
        cat <<EOT > "$slurm_script"
#!/bin/bash
#SBATCH --job-name=process_${bam_filename}_${gene}
#SBATCH --output=${DIR_slurm_scripts}/${bam_filename}_${gene}_process_bam.out
#SBATCH --error=${DIR_slurm_scripts}/${bam_filename}_${gene}_process_bam.err
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

source /home/amunzur/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

# Subset the BAM file to the regions in the BED file
samtools view -b "$bam_file" "$chrom:$start-$end" -o "$output_bam"

# Index the subset BAM file
samtools index "$output_bam"

# Run Picard CollectInsertSizeMetrics

echo "Processed $bam_file for $gene and generated insert size metrics."
EOT

        # Optionally, submit the SLURM script
        # sbatch "$slurm_script"
        
        echo "Generated SLURM script for $bam_file for gene $gene."
    done
done < "$PATH_bed"
