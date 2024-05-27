#!/bin/bash
#SBATCH --job-name=fp
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH -n 3 # number of cores
#SBATCH --mem 200000 # memory pool for all cores
#SBATCH -t 5:00:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --output=/groups/wyattgrp/log/%j.out
#SBATCH --error=/groups/wyattgrp/log/%j.out

fr=$(free -hm);
echo "total memory:";
echo "$fr";
datenow=$(date +"%D %H:%M");
printf "starttime: $datenow\n"

## need for cp pipe
bam="$1"
basen="$2"
outputdir="$3"

ref="/groups/wyattgrp/reference/hg38/hg38.fa";

#targetbed="/groups/wyattgrp/reference/target/72gene_panel/grch38genes4intersect.bed";
#gcbed="/groups/wyattgrp/reference/grch38/grch38.gc50Base.bed";
#mindepth=40
#cores=

# printf "Arguments:\n";
# printf "read1 tumor = $read1_t\n";
# printf "read2 tumor = $read2_t\n";
# printf "read1 normal = $read1_n\n";
# printf "read2 normal = $read2_n\n";
# printf "basen = $basen\n";
# printf "outputdir = $outputdir\n";
# printf "reference = $ref\n";
# printf "background snv = $bgsnv\n";
# printf "background indel = $bgindel\n";

source activate pysamstats;
shopt -s extglob;
mkdir -p $outputdir/bam;
mkdir -p $outputdir/tmp;
mkdir -p $outputdir/bam/tmp;
###########################
#tumor alignment to vstat #
###########################

#run pysamstats

pysamstats --type variation -f $ref --no-dup --min-baseq=20  $bam | grep -v ^KI270 | grep -v ^GL0002 > $outputdir/bam/${basen}.tmp.vstat;

#format pysamstats
datenow=$(date +"%D %H:%M");
echo "start adding sample name and linux formatting on $basen tumor at $datenow";
paste <(cat $outputdir/bam/${basen}.tmp.vstat | tr -d '\r' | head -n1) <(echo sample) > $outputdir/bam/${basen}.vstat;
cat $outputdir/bam/${basen}.tmp.vstat | tr -d '\r' | tail -n +2 | awk -F $'\t' -v basen=$basen '{print $0,basen}' FS="\t" OFS="\t" >> $outputdir/bam/${basen}.vstat;
datenow=$(date +"%D %H:%M");
echo "finished adding sample name and linux formatting on $basen tumor at $datenow";
