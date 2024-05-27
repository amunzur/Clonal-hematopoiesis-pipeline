#!/bin/bash
#SBATCH --job-name=ErrDist
#SBATCH -p debug,express,normal,big-mem,long
#SBATCH -n 40 # number of cores
#SBATCH --mem 500000 # memory pool for all cores
#SBATCH -t 23:59:00 # time (D-HH:MM or HH:MM:SS)
#SBATCH --export=all
#SBATCH --chdir=/groups/wyattgrp/log
#SBATCH --output=/groups/wyattgrp/log/%j.log
#SBATCH --error=/groups/wyattgrp/log/%j.log
#printf "SLURM_JOB_ID=$SLURM_JOB_ID\n";

conda_profile_path="/home/$USER/mambaforge/etc/profile.d/conda.sh";
conda_env="samspider";
threads=20;
vstat_dir="/groups/wyattgrp/users/jbacon/err_rates/prostate_exons/NovaSeq/vstat";
output_dir="/groups/wyattgrp/users/jbacon/err_rates/prostate_exons/NovaSeq/error_rate/";
compute_error_script="/groups/wyattgrp/users/jbacon/elie_bgerrdata/scripts/compute_error.py";

source $conda_profile_path;
conda activate $conda_env;
shopt -s extglob;

mkdir -p $output_dir;
tmp_dir=$(mktemp -d -p ${output_dir});
cd ${tmp_dir};

chromchr=$(printf "$(seq 1 22)\nX\nY\nM\n" | sed 's/^/chr/g');

##########################################################
# subset by chromosome 
##########################################################

printf "subset by chromosome $(date +"%D %H:%M")\n"
echo "$chromchr" | parallel --env tmp_dir --env vstat_dir -j $threads "cat $vstat_dir/*.vstat | grep {}[[:space:]] | cut -f 1-4,6,8,10,12,14,16,18,20,22,24 | awk '\$4 >= 1' | sort --temporary-directory=$tmp_dir -k2,2n >> $tmp_dir/{}.tsv";

##########################################################
# split chromosome into smaller files with suffix _split
##########################################################

printf "split into n files $(date +"%D %H:%M")\n"
echo "$chromchr" | parallel --env tmp_dir -j $threads "split --suffix-length=4 --numeric-suffixes=1000 --lines=1000000 --verbose $tmp_dir/{}.tsv $tmp_dir/{}_split";

# Make sure that positions arent split between files
printf "check split files for position splits $(date +"%D %H:%M")\n"
for i in $(echo "$chromchr"); do
	chromosome=$(echo $i);
	for j in $(ls | grep ^${chromosome}_split | sort -n); do
		number=$(echo $j | sed "s/${i}_split//g");
		number_plus1=$(echo "$number + 1" | bc -l);
		onenottwo=$(echo ${#number_plus1}); #when remove numberplusone_if there is still a number left. then a zero was added to front. This happens at splits of large files. ie exome or greater. Need to put in if here for this occurence.
		if [ "$onenottwo" -eq "1" ]; then
			fileplusone=${chromosome}_split0${number_plus1};
		else
			fileplusone=${chromosome}_split${number_plus1};
		fi;

		if [ -f $fileplusone ]; then
			datenow=$(date +"%D %H:%M");
			printf "\nchromosome=$chromosome \nj=$j \nnumber=$number \nnumber_plus1=$number_plus1 \nonenottwo=$onenottwo \nfileplusone=$fileplusone\n";
			value_needs_adding=$(tail -n1 $j | cut -f 2);
			printf "value_needs_adding = $value_needs_adding\n";
			printf "head -n 100 $fileplusone | grep ^$chromosome[[:space:]]$value_needs_adding >> $j\n";
			head -n 100 $fileplusone | grep ^$chromosome[[:space:]]$value_needs_adding >> $j;
			#printf "sed -i.bak '/^$chromosome[[:space:]]$value_needs_adding[[:space:]]/d' $fileplusone\n";
			#sed -i.bak '/^$chromosome[[:space:]]$value_needs_adding[[:space:]]/d' $fileplusone;
			printf "grep -v ^$chromosome[[:space:]]$value_needs_adding[[:space:]] $fileplusone > ${fileplusone}_fix\n";
			grep -v ^$chromosome[[:space:]]$value_needs_adding[[:space:]] $fileplusone > ${fileplusone}_fix;
			rm $fileplusone;
			mv ${fileplusone}_fix $fileplusone;
		else
			printf "$fileplusone does not exist; $j is the last file for $chromosome\n";
		fi;
	done;
done;
datenow=$(date +"%D %H:%M");
printf "finished repair after split $datenow\n";

##########################################################
# compute lambda 
##########################################################

printf "Compute mean errors $(date +"%D %H:%M")\n"
ls *_split* | parallel --env conda_profile_path --env conda_env --env compute_error_script -j $threads "source $conda_profile_path; conda activate $conda_env; python $compute_error_script -in {} -out {}_errors"; 

##########################################################
# make final file 
##########################################################

#faster than doing in python
printf "final processing $(date +"%D %H:%M")\n"
head -n1 $(ls *_errors | head -n1) > ../error_rates.tsv;
tail -n +2 *_errors | grep -v "[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0[[:space:]]0\.0" | sort -k1,1V -k2,2n >> ../error_rates.tsv;

cd ../;
#cp /groups/wyattgrp/log/$SLURM_JOB_ID.log ${output_dir}/$SLURM_JOB_ID.log;
rm -r ${tmp_dir};
final_file=$(readlink -ve error_rates.tsv)
printf "Finished.\nFinal error rate file $final_file completed at $(date +"%D %H:%M")\n"
