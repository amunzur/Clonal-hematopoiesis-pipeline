script_dir="/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/tnvstats";
configfile="/groups/wyattgrp/users/amunzur/pipeline/workflow/scripts/tnvstats/config.txt";
backgrounderror="/groups/wyattgrp/eritch/projects/gene_panel_pipeline/supp_files/genepanel_allchr_errorrates.txt"

source ${configfile};
source ${conda_profile_path};
conda activate ${conda_pysam_env};
while read name_tumor path_tumor name_normal path_normal
do
    echo $name_tumor
    echo $name_normal
    outputdir="/groups/wyattgrp/users/amunzur/pipeline/results/tnvstats/${name_tumor}";

    tbam=$(readlink -ve "${path_tumor}")
    nbam=$(readlink -ve "${path_normal}")

    printf "bash ${script_dir}/make_tnvstat_file.bash ${tbam} ${nbam} ${name_tumor} ${name_normal} ${outputdir} ${configfile}\n";
    sbatch ${script_dir}/make_tnvstat_file.bash ${tbam} ${nbam} ${name_tumor} ${name_normal} ${outputdir} ${configfile};

done < /groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/tnvstats_sample_list.txt

