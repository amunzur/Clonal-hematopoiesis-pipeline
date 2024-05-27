rule make_vstat:
    input:
        bam=DIR_results + "/data/bam/SSCS1_filtered/{wildcard}.bam",
        PATH_hg38=PATH_hg38,
    output:
        DIR_results + "/data/bg/vstat/{wildcard}.tmp.vstat",
    threads: 12
    params:
        min_baseq=20,
    conda:
        "../envs/pysamstats.yaml"
    shell:
        """
            pysamstats --type variation -f {input.PATH_hg38} --no-dup --min-baseq={params.min_baseq}  {input.bam} | grep -v ^KI270 | grep -v ^GL0002 > {output}
        """


rule format_vstats:
    input:
        DIR_results + "/data/bg/vstat/{wildcard}.tmp.vstat",
    output:
        DIR_results + "/data/bg/vstat/{wildcard}.vstat",
    params:
        basen="{wildcard}"
    shell:
        """
        paste <(cat {input} | tr -d '\r' | head -n1) <(echo sample) > {output};
        cat {input} | tr -d '\r' | tail -n +2 | awk -F $'\t' -v basen={params} '{print $0,basen}' FS="\t" OFS="\t" >> {output};
        """

# rule make_bg_file:
#     input:
#         DIR_results + "/data/bg/vstat/{wildcard}.vstat",
#     output:
#         dir(DIR_results + "/data/bg/bg")
#     params:
#         script="/groups/wyattgrp/scripts/make_background_error2/make_background_error.bash"
#     shell:
#         """
#             CHROM="$(seq 1 22) X Y MT";
#             for i in $(echo $CHROM);do
#                 {params.script} ${vstatdir} {output} $i;
#             done
#         """




# STEP 2:
# outdir="/groups/wyattgrp/users/amunzur/chip_project/data/background";
# vstatdir="/groups/wyattgrp/users/amunzur/chip_project/data/vstats";
# ;
# mkdir -p ${outdir};
