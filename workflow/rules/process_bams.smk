rule indel_realignment:
    input:
        MAPPED_bam=DIR_bams + "/mBAM/{wildcard}.bam",
        PATH_bed=PATH_bed,
        PATH_hg38=PATH_hg38,
    params:
        min_mapping_quality=20,
    output:
        DIR_bams + "/abra2/{wildcard}.bam",
    threads: 12
    run:
        shell(
            "java -Xms64G -jar /home/amunzur/anaconda3/envs/snakemake/share/abra2-2.24-1/abra2.jar \
        --in {input.SORTED_bam} --out {output} --ref {input.PATH_hg38} --threads {threads} --index \
        --no-edge-ci --nosort --mmr=0.1 --cons --targets {input.PATH_bed} --tmpdir /groups/wyattgrp/users/amunzur/COMPOST_BIN > \
        '/groups/wyattgrp/users/amunzur/pipeline/results/logs_slurm/indel_realignment/{wildcards.wildcard}'"
        )


# tee sends the stdout (the fixmate bam to multiple commands to avoid writing intermediate files)
rule fixmate_and_recalibrate_bases:
    input:
        DIR_bams + "/abra2/{wildcard}.bam",
    output:
        base_scores=DIR_recalibrated_base_scores + "/{wildcard}.table",
        uncollapsed_BAM=DIR_bams + "/uncollapsed_BAM/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        PATH_known_indels=PATH_known_indels,
        PATH_gold_std_indels=PATH_gold_std_indels,
        PATH_SNP_db=PATH_SNP_db,
    threads: 12
    run:
        shell(
            "picard -Xmx40g FixMateInformation I={input.ABRA2_bam} O=/dev/stdout IGNORE_MISSING_MATES=true SORT_ORDER=coordinate | \
        tee >(gatk BaseRecalibrator -I /dev/stdin -R {params.PATH_hg38} --known-sites {params.PATH_known_indels} --known-sites {params.PATH_gold_std_indels} --known-sites {params.PATH_SNP_db} -O {output.base_scores}) \
        tee >(gatk ApplyBQSR --input /dev/stdin --output {output.uncollapsed_BAM} --bqsr-recal-file {output.base_scores})"
        )
