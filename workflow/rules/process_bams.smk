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
    shell:
        "java -Xms64G -jar /home/amunzur/anaconda3/envs/snakemake/share/abra2-2.24-1/abra2.jar \
        --in {input.SORTED_bam} --out {output} --ref {input.PATH_hg38} --threads {threads} --index \
		--no-edge-ci --nosort --mmr=0.1 --cons --targets {input.PATH_bed} --targets /groups/wyattgrp/users/jbacon/reference/CHIP_capturetargets.bed --tmpdir /groups/wyattgrp/users/amunzur/COMPOST_BIN > \
        '/groups/wyattgrp/users/amunzur/pipeline/results/logs_slurm/indel_realignment/{wildcards.wildcard}'"


rule PICARD_fixmate:
    input:
        ABRA2_bam=DIR_bams + "/abra2/{wildcard}.bam",
        SORTED_bam_index=DIR_bams + "/sorted/{wildcard}.bam.bai",  # to make sure this doesnt run before we index the aligned bams
    output:
        DIR_bams + "/fixmate/{wildcard}.bam",
    threads: 12
    shell:
        "picard -Xmx40g FixMateInformation \
        I={input.ABRA2_bam} \
        O={output}"



# Add read groups after removing duplicates
rule add_read_groups_PICARD:
    input:
        DIR_bams + "/markdup/{wildcard}.bam",
    params:
        rglb="library",
        rgpl="ILLUMINA",
        rgpu="unit",
        rgsm="sample",
    output:
        DIR_bams + "/readGroup/{wildcard}.bam",
    threads: 12
    shell:
        "picard -Xmx40g AddOrReplaceReadGroups I={input} O={output} RGID=1 RGLB={params.rglb} RGPL={params.rgpl} RGPU={params.rgpu} RGSM={params.rgsm}"
