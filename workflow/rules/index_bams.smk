rule index_bams:
    input:
        DIR_results + "/data/bam/{consensus_type}_filtered/{wildcard}.bam"
    output:
        DIR_results + "/data/bam/{consensus_type}_filtered/{wildcard}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"