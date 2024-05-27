rule make_cnvkit_bams:
    input:
        DIR_bams + "/SSCS2_final/{wildcard}.bam",
    output:
        DIR_bams + "/cnvkit_bams/{wildcard}.bam",
    shell:
        """
        samtools view -h -b -q 60 -F 2820 {input} > {output}
        """

rule run_cnvkit_coverage_TARGET:
    input:
        bam = DIR_bams + "/cnvkit_bams/{wildcard}.bam",
        chip_probes_bed = chip_probes_bed,
    output:
        DIR_results + "/cnvkit/coverage/SSCS2/{wildcard}.targetcoverage.cnn",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        cnvkit.py coverage {input.bam} {input.chip_probes_bed} -o {output}
        """

# # DO THIS LATER
# rule generate_cfDNA_reference:
#     input:
#         sample_target=DIR_results + "/cnvkit/coverage/SSCS2/{wildcard}.targetcoverage.cnn"
#         sample_antitarget=DIR_results + "/cnvkit/coverage/SSCS2/{wildcard}.antitargetcoverage.cnn"
#     output:
#         DIR_results + "/references/VIP_cfDNA_reference.cnn",
#     conda:
#         "../envs/cnvkit.yaml"
#     shell:
#         """
#         cnvkit.py reference "${DIR_results}/cnvkit/coverage/SSCS2/"*VIP*cfDNA*coverage.cnn -f ${PATH_hg38} -o "${DIR_results}/cnvkit/references/VIP_cfDNA_reference.cnn" &
#         """

# rule generate_EBC_reference:
#     input:
#         sample_target=DIR_results + "/cnvkit/coverage/SSCS2/{wildcard}.targetcoverage.cnn"
#         sample_antitarget=DIR_results + "/cnvkit/coverage/SSCS2/{wildcard}.antitargetcoverage.cnn"
#     output:
#         DIR_results + "/references/VIP_WBC_reference.cnn",
#     conda:
#         "../envs/cnvkit.yaml"
#     shell:
#         """
#         cnvkit.py reference "${DIR_results}/cnvkit/coverage/SSCS2/"*VIP*WBC*coverage.cnn -f ${PATH_hg38} -o "${DIR_results}/cnvkit/references/VIP_WBC_reference.cnn" &
#         """

# STEP 3. CONSTRUCT A REFERENCE USING THE cfDNA VIP NORMALS


# For each cfDNA sample use the reference to correct for regional coverage and GC bias
rule fx_gc_bias_and_regional_coverage_cfDNA:
    input:
        sample_target=DIR_results + "/cnvkit/coverage/SSCS2/{wildcard}.targetcoverage.cnn",
        dummy_antitarget="/groups/wyattgrp/users/amunzur/pipeline/resources/cnvkit/dummy_antitarget.cnn",
        reference=DIR_results + "/cnvkit/references/VIP_cfDNA_reference.cnn"
    output:
        DIR_results + "/cnvkit/coverage_fixed/cfDNA/{wildcard}.cnr",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        cnvkit.py fix {input.sample_target} {input.dummy_antitarget} {input.reference} -o {output}
        """

# For each cfDNA sample use the reference to correct for regional coverage and GC bias
rule fx_gc_bias_and_regional_coverage_WBC:
    input:
        sample_target=DIR_results + "/cnvkit/coverage/SSCS2/{wildcard}.targetcoverage.cnn",
        dummy_antitarget="/groups/wyattgrp/users/amunzur/pipeline/resources/cnvkit/dummy_antitarget.cnn",
        reference=DIR_results + "/cnvkit/references/VIP_WBC_reference.cnn"
    output:
        DIR_results + "/cnvkit/coverage_fixed/WBC/{wildcard}.cnr",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        cnvkit.py fix {input.sample_target} {input.dummy_antitarget} {input.reference} -o {output}
        """

# For each cfDNA sample use the reference to correct for regional coverage and GC bias
rule segment:
    input:
        DIR_results + "/cnvkit/coverage_fixed/{sample_type}/{wildcard}.cnr"
    output:
        DIR_results + "/cnvkit/segmented/{sample_type}/{wildcard}.cns",
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        cnvkit.py segment {input} -o {output}
        """