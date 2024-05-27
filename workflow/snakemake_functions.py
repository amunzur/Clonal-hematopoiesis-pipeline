# def get_cfDNA_vcf_SSCS3(wildcard):
#     cfDNA_vcf =  wildcard.replace("gDNA", "cfDNA") + ".vcf"
#     PATH_complete = "results/variant_calling/Vardict/SSCS3/" + cfDNA_vcf
#     return(PATH_complete)

# def get_cfDNA_bam_SSCS1(wildcard):
#     cfDNA_bam = wildcard.replace("gDNA", "cfDNA") + ".bam"
#     PATH_complete = "results/data/bam/SSCS1_clipped/" + cfDNA_bam
#     return(PATH_complete)

# def get_cfDNA_bam_DCS(wildcard):
#     cfDNA_bam = wildcard.replace("gDNA", "cfDNA") + ".bam"
#     PATH_complete = "results/data/bam/DCS_filtered/" + cfDNA_bam
    # return(PATH_complete)

def get_wbc(wildcards, consensus, bam_type):
    DIR_bams=f"/groups/wyattgrp/users/amunzur/pipeline/results/data/bam/{consensus}_{bam_type}"
    sample_list_path = f"/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/paired_samples.tsv"
    samples_df = pd.read_csv(sample_list_path, sep = "\t")
    wbc_name = samples_df[samples_df["cfDNA"] == wildcards]["WBC"].iloc[0]
    wbc_bam_path = DIR_bams + "/" + wbc_name + ".bam"
    wbc_bai_path = wbc_bam_path + ".bai"    
    return([wbc_bam_path, wbc_bai_path])

def get_wbc_name(wildcards):
    sample_list_path = f"/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/paired_samples.tsv"
    with open(sample_list_path, "r") as f:
        for line in f:
            wbc_sample, cfDNA_sample = line.strip().split("\t")
            if cfDNA_sample == wildcards:
                return(wbc_sample) 

def get_samples_coverage(keyword, sample_type):
    """
    Returns the name of the VIP normals.
    """
    DIR_coverage = "/groups/wyattgrp/users/amunzur/pipeline/results/cnvkit/coverage/SSCS2"
    files = [os.path.join(DIR_coverage, file) for file in os.listdir(DIR_coverage) if keyword in file and sample_type in file and ".targetcoverage.cnn" in file]
    
    print("Matching files:")
    for file in files:
        print(file)
    
    return files

def get_relevant_bams_abra2(wildcards):
    """
    Given a wildcard, which is a full sample including the date, return the cfDNA and WBC pair from that patient at that time point. 
    This function is used to run ABRA2 jointly on all time matched samples from the sample patient.
    To pull the sample names, the function uses a tsv file where all samples are noted. They don't need to be ordered in any way. Single column file.
    """
    sample_list_path = "/groups/wyattgrp/users/amunzur/pipeline/resources/sample_lists/sample_list_bladder.tsv"
    all_samples = pd.read_csv(sample_list_path, sep = "\t")
    all_samples["patient_names"] = all_samples["sample_names"].str.split("_").str.get(0)
    all_samples["timepoint"] = all_samples["sample_names"].str.split("-").str.get(-1)
    patient = wildcards.split("_")[0]
    timepoint = wildcards.split("-")[-1]
    relevant_samples = all_samples[(all_samples["patient_names"] == patient) & (all_samples["timepoint"] == timepoint)]["sample_names"].tolist() # all relevant samples from patient
    # Ensure cfDNA is given first, then the WBC is given.
    sorted_samples = sorted(relevant_samples, key=lambda x: ('cfDNA' not in x, x))
    # now add the paths to the bam files
    DIR_bams = "/groups/wyattgrp/users/amunzur/pipeline/results/data/bam/"
    return [DIR_bams + x + ".bam" for x in sorted_samples]