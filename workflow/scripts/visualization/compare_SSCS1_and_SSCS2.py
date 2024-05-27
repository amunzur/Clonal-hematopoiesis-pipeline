import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


def plot_comparison(df, vaf_column, depth_column, vaf_ylabel, depth_ylabel, title, filename):
    fig, ax1 = plt.subplots(figsize=(9, 5))
    # Plot VAF as bars on the primary y-axis
    ax1.bar(df.index, df[vaf_column], color='blue')
    ax1.set_xlabel('Sample specific mutations')
    ax1.set_ylabel(vaf_ylabel, color='blue')
    ax1.set_title(title)
    ax1.set_xticks([])
    ax1.tick_params(axis='y', labelsize=12)
    # Create a secondary y-axis for depth
    ax2 = ax1.twinx()
    ax2.plot(df.index, df[depth_column], 'ro', label='Depth', markersize = 3)  # Use dots for depth
    ax2.set_ylabel(depth_ylabel, color='red')
    ax2.tick_params(axis='y', labelsize=12)
    # Remove top and right spines
    ax1 = remove_frame(ax1, ["top", "bottom"])
    ax2 = remove_frame(ax2, ["top", "bottom"])
    # ylims
    ax1.set_ylim((-3, 3))
    ax2.set_ylim((-1500, 1500))
    fig.tight_layout()
    fig.savefig(filename)

def plot_gene_counts(df, title, save_path):
    gene_counts = df['Gene'].value_counts()
    # Create a bar chart using ax
    fig, ax = plt.subplots(figsize=(10, 6))
    gene_counts.plot(kind='bar', color='skyblue', ax=ax)
    ax.set_title(title)
    ax.set_xlabel('Genes')
    ax.set_ylabel('Count')
    ax.tick_params(axis='x', rotation=45, ha='right')
    fig.tight_layout()
    fig.savefig(save_path)

def remove_frame(ax, spines_to_remove):
    for spine in spines_to_remove:
        ax.spines[spine].set_visible(False)
    return ax

def plot_gene_vaf(df, figuresize, title, save_path, annotate_protein, annotate_depth):
    # Assuming 'Gene' is the column containing gene names, 'VAF_tSSCS1' contains the y-axis values,
    # and 'Protein_annotation' contains the annotation values
    gene_vaf = df[['Gene', 'VAF_tSSCS1', 'Protein_annotation', 'Depth_tSSCS1']].copy()
    gene_vaf["VAF_tSSCS1"] = gene_vaf["VAF_tSSCS1"]*100
    # Create a scatter plot using ax
    fig, ax = plt.subplots(figsize=figuresize)
    ax.scatter(gene_vaf['Gene'], gene_vaf['VAF_tSSCS1'], color='red', s = 4)
    ax.set_title(title)
    ax.set_xlabel('Genes')
    ax.set_ylabel('SSCS1 VAF (%)')
    ax.tick_params(axis='x', rotation=90)  # Adjust the rotation angle as needed
    ax.tick_params(axis='x')
    ax = remove_frame(ax, ["top", "right"])
    # Annotate each dot with the corresponding 'Protein_annotation' value
    if annotate_protein:
        for i, txt in enumerate(gene_vaf['Protein_annotation']):
            ax.annotate(txt, (gene_vaf['Gene'].iloc[i], gene_vaf['VAF_tSSCS1'].iloc[i]), ha='left', va='center', fontsize=8, color='black')
    elif annotate_depth: 
        for i, txt in enumerate(gene_vaf['Depth_tSSCS1']):
            ax.annotate(txt, (gene_vaf['Gene'].iloc[i], gene_vaf['VAF_tSSCS1'].iloc[i]), ha='left', va='center', fontsize=8, color='black')
    fig.tight_layout()
    fig.savefig(save_path)


type = "somatic"
variant_callers = ["Mutect2", "Vardict"]

sscs1_df = pd.read_csv(f"/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/{type}_curated_SSCS1.csv")[["Sample_name_t", "Diagnosis", "Gene", "Protein_annotation", "VAF_t", "Depth_t"]]
sscs2_df = pd.read_csv(f"/groups/wyattgrp/users/amunzur/pipeline/results/variant_calling/{type}_curated_SSCS2.csv")[["Sample_name_t", "Diagnosis", "Gene", "Protein_annotation", "VAF_t", "Depth_t"]]

combined = sscs1_df.merge(sscs2_df, how = "outer", on = ["Sample_name_t", "Diagnosis", "Gene", "Protein_annotation"], indicator = True, suffixes = ["SSCS1", "SSCS2"])
kidney = combined[combined["Diagnosis"] == "Kidney"].reset_index(drop = True)
bladder = combined[combined["Diagnosis"] == "Bladder"].reset_index(drop = True)

kidney["_merge"].value_counts()
bladder["_merge"].value_counts()

# plot to compare the vafs of mutations found in both
both_kidney = kidney[kidney["_merge"] == "both"]
both_kidney["vaf_diff"] = both_kidney['VAF_tSSCS1']*100 - both_kidney['VAF_tSSCS2']*100
both_kidney["depth_diff"] = both_kidney['Depth_tSSCS1'] - both_kidney['Depth_tSSCS2']
both_kidney = both_kidney.sort_values(by = "vaf_diff").reset_index(drop = True)

both_bladder = bladder[bladder["_merge"] == "both"]
both_bladder["vaf_diff"] = both_bladder['VAF_tSSCS1']*100 - both_bladder['VAF_tSSCS2']*100
both_bladder["depth_diff"] = both_bladder['Depth_tSSCS1'] - both_bladder['Depth_tSSCS2']
both_bladder = both_bladder.sort_values(by = "vaf_diff").reset_index(drop = True)

# Call the function for VAF and depth comparison on the same ax
plot_comparison(
    df=both_bladder,
    vaf_column="vaf_diff",
    depth_column="depth_diff",
    vaf_ylabel="SSCS1 - SSCS2 VAF (%)",
    depth_ylabel="SSCS1 depth - SSCS2 depth",
    title="Comparing the VAF and depth of mutations called by both SSCS1 and SSCS2 [Bladder, somatic]",
    filename="/groups/wyattgrp/users/amunzur/pipeline/results/figures/sscs1_vs_sscs2_comparison/vaf_depth_comparison_BLADDER.png"
)

plot_comparison(
    df=both_kidney,
    vaf_column="vaf_diff",
    depth_column="depth_diff",
    vaf_ylabel="SSCS1 - SSCS2 VAF (%)",
    depth_ylabel="SSCS1 depth - SSCS2 depth",
    title="Comparing the VAF and depth of mutations called by both SSCS1 and SSCS2 [Kidney, somatic]",
    filename="/groups/wyattgrp/users/amunzur/pipeline/results/figures/sscs1_vs_sscs2_comparison/vaf_depth_comparison_KIDNEY.png"
)

kidney_sscs1_only = kidney[kidney["_merge"] == "left_only"]
plot_gene_vaf(df=kidney_sscs1_only, 
              figuresize=(5,3),
              title="SSCS1 only mutations [Kidney]", 
              save_path="/groups/wyattgrp/users/amunzur/pipeline/results/figures/sscs1_vs_sscs2_comparison/kidney_sscs1_only_mutations.png", 
              annotate_protein=True, 
              annotate_depth=False)

bladder_sscs1_only = bladder[bladder["_merge"] == "left_only"]
plot_gene_vaf(df=bladder_sscs1_only, 
              figuresize=(10,6), 
              title="SSCS1 only mutations [Bladder]", 
              save_path="/groups/wyattgrp/users/amunzur/pipeline/results/figures/sscs1_vs_sscs2_comparison/bladder_sscs1_only_mutations.png", 
              annotate_protein=False,
              annotate_depth=True)




