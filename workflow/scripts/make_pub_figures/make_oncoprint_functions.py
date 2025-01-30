def remove_frame(ax, list_to_remove):
    for i in list_to_remove:
        ax.spines[i].set_visible(False)

def enumerate_col(df, old_col, new_col):
    df[new_col] = df[old_col].map({gene: number + 1 for number, gene in enumerate(df[old_col].unique())})
    return(df)

def plot_mut_counts_per_gene(mut_counts_df, vaf_df, ax_mut_counts, mut_dict, bar_height, diagnosis = "Both", logscale = False, vaf_col = "VAF_t"):
    ax_secondary = ax_mut_counts.twiny()
    values = mut_counts_df.values.T
    bottom = np.zeros_like(mut_counts_df.index)
    
    for i, consequence in enumerate(mut_counts_df.columns):
        color = mut_dict[consequence]
        ax_mut_counts.barh(mut_counts_df.index, values[i], left=bottom, color=color, label=consequence, edgecolor="none", align = "edge", height=bar_height)
        bottom += values[i]
    
    remove_frame(ax_mut_counts, ["top", "right"])
    if diagnosis == "Kidney":
        ax_mut_counts.set_xlim((0, 150))
        ax_mut_counts.set_xticks([])
        ax_mut_counts.set_xticks([0, 50, 100, 150])
        ax_mut_counts.set_xticklabels(["0", "50", "100", "150"], fontsize = 6)
    elif diagnosis == "Both":
        ax_mut_counts.set_xlim((0, 250))
        ax_mut_counts.set_xticks([])
        ax_mut_counts.set_xticks([0, 100, 200])
        ax_mut_counts.set_xticklabels(["0", "100", "200"], fontsize = 6)
    else:
        ax_mut_counts.set_xlim((0, 100))
        ax_mut_counts.set_xticks([])
        ax_mut_counts.set_xticks([0, 50, 100])
        ax_mut_counts.set_xticklabels(["0", "50", "100"], fontsize = 5)
    plt.setp(ax_mut_counts.get_yticklabels(), visible=False)
    ax_mut_counts.set_xlabel("Mutation counts\nper gene", fontsize = 6)
    ax_mut_counts.tick_params(axis="x", direction="out", bottom=True, top = True, colors='k')
    ax_mut_counts.tick_params(axis='y', which='both', left=False, right = False)
    
    # Plot the vafs of mutations
    if logscale:
        vaf_df[vaf_col] = np.log10(vaf_df[vaf_col])
        ax_secondary.set_xlim((np.log10(0.25), np.log10(100)))
        ax_secondary.set_xticks([])
        ax_secondary.set_xticks([np.log10(0.25), np.log10(10), np.log10(100)])
        ax_secondary.set_xticklabels(["0.25", "10", "100"], fontsize = 5)
    else:
        ax_secondary.set_xlim((0, 100))
        ax_secondary.set_xticks([])
        ax_secondary.set_xticks([0, 50, 100])
        ax_secondary.set_xticklabels(["0", "50", "100"], fontsize = 5)
    vaf_df = vaf_df.reset_index()
    for i, row in vaf_df.iterrows():
        jitter = np.random.uniform(-0.08, 0.08)
        ax_secondary.scatter(row[vaf_col], row["Order on oncoprint"] +bar_height/2 + jitter, s = 0.5, color = "black", alpha = 0.5)
    
    # Indicate the median VAF per gene on the plot.
    medians_df = vaf_df.reset_index().groupby("Order on oncoprint")[vaf_col].median().reset_index()
    ax_secondary.scatter(medians_df[vaf_col], medians_df["Order on oncoprint"]+bar_height/2, linewidths = 0, marker = "<", s = 10, color = "red", zorder = 99)
    
    # max_vaf = vaf_df["VAF_n"].max()
    ax_secondary.tick_params(axis="x", direction="out", which="both", bottom=True, top = True, colors='k')
    plt.setp(ax_secondary.get_yticklabels(), visible=False)
    ax_secondary.tick_params(axis='y', which='both', left=False, right = False)
    
    ax_secondary.set_xlabel("VAF%", fontsize = 6)
    remove_frame(ax_secondary, ["right"])
    
    # Set the spine thickness
    ax_secondary.tick_params(width=0.5)
    ax_secondary.spines[["top", "right", "bottom", "left"]].set_linewidth(0.5)
    ax_mut_counts.tick_params(width=0.5)
    ax_mut_counts.spines[["top", "right", "bottom", "left"]].set_linewidth(0.5)
    return [ax_mut_counts, ax_secondary]

def plot_mut_counts_per_patient(ax_mut_counts_patient, pt_counts, bar_width):
    # 3. PLOTTING THE MUTATION COUNTS PER PATIENT, ON TOP OF THE ONCOPRINT
    ax_mut_counts_patient.bar(x = pt_counts["Samples_enumerated"], height = pt_counts["Count"], width=bar_width, capstyle='butt', color = "black", edgecolor = "none", linewidth = 0.25, zorder = 10)
    plt.setp(ax_mut_counts_patient.get_xticklabels(), visible=False)
    ax_mut_counts_patient.tick_params(axis='x', which='both', left=False)
    remove_frame(ax_mut_counts_patient, ["top", "right"])
    ax_mut_counts_patient.set_ylabel("Mutation\ncounts", rotation = 0, labelpad = 5, fontsize = 6)
    ax_mut_counts_patient.set_ylim(0, 15)
    ax_mut_counts_patient.set_yticks([0, 5, 10])
    ax_mut_counts_patient.set_yticklabels(["0", "5", "10"])
    ax_mut_counts.tick_params(axis='y', which='both', left=False, right = False)
    ax_mut_counts_patient.spines[["left", "bottom"]].set_linewidth(0.5)
    return ax_mut_counts_patient

def plot_mut_counts_per_patient_CH_and_ctDNA_combined(ax_mut_counts_patient, ctdna_pt_counts, ch_pt_counts, bar_width):
    """
    Plots the patient counts for ctDNA and CH on top of the barchart.
    """
    combined_counts = ctdna_pt_counts.merge(ch_pt_counts, how = "outer", on = ["Patient_id", "Samples_enumerated"], suffixes = ["_ctdna", "_ch"]).fillna(0)
    for i, row in combined_counts.iterrows():
        ch_count = row["Count_ch"]
        ctdna_count = row["Count_ctdna"]
        if ctdna_count>0:
            ax_mut_counts_patient.bar(x = row["Samples_enumerated"], height = ctdna_count, width=bar_width, capstyle='butt', color = "forestgreen", edgecolor = "none", linewidth = 0.25, zorder = 10)
        if ch_count>0:
            ax_mut_counts_patient.bar(x = row["Samples_enumerated"], height = ch_count, bottom = ctdna_count, width=bar_width, capstyle='butt', color = "red", edgecolor = "none", linewidth = 0.25, zorder = 10)
    
    plt.setp(ax_mut_counts_patient.get_xticklabels(), visible=False)
    ax_mut_counts_patient.tick_params(axis='x', which='both', left=False)
    remove_frame(ax_mut_counts_patient, ["top", "right"])
    ax_mut_counts_patient.set_ylabel("Mutation\ncounts", rotation = 0, labelpad = 5, fontsize = 6)
    ax_mut_counts_patient.set_ylim(0, 20)
    ax_mut_counts_patient.set_yticks([0, 10, 20])
    ax_mut_counts_patient.set_yticklabels(["0", "10", "20"], fontsize = 6)
    ax_mut_counts.tick_params(axis='y', which='both', left=False, right = False)
    ax_mut_counts_patient.spines["top"].set_linewidth(0.5)
    ax_mut_counts_patient.spines["right"].set_linewidth(0.5)
    ax_mut_counts_patient.tick_params(width=0.5)
    return ax_mut_counts_patient

def plot_ctdna_fractions(ax_tc, ctdna_fr, bar_width):
    # 4. PLOT THE CTDNA FRACTIONS
    ax_tc.bar(x = ctdna_fr["Samples_enumerated"], height = ctdna_fr["ctDNA fraction (%)"], width=bar_width, capstyle='butt', color = "black", edgecolor = "none", linewidth = 0.25, zorder = 10)
    ax_tc.set_ylabel("ctDNA %", rotation = 0, labelpad = 5)
    
    plt.setp(ax_tc.get_xticklabels(), visible=False)
    ax_tc.tick_params(axis='x', which='both', bottom=False)
    ax_tc.tick_params(axis='y', which='both', left=True)
    
    remove_frame(ax_tc, ["top", "right"])
    ax_tc.spines['bottom'].set_visible(True)
    
    ax_tc.set_ylim((0, 100))
    ax_tc.set_yticks([0, 50, 100])
    ax_tc.set_yticklabels(["0", "50", "100"])
    return(ax_tc)

def plot_chip_fractions(ax_chip_fr, chip_fr, bar_width):
    # 5. PLOT THE CHIP FRACTIONS
    ax_chip_fr.bar(x = chip_fr["Samples_enumerated"], height = chip_fr["CHIP_fr"], width=bar_width, capstyle='butt', color = "black", edgecolor = "none", linewidth = 0.25, zorder = 10)
    ax_chip_fr.set_ylabel("CH %", rotation = 0, labelpad = 5)
    plt.setp(ax_chip_fr.get_xticklabels(), visible=False)
    ax_chip_fr.tick_params(axis='x', which='both', bottom=False)
    ax_chip_fr.set_ylim((0, 100))
    ax_chip_fr.set_yticks([0, 50])
    ax_chip_fr.set_yticklabels(["0", "50"])
    ax_chip_fr.tick_params(axis='y', which='both', left=True, right = False)
    remove_frame(ax_chip_fr, ["top", "right", "bottom"])
    return ax_chip_fr

def add_legend(mut_dict, ax, age_df, fig, ct_and_ctdna_combined = False):
    # Create the first legend for mutations and position it to the left
    legend_handles = [Patch(color=color, label=effect) for effect, color in mut_dict.items()]
    legend_handles.append(Line2D([0], [0], marker='D', color='black', markersize=5, label='4 or more mutations'))
    
    mut_legend = ax.legend(handles=legend_handles, frameon=False, title = "Mutations", loc='upper left', bbox_to_anchor=(0.05, 1))
    ax.add_artist(mut_legend)  # Manually add the first legend to avoid being overwritten
    
    # Legend for diagnosis
    diagnosis_handles = [Patch(color=color, label=diagnosis) for diagnosis, color in {"mUC": "deepskyblue", "mRCC": "orangered"}.items()]
    diagnosis_legend = ax.legend(handles=diagnosis_handles, frameon=False, title = "Diagnosis", loc='upper left', bbox_to_anchor=(0.3, 1))
    ax.add_artist(diagnosis_legend)  # Manually add the first legend to avoid being overwritten
    
    # Add the colorbar for age right beside the first legend
    norm = mcolors.Normalize(vmin=age_df['Age'].min(), vmax=age_df['Age'].max())
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.Greys), ax=ax, orientation='horizontal', pad=0.15)
    cbar.set_label('Age')
    tick_values = np.linspace(age_df['Age'].min(), age_df['Age'].max(), num=5)
    cbar.set_ticks(tick_values)
    cbar.set_ticklabels([f'{int(tick)}' for tick in tick_values])
    
    # Create the second legend for biological sex and position it beside the first legend
    sex_dict = {"Male": "yellow", "Female": "slateblue"}
    sex_legend_handles = [Patch(color=color, label=sex) for sex, color in sex_dict.items()]
    sex_legend = ax.legend(handles=sex_legend_handles, frameon=False, title="Biological sex", loc='upper left', bbox_to_anchor=(0.2, 1))
    ax.add_artist(sex_legend)
    
    mutation_source_dict = {"CH": "whitesmoke", "ctDNA": "gainsboro"}
    if ct_and_ctdna_combined: 
        legend_handles = [Patch(color=color, label=effect) for effect, color in mutation_source_dict.items()]
        legend = ax.legend(handles=legend_handles, frameon=False, title = "Mutation source", loc='upper left', bbox_to_anchor=(0.7, 1))
        ax.add_artist(legend)  # Manually add the first legend to avoid being overwritten
        
        mydict = {"CH": "gray", "ctDNA": "black"}
        legend_handles = [Patch(color=color, label=effect) for effect, color in mutation_source_dict.items()]
        legend2 = ax.legend(handles=legend_handles, frameon=False, title = "", loc='upper left', bbox_to_anchor=(0.85, 0.2))
        ax.add_artist(legend)  # Manually add the first legend to avoid being overwritten
    
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax.tick_params(axis='y', which='both', left=False, labelleft=False)
    return ax

def plot_oncoprint(ax_main, muts_df, bar_height, bar_width):
    """
    If plot_mutations_as_squares is set to False, muts will be plotted as colors.
    """
    pts = muts_df["Samples_enumerated"].unique()
    genes = muts_df["Order on oncoprint"].unique()
    
    for i in genes: 
        for j in pts: 
            ax_main.bar(x = j, bottom = i, height = bar_height, width=bar_width, capstyle='butt', color = "whitesmoke", edgecolor = "none", linewidth = 0.25, zorder = 10)
    
    # 1.2 Plot the individual mutations as scatter on top of the grid
    size_scatter = 6
    for index, row in muts_df.iterrows():
        if row["Total_mutations_nonsilent_per_sample_per_gene"] == 1: # if the patient has only 1 mutation in a given gene
            ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s",  lw = 0, zorder = 9999)
        elif row["Total_mutations_nonsilent_per_sample_per_gene"] == 2:
            if row["Mutation_order"] == 1: # first mutation in the pair
                ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] - 0.12 + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s", lw = 0, zorder = 9999)
            elif row["Mutation_order"] == 2: # second mutation in the pair
                ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] + 0.12 + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s", lw = 0, zorder = 9999)
        if row["Total_mutations_nonsilent_per_sample_per_gene"] == 3:
            if row["Mutation_order"] == 1: # first mutation in the triplet
                ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] - 0.24 + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
            elif row["Mutation_order"] == 2: # second mutation in the triplet
                ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
            elif row["Mutation_order"] == 3: # third mutation in the triplet
                ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + 0.24 + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
        elif row["Total_mutations_nonsilent_per_sample_per_gene"] > 3:
            ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + bar_height/2, marker="d", lw=1, facecolor="black", edgecolors="none", s = size_scatter, zorder = 9999)
    return ax_main

def beautify_OP_ax(ax_main, gene_order_df, DF_Samples_enumerated, xlabel, fontsize = 10):
    
    # y tick positions
    gene_order_df["ytick position"] = (gene_order_df["CH position"] + gene_order_df["ctDNA position"])/2
    gene_order_df['ytick position'] = gene_order_df['ytick position'].fillna(gene_order_df['CH position'])
    gene_order_df['ytick position'] = gene_order_df['ytick position'].fillna(gene_order_df['ctDNA position'])
    
    ax_main.set_yticks(gene_order_df['ytick position'] + bar_height/2)
    ax_main.set_yticklabels([f'$\it{{{gene}}}$' for gene in gene_order_df["Gene"]], fontsize=8)
    ax_main.set_ylim((gene_order_df["ytick position"].min() - bar_height*2, gene_order_df["ytick position"].max() + bar_height*2))
    
    ax_main.tick_params(axis='y', which='both', left=False)
    ax_main.tick_params(axis='y', which='both', pad=-2) # brings the y ticklabels closer to the plot
    
    ax_main.set_xticks(DF_Samples_enumerated["Samples_enumerated"])
    ax_main.set_xticklabels(DF_Samples_enumerated["Patient_id"], rotation = 90, fontsize = 6)
    ax_main.set_xlim((DF_Samples_enumerated["Samples_enumerated"].min() - 1, DF_Samples_enumerated["Samples_enumerated"].max() + 1))
    
    ax_main.tick_params(axis='x', which='both', bottom=False)
    ax_main.tick_params(axis='x', which='both', pad=-2) # brings the y ticklabels closer to the plot
    ax_main.set_xlabel(xlabel, fontsize = 10)
    
    ax_main.set_xticklabels("") # Decided not to show the patient ids
    
    remove_frame(ax_main, ["top", "bottom", "right", "left"])
    return ax_main

def calculate_mutation_burden_per_gene(muts_df):
    """
    Calculates the total number of nonsilent mutations per gene.
    Decides the mutation order per gene.
    """
    muts_df = muts_df[["Patient_id", "Sample_name_t", "Timepoint", "Gene", "Consequence", "VAF_t", "Chrom"]]
    muts_df["Consequence"] = muts_df["Consequence"].str.lower().replace({"startloss": "missense", "frameshift deletion": "frameshift indel", "frameshift insertion": "frameshift indel", "non-frameshift deletion": "nonframeshift indel", "non-frameshift insertion": "nonframeshift indel"})
    muts_df = muts_df.replace("U2AF1\\x3bU2AF1L5", "U2AF1")
    muts_df["Total_mutations_nonsilent_per_sample_per_gene"] = muts_df.groupby(['Sample_name_t', 'Gene'])['Gene'].transform('count')
    muts_df["Mutation_order"] = muts_df["Mutation_order"] = muts_df.groupby(['Sample_name_t', 'Gene']).cumcount() + 1
    
    return(muts_df)

def assign_mutations_to_colors(muts_df, mut_dict):
    """
    Maps mut types to colors.
    """
    set1_palette = sns.color_palette("Set1", n_colors=len(muts_df["Consequence"].unique()))
    consequence_colors = dict(zip(sorted(muts_df["Consequence"].unique()), set1_palette))
    mut_dict.update((consequence, color) for consequence, color in consequence_colors.items() if consequence not in mut_dict)    
    muts_df["Mutation_color"] = muts_df["Consequence"].map(mut_dict)
    
    return(muts_df)

def calculate_ch_fraction(muts_df):
    """
    Calculates the CHIP fraction, which is just highest VAF times two.
    """
    chip_fr = muts_df.copy()[["Patient_id", "VAF_t", "Chrom"]]
    chip_fr = chip_fr.loc[chip_fr.groupby('Patient_id')["VAF_t"].idxmax()].reset_index(drop = True)
    chip_fr["CHIP_fr"] = chip_fr.apply(lambda row: row["VAF_t"] * 2 if row["Chrom"] != "chrX" else row["VAF_t"], axis=1)
    
    return(chip_fr)

def calculate_patient_counts(muts_df):
    """
    Calculates the mutation burden per patient in the muts_df. For the bar chart on top of the df.
    """
    pt_counts = muts_df["Patient_id"].value_counts().reset_index().rename(columns = {"index": "Patient_id", "Patient_id": "Count"})
    return(pt_counts)

def get_mutation_counts_per_gene(muts_df):
    """
    This is for the sideways stacked bar chart, showing the number of mutations per gene.
    """
    # Get genes with multiple mutations
    df = muts_df[["Patient_id", "Gene", "Consequence"]]
    df["Consequence"] = df["Consequence"].str.lower()
    counts_df = df.groupby(["Patient_id", "Gene"]).count().reset_index().rename(columns = {"Consequence": "Counts"})
    df = df.merge(counts_df)
    
    # Annotate gene-patient pairs with multiple alterations
    df.loc[df["Counts"] > 1, "Consequence"] = "multiple"

    # Get mut counts per gene
    mut_counts_df = df[["Gene", "Consequence"]].pivot_table(index='Gene', columns='Consequence', aggfunc=len, fill_value=0).reset_index()
    
    # Generate the VAF df
    if "Status" in muts_df.columns:
        vaf_col = "VAF_t"
    else:
        vaf_col = "VAF_n"
    
    vaf_df = muts_df[["Gene", vaf_col]]
    return(mut_counts_df, vaf_df)

def return_ctDNA_fractions(PATH_sample_information, PATH_mutation_ctfractions, ctdna_muts, ch_muts):
    """
    Returns the ctDNA fractions per patient.
    """
    ctdna_df = pd.read_csv(PATH_mutation_ctfractions)[["Patient_id", "Date_collected", "Mutation_ctDNA_fraction"]]
    ctdna_df["ctDNA fraction (%)"] = ctdna_df["Mutation_ctDNA_fraction"] * 100
    ctdna_df['Date_collected'] = pd.to_datetime(ctdna_df['Date_collected'])
    ctdna_df["Patient_id"] = ctdna_df["Patient_id"].str.replace("GU-", "")
    
    # Subset to baseline samples
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    sample_info["Date_collected"] = pd.to_datetime(sample_info["Date_collected"], format = '%Y%b%d')
    sample_info = sample_info[sample_info["Timepoint"] == "Baseline"]
    ctdna_df["Date_collected"] = pd.to_datetime(ctdna_df["Date_collected"])
    ctdna_df = ctdna_df.merge(sample_info, on = ["Patient_id", "Date_collected"], how = "inner")
    
    # Subset to patients included in the figure
    pts_list = pd.concat([ctdna_muts, ch_muts]).reset_index()["Patient_id"].unique().tolist()
    ctdna_df = ctdna_df[ctdna_df["Patient_id"].isin(pts_list)].reset_index(drop = True)
    ctdna_df = ctdna_df.sort_values(by = "ctDNA fraction (%)", ascending = False)
    
    ctdna_df = ctdna_df.reset_index(drop = True)
    
    return(ctdna_df)

def get_sex_and_age(PATH_clinical_data):
    clin_df = pd.read_csv(PATH_clinical_data)
    
    if "Age at blood draw" in clin_df.columns:
        age_col = "Age at blood draw"
    else:
        age_col = "Age at GUBB draw"
    
    if "First sample?" in clin_df.columns:
        clin_df = clin_df[clin_df["First sample?"] == True]
    
    age_df = clin_df[["Patient_id", age_col]].rename(columns = {age_col: "Age"})
    sex_df = clin_df[["Patient_id", "Sex"]].drop_duplicates()
    
    # Setting color for age and sex
    norm = mcolors.Normalize(vmin=age_df['Age'].min(), vmax=age_df['Age'].max())
    cmap = plt.cm.Greys
    age_df['age_color'] = age_df['Age'].apply(lambda x: mcolors.to_hex(cmap(norm(x))))    
    sex_df["sex_color"] = sex_df["Sex"].map({"Male": "yellow", "Female": "slateblue"})
    
    return(age_df, sex_df)

def prepare_datasets_for_plotting(main_muts_df, mut_dict, PATH_gene_categories, PATH_mutation_ctfractions, PATH_bladder_clinical_data, PATH_kidney_clinical_data):
    muts_df = calculate_mutation_burden_per_gene(muts_df)
    muts_df = assign_mutations_to_colors(muts_df, mut_dict)
    chip_fr = calculate_ch_fraction(muts_df)
    pt_counts = calculate_patient_counts(muts_df)
    muts_df = get_mutation_counts_per_gene(muts_df)
    ctdna_df = return_ctDNA_fractions(PATH_mutation_ctfractions, muts_df)
    age_df, sex_df = get_sex_and_age(PATH_clinical_data)

    
    
    # Helps labelling the y axis ticks (gene names). 
    # Can be a path or a dict.
    if isinstance(PATH_gene_categories, str):
        gene_order_df = pd.read_csv(PATH_gene_categories, sep = "\t").rename(columns = {"Panel genes": "Gene"}).sort_values(by = "Order on oncoprint")
    elif isinstance(PATH_gene_categories, dict):
        gene_order_df = pd.DataFrame.from_dict(PATH_gene_categories, orient='index', columns=['Order on oncoprint']).reset_index().rename(columns = {"index": "Gene"})
    
    gene_order_df = gene_order_df[~pd.isnull(gene_order_df["Order on oncoprint"])]
    genes_present_in_df = muts_df["Gene"].unique()
    gene_order_df = gene_order_df[gene_order_df["Gene"].isin(genes_present_in_df)]
    # gene_order_df["Order on oncoprint"] = range(len(gene_order_df), 0, -1)
    muts_df = muts_df.merge(gene_order_df, how = "left")
    muts_df = muts_df[~pd.isnull(muts_df["Order on oncoprint"])]
    
    # DATASET 3: MUTATION COUNTS FOR EACH PATIENT, ON TOP OF THE ONCOPRINT
    # pt_counts = muts_df["Patient_id"].value_counts().reset_index().rename(columns = {"index": "Patient_id", "Patient_id": "Count"})
    # # pt_counts = enumerate_col(pt_counts, "Patient_id", "Samples_enumerated")
    # pt_counts = pt_counts.merge(DF_Samples_enumerated, how = "left")
    # DF_Samples_enumerated = pt_counts[["Patient_id", "Samples_enumerated"]].drop_duplicates().reset_index(drop = True).sort_values(by = "Samples_enumerated")
    # muts_df = muts_df.merge(DF_Samples_enumerated, how = "left")
    
    muts_df_pivot = muts_df[["Patient_id", "Gene"]].value_counts().reset_index().pivot(columns = "Gene", values = 0, index = "Patient_id").fillna(0)
    # Sort by the first column
    muts_df_pivot = muts_df_pivot.sort_values(by=muts_df_pivot.columns[0], ascending=False) 
    muts_df_pivot = muts_df_pivot.reset_index()
    # gene_order = ['DNMT3A', 'TET2', 'ASXL1', 'ATM', 'PPM1D', 'TP53', 'JAK2', 'CHEK2', 'BRCA2']
    gene_order = muts_df.sort_values(by="Order on oncoprint")["Gene"].unique()[::-1].tolist()
    muts_df_pivot = muts_df_pivot.sort_values(by=gene_order, ascending=False).reset_index(drop=True)
    muts_df_pivot = enumerate_col(muts_df_pivot, "Patient_id", "Samples_enumerated")
    DF_Samples_enumerated = muts_df_pivot[["Patient_id", "Samples_enumerated"]].drop_duplicates().reset_index(drop = True).sort_values(by = "Samples_enumerated")
    muts_df = muts_df.merge(DF_Samples_enumerated, how = "left")
    
    
    # DATASET 5 - CHIP FRACTIONS PER PATIENT
    # chip_fr = muts_df.copy()[["Patient_id", vaf_colname, "Chrom"]]
    # chip_fr = chip_fr.loc[chip_fr.groupby('Patient_id')[vaf_colname].idxmax()].reset_index(drop = True)
    # chip_fr["CHIP_fr"] = chip_fr.apply(lambda row: row[vaf_colname] * 2 if row["Chrom"] != "chrX" else row[vaf_colname], axis=1)
    # chip_fr = chip_fr.merge(DF_Samples_enumerated, how = "left")
    
    # DATASET 2: FOR THE SIDEWAYS STACKED BARCHART FOR MUTATION COUNTS
    # mut_counts_df = muts_df[["Gene", "Consequence"]].pivot_table(index='Gene', columns='Consequence', aggfunc=len, fill_value=0).reset_index()
    # mut_counts_df = gene_order_df.merge(mut_counts_df, how = "left").fillna(0)
    # mut_counts_df = mut_counts_df.sort_values(by="Order on oncoprint", inplace=False)
    # mut_counts_df["Order on oncoprint"] = mut_counts_df["Order on oncoprint"] + bar_height/2
    # mut_counts_df.set_index(["Order on oncoprint"], inplace=True)
    # del mut_counts_df["Gene"]
    
    # DATASET 2.5: THE VAFS TO SHOW ON TOP OF THE SIDEWAYS STACKED BAR CHART
    # vaf_df = muts_df[["Gene", vaf_colname]]
    # vaf_df = gene_order_df.merge(vaf_df, how = "left").fillna(0)
    # vaf_df = vaf_df.sort_values(by="Order on oncoprint", inplace=False)
    # vaf_df["Order on oncoprint"] = vaf_df["Order on oncoprint"] + bar_height/2
    # vaf_df.set_index(["Order on oncoprint"], inplace=True)
    # del vaf_df["Gene"]  
    
    # DATASET 4: CTDNA FRACTIONS
    # ctdna_df = pd.read_csv(PATH_mutation_ctfractions)[["Patient_id", "Date_collected", "Mutation_ctDNA_fraction"]]
    # ctdna_df["ctDNA fraction (%)"] = ctdna_df["Mutation_ctDNA_fraction"] * 100
    # ctdna_df['Date_collected'] = pd.to_datetime(ctdna_df['Date_collected'])
    # ctdna_df["Patient_id"] = ctdna_df["Patient_id"].str.replace("GU-", "")
    # pts = DF_Samples_enumerated.merge(main_muts_df[["Patient_id", "Date_collected"]].drop_duplicates(), how = "left")
    # pts['Date_collected'] = pd.to_datetime(pts['Date_collected'], format='%Y%b%d')
    # pts_ctdna = pts.merge(ctdna_df, how = "left", on = ["Patient_id", "Date_collected"])
    
    # pts_ctdna_NO_nan = pts_ctdna[~pd.isnull(pts_ctdna["ctDNA fraction (%)"])]
    # pts_ctdna_nan = pts_ctdna[pd.isnull(pts_ctdna["ctDNA fraction (%)"])]
    
    # Make the age and sex dfs
    bladder_age = pd.read_csv(PATH_bladder_clinical_data)[["Patient_id", "Age at blood draw"]].rename(columns = {"Age at blood draw": "Age"})
    kidney_age = pd.read_csv(PATH_kidney_clinical_data)[["Patient_id", "Age at GUBB draw"]].rename(columns = {"Age at GUBB draw": "Age"})
    combined_age = pd.concat([bladder_age, kidney_age]).reset_index(drop = True)
    diagnosis = main_muts_df["Diagnosis"].unique()[0]
    if diagnosis == "Bladder": 
        age_df = bladder_age.copy()
        sex_df = pd.read_csv(PATH_bladder_clinical_data)[["Patient_id", "Sex"]].drop_duplicates()
    else: 
        age_df = kidney_age.copy()
        sex_df = pd.read_csv(PATH_kidney_clinical_data)[["Patient_id", "Sex"]]
    age_df = age_df.merge(DF_Samples_enumerated, how = "inner", on = "Patient_id")
    sex_df = sex_df.merge(DF_Samples_enumerated, how = "inner")
    # Setting color for age and sex
    norm = mcolors.Normalize(vmin=combined_age['Age'].min(), vmax=combined_age['Age'].max())
    cmap = plt.cm.Greys
    age_df['age_color'] = age_df['Age'].apply(lambda x: mcolors.to_hex(cmap(norm(x))))    
    sex_df["sex_color"] = sex_df["Sex"].map({"Male": "yellow", "Female": "deepskyblue"})    
    return [muts_df, gene_order_df, chip_fr, mut_counts_df, vaf_df, pt_counts, pts_ctdna_NO_nan, pts_ctdna_nan, DF_Samples_enumerated, mut_dict, sex_df, age_df, combined_age]

def plot_sex_and_age(ax_sex, ax_age, sex_df, age_df, bar_width, bar_height): 
    # ax_sex.set_ylim((0, 1))
    # ax_age.set_ylim((0, 1))
    ax_sex.bar(sex_df["Samples_enumerated"], bar_height, color=sex_df["sex_color"], edgecolor="none", width = bar_width)
    ax_age.bar(age_df["Samples_enumerated"], bar_height, color=age_df["age_color"], edgecolor="none", width = bar_width)
    for ax in [ax_sex, ax_age]:
        remove_frame(ax, ["top", "right", "bottom", "left"])
        ax.tick_params(axis='y', which='both', left=False, right = False, labelleft = False, labelright = False)
        ax.tick_params(axis='x', which='both', top=False, bottom = False, labelbottom = False, labeltop = False)
    return [ax_sex, ax_age]

def plot_sex_age_diagnosis(ax_sex, ax_age, ax_diagnosis, sex_df, age_df, diagnosis_df, bar_width, bar_height): 
    # ax_sex.set_ylim((0, 1))
    # ax_age.set_ylim((0, 1))
    ax_sex.bar(sex_df["Samples_enumerated"], bar_height, color=sex_df["sex_color"], edgecolor="none", width = bar_width)
    ax_age.bar(age_df["Samples_enumerated"], bar_height, color=age_df["age_color"], edgecolor="none", width = bar_width)
    ax_diagnosis.bar(diagnosis_df["Samples_enumerated"], bar_height, color=diagnosis_df["diagnosis_color"], edgecolor="none", width = bar_width)
    for ax in [ax_sex, ax_age, ax_diagnosis]:
        remove_frame(ax, ["top", "right", "bottom", "left"])
        ax.tick_params(axis='y', which='both', left=False, right = False, labelleft = False, labelright = False)
        ax.tick_params(axis='x', which='both', top=False, bottom = False, labelbottom = False, labeltop = False)
    return [ax_sex, ax_age, ax_diagnosis]

def assign_mutation_values(genes, muts_ch, muts_ctDNA, bar_height=1, omit_missing_row = True):
    """
    Assigns non-overlapping values to CH and ctDNA columns for genes based on their mutations,
    with configurable spacing between and within genes.
    
    Parameters:
    genes (list): List of gene names to check.
    muts_ch (pd.DataFrame): DataFrame containing CH mutation data with a 'Gene' column.
    muts_ctDNA (pd.DataFrame): DataFrame containing ctDNA mutation data with a 'Gene' column.
    bar_height (float): The height of the bars to adjust the buffer values.
    omit_missing_row: If a gene has no CH or ctDNA mutations omits that row entirely.
    
    Returns:
    pd.DataFrame: A DataFrame with columns 'Gene', 'CH position', and 'ctDNA position' with assigned values.
    """
    
    # Calculate buffer values based on the provided bar height
    buffer_within_gene = 0.5
    buffer_between_genes = 0.7
    
    # Initialize empty lists for CH and ctDNA values
    ch_values = []
    ctdna_values = []
    
    # Start with the highest possible value
    current_value = (len(genes) - 1) * 2 * buffer_between_genes + buffer_between_genes
    
    if omit_missing_row:
        for gene in genes:
            ch_values.append(current_value)
            current_value -= buffer_within_gene  # Apply buffer within the same gene
            if gene in ["TP53", "ERBB2"]:
                ctdna_values.append(current_value)
                current_value -= buffer_between_genes*2  # Apply buffer after assigning ctDNA value
            else:
                ctdna_values.append(current_value)
                current_value -= buffer_between_genes  # Apply buffer after assigning ctDNA value
    else:
        for gene in genes:
            # Check if the gene has a CH mutation
            has_ch = not muts_ch[muts_ch['Gene'] == gene].empty
            has_ctdna = not muts_ctDNA[muts_ctDNA['Gene'] == gene].empty
            
            if has_ch:
                ch_values.append(current_value)
                current_value -= buffer_within_gene  # Apply buffer within the same gene
            else:
                ch_values.append(None)
            
            if has_ctdna:
                ctdna_values.append(current_value)
                current_value -= buffer_between_genes  # Apply buffer after assigning ctDNA value
            else:
                ctdna_values.append(None)
    
    return pd.DataFrame({
        'Gene': genes,
        'CH position': ch_values,
        'ctDNA position': ctdna_values
    })

def plot_oncoprint_combined(ax_main, gene_order_df, muts_ctDNA, muts_ch, bar_height, bar_width):
    
    pts = pd.concat([muts_ctDNA, muts_ch])["Samples_enumerated"].unique()
    
    # Layout the grids
    for i, row in gene_order_df.iterrows():
        for j in pts:
            CH_position = row["CH position"]
            ctDNA_position = row["ctDNA position"]
            
            if CH_position is not None:
                ax_main.bar(x = j, bottom = CH_position, height = bar_height, width=bar_width, capstyle='butt', color = "whitesmoke", edgecolor = "none", linewidth = 0.25, zorder = 10)
            
            if ctDNA_position is not None:
                ax_main.bar(x = j, bottom = ctDNA_position, height = bar_height, width=bar_width, capstyle='butt', color = "gainsboro", edgecolor = "none", linewidth = 0.25, zorder = 10)
    
    # 1.2 Plot the individual mutations as scatter on top of the grid
    size_scatter = 4
    for muts_df in [muts_ctDNA, muts_ch]:
        for index, row in muts_df.iterrows():
            if row["Total_mutations_nonsilent_per_sample_per_gene"] == 1: # if the patient has only 1 mutation in a given gene
                ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s",  lw = 0, zorder = 9999)
            elif row["Total_mutations_nonsilent_per_sample_per_gene"] == 2:
                if row["Mutation_order"] == 1: # first mutation in the pair
                    ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] - 0.12 + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s", lw = 0, zorder = 9999)
                elif row["Mutation_order"] == 2: # second mutation in the pair
                    ax_main.scatter(x = row["Samples_enumerated"], y = row["Order on oncoprint"] + 0.12 + bar_height/2 , color = row["Mutation_color"], s = size_scatter, marker = "s", lw = 0, zorder = 9999)
            if row["Total_mutations_nonsilent_per_sample_per_gene"] == 3:
                if row["Mutation_order"] == 1: # first mutation in the triplet
                    ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] - 0.24 + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
                elif row["Mutation_order"] == 2: # second mutation in the triplet
                    ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
                elif row["Mutation_order"] == 3: # third mutation in the triplet
                    ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + 0.24 + bar_height/2, color=row["Mutation_color"], s=size_scatter, marker="s", lw=0, zorder=9999)
            elif row["Total_mutations_nonsilent_per_sample_per_gene"] > 3:
                ax_main.scatter(x=row["Samples_enumerated"], y=row["Order on oncoprint"] + bar_height/2, marker="d", lw=1, facecolor="black", edgecolors="none", s = size_scatter, zorder = 9999)
    return ax_main

def generate_plotting_df(muts_df, mut_dict, gene_order_df):
    if "Status" in muts_df.columns:
        mut_type = "ctDNA"
    else:
        mut_type = "CH"
    
    df = muts_df[["Patient_id", "Gene", "Consequence"]]
    df["Consequence"] = df["Consequence"].str.lower()
    counts_df = df.groupby(["Patient_id", "Gene"]).count().reset_index().rename(columns = {"Consequence": "Counts"})
    df = df.merge(counts_df)
    
    # Annotate gene-patient pairs with multiple alterations
    df.loc[df["Counts"] > 1, "Consequence"] = "multiple"
    df["Mutation_color"] = df["Consequence"].map(mut_dict)
    
    # Add the gene orders
    df = df.merge(gene_order_df, how = "left")
    if mut_type == "ctDNA":
        del df["CH position"]
    else:
        del df["ctDNA position"]
    
    # Add samples enumerated
    # df = df.merge(samples_enumerated)
    
    df.columns = ["Patient_id", "Gene", "Consequence", "Counts", "mut color", "gene_position"]    
    return(df)

def plot_oncoprint_compressed(ax_main, gene_order_df, muts_ctDNA, muts_ch, bar_height, bar_width):
    
    pts = pd.concat([muts_ctDNA, muts_ch])["Samples_enumerated"].unique()
    
    # Layout the grids
    for i, row in gene_order_df.iterrows():
        for j in pts:
            CH_position = row["CH position"]
            ctDNA_position = row["ctDNA position"]
            
            if CH_position is not None:
                ax_main.bar(x = j, bottom = CH_position, height = bar_height, width=bar_width, capstyle='butt', color = "whitesmoke", edgecolor = "none", linewidth = 0.25, zorder = 10)
            
            if ctDNA_position is not None:
                ax_main.bar(x = j, bottom = ctDNA_position, height = bar_height, width=bar_width, capstyle='butt', color = "gainsboro", edgecolor = "none", linewidth = 0.25, zorder = 10)
    
    # 1.2 Plot the mutations in genes
    size_scatter = 4
    for muts_df in [muts_ctDNA, muts_ch]:
        for index, row in muts_df.iterrows():
            ax_main.bar(x = row["Samples_enumerated"], bottom = row["gene_position"], height = bar_height, width=bar_width, capstyle='butt', color = row["mut color"], edgecolor = "none", linewidth = 0.25, zorder = 10)
    return ax_main

def get_diagnosis(df, PATH_sample_information):
    """
    Returns the diagnosis of patients
    """
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    sample_info = sample_info[["Patient_id", "Diagnosis"]]
    
    # Subset to patients in the df
    sample_info = sample_info[sample_info["Patient_id"].isin(df["Patient_id"])]
    
    return(sample_info)

def plot_highest_VAF_per_patient(ax_max_vaf, muts_df, samples_enumerated, bar_width):
    if "Status" in muts_df.columns:
        vaf_col = "VAF_t"
    else:
        vaf_col = "VAF_n"
    
    max_vaf_df = muts_df.groupby("Patient_id")[vaf_col].max().reset_index()
    max_vaf_df = max_vaf_df.merge(samples_enumerated)
    ax_max_vaf.bar(x = max_vaf_df["Samples_enumerated"], height = max_vaf_df[vaf_col], width=bar_width, capstyle='butt', color = "black", edgecolor = "none", linewidth = 0.25, zorder = 10)
    
    # AES
    remove_frame(ax_max_vaf, ["top", "right"])
    ax_max_vaf.spines[["left", "bottom"]].set_linewidth(0.5)
    ax_max_vaf.tick_params(width=0.5)
    ax_max_vaf.tick_params(axis='x', bottom=False, labelbottom=False)
    if vaf_col == "VAF_t":
        ax_max_vaf.set_ylim((0, 100))
        ax_max_vaf.set_yticks([0, 50, 100])
        ax_max_vaf.set_yticklabels(["0", "50", "100"], fontsize = 5)
        ax_max_vaf.set_ylabel("Max ctDNA\nVAF%", rotation = 0, fontsize = 6)
    else:
        ax_max_vaf.set_ylim((0, 50))
        ax_max_vaf.set_yticks([0, 25, 50])
        ax_max_vaf.set_yticklabels(["0", "25", "50"], fontsize = 5)
        ax_max_vaf.set_ylabel("Max CH\nVAF%", rotation = 0, fontsize = 6)
    return(ax_max_vaf)

def add_gene_group_annotations(ax_main, bar_height): 
    """
    Adds vertical lines helping group genes
    """
    # Sorry this part had to be hard coded!
    
    ch_start = 21.4
    ch_end = 31.5+bar_height
    ch_text_pos_center = (ch_end-ch_start)/2+ch_start
    
    muc_start = 9.9
    muc_end = 20+bar_height
    muc_text_pos_center = (muc_end-muc_start)/2+muc_start
    
    mrcc_start = 3.2
    mrcc_end = 8.5+bar_height
    mrcc_text_pos_center = (mrcc_end-mrcc_start)/2+mrcc_start
    
    ax_main.vlines(x=-0.2, ymin=ch_start, ymax=ch_end, color='black', linewidth=0.5) # CH genes
    ax_main.vlines(x=-0.2, ymin=muc_start, ymax=muc_end, color='black', linewidth=0.5) # Bladder genes
    ax_main.vlines(x=-0.2, ymin=mrcc_start, ymax=mrcc_end, color='black', linewidth=0.5) # Kidney genes
    
    ax_main.text(-0.1, ch_text_pos_center, 'CH genes', fontsize=6, color='black', ha='center', va='center', rotation = 90, rotation_mode='anchor', transform=ax_main.transAxes)
    ax_main.text(-0.1, muc_text_pos_center, 'mUC genes', fontsize=6, color='black', ha='center', va='center', rotation = 90, rotation_mode='anchor', transform=ax_main.transAxes)
    ax_main.text(-0.1, mrcc_text_pos_center, 'mRCC genes', fontsize=6, color='black', ha='center', va='center', rotation = 90, rotation_mode='anchor', transform=ax_main.transAxes)
    
    return(ax_main)
