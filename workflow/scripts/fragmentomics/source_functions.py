# from itertools import chain
# import ast
# import random
# from scipy.stats import gaussian_kde


def get_fragments_with_deletion(bam_file, chromosome, position):
    """
    Returns the fragment lengths of the reads that have a deletion event at a specific genomic position.
    """
    read_info_with_deletion = []  # List to store read information with deletion
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate over reads in the BAM file
        for read in bam.fetch(chromosome, position - 1, position):
            ref_pos = read.reference_start  # Initialize reference position
            fragment_id = read.query_name.split('/')[0]  # Extract fragment ID from read name
            fragment_length = None  # Placeholder for fragment length
            # Iterate over the CIGAR tuples to find deletions
            for op, length in read.cigartuples:
                if op == 2:  # Deletion
                    # Check if the specified position falls within this deletion
                    if ref_pos <= position - 1 < ref_pos + length:
                        # Get the fragment length if not already obtained
                        if fragment_length is None:
                            if read.is_paired:
                                fragment_length = abs(read.template_length)
                        # Add read information to the list
                        read_info_with_deletion.append({
                            'read_id': read.query_name,
                            'fragment_id': fragment_id,
                            'fragment_length': fragment_length
                        })
                        break  # Stop checking further CIGAR tuples for this read
                # Update the reference position
                if op in {0, 2, 3, 7, 8}:  # M, D, N, =, X
                    ref_pos += length
                elif op == 1 or op == 4:  # I or S
                    continue
    fragment_lengths = pd.DataFrame(read_info_with_deletion).drop_duplicates("fragment_id")["fragment_length"].tolist()
    return(fragment_lengths)

def get_fragments_with_snv(bam_file, chromosome, position, ref_base, alt_base):
    """
    Returns the fragment lengths of the reads that have an SNV event (either ref_base or alt_base)
    at a specific genomic position.
    """
    read_info_with_snv = []  # List to store read information with SNV
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate over reads in the BAM file
        for read in bam.fetch(chromosome, position - 1, position):
            ref_pos = read.reference_start  # Initialize reference position
            read_pos = 0  # Initialize read position
            fragment_id = read.query_name.split('/')[0]  # Extract fragment ID from read name
            fragment_length = None  # Placeholder for fragment length
            # Iterate over the CIGAR tuples to find matches and mismatches
            for op, length in read.cigartuples:
                if op == 0 or op == 7 or op == 8:  # Match or Mismatch
                    for i in range(length):
                        if ref_pos == position - 1:
                            read_base = read.query_sequence[read_pos]
                            if read_base == alt_base:
                                # Get the fragment length if not already obtained
                                if fragment_length is None:
                                    if read.is_paired:
                                        fragment_length = abs(read.template_length)
                                # Add read information to the list
                                read_info_with_snv.append({
                                    'read_id': read.query_name,
                                    'fragment_id': fragment_id,
                                    'fragment_length': fragment_length
                                })
                                break  # Stop checking further CIGAR tuples for this read
                        ref_pos += 1
                        read_pos += 1
                    if ref_pos > position - 1:
                        break  # If we've passed the position, no need to check further
                elif op == 2 or op == 3:  # Deletion or N (skip region in reference)
                    ref_pos += length
                elif op == 1:  # Insertion
                    read_pos += length
                elif op == 4 or op == 5:  # Soft clipping or hard clipping
                    read_pos += length
    fragment_lengths = pd.DataFrame(read_info_with_snv).drop_duplicates("fragment_id")["fragment_length"].tolist()
    return fragment_lengths

def get_fragments_with_insertion(bam_file, chromosome, position):
    """
    Returns the fragment lengths of the reads that have an insertion event at a specific genomic position.
    """
    read_info_with_insertion = []  # List to store read information with insertion
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate over reads in the BAM file
        for read in bam.fetch(chromosome, position - 1, position):
            ref_pos = read.reference_start  # Initialize reference position
            read_pos = 0  # Initialize read position
            fragment_id = read.query_name.split('/')[0]  # Extract fragment ID from read name
            fragment_length = None  # Placeholder for fragment length
            # Iterate over the CIGAR tuples to find insertions
            for op, length in read.cigartuples:
                if op == 1:  # Insertion
                    # Get the fragment length if not already obtained
                    if fragment_length is None:
                        if read.is_paired:
                            fragment_length = abs(read.template_length)
                    # Add read information to the list
                    read_info_with_insertion.append({
                        'read_id': read.query_name,
                        'fragment_id': fragment_id,
                        'fragment_length': fragment_length
                    })
                    break  # Stop checking further CIGAR tuples for this read
                # Update the reference position and read position
                if op in {0, 7, 8}:  # M, =, X
                    ref_pos += length
                    read_pos += length
                elif op == 2 or op == 3:  # D, N
                    ref_pos += length
                elif op == 1:  # I
                    read_pos += length
                elif op == 4 or op == 5:  # S, H
                    read_pos += length
    fragment_lengths = pd.DataFrame(read_info_with_insertion).drop_duplicates("fragment_id")["fragment_length"].tolist()
    return fragment_lengths

def get_fragments_without_variant(bam_file, chromosome, position, variant_type, ref_base=None, alt_base=None):
    """
    Returns the fragment lengths of the reads that do not have a variant (either SNV, deletion, or insertion)
    at a specific genomic position.
    """
    read_info_without_variant = []  # List to store read information without variant
    wild_type_count = 0  # Counter
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate over reads in the BAM file
        for read in bam.fetch(chromosome, position - 1, position):
            ref_pos = read.reference_start  # Initialize reference position
            read_pos = 0  # Initialize read position
            fragment_id = read.query_name.split('/')[0]  # Extract fragment ID from read name
            fragment_length = None  # Placeholder for fragment length
            has_variant = False  # Flag to check if read has variant
            # Iterate over the CIGAR tuples to find matches and mismatches
            for op, length in read.cigartuples:
                if op == 0 or op == 7 or op == 8:  # Match or Mismatch
                    for i in range(length):
                        if ref_pos == position - 1:
                            read_base = read.query_sequence[read_pos]
                            if variant_type.lower() == "snv":
                                if read_base == alt_base or read_base != ref_base:
                                    has_variant = True
                                    break
                                break
                        ref_pos += 1
                        read_pos += 1
                    if ref_pos > position - 1:
                        break  # If we've passed the position, no need to check further
                elif op == 2 or op == 3:  # Deletion or N (skip region in reference)
                    if ref_pos <= position - 1 < ref_pos + length:
                        if variant_type.lower() == "deletion":
                            has_variant = True
                            break
                    ref_pos += length
                elif op == 1:  # Insertion
                    # Check if the insertion spans the specified position
                    if ref_pos <= position - 1 < ref_pos + length:
                        if variant_type.lower() == "insertion":
                            has_variant = True
                            break
                    read_pos += length  # Move to the next position in the read
                elif op == 4 or op == 5:  # Soft clipping or hard clipping
                    read_pos += length
            if not has_variant:
                wild_type_count += 1  
                # Get the fragment length if not already obtained
                if fragment_length is None:
                    if read.is_paired:
                        fragment_length = abs(read.template_length)
                # Add read information to the list
                read_info_without_variant.append({
                    'read_id': read.query_name,
                    'fragment_id': fragment_id,
                    'fragment_length': fragment_length
                })
    fragment_lengths = pd.DataFrame(read_info_without_variant).drop_duplicates("fragment_id")["fragment_length"].tolist()
    return fragment_lengths

def make_density_plot(values, title, color, ax, bw_adjust=0.05, annotate_median = True, alpha = 0.75):
    """
    Given a list of values make histogram.
    """
    sns.kdeplot(values, ax=ax, color=color, fill=True, alpha=alpha, bw_adjust=bw_adjust, linewidth = 3)
    if annotate_median:
        ax.text(median_value+10, ax.get_ylim()[1] * 0.9, f'Median: {int(median_value)}', color='black', ha='left')
        median_value = np.median(values)
        ax.axvline(median_value, color='black', linestyle='dashed', linewidth=1)
    ax.set_xlabel('Fragment length')
    ax.set_ylabel('Frequency')
    ax.set_title(title)
    ax.set_xlim((0, 400))
    ax.set_xticks([0, 100, 200, 300, 400])
    ax.set_xticklabels(["0", "100", "200", "300", "400"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return(ax)

def make_histogram(values, title, color, color_kde, ax, annotate_median = True, alpha = 0.75, bw_adjust= 0.05, show_hist = True):
    """
    Given a list of values make histogram.
    """
    if show_hist:
        ax.hist(values, bins=250, alpha=0.75, color = color)
        kde_ax = ax.twinx()  # Create a twin axes that shares the same x-axis
        ax.spines["top"].set_visible(False)
        kde_ax.spines["top"].set_visible(False)
    else:
        kde_ax = ax
        kde_ax.spines[["top", "right"]].set_visible(False)
    sns.kdeplot(values, ax=kde_ax, color=color_kde, fill=False, bw_adjust=bw_adjust, linewidth=1, zorder=10)
    kde_ax.set_ylabel('Density')  # Label the secondary y-axis
    if annotate_median:
        median_value = np.median(values)
        ax.axvline(median_value, color='black', linestyle='dashed', linewidth=1)
        ax.text(median_value+10, ax.get_ylim()[1] * 0.9, f'Median: {int(median_value)}', color='black', ha='left')
    ax.set_xlabel('Fragment length')
    ax.set_ylabel('Frequency')
    ax.set_title(title)
    ax.set_xlim((0, 400))
    ax.set_xticks([0, 100, 200, 300, 400])
    ax.set_xticklabels(["0", "100", "200", "300", "400"])
    return(ax)

def MWU_mutated_to_WT_fragments(path_to_df, mutated_col_name, WT_col_name): 
    """
    Runs MWU and corrects for multiple testing when comparing the fragments carrying mutated fragments to WT fragments.
    Returns an array of corrected p values. 
    """
    df = pd.read_csv(path_to_df)
    p_values = []   
    # For each ctDNA mutation, run MWU to test if the medians are different
    for i, row in df.iterrows():
        mutant_fr = ast.literal_eval(row[mutated_col_name])
        WT_fr = ast.literal_eval(row[WT_col_name])
        statistic, p_value = mannwhitneyu(mutant_fr, WT_fr)
        p_values.append(p_value)
    
    # Convert p-values to a numpy array
    p_values = np.array(p_values)   
    _, p_values_corrected_bh, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh') # correct for multiple testing
    
    # Add original and corrected p-values to the DataFrame
    df['p_value'] = p_values
    df['p_value_corrected_bh'] = p_values_corrected_bh
    
    # Determine the number of significant p-values
    significant_pvals_count = np.sum(p_values_corrected_bh < 0.05)
    print(f"Number of significant p-values after correction: {significant_pvals_count}") 
    return(df)  

def violin_box(data):
    df_plot = pd.DataFrame(data)
    
    # Map groups to numerical labels
    group_mapping = {group: i for i, group in enumerate(df_plot['Group'].unique())}
    df_plot['Group_num'] = df_plot['Group'].map(group_mapping)
    
    # Create the plot using the ax method
    fig, ax = plt.subplots(figsize=(10, 6)) 
    groups = df_plot['Group'].unique()
    colors = ['blue', 'green', 'red', 'orange']
    violin_width = 0.8  # Adjust the width of violins as needed
    
    violins = ax.violinplot(
        [df_plot[df_plot['Group_num'] == group_mapping[group]]['Fragment Length'] for group in groups], 
        positions=[group_mapping[group] for group in groups],
        widths=violin_width,
        showmedians=True
    )
    
    # Customize the violins' appearance
    for pc in violins['bodies']:
        pc.set_facecolor('grey')
        pc.set_alpha(1.0)  # Normal opacity
    
    # Annotate medians on the violins
    for group, color in zip(groups, colors):
        median_val = np.median(df_plot[df_plot['Group'] == group]['Fragment Length'])
        ax.text(group_mapping[group], median_val, f'{median_val:.2f}', color='black', ha='center', va='bottom')
    
    # Adding labels and title
    ax.set_xlabel('Group')
    ax.set_ylabel('Fragment Length')
    ax.set_title('Fragment Lengths by Group')
    ax.set_xticks(list(group_mapping.values()))
    ax.set_xticklabels(list(group_mapping.keys()))
    ax.legend()
    
    fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/fragmentomics/violin.png")

def create_box_plot_with_jitter_grouped(data, save_path):
    """
    # Function to create box plots with jittered points and median values
    """
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(6, 3))
    
    # Define fragment types and diagnoses
    fragment_types = data['Fragment Type'].unique()
    diagnoses = data['Diagnosis'].unique()
    
    # Colors for the boxes and jitter points
    colors = {"Bladder": (143/255, 215/255, 239/255, 255/255), "Kidney": (239/255, 169/255, 143/255, 255/255)}
    jitter_colors = {"Bladder": "deepskyblue", "Kidney": "orangered"}
    
    # Create box plots
    box_data = []
    positions = []
    color_list = []
    
    for i, fragment_type in enumerate(fragment_types):
        for j, diagnosis in enumerate(diagnoses):
            subset = data[(data['Fragment Type'] == fragment_type) & (data['Diagnosis'] == diagnosis)]
            box_data.append(subset['Fragment Count'].values)
            positions.append(i + (j * 0.4) - 0.125)  # Adjust positions for dodge effect
            color_list.append(colors[diagnosis])
    
    boxprops = dict(linestyle='-', linewidth=0.5, color='black')
    whiskerprops = dict(color='black', linewidth=0.5)
    capprops = dict(color='black', linewidth=0.5)
    medianprops = dict(color='black', linewidth=0.5)
    
    bp = ax.boxplot(box_data, positions=positions, widths=0.4, patch_artist=True, boxprops=boxprops, whiskerprops=whiskerprops, capprops=capprops, medianprops=medianprops, showfliers=False)
    
    for patch, color in zip(bp['boxes'], color_list):
        patch.set_facecolor(color)
    
    # Add jittered strip plots
    for i, fragment_type in enumerate(fragment_types):
        for j, diagnosis in enumerate(diagnoses):
            subset = data[(data['Fragment Type'] == fragment_type) & (data['Diagnosis'] == diagnosis)]
            jitter_x = np.random.normal(loc=i + (j * 0.4) - 0.125, scale=0.02, size=subset.shape[0])
            ax.scatter(jitter_x, subset['Fragment Count'], alpha=0.7, color=jitter_colors[diagnosis], s=10, zorder = 10)
    
    # Add median values
    for i, fragment_type in enumerate(fragment_types):
        for j, diagnosis in enumerate(diagnoses):
            median_val = data[(data['Fragment Type'] == fragment_type) & (data['Diagnosis'] == diagnosis)]['Fragment Count'].median()
            pos = i + (j * 0.4) - 0.125
            ax.text(pos, median_val, f"{median_val:.1f}", ha='center', va='bottom', color='black', zorder = 100, fontsize = 7)
    
    # Set plot labels and title
    ax.set_xticks(range(len(fragment_types)))
    ax.set_xticklabels(fragment_types)
    ax.set_xlabel("Fragment Type")
    ax.set_ylabel("Fragment Count")
    ax.set_title("")
    
    # Remove legend
    ax.legend().set_visible(False)
    
    # Adjust layout and save plot
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)

def create_box_plot_with_jitter(data, diagnosis, ax):
    """
    Function to create box plots with jittered points and median values for a specific diagnosis. Kidney and bladder not together.
    """
    
    # Filter data for the specified diagnosis
    filtered_data = data[data['Diagnosis'] == diagnosis]
    
    # Define fragment types for the filtered data
    fragment_types = filtered_data['Fragment Type'].unique()
    
    # Colors for the boxes and jitter points
    colors = {"Bladder": (143/255, 215/255, 239/255, 255/255), "Kidney": (239/255, 169/255, 143/255, 255/255)}
    jitter_color = {"Bladder": "deepskyblue", "Kidney": "orangered"}
    
    # Create box plots
    box_data = []
    positions = []
    color_list = []
    
    for i, fragment_type in enumerate(fragment_types):
        subset = filtered_data[filtered_data['Fragment Type'] == fragment_type]
        box_data.append(subset['Fragment Count'].values)
        positions.append(i)
        color_list.append(colors[diagnosis])
    
    bp = ax.boxplot(box_data, positions=positions, widths=0.6, patch_artist=True, boxprops=dict(linestyle='-', linewidth=0.5, color='black'), whiskerprops=dict(color='black', linewidth=0.5), capprops=dict(color='black', linewidth=0.5), medianprops=dict(color='black', linewidth=0.5), showfliers=False)
    
    for patch, color in zip(bp['boxes'], color_list):
        patch.set_facecolor(color)
    
    # Add jittered strip plots
    for i, fragment_type in enumerate(fragment_types):
        subset = filtered_data[filtered_data['Fragment Type'] == fragment_type]
        jitter_x = np.random.normal(loc=i, scale=0.1, size=subset.shape[0])
        ax.scatter(jitter_x, subset['Fragment Count'], alpha=0.7, color=jitter_color[diagnosis], s=5, zorder=10)
        # # Update x tick labels to include the number of dots (n) in the jitter
        # ax.set_xticks(range(len(fragment_types)))
        # ax.set_xticklabels([f"{fragment_type}\nn={subset.shape[0]}" for fragment_type in fragment_types], fontsize=8, rotation=0, ha='center')
    
    # Add median values
    for i, fragment_type in enumerate(fragment_types):
        median_val = filtered_data[filtered_data['Fragment Type'] == fragment_type]['Fragment Count'].median()
        formatted_median = f"{int(median_val)}" if median_val.is_integer() else f"{median_val:.1f}"
        ax.text(i, median_val, formatted_median, ha='center', va='bottom', color='black', zorder=100, fontsize=7)
    
    # Mann-Whitney U test and p-value display
    # ctDNA
    group1 = filtered_data[filtered_data['Fragment Type'] == "ctDNA WT"]['Fragment Count']
    group2 = filtered_data[filtered_data['Fragment Type'] == "ctDNA mutant"]['Fragment Count']
    statistic, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
    ax.text(0.25, 0.95, f"p-value:\n{p_value:.2e}", ha='center', va='center', transform=ax.transAxes, fontsize=8)
    # CH
    group1 = filtered_data[filtered_data['Fragment Type'] == "CH WT"]['Fragment Count']
    group2 = filtered_data[filtered_data['Fragment Type'] == "CH mutant"]['Fragment Count']
    statistic, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
    ax.text(0.75, 0.95, f"p-value:\n{p_value:.2f}", ha='center', va='center', transform=ax.transAxes, fontsize=8)
    
    # Set plot labels and title
    ax.set_xticks(range(len(fragment_types)))
    ax.set_xticklabels(fragment_types, fontsize = 8)
    ax.set_xlabel("Fragment Type")
    ax.set_ylabel("Fragment Count")
    ax.set_ylim((95, 355))
    ax.set_title(diagnosis)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    # Remove legend
    ax.legend().set_visible(False)
    return(ax)

def scatter_n_muts_and_fragment_length(melted_df, muts_df, type_of_muts, ax):
    """
    Is there a relationship between the number of mutations a person may have and the median mutated fragment length?
    melted_df: Df that tells the median fragment length per patient
    muts_df: mutation calls. could be CH or ctDNA muts.
    """
    if type_of_muts == "CH": 
        filtered = melted_df[melted_df["Fragment Type"] == "CH mutant"]
    else: 
        filtered = melted_df[melted_df["Fragment Type"] == "ctDNA mutant"]
    # combined
    n_muts = muts_df["Patient_id"].value_counts().reset_index()
    n_muts.columns = ["Patient_id", "n_muts"]
    combined = filtered.merge(n_muts, how = "inner", on = "Patient_id")
    
    # Plotting
    ax = ax.scatter(combined['n_muts'], combined['Fragment Count'], 
                     alpha=0.7, c=combined['Diagnosis'].apply(lambda x: {'Kidney': 'orangered', 'Bladder': 'deepskyblue'}[x]), label=combined['Diagnosis'])
    return(ax)

def generate_histograms(df, diagnosis, title, fill_color, trace_color, save_path):
    if "CH fragments" in df.columns: 
        col_to_plot = "CH fragments"
    elif "ctDNA fragments" in df.columns:
        col_to_plot = "ctDNA fragments"
    
    if diagnosis == "both": 
        ch_del = list(chain.from_iterable(df[df["Type"] == "Deletion"][col_to_plot]))
        ch_ins = list(chain.from_iterable(df[df["Type"] == "Insertion"][col_to_plot]))
        ch_snv = list(chain.from_iterable(df[df["Type"] == "SNV"][col_to_plot]))
    else:
        ch_del = list(chain.from_iterable(df[(df["Type"] == "Deletion") & (df["Diagnosis"] == diagnosis)][col_to_plot]))
        ch_ins = list(chain.from_iterable(df[(df["Type"] == "Insertion") & (df["Diagnosis"] == diagnosis)][col_to_plot]))
        ch_snv = list(chain.from_iterable(df[(df["Type"] == "SNV") & (df["Diagnosis"] == diagnosis)][col_to_plot]))
    all_ch = ch_del + ch_ins + ch_snv
    
    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax1 = make_histogram(ch_del, "Deletions", color=fill_color, color_kde=trace_color, ax=ax1)
    ax2 = fig.add_subplot(gs[0, 1])
    ax2 = make_histogram(ch_ins, "Insertions", color=fill_color, color_kde=trace_color, ax=ax2)
    ax3 = fig.add_subplot(gs[1, 0])
    ax3 = make_histogram(ch_snv, "SNVs", color=fill_color, color_kde=trace_color, ax=ax3)
    ax4 = fig.add_subplot(gs[1, 1])
    ax4 = make_histogram(all_ch, "All events", color=fill_color, color_kde=trace_color, ax=ax4)
    
    fig.suptitle(title)
    gs.tight_layout(fig)
    fig.savefig(save_path)
    plt.close(fig)

# Scatter visualizations 
# DIR_output = "/groups/wyattgrp/users/amunzur/pipeline/results/fragmentomics"

# p_value_list = []
# for i, row in df_ctDNA.iterrows():
#     ctDNA_fragments = ast.literal_eval(row["ctDNA fragments"])
#     wt_fragments = ast.literal_eval(row["WT fragments"])
#     statistic, p_value = mannwhitneyu(ctDNA_fragments, wt_fragments)
#     p_value_list.append(p_value)

# _, corrected_p_values_ctDNA, _, _ = multipletests(p_value_list, method='fdr_bh')

# len(corrected_p_values_ctDNA[corrected_p_values_ctDNA<0.05])

def compare_mutated_to_WT(df_path, genes = None):
    """
    Given a muts_df, with CH and WT fragments added as two additional columns
    """
    df = pd.read_csv(df_path)
    
    if genes is not None:
        if isinstance(genes, str): 
            df = df[df["Gene"] == genes]
        elif isinstance(genes, list): 
            df = df[df["Gene"].isin(genes)]
    
    if "ctDNA fragments" in df.columns:
        col_to_use = "ctDNA fragments"
    else: 
        col_to_use = "CH fragments"
    
    # Go through each mut and compare the WT to mutated
    p_value_list = []
    for i, row in df.iterrows():
        mutated_fragments = ast.literal_eval(row[col_to_use])
        wt_fragments = ast.literal_eval(row["WT fragments"])
        statistic, p_value = mannwhitneyu(mutated_fragments, wt_fragments)
        p_value_list.append(p_value)
    
    _, corrected_p_values, _, _ = multipletests(p_value_list, method='fdr_bh')
    
    # add to the df
    df["MWU p value"] = p_value_list
    df["corrected MWU p value"] = corrected_p_values
    
    n_sig = len(corrected_p_values[corrected_p_values<0.05])
    
    print(f"There are {n_sig} significant p values.")
    return(df)

def plot_mutated_median_vs_wt_median(df, ax, title, set_y_label = True, add_legend = True, add_pie = True):
    if "ctDNA fragments" in df.columns:
        col_to_use = "ctDNA fragments"
    else: 
        col_to_use = "CH fragments"
    
    color_dict = {True: "cornflowerblue", False: "black"}
    df["sig"] = df["corrected MWU p value"] < 0.05
    df["color"] = df["sig"].map(color_dict)
    
    for i, row in df.iterrows():
        mutated_fragments = ast.literal_eval(row[col_to_use])
        wt_fragments = ast.literal_eval(row["WT fragments"])
        
        mutated_median = np.median(mutated_fragments)
        wt_median = np.median(wt_fragments)
        color = row["color"]
        
        ax.scatter(mutated_median, wt_median, alpha=0.5, color=color, s=1, edgecolor = None)
    
    # set ax lims
    ax_max = max(max(ax.get_xlim()), max(ax.get_ylim()))
    ax_min = min(max(ax.get_xlim()), min(ax.get_ylim())) 
    ax.set_xlim((ax_min-20, ax_max))
    ax.set_ylim((ax_min-10, ax_max))
    
    # Aes
    if set_y_label: 
        ax.set_ylabel('Median wildtype fragment length', fontsize = 9)
    ax.set_title(title, loc="center")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Median mutated fragment length")
    
    if add_legend:
        legend_colors = ["cornflowerblue", "black"]
        legend_labels = ["p < 0.05", "p > 0.05"]
        legend_handles = [plt.Line2D([0], [0], color=color, marker='o', label=label, linestyle='None', markersize=5) for color, label in zip(legend_colors, legend_labels)]
        ax.legend(handles=legend_handles, loc="upper right", frameon=False, handletextpad=0.2)
    
    if add_pie:
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        ax_inset = inset_axes(ax, width="30%", height="30%", loc='upper right')
        n_sig_p = df[df["sig"] == True].shape[0]
        n_nonsig_p = df[df["sig"] == False].shape[0]
        
        # Data for pie chart
        # labels = ['', '']
        sizes = [n_sig_p, n_nonsig_p]
        colors = ['cornflowerblue', 'black']
        explode = (0.1, 0)  # explode 1st slice
        
        # Plot pie chart on inset axis
        ax_inset.pie(sizes, explode=explode, colors=colors, shadow=False, startangle=140, autopct='%1.1f%%')
        ax_inset.axis('equal')
        ax_inset.spines[["top", "right", "bottom", "left"]].set_visible(False)
        return(ax, ax_inset)
    else: 
        return(ax)

def plot_gs(path_save, fig_title = None):
    fig = plt.figure(figsize=(6, 2.5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
    ax1 = fig.add_subplot(gs[0])
    ax1 = plot_mutated_median_vs_wt_median(df_ctDNA, ax=ax1, title = "ctDNA mutations", add_legend = False)
    ax2 = fig.add_subplot(gs[1])
    ax2 = plot_mutated_median_vs_wt_median(df_CH, ax=ax2, title = "CH mutations", set_y_label = False)
    
    fig.text(0.5, 0.02, 'Mutated fragments', ha='center', fontsize = 9)
    gs.tight_layout(fig)
    
    if fig_title is not None:
        fig.suptitle(fig_title, fontsize = 9)
    
    fig.savefig(path_save)

def plot_patient_jitter_box(patient_df, save_path):
    # Set up the figure and axes
    fig, ax = plt.subplots(figsize=(10, 8))
    # Create a box plot with jitter using Seaborn
    sns.boxplot(data=patient_df, x="Patient_id", y="Median WT ctDNA fragments", hue="Median ctDNA fragments", dodge=True, linewidth=1.5, ax=ax)
    sns.stripplot(data=patient_df, x="Patient_id", y="Median WT ctDNA fragments", hue="Median ctDNA fragments", dodge=True, size=8, alpha=0.5, ax=ax)
    # Set labels and title
    ax.set_xlabel("Patient ID")
    ax.set_ylabel("Median WT ctDNA Fragments")
    ax.set_title("Patient Fragmentomics Analysis")
    # Add legend
    ax.legend(title="Median ctDNA Fragments", loc="upper right")
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha="right")
    # Add vertical lines for median values of each group
    medians = patient_df[["Median WT ctDNA fragments", "Median ctDNA fragments", "Median WT CH fragments", "Median CH fragments"]].median()
    for median_val in medians:
        ax.axhline(y=median_val, color="r", linestyle="--", linewidth=2)
    # Save the plot
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

# def make_kde_fragmentomics(path_ch_fr, path_ctDNA_fr, ax1, ax2):
#     """
#     Shows the normal, ctDNA and CH in the same panel.
#     """
#     # path_to_df = os.path.join(DIR_fragment_counts, "CH.csv")
#     df_ch = MWU_mutated_to_WT_fragments(path_ch_fr, mutated_col_name = "CH fragments", WT_col_name = "WT fragments") # 18 sig
#     df_ch['CH fragments'] = df_ch['CH fragments'].apply(ast.literal_eval)
#     df_ch['WT fragments'] = df_ch['WT fragments'].apply(ast.literal_eval)
    
#     # ctDNA vs WT fragments
#     # path_to_df = os.path.join(DIR_fragment_counts, "ctDNA.csv")
#     df_ctDNA = MWU_mutated_to_WT_fragments(path_ctDNA_fr, mutated_col_name = "ctDNA fragments", WT_col_name = "WT fragments") # 189 sig
#     df_ctDNA['ctDNA fragments'] = df_ctDNA['ctDNA fragments'].apply(ast.literal_eval)
#     df_ctDNA['WT fragments'] = df_ctDNA['WT fragments'].apply(ast.literal_eval)
    
#     CH_mutant_fragments = list(chain.from_iterable(df_ch['CH fragments']))
#     CH_WT_fragments = list(chain.from_iterable(df_ch['WT fragments']))
#     ctDNA_mutant_fragments = list(chain.from_iterable(df_ctDNA['ctDNA fragments']))
#     ctDNA_WT_fragments = list(chain.from_iterable(df_ctDNA['WT fragments']))
    
#     # generate kde values and plot
#     kde1 = gaussian_kde(CH_mutant_fragments, bw_method=0.01)
#     kde2 = gaussian_kde(CH_WT_fragments, bw_method=0.01)
#     kde3 = gaussian_kde(ctDNA_mutant_fragments, bw_method=0.01)
#     kde4 = gaussian_kde(ctDNA_WT_fragments, bw_method=0.01)

#     x_range = np.linspace(0, 300, 1000)

#     kde1_values = kde1(x_range)
#     kde2_values = kde2(x_range)
#     kde3_values = kde3(x_range) 
#     kde4_values = kde4(x_range)


# ks_statistic, p_value = ks_2samp(kde1_values, kde2_values)
# print("KS statistic:", ks_statistic)
# print("p-value:", p_value)


#     # plotting
#     ax1.plot(x_range, kde1_values, label='Distribution 1')
#     ax1.plot(x_range, kde2_values, label='Distribution 2')

#     ax2.plot(x_range, kde3_values, label='Distribution 1')
#     ax2.plot(x_range, kde4_values, label='Distribution 2')

#     fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/pub_figures/test.png")

    
#     plt.plot(x_range, kde3_values, label='Distribution 3')
#     plt.fill_between(x_range, kde1_values, alpha=0.3)
#     plt.fill_between(x_range, kde2_values, alpha=0.3)
#     plt.fill_between(x_range, kde3_values, alpha=0.3)


# ks_2samp(kde3_values, kde4_values)



    
#     ax1 = make_histogram(CH_WT_fragments, title = "", color = "black", color_kde = "black", ax = ax1, annotate_median = False, alpha = 0.75, bw_adjust= 0.07, show_hist = False)
#     ax1 = make_histogram(CH_mutant_fragments, title = "", color = "lightcoral", color_kde = "lightcoral", ax = ax1, annotate_median = False, alpha = 0.75, bw_adjust= 0.07, show_hist = False)
#     ax2 = make_histogram(ctDNA_WT_fragments, title = "", color = "black", color_kde = "black", ax = ax2, annotate_median = False, alpha = 0.75, bw_adjust= 0.07, show_hist = False)
#     ax2 = make_histogram(ctDNA_mutant_fragments, title = "", color = "lightcoral", color_kde = "lightcoral", ax = ax2, annotate_median = False, alpha = 0.75, bw_adjust= 0.07, show_hist = False)
    
#     ax1.set_title("CH")
#     ax2.set_title("ctDNA")
#     # ax.set_xlim((100, 200))
#     # legend_colors = ["brown", "mediumblue", "black"]
#     # legend_labels = ["CH", "ctDNA", "Control"]
#     # legend_handles = [plt.Line2D([0], [0], color=color, label=label, linewidth=2) for color, label in zip(legend_colors, legend_labels)]
#     # ax.legend(handles=legend_handles, loc="upper right", frameon=False)
    
#     return(ax1, ax2)


