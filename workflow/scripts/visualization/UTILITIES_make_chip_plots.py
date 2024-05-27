def annotate_CH_status_of_cohort(PATH_sample_info, df_ch): 
    all_pts = pd.read_csv(PATH_sample_info, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
    ch_df["CH_status"] = "CH+"
    combined = all_pts.merge(ch_df[["Patient_id", "CH_status"]], how="left")
    combined["CH_status"] = combined["CH_status"].fillna("CH-")
    combined = combined[["Patient_id", "CH_status"]]
    return combined

def plot_change_in_vaf(df, figure_dir):
    """
    Plots the change in VAF going from baseline to progression.
    """
    df = df.sort_values(by = ["Gene", "VAF difference"]).reset_index(drop = True)
    fig, ax = plt.subplots(figsize = (8, 4))
    ax.bar(df.index, df['VAF difference'])
    ax.set_xlabel('Individual mutations')
    ax.set_ylabel('Progression VAF - Baseline VAF')
    ax.set_title('VAF difference in progression and baseline samples - CRPC 2022 - TheraP')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks(df.index)
    ax.set_xticklabels(df["Gene"], rotation=90)
    fig.tight_layout()
    fig.savefig(os.path.join(figure_dir, "vaf_difference.png"))

def plot_age_plots(baseline_df, PATH_sample_information, figure_dir, path_clinical, fontsize, ext = "png"):
    """
    Makes a figure with two subplots stacked on top of each other.
    Top one is a line plot showing the percentage of the cohort that is CH positive for each age bin.
    Bottom is a histogram showing the number of patients we have in each bin.
    """
    # Add mut status to the samples
    diagnosis = baseline_df["Diagnosis"].unique()[0]
    timepoint = baseline_df["Timepoint"].unique()[0]
    if "Status_n" in baseline_df.columns:
        status_column_name = "Status_n"
    elif "Status" in baseline_df.columns:
        status_column_name = "Status"
    else:
        raise ValueError("Neither 'Status_n' nor 'Status' column found in DataFrame.")
    status = baseline_df[status_column_name].unique()[0].replace("IP", "").replace("SOMATIC", "ctDNA")
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    sample_info = sample_info[(sample_info["Diagnosis"] == diagnosis) & (sample_info["Timepoint"] == timepoint)].reset_index(drop = True)
    positive_pts = baseline_df[["Patient_id", "Diagnosis"]].drop_duplicates().reset_index(drop=True).assign(Status="Positive")
    all_pts = sample_info.merge(positive_pts, how = "left")
    all_pts["Status"] = all_pts["Status"].fillna("Negative")
    all_pts['Date_collected'] = pd.to_datetime(all_pts['Date_collected'], format='%Y%b%d')
    
    # clinical data
    clin_df = pd.read_csv(path_clinical)[["Patient_id", "Date of birth"]].drop_duplicates()
    if diagnosis == "Bladder":
        clin_df['Date of birth'] = pd.to_datetime(clin_df['Date of birth'])
    else:
        clin_df['Date of birth'] = pd.to_datetime(clin_df['Date of birth'], format='%d-%b-%y', yearfirst=True)
    # Correct dates where the parsed year is in the future (greater than the current year)
    current_year = pd.Timestamp.now().year
    corrected_dates = clin_df['Date of birth'].apply(lambda x: x - pd.DateOffset(years=100) if x.year > current_year else x)
    clin_df['Date of birth'] = corrected_dates
    
    all_pts = all_pts.merge(clin_df, how = "left", on = "Patient_id")
    all_pts['Age at blood draw'] = (all_pts['Date_collected'] - all_pts['Date of birth']).astype('<m8[Y]')
    all_pts = all_pts[~pd.isnull(all_pts["Age at blood draw"])]
    
    # Bin patients' age
    if diagnosis == "Bladder":
        bins = [20, 40, 50, 60, 70, 80, 90]
        labels = ['20-40', '40-50', '50-60', '60-70', '70-80', '80-90']
    else: 
        bins = [20, 40, 50, 60, 70, 80, 90]
        labels = ['20-40', '40-50', '50-60', '60-70', '70-80', '80-90']
    all_pts['Age_bin'] = pd.cut(all_pts['Age at blood draw'], bins=bins, labels=labels, right=False)
    
    # Calculate percentage of positive patients in each age group
    age_groups = all_pts.groupby('Age_bin')['Status'].value_counts(normalize=True).unstack().fillna(0)
    age_groups['Positive_percentage'] = age_groups['Positive'] * 100
    
    sorted_labels = sorted(labels, key=lambda x: int(x.split('-')[0]))
    patient_counts = all_pts['Age_bin'].value_counts().reindex(sorted_labels).fillna(0)
    
    # Create figure and gridspec
    fig = plt.figure(figsize=(4.8, 3))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    
    # Plot dot plot in the first subplot
    ax1 = plt.subplot(gs[0])
    ax1.plot(age_groups.index, age_groups['Positive_percentage'], marker='o', color='black')
    ax1.set_xlabel('', fontsize=fontsize)
    ax1.set_ylabel(f'Percentage of\n{status}+ patients', fontsize=fontsize)
    ax1.set_xticklabels(labels, fontsize=fontsize)
    ax1.set_title(f'{diagnosis} - Percentage of {status}+ patients in each age group')
    ax1.grid(True)
    ax1.tick_params(axis='x')
    ax1.set_yticks([0, 20, 40, 60, 80, 100])
    ax1.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=fontsize)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    
    # Plot histogram in the second subplot
    ax2 = plt.subplot(gs[1], sharex=ax1)
    patient_counts = all_pts['Age_bin'].value_counts().reindex(labels).fillna(0)
    ax2.bar(labels, patient_counts, color="black", align='center')
    ax2.set_xlabel('Age at baseline blood draw', fontsize=fontsize)
    ax2.set_ylabel('Number of\npatients', fontsize=fontsize)
    ax2.set_xticklabels(labels, fontsize=fontsize)
    ax2.set_yticks([0, 10, 20, 30, 40, 50])
    ax2.set_yticklabels(["0", "10", "20", "30", "40", "50"], fontsize=fontsize)
    ax2.grid(False)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    for i, count in enumerate(patient_counts):
        ax2.text(i, count, str(int(count)), ha='center', va='bottom', fontsize=fontsize, color = "black")
    
    # Save the figure
    figure_path = os.path.join(figure_dir, f'{diagnosis}_{status}_{timepoint}_age_plots.{ext}')
    fig.tight_layout()
    fig.savefig(figure_path)

def plot_age_plots_stacked(df, age_df, PATH_sample_information, figure_dir, fontsize, color_dict, ext = "png"):
    """
    SIMILAR TO THE plot_age_plots FUNCTION - ONLY DIFFERENCE IS IT PLOTS BLADDER AND KIDNEY TOGETHER.
    Makes a figure with two subplots stacked on top of each other.
    Top one is a line plot showing the percentage of the cohort that is CH positive for each age bin.
    Bottom is a histogram showing the number of patients we have in each bin.
    """
    timepoint = df["Timepoint"].unique()[0]
    if "Status_n" in df.columns:
        status_column_name = "Status_n"
    elif "Status" in df.columns:
        status_column_name = "Status"
    else:
        raise ValueError("Neither 'Status_n' nor 'Status' column found in DataFrame.")
    status = df[status_column_name].unique()[0].replace("IP", "").replace("SOMATIC", "ctDNA")
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    sample_info = sample_info[(sample_info["Timepoint"] == timepoint)].reset_index(drop = True)
    positive_pts = df[["Patient_id", "Diagnosis"]].drop_duplicates().reset_index(drop=True).assign(Status="Positive")
    all_pts = sample_info.merge(positive_pts, how = "left")
    all_pts["Status"] = all_pts["Status"].fillna("Negative")
    all_pts['Date_collected'] = pd.to_datetime(all_pts['Date_collected'], format='%Y%b%d')
    all_pts = all_pts.merge(age_df, how = "left")
    all_pts.dropna(inplace = True)
    
    # Bin patients' age
    bins = [20, 40, 50, 60, 70, 80, 90]
    labels = ['20-40', '40-50', '50-60', '60-70', '70-80', '80-90']
    all_pts['Age_bin'] = pd.cut(all_pts['Age'], bins=bins, labels=labels, right=False)
    bladder = all_pts[all_pts["Diagnosis"] == "Bladder"]
    kidney = all_pts[all_pts["Diagnosis"] == "Kidney"]
    
    # Calculate percentage of positive patients in each age group
    age_groups_bladder = bladder.groupby('Age_bin')['Status'].value_counts(normalize=True).unstack().fillna(0)
    age_groups_bladder['Positive_percentage'] = age_groups_bladder['Positive'] * 100
    age_groups_kidney = kidney.groupby('Age_bin')['Status'].value_counts(normalize=True).unstack().fillna(0)
    age_groups_kidney['Positive_percentage'] = age_groups_kidney['Positive'] * 100
    
    sorted_labels = sorted(labels, key=lambda x: int(x.split('-')[0]))
    patient_counts_bladder = bladder['Age_bin'].value_counts().reindex(sorted_labels).fillna(0)
    patient_counts_kidney = kidney['Age_bin'].value_counts().reindex(sorted_labels).fillna(0)
       
    # Create figure and gridspec
    fig = plt.figure(figsize=(3, 2.7))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    
    # Plot dot plot in the first subplot
    ax1 = plt.subplot(gs[0])
    ax1.plot(age_groups_bladder.index, age_groups_bladder['Positive_percentage'], marker='o', color=color_dict["Bladder"])
    ax1.plot(age_groups_kidney.index, age_groups_kidney['Positive_percentage'], marker='o', color=color_dict["Kidney"])        
    ax1.set_xlabel('', fontsize=fontsize)
    ax1.set_ylabel(f'Percentage of\n{status}+ patients', fontsize=fontsize)
    ax1.set_xticklabels(labels, fontsize=fontsize-2)
    # ax1.set_title(f'{diagnosis} - Percentage of {status}+ patients in each age group')
    # ax2.xaxis.grid(False)
    ax1.yaxis.grid(True, linestyle='--', linewidth=0.5)
    ax1.tick_params(axis='x')
    ax1.set_yticks([0, 20, 40, 60, 80, 100])
    ax1.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=fontsize)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    
    # Plot histogram in the second subplot
    ax2 = plt.subplot(gs[1], sharex=ax1)
    patient_counts = all_pts['Age_bin'].value_counts().reindex(labels).fillna(0)
    ax2.bar(labels, patient_counts_bladder, color=color_dict["Bladder"], label='Bladder', edgecolor = None, linewidth=0.0)
    ax2.bar(labels, patient_counts_kidney, bottom=patient_counts_bladder, color=color_dict["Kidney"], label='Kidney', edgecolor = None, linewidth=0.0)
    ax2.set_xlabel('Age at baseline blood draw', fontsize=fontsize)
    ax2.set_ylabel('Number of\npatients', fontsize=fontsize)
    ax2.set_xticklabels(labels, fontsize=fontsize-2)
    ax2.set_yticks([0, 20, 40, 60, 80, 100])
    ax2.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=fontsize)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    for i, count in enumerate(patient_counts):
        ax2.text(i, count, str(int(count)), ha='center', va='bottom', fontsize=fontsize, color = "black")
    
    # Save the figure
    figure_path = os.path.join(figure_dir, f"age_and_ch_status_plots.{ext}")
    fig.tight_layout()
    fig.savefig(figure_path)






def generate_patient_profiles(baseline_df, prog_df, figure_dir, ext):
    """
    Plot all muts found in all samples of the patient in the same panel.
    """
    # Merge baseline and progression DataFrames
    merged_df = pd.concat([base_bladder_chip, prog_bladder_chip])
    
    # Group by patient information
    grouped_by_patient = merged_df.groupby(["Patient_id", "Diagnosis"])
    # Define colors for different sample types and progression samples within a patient
    sample_type_colors = {'Baseline': '#8ac926'}
    progression_sample_colors = {}
    for (patient_id, diagnosis), group in grouped_by_patient:
        print(patient_id)
        if group.shape[0] < 10: 
            fig_width = 11
            rotation = 0
        elif group.shape[0] >= 10: 
            fig_width = 15
            rotation = 45
        fig, ax = plt.subplots(figsize=(fig_width, 4))
        progression_sample_colors[patient_id] = {date: f'C{i}' for i, date in enumerate(group[group['Timepoint'] == 'During treatment']['Date_collected'].unique())}
        unique_mutations = group[['Chrom', 'Position', 'Gene', 'Protein_annotation', 'Consequence']].drop_duplicates()
        bar_positions = np.arange(len(unique_mutations))
        bar_width = 0.30
        x_ticks = []
        # patient_unique_mutations = []
        for i, (sample_type, date_collected) in enumerate(group.groupby(['Timepoint', 'Date_collected']).groups.keys()):
            sample_group = group[(group['Timepoint'] == sample_type) & (group['Date_collected'] == date_collected)].sort_values(by=['VAF_n'])
            x_values = bar_positions + i * bar_width
            sample_group = sample_group.drop_duplicates()
            sample_group = sample_group.reset_index().set_index(['Chrom', 'Position', 'Gene', 'Protein_annotation', 'Consequence']).reindex(unique_mutations)
            y_values = sample_group['VAF_n'].values
            if sample_type == 'Baseline':
                color = sample_type_colors[sample_type]
            else:
                color = progression_sample_colors[patient_id][date_collected]
            ax.bar(x_values, y_values, bar_width, label=f'{patient_id} - {sample_type} - {date_collected}', color=color)
            # Append unique mutations for this patient and sample
            # patient_unique_mutations.extend(unique_mutations['Protein_annotation'].tolist())
        # Add x ticks for each mutation
        x_ticks.extend(np.arange(len(unique_mutations)) + i * bar_width / 2)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(unique_mutations['Gene'] + "\n" + unique_mutations["Protein_annotation"], rotation=rotation)  # Rotate labels for better visibility
        ax.set_xlabel('')
        ax.set_ylabel('WBC VAF%')
        ax.set_title(f'All mutations identified in Patient {patient_id} - {diagnosis}')
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(figure_dir, "Patient_profiles_bar_graph", f'{diagnosis}_{patient_id}.{ext}'))

def add_vaf_0_for_missing_mutations(df, PATH_sample_information):
    sample_info_main = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    results_df = []
    for i, group in df.groupby(["Patient_id", "Gene", "Chrom", "Position", "Protein_annotation", "Consequence", "Variant_type"]):
        patient_id = group["Patient_id"].unique()[0]
        sample_info = sample_info_main[sample_info_main["Patient_id"] == patient_id]
        merged = sample_info.merge(group, how = "outer")
        merged['VAF_t'] = merged['VAF_t'].fillna(0)
        merged['VAF_n'] = merged['VAF_n'].fillna(0)
        merged = merged.ffill()
        results_df.append(merged)
    concatted = pd.concat(results_df).reset_index(drop = True)
    return(concatted)

def plot_patient_timeline(ax_clinical_timeline, ax_clinical_timeline_legend, PATH_treatment_landscape, patient_id, sample_info_main):
    # Some modifications in treatment df for plotting
    treatment_main = pd.read_csv(PATH_treatment_landscape)
    treatment_main = treatment_main[treatment_main["Drug"] != "Palliative / best supportive care"].reset_index(drop = True)
    treatment_main = treatment_main[treatment_main["Treatment"] != "Other (describe)"]
    treatment_main = treatment_main[~pd.isnull(treatment_main["Treatment"])]
    treatment_main['Date start'] = pd.to_datetime(treatment_main['Date start'])
    treatment_main['Date discontinuation'] = pd.to_datetime(treatment_main['Date discontinuation'])
        
    # Assign colors to to each drug for plotting
    unique_values = treatment_main['Drug'].dropna().unique()
    colors = plt.cm.get_cmap('tab20')(np.linspace(0, 1, len(unique_values))).tolist()  # Using the Tab20 colormap for distinct colors
    color_map = dict(zip(unique_values, colors)) # Create a dictionary to map colors to unique values
    color_map_tuples = {key: tuple(value) for key, value in color_map.items()}
    treatment_main['Drug color'] = treatment_main['Drug'].map(color_map_tuples) # Map colors to the dataframe
    
    # Subset the treatment dfs to the patient   
    pt_treatment = treatment_main[treatment_main["Patient_id"] == patient_id].reset_index(drop=True)        
    legend_labels = []
    for idx, row in pt_treatment.iterrows():
        if pd.isnull(row['Date discontinuation']) or pd.isna(row['Date discontinuation']):
            dummy_duration = pd.Timedelta(days=10)
            ax_clinical_timeline.barh(idx, dummy_duration, height=0.3, left=row['Date start'], color=row["Drug color"], alpha=0.8, label=row['Drug'])
            ax_clinical_timeline.text(row['Date start'], idx, 'End date unknown', va='center', ha='left', color='black', fontsize=10)
            ax_clinical_timeline.text(row['Date start'], idx-0.05, row["Treatment line"], va='center', ha='left', color='black', fontsize=10) 
        else:
            duration = row['Date discontinuation'] - row['Date start']
            ax_clinical_timeline.barh(idx, duration.days, height=0.3, left=row['Date start'], color=row["Drug color"], alpha=0.8, label=row['Drug'])
            ax_clinical_timeline.text(row['Date start']+pd.Timedelta(days=5), idx, row["Treatment line"], va='center', ha='left', color='black', fontsize=10) 
        if row['Drug'] not in legend_labels:
            legend_labels.append(row['Drug'])
    # Plot when the GUBB blood samples are taken
    gubb_blood = sample_info_main[sample_info_main["Patient_id"] == patient_id].reset_index(drop = True)
    gubb_blood["Date_collected"] = pd.to_datetime(gubb_blood["Date_collected"], format='%Y%b%d')
    for date in gubb_blood["Date_collected"]:
        ax_clinical_timeline.scatter(date, 0.8, marker='v', color='red', s=100)
        ax_clinical_timeline.axvline(date, color='red', linestyle='--')
    # Customize the plot
    ax_clinical_timeline.set_yticks([])
    ax_clinical_timeline.set_xlabel('')
    ax_clinical_timeline.set_ylabel('')
    ax_clinical_timeline.spines["top"].set_visible(False)
    ax_clinical_timeline.spines["right"].set_visible(False)
    ax_clinical_timeline.set_title(f"Clinical timeline")
    handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in pt_treatment['Drug color'].unique()]
    legend_handles = handles + [plt.Line2D([0], [0], marker='v', color='red', markersize=10, linestyle='None', label='GUBB blood draw')]
    legend_labels.append('GUBB blood draw')
    ax_clinical_timeline_legend.legend(handles=legend_handles, labels=legend_labels, loc='best')
    
    ax_clinical_timeline_legend.spines[["top", "right", "bottom", "left"]].set_visible(False)
    ax_clinical_timeline_legend.set_xticks([])
    ax_clinical_timeline_legend.set_xticklabels([])
    ax_clinical_timeline_legend.set_yticks([])
    ax_clinical_timeline_legend.set_yticklabels([])
    
    # Format x-axis ticks
    interval = 6
    ax_clinical_timeline.xaxis.set_major_locator(mdates.MonthLocator(interval=6))   
    ax_clinical_timeline.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
    
    # Calculate the number of months elapsed from the minimum date
    min_date_pt = pt_treatment["Date start"].min()
    def months_elapsed(date):
        return (date - min_date_pt).days // 30
    
    # Modifications on the x ticks
    n_xticks = len(ax_clinical_timeline.get_xticks())
    if n_xticks > 10:
        if n_xticks <= 15:
            interval = 10
        elif n_xticks <= 20:
            interval = 20
        else:
            interval = 30
    # Set the locator and formatter for x-axis ticks
    ax_clinical_timeline.xaxis.set_major_locator(mdates.MonthLocator(interval=interval)) 
    ax_clinical_timeline.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))  
    # ax_clinical_timeline.xaxis.set_major_formatter(mdates.FuncFormatter(lambda x, _: f'{months_elapsed(pd.Timestamp(x))}\nm'))
    
    # for tick in ax.get_xticklabels():
    #     tick.set_rotation(45)    
    # Add buffer around x-axis
    start_dates = pd.concat([pt_treatment['Date start'], gubb_blood["Date_collected"]]).dropna()
    end_dates = pd.concat([pt_treatment['Date discontinuation'], gubb_blood["Date_collected"]]).dropna()
    if not start_dates.empty and not end_dates.empty:
        min_date = min(start_dates) - pd.DateOffset(months=6)
        max_date = max(end_dates) + pd.DateOffset(months=6)
        ax_clinical_timeline.set_xlim(min_date, max_date)
    return[ax_clinical_timeline, ax_clinical_timeline_legend]

# def annotate_patient_samples():

def CALL_generate_patient_profiles_dot_plot(baseline_chip, prog_chip, baseline_somatic, prog_somatic, figure_dir, PATH_sample_information, PATH_treatment_landscape, show_dates_on_x_axis = False, showtimeline = True, ext = "png"):    
    sample_info_main = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    merged_df_chip = pd.concat([baseline_chip, prog_chip]).reset_index(drop=True).assign(Variant_type = "CHIP")
    merged_df_somatic = pd.concat([baseline_somatic, prog_somatic]).reset_index(drop=True).assign(Variant_type = "Somatic")
    # merged_df_somatic = add_vaf_0_for_missing_mutations(merged_df_somatic, PATH_sample_information)
    merged = pd.concat([merged_df_chip, merged_df_somatic]).reset_index(drop=True)
    for (patient_id, diagnosis), group in merged.groupby(["Patient_id", "Diagnosis"]):
        if any(group["Patient_id"].unique()[0] == str(pid) for pid in sample_info_main["Patient_id"].astype(str).tolist()):
            print(patient_id, diagnosis)
            sample_info = sample_info_main[sample_info_main["Patient_id"] == patient_id]
            if sample_info["Diagnosis"].unique()[0] == "Bladder":
                showtimeline = True
            else:
                showtimeline = False
            if showtimeline:
                fig = plt.figure(figsize=(14, 10))
                gs = fig.add_gridspec(3, 3, width_ratios=[10, 5, 15], height_ratios = [1, 1, 0.1], wspace=0.4, hspace=0.5)
                ax_chip = fig.add_subplot(gs[0, 0])
                ax_somatic = fig.add_subplot(gs[1, 0], sharex = ax_chip)
                ax_treatment = fig.add_subplot(gs[2, 0], sharex = ax_chip)
                ax_clinical_timeline = fig.add_subplot(gs[0, 2])
                ax_clinical_timeline_legend = fig.add_subplot(gs[1, 2])
                [ax_clinical_timeline, ax_clinical_timeline_legend] = plot_patient_timeline(ax_clinical_timeline, ax_clinical_timeline_legend, PATH_treatment_landscape, patient_id, sample_info_main)
                ax_treatment = plot_treatment_landscape(ax_treatment, ax_chip, patient_id, PATH_treatment_landscape, sample_info)                
            else:
                fig = plt.figure(figsize=(8, 8))
                gs = fig.add_gridspec(3, 3, width_ratios=[10, 5], height_ratios = [1, 1, 0.1], wspace=0.4, hspace=0.5)
                ax_chip = fig.add_subplot(gs[0, 0])
                ax_somatic = fig.add_subplot(gs[1, 0], sharex = ax_chip)                    
            plt.subplots(layout="constrained")
            # Prepare to plot
            patient_chip_df = group[group["Variant_type"] == "CHIP"]
            patient_somatic_df = group[group["Variant_type"] == "Somatic"]
            if patient_id == "20-383":
                condition = (patient_somatic_df['Position'] == 1805767) & ((patient_somatic_df['Ref'] == 'G') & ((patient_somatic_df['Alt'] == 'T') | (patient_somatic_df['Alt'] == 'A')))
                patient_somatic_df = patient_somatic_df[~condition].reset_index(drop = True)
            # Plotting
            ax_chip = generate_patient_profiles_dot_plot(patient_chip_df, figure_dir, PATH_sample_information, show_dates_on_x_axis, ax_chip, "tab10", patient_id, diagnosis, "Clonal hematopoiesis", title = "CH mutations")
            ax_somatic = generate_patient_profiles_dot_plot(patient_somatic_df, figure_dir, PATH_sample_information, show_dates_on_x_axis, ax_somatic, "Accent", patient_id, diagnosis, "ctDNA", title = "ctDNA mutations")
            # SAVE            
            fig.suptitle(f'Somatic mutations in {patient_id} - {diagnosis.upper()}', fontsize=15)
            fig.savefig(os.path.join(figure_dir, "Patient_profiles_line_graph", f'{diagnosis}_{patient_id}_dot_plot.{ext}'))
        else: 
            print(f"Patient {patient_id} not in the sample sheet. Skipping.")  

def generate_patient_profiles_dot_plot(patient_df, figure_dir, PATH_sample_information, show_dates_on_x_axis, ax_main, color_map_name, patient_id, diagnosis, variant_type, title, show_legend = True, add_figure_title = True, spine_to_hide = ["top", "right"], add_somatic_to_legend = False, gene=None):    
    """
    Plot all mutations found in all samples of the patient in the same panel.
    Mutations are represented as dots and connected by lines, showing the trajectory of a mutation through the samples taken from the patient.
    """
    unique_mutations = patient_df[['Chrom', 'Position', 'Gene', 'Protein_annotation', 'Consequence']].drop_duplicates()
    def colormap_exists(name):
        try:
            cm.get_cmap(name)
            return True
        except ValueError:
            return False
    if colormap_exists(color_map_name): 
        unique_mutations_colors = {mutation: get_cmap(color_map_name)(i) for i, mutation in enumerate(unique_mutations.itertuples(index=False, name=None))}
    else: 
        unique_mutations_colors = {mutation: color_map_name  for i, mutation in enumerate(unique_mutations.itertuples(index=False, name=None))}
    x_ticks = []
    x_values_list = []
    y_values_list = []
    legend_elements = []  # List to store legend elements
    
    ####################################################
    # X axis tick labelling with the dates
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    sample_info = sample_info[sample_info["Patient_id"] == patient_id]
    sample_info["Date_collected"] = pd.to_datetime(sample_info["Date_collected"], format='%Y%b%d')
    
    baseline_date = sample_info[sample_info["Timepoint"] == "Baseline"]["Date_collected"].iloc[0]
    sample_info["months_from_baseline"] = (sample_info["Date_collected"] - baseline_date).dt.days // 30
    ax_main.set_xticks(list(sample_info["Date_collected"])) # number of samples
    if not show_dates_on_x_axis:
        def generate_timepoint_labels(months_from_baseline):
            if months_from_baseline == 0:
                return "Base"
            else:
                return f"{months_from_baseline}\nm"
        tick_labels = sample_info["months_from_baseline"].apply(generate_timepoint_labels)
        ax_main.set_xticklabels(tick_labels)
    ####################################################
    # Y axis labelling
    if variant_type == "ctDNA": 
        ax_main.set_ylabel('ctDNA VAF%')
    elif variant_type == "Clonal hematopoiesis":  
        ax_main.set_ylabel('WBC VAF%')
        col_to_use_for_chip = "VAF_n"
    ####################################################
    for i, (sample_type, date_collected) in enumerate(patient_df.groupby(['Timepoint', 'Date_collected']).groups.keys()):
        sample_group = patient_df[(patient_df['Timepoint'] == sample_type) & (patient_df['Date_collected'] == date_collected)].sort_values(by=['Chrom', 'Position', 'Gene', 'Protein_annotation', 'Consequence'])
        sample_group = sample_group.drop_duplicates()
        sample_group = sample_group.reset_index().set_index(['Chrom', 'Position', 'Gene', 'Protein_annotation', 'Consequence']).reindex(unique_mutations)
        if variant_type == "ctDNA": 
            y_values = sample_group['VAF_t'].values
        else: 
            y_values = sample_group[col_to_use_for_chip].values
        y_values_list.append(y_values)
        x_values_list.append(np.repeat(pd.to_datetime(date_collected, format = '%Y%b%d'), len(y_values)).tolist())
    for j, mutation in enumerate(unique_mutations.iterrows()):
        zorder = 1
        mutation_key = tuple(mutation[1].to_numpy())
        if mutation_key not in legend_elements:
            legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=unique_mutations_colors[mutation_key], markersize=8, label=str(mutation[1]["Gene"] + " " + mutation[1]["Protein_annotation"])))
            mutation_color = unique_mutations_colors.get(mutation_key, 'black')
            if gene is not None and mutation[1]["Gene"] == gene:
                mutation_color = 'red'
                zorder = 10000
            elif gene is not None and mutation[1]["Gene"] != gene:
                mutation_color = 'grey'
            if len(y_values_list) == 1:
                ax_main.plot(baseline_date, [y_values_list[0][j]][0], marker='o', color=mutation_color, markersize=8, linestyle='-', linewidth=1, zorder = zorder)
            elif len(y_values_list) > 1:
                x_values = [vals[j] for vals in x_values_list]
                y_values = [vals[j] for vals in y_values_list]
                ax_main.plot(x_values, y_values, marker='o', color=mutation_color, markersize=8, linestyle='-', linewidth=1, zorder = zorder)
    # Legend plot
    # if color_mutations: 
    #     ax_main.legend(handles=legend_elements, title="Mutations", loc='center left', bbox_to_anchor=(1, 0.5))
    #     ax_legend.axis('off')  # Turn off the axis for the legend subplot
    #     if len(x_ticks) == 1:
    #         ax.set_xlim((-1, 1))
    if show_legend:
        if add_somatic_to_legend:
            legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor="lightgrey", markersize=8, label="ctDNA alterations"))
        ax_main.legend(handles=legend_elements, title="Mutations", loc='center left', bbox_to_anchor=(1.2, 0.5), frameon=False)
    if len(x_ticks) == 1:
        ax.set_xlim((-1, 1))
    # ax_main.set_xticks(range(len(x_ticks)))
    # ax_main.set_xticklabels(x_ticks)
    ax_main.set_xlabel('')
    ax_main.set_title(title)
    for spine in spine_to_hide: 
        ax_main.spines[spine_to_hide].set_visible(False)
    return ax_main

def plot_treatment_landscape(ax_treatment, ax_chip, patient_id, PATH_treatment_landscape, sample_info):
    # Some modifications in treatment df for plotting
    treatment_main = pd.read_csv(PATH_treatment_landscape)
    treatment_main = treatment_main[treatment_main["Drug"] != "Palliative / best supportive care"].reset_index(drop = True)
    treatment_main = treatment_main[treatment_main["Treatment"] != "Other (describe)"]
    treatment_main = treatment_main[~pd.isnull(treatment_main["Treatment"])]
    treatment_main['Date start'] = pd.to_datetime(treatment_main['Date start'])
    treatment_main['Date discontinuation'] = pd.to_datetime(treatment_main['Date discontinuation'])
        
    # Assign colors to to each drug for plotting
    unique_values = treatment_main['Drug'].dropna().unique()
    colors = plt.cm.get_cmap('tab20')(np.linspace(0, 1, len(unique_values))).tolist()  # Using the Tab20 colormap for distinct colors
    color_map = dict(zip(unique_values, colors)) # Create a dictionary to map colors to unique values
    color_map_tuples = {key: tuple(value) for key, value in color_map.items()}
    treatment_main['Drug color'] = treatment_main['Drug'].map(color_map_tuples) # Map colors to the dataframe
    
    # Subset the treatment dfs to the patient   
    pt_treatment = treatment_main[treatment_main["Patient_id"] == patient_id].reset_index(drop=True)
    pt_treatment['Date start'] = pd.to_datetime(pt_treatment['Date start'])
    pt_treatment['Date discontinuation'] = pd.to_datetime(pt_treatment['Date discontinuation'])
    
    # Determine boundaries of interest
    sample_info.loc[:, "Date_collected"] = pd.to_datetime(sample_info["Date_collected"], format='%Y%b%d')
    min_date = sample_info["Date_collected"].min()
    max_date = sample_info["Date_collected"].max()
    
    # Get the treatments the pt got within the timeframe we are interested in
    # filtered_df = pt_treatment[(pt_treatment['Date start'] >= min_date) & (pt_treatment['Date discontinuation'] <= max_date)]
    
    # Calculate the width of the bars based on the time duration
    bar_width = 0.3  # Adjust this value as needed to align with ax_chip
    
    legend_labels = []
    for idx, row in pt_treatment.iterrows():
        if not pd.isnull(row['Date discontinuation']) and not pd.isna(row['Date discontinuation']):
            duration = row['Date discontinuation'] - row['Date start']
            # Calculate the left position of the bar based on the normalized date
            left_pos = row['Date start']
            ax_treatment.barh(0, duration.days, height=1, left=left_pos, color=row["Drug color"], alpha=0.8, label=row['Drug'], align='edge')
            if row['Drug'] not in legend_labels:
                legend_labels.append(row['Drug'])
    
    # Set the same x-axis limits and ticks as ax_chip
    # ax_treatment.set_xlim(ax_chip.get_xlim())
    ax_treatment.tick_params(axis='x', bottom=False, labelbottom=False)
        
    # Hide the y-axis ticks and label since they are not needed
    ax_treatment.set_yticks([])
    ax_treatment.set_ylabel('')
    
    ax_treatment.spines[["top", "right", "left", "bottom"]].set_visible(False)
    
    return ax_treatment

def plot_vaf(df, figure_dir, figure_name, select_genes=None):
    """
    Plots the VAF of mutations at baseline as a bar chart. Sorts by gene and VAF.
    """
    df = df.sort_values(by=["Gene", "VAF_t"]).reset_index(drop=True)
    diagnosis = df["Diagnosis"].unique()[0].lower()
    timepoint = df["Timepoint"].unique()[0].lower()
    if select_genes is None:
        select_genes = ["DNMT3A", "TET2", "ASXL1", "JAK2", "BRCA1", "BRCA2", "ATM", "TP53", "PPM1D", "CHEK2"]
    
    fig, axes = plt.subplots(2, 5, figsize=(15, 6), sharex=False)
    
    for i, gene in enumerate(select_genes):
        gene_df = df[df["Gene"] == gene].reset_index(drop = True)
        row, col = divmod(i, 5)  # Determine the row and column for subplot
        axes[row, col].bar(gene_df.index, gene_df['VAF_t'], label=gene, color = "black")
        axes[row, col].set_title(f"{gene}")
        
        # Set labels and title for each subplot
        axes[row, col].set_xticks([])
        axes[row, col].set_xticklabels([])
        axes[row, col].set_ylabel('VAF%')
        axes[row, col].spines["top"].set_visible(False)
        axes[row, col].spines["right"].set_visible(False)
        
        # Add horizontal dashed line through the median VAF
        median_vaf = gene_df['VAF_t'].median()
        axes[row, col].axhline(median_vaf, linestyle='--', color='red', label=f'Median VAF ({median_vaf:.2f})')       
    fig.suptitle(f'VAF of CH events in select genes {diagnosis} - {timepoint}')
    # Adjust layout and save the figure
    fig.tight_layout()
    fig.savefig(os.path.join(figure_dir, figure_name))

# Example usage
# plot_vaf(your_dataframe, '/path/to/figure_directory', 'your_figure_name.png')

def print_perc_CHpositive_baseline(base, PATH_sample_info, gene = None):
    """
    Prints some stats about the cohort.
    """
    if gene != None:
        base = base[base["Gene"] == gene]
        message = f"Percentage of the cohort CH positive in {gene} at baseline is " 
    else: 
        message = "Percentage of the cohort CH positive at in any gene at baseline is " 
    baseline_info = pd.read_csv(PATH_sample_info, "\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"]).query('Timepoint == "Baseline"').reset_index(drop = True)
    baseline_pts = baseline_info.assign(Type="all_pts")[["Patient_id", "Type"]]
    base = base.assign(Type="Baseline_CH_positive_pts")[["Patient_id", "Type"]]
    combined = baseline_pts.merge(base, how = "left", on = "Patient_id", indicator = True)
    x = combined[combined["_merge"] == "both"].shape[0]
    value1 = round((x/baseline_pts.shape[0])*100, 1)
    return(f"{message}{value1}")

def per_gene_counts(df): 
    """
    Plots a bar chart showing how mnay 
    """
    # Count the occurrences of each combination of 'Gene' and 'Consequence'
    counts = df.groupby(['Gene', 'Consequence']).size().reset_index(name='Count')
    # Pivot the DataFrame for plotting
    pivot_df = counts.pivot(index='Gene', columns='Consequence', values='Count').fillna(0)
    # Order rows (genes) by the sum of counts in descending order
    row_order = pivot_df.sum(axis=1).sort_values(ascending=False).index
    pivot_df = pivot_df.loc[row_order]
    sns.set(style="white")
    ax = pivot_df.plot(kind='bar', stacked=True, colormap='viridis', figsize=(10, 6))
    # Adding labels and title
    ax.set_xlabel("Gene")
    ax.set_ylabel("Count")
    ax.set_title("Mutations effects")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 0)
    plt.tight_layout()
    plt.savefig('/groups/wyattgrp/users/amunzur/therapy_chip/results/figures/mutation_effects_and_genes.png')

def plot_gene_counts(df, figure_dir, filename, figure_title, unique = False, fontsize = 10):
    """
    Plots the number of mutations each individual gene has in the cohort as a BAR CHART. If unique, non unique mutations in the same gene are counted only once.
    Colored by the type of mutation (like missense, nonsense etc)
    """
    df["Consequence"] = df["Consequence"].replace({"Frameshift Deletion": "Frameshift deletion", "Startloss": "Missense"})
    if unique:
        df = df.drop_duplicates(subset=["Patient_id", "Gene"]).reset_index(drop = True)
    counts = df.groupby(['Gene', 'Consequence']).size().reset_index(name='Count')
    # Pivot the DataFrame for plotting
    pivot_df = counts.pivot(index='Gene', columns='Consequence', values='Count').fillna(0)
    # Order rows (genes) by the sum of counts in descending order
    row_order = pivot_df.sum(axis=1).sort_values(ascending=False).index
    pivot_df = pivot_df.loc[row_order]
    sns.set(style="white")
    fig, ax = plt.subplots(figsize=(7.5, 3.5))
    pivot_df.plot(kind='bar', stacked=True, colormap='viridis', ax=ax, width=0.9)
    ax.set_xlabel('Genes with at least one CH event', fontsize=fontsize, color='k')
    ax.set_ylabel('Number of mutations', fontsize=fontsize, color='k')
    ax.set_title(figure_title, fontsize=fontsize, color='k')  # Increase title fontsize by 2 for emphasis
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(title="", frameon=False, fontsize=fontsize)
    ax.tick_params(axis="both", direction="out", which="both", left=True, bottom=False)
    # Set spines, ticks, tick labels, and legend color to black
    ax.spines["top"].set_color('k')
    ax.spines["right"].set_color('k')
    ax.tick_params(axis="both", direction="out", which="both", left=True, bottom=False, colors='k')
    ax.legend(title="", frameon=False, fontsize=fontsize-3, title_fontsize=fontsize-3)
    for patch in ax.patches:
        patch.set_linewidth(0)
    fig.tight_layout()
    fig.savefig(os.path.join(figure_dir, filename))

def plot_gene_counts_mirrored(df, figure_dir, filename, kidney_clinical, bladder_clinical):
    """
    Plots the number of mutations each individual gene has in the cohort as a BAR CHART. If unique, non unique mutations in the same gene are counted only once.
    Plots kidney and bladder separately. Mirrors. 
    """
    # sex info
    kidney_sex = kidney_clinical[["Patient_id", "Sex"]].drop_duplicates()
    bladder_sex = bladder_clinical[["Patient_id", "Sex"]].drop_duplicates()
    sex = pd.concat([kidney_sex, bladder_sex]).reset_index(drop = True)
    
    df = df.drop_duplicates(subset=["Patient_id", "Gene"]).reset_index(drop = True)
    genes_to_include = ["TERT", "U2AF1", "CREBBP", "BRCC3", "BRCA2", "SRSF2", "SH2B3", "NF1", "JAK2", "CBL", "SF3B1", "RUNX1", "RAD21", "KMT2D", "GNAS", "STAG2", "TP53", "CHEK2", "ATM", "ASXL1", "PPM1D", "TET2", "DNMT3A"]
    df = df[df["Gene"].isin(genes_to_include)].reset_index(drop = True)
    df = df.merge(sex)
    counts = df.groupby(['Gene', 'Diagnosis', 'Sex']).size().reset_index(name='Count')
    bladder_df = counts[counts["Diagnosis"] == "Bladder"].pivot(index='Gene', columns='Sex', values='Count').fillna(0)
    bladder_total = bladder_df["Female"] + bladder_df["Male"]
    bladder_total = bladder_total.reset_index().sort_values(by = 0, ascending = True)
    
    kidney_df = counts[counts["Diagnosis"] == "Kidney"].pivot(index='Gene', columns='Sex', values='Count').fillna(0)
    kidney_df = kidney_df.reindex(index=bladder_total["Gene"])
    bladder_df = bladder_df.reindex(index=bladder_total["Gene"])
    kidney_df["Male"] = kidney_df["Male"]*-1
    kidney_df["Female"] = kidney_df["Female"]*-1
    
    kidney_m = kidney_df["Male"].reset_index()
    kidney_f = kidney_df["Female"].reset_index()
    bladder_m = bladder_df["Male"].reset_index()
    bladder_f = bladder_df["Female"].reset_index()
    
    # need to melt
    # kidney_df = kidney_df.reset_index().melt(id_vars='Gene', var_name='Sex', value_name='Count')
    # bladder_df = bladder_df.reset_index().melt(id_vars='Gene', var_name='Sex', value_name='Count')
    
    # Plotting
    fig = plt.figure(figsize=(5, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], wspace = 0)
    ax1 = plt.subplot(gs[0]) # kidney
    ax2 = plt.subplot(gs[1], sharey = ax1) # bladder
    ax1_twin = ax1.twiny()
    ax2_twin = ax2.twiny()
    
    # Plotting
    ax1_twin.barh(kidney_m["Gene"], kidney_m["Male"], color='black', edgecolor = "None")
    ax1_twin.barh(kidney_f["Gene"], kidney_f["Female"], color='dimgray', left = kidney_m["Male"].fillna(0), edgecolor = "None")
    ax2_twin.barh(bladder_m["Gene"], bladder_m["Male"], color='black', edgecolor = "None")
    ax2_twin.barh(bladder_f["Gene"], bladder_f["Female"], color='dimgray', left = bladder_m["Male"], edgecolor = "None")
    
    # sns.barplot(x='Count', y='Gene', label = "Sex", data=kidney_df, ax=ax1_twin, orient='h')
    # sns.barplot(x='Count', y='Gene', label = "Sex", data=bladder_df, ax=ax2_twin, orient='h')
    
    ax1.spines[["left", "bottom"]].set_visible(False)
    ax2.spines[["bottom", "right"]].set_visible(False)
    ax1_twin.spines[["left", "bottom"]].set_visible(False)
    ax2_twin.spines[["bottom", "right"]].set_visible(False)
    
    ax1_twin.set_ylabel("")
    ax2_twin.set_ylabel("")
    ax1_twin.set_xlabel("Kidney")
    ax2_twin.set_xlabel("Bladder")
    ax1_twin.set_xticks([0, -20, -40, -60, -80])
    ax1_twin.set_xticklabels(["0", "20", "40", "60", "80"])
    ax2_twin.set_xticks([0, 20, 40, 60, 80])
    ax2_twin.set_xticklabels(["0", "20", "40", "60", "80"])  
    
    ax1_twin.set_ylim((-1, 41.5))
    ax2_twin.set_ylim((-1, 41.5))
    ax1.set_ylim((-1, 41.5))
    ax2.set_ylim((-1, 41.5))
    
    # set the y tick labels
    gene_names = bladder_f["Gene"].tolist()
    ax2.set_yticklabels(gene_names, ha='center')
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax2.yaxis.set_tick_params(pad=30)   
    ax1.set_xticks([])
    ax2.set_xticks([])
    # Add annotations on top of the bars
    # def annotate_bars(ax, offset):
    #     # Iterate through patches in pairs
    #     for i in range(0, len(ax.patches), 2):
    #         male_rect = ax.patches[i]
    #         female_rect = ax.patches[i + 1]
    #         male_width = int(male_rect.get_x())
    #         female_width = int(female_rect.get_width()) - male_width
    #         total = male_width + female_width
    #         ax.annotate(f'{total}', 
    #                     xy=(male_rect.get_x() + male_rect.get_width(), 
    #                         male_rect.get_y() + male_rect.get_height() / 2),
    #                     xytext=(offset, 0),  # 3 points horizontal offset
    #                     textcoords="offset points",
    #                     ha='left', va='center', fontsize=10)
    
    # annotate_bars(ax1_twin, -12)
    # annotate_bars(ax2_twin, 3)
    
    # add legend
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=10, linestyle='') for color, label in zip(["black", "dimgray"], ["Male", "Female"])]
    ax2_twin.legend(handles=legend_handles, loc="lower right", frameon=False)
    
    gs.tight_layout(fig)
    fig.savefig(os.path.join(figure_dir, filename))

def plot_gene_counts_grouped(all_vars_chip, figure_dir, filename, kidney_clinical, bladder_clinical, gene_list = None):
    """
    Plots the number of mutations each individual gene has in the cohort as a BAR CHART. If unique, non unique mutations in the same gene are counted only once.
    Plots kidney and bladder separatelyin grouped bar charts. 
    """
    df = all_vars_chip[all_vars_chip["Gene"].isin(gene_list)].reset_index(drop = True)
    # df["VAF_n"] = df["VAF_n"]
    #
    # subset
    bladder_df = df[df["Diagnosis"] == "Bladder"].reset_index(drop = True)["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "Counts_bladder"}).reset_index()
    kidney_df = df[df["Diagnosis"] == "Kidney"].reset_index(drop = True)["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "Counts_kidney"})
    kidney_df = kidney_df.merge(bladder_df[["Gene", "index"]], how = "inner") # gene order will be determined based on bladder df
    combined = bladder_df.merge(kidney_df, how = "inner")
    # vaf info
    kidney_vaf = df[df["Diagnosis"] == "Kidney"][["Gene", "VAF_n"]].reset_index(drop = True).merge(bladder_df[["Gene", "index"]])
    bladder_vaf = df[df["Diagnosis"] == "Bladder"][["Gene", "VAF_n"]].reset_index(drop = True).merge(bladder_df[["Gene", "index"]])
    #
    kidney_color = "orangered"
    bladder_color = "deepskyblue"
    #
    # plotting
    fig, ax = plt.subplots(figsize = (5.6, 2.8))
    ax2 = ax.twinx()
    for i, row in combined.iterrows():
        ax.bar(row["index"]+0.2, row["Counts_kidney"], color=kidney_color, width = 0.4)
        ax.bar(row["index"]-0.2, row["Counts_bladder"], color=bladder_color, width = 0.4)
    # plot the vafs for kidney
    for i, row in kidney_vaf.iterrows():
        jitter = np.random.normal(-0.05, 0.05, 1)
        ax2.scatter(row["index"]+0.2+jitter, row["VAF_n"], color="black", s = 2)
    # plot the vafs for bladder
    for i, row in bladder_vaf.iterrows():
        jitter = np.random.normal(-0.05, 0.05, 1)
        ax2.scatter(row["index"]-0.2+jitter, row["VAF_n"], color="black", s = 2)
    # aesthetics
    ax.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    # x ticks and labels
    tick_df = bladder_df[["index", "Gene"]].sort_values(by = "index")
    ax.set_xticks(tick_df["index"].tolist())
    ax.set_xticklabels(tick_df["Gene"].tolist(), rotation = 90)
    ax.tick_params(axis = "x", which='both', length=0)
    ax.set_ylabel("Count")
    ax.set_yticks([0, 25, 50, 75, 100, 125, 150])
    ax.set_yticklabels(["0", "25", "50", "75", "100", "125", "150"])
    ax.set_xlim((-0.6, 14.6))
    ax2.set_ylabel("WBC VAF%")
    ax2.set_ylim((0, 50))
    ax2.set_yticks([0, 10, 20, 30, 40, 50])
    ax2.set_yticklabels(["0", "10", "20", "30", "40", "50"])
    ax2.tick_params(axis='x', pad=1)
    fig.tight_layout()
    fig.savefig(os.path.join(figure_dir, filename))








    # sex info
    kidney_sex = kidney_clinical[["Patient_id", "Sex"]].drop_duplicates()
    bladder_sex = bladder_clinical[["Patient_id", "Sex"]].drop_duplicates()
    sex = pd.concat([kidney_sex, bladder_sex]).reset_index(drop = True)
    
    df = df.drop_duplicates(subset=["Patient_id", "Gene"]).reset_index(drop = True)
    genes_to_include = ["TERT", "U2AF1", "CREBBP", "BRCC3", "BRCA2", "SRSF2", "SH2B3", "NF1", "JAK2", "CBL", "SF3B1", "RUNX1", "RAD21", "KMT2D", "GNAS", "STAG2", "TP53", "CHEK2", "ATM", "ASXL1", "PPM1D", "TET2", "DNMT3A"]
    df = df[df["Gene"].isin(genes_to_include)].reset_index(drop = True)
    df = df.merge(sex)
    counts = df.groupby(['Gene', 'Diagnosis', 'Sex']).size().reset_index(name='Count')
    bladder_df = counts[counts["Diagnosis"] == "Bladder"].pivot(index='Gene', columns='Sex', values='Count').fillna(0)
    bladder_total = bladder_df["Female"] + bladder_df["Male"]
    bladder_total = bladder_total.reset_index().sort_values(by = 0, ascending = True)
    
    kidney_df = counts[counts["Diagnosis"] == "Kidney"].pivot(index='Gene', columns='Sex', values='Count').fillna(0)
    kidney_df = kidney_df.reindex(index=bladder_total["Gene"])
    bladder_df = bladder_df.reindex(index=bladder_total["Gene"])
    kidney_df["Male"] = kidney_df["Male"]*-1
    kidney_df["Female"] = kidney_df["Female"]*-1
    
    kidney_m = kidney_df["Male"].reset_index()
    kidney_f = kidney_df["Female"].reset_index()
    bladder_m = bladder_df["Male"].reset_index()
    bladder_f = bladder_df["Female"].reset_index()
    
    # need to melt
    # kidney_df = kidney_df.reset_index().melt(id_vars='Gene', var_name='Sex', value_name='Count')
    # bladder_df = bladder_df.reset_index().melt(id_vars='Gene', var_name='Sex', value_name='Count')
    
    # Plotting
    fig = plt.figure(figsize=(5, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], wspace = 0)
    ax1 = plt.subplot(gs[0]) # kidney
    ax2 = plt.subplot(gs[1], sharey = ax1) # bladder
    ax1_twin = ax1.twiny()
    ax2_twin = ax2.twiny()
    
    # Plotting
    ax1_twin.barh(kidney_m["Gene"], kidney_m["Male"], color='black', edgecolor = "None")
    ax1_twin.barh(kidney_f["Gene"], kidney_f["Female"], color='dimgray', left = kidney_m["Male"].fillna(0), edgecolor = "None")
    ax2_twin.barh(bladder_m["Gene"], bladder_m["Male"], color='black', edgecolor = "None")
    ax2_twin.barh(bladder_f["Gene"], bladder_f["Female"], color='dimgray', left = bladder_m["Male"], edgecolor = "None")
    
    # sns.barplot(x='Count', y='Gene', label = "Sex", data=kidney_df, ax=ax1_twin, orient='h')
    # sns.barplot(x='Count', y='Gene', label = "Sex", data=bladder_df, ax=ax2_twin, orient='h')
    
    ax1.spines[["left", "bottom"]].set_visible(False)
    ax2.spines[["bottom", "right"]].set_visible(False)
    ax1_twin.spines[["left", "bottom"]].set_visible(False)
    ax2_twin.spines[["bottom", "right"]].set_visible(False)
    
    ax1_twin.set_ylabel("")
    ax2_twin.set_ylabel("")
    ax1_twin.set_xlabel("Kidney")
    ax2_twin.set_xlabel("Bladder")
    ax1_twin.set_xticks([0, -20, -40, -60, -80])
    ax1_twin.set_xticklabels(["0", "20", "40", "60", "80"])
    ax2_twin.set_xticks([0, 20, 40, 60, 80])
    ax2_twin.set_xticklabels(["0", "20", "40", "60", "80"])  
    
    ax1_twin.set_ylim((-1, 41.5))
    ax2_twin.set_ylim((-1, 41.5))
    ax1.set_ylim((-1, 41.5))
    ax2.set_ylim((-1, 41.5))
    
    # set the y tick labels
    gene_names = bladder_f["Gene"].tolist()
    ax2.set_yticklabels(gene_names, ha='center')
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax2.yaxis.set_tick_params(pad=30)   
    ax1.set_xticks([])
    ax2.set_xticks([])
    # Add annotations on top of the bars
    # def annotate_bars(ax, offset):
    #     # Iterate through patches in pairs
    #     for i in range(0, len(ax.patches), 2):
    #         male_rect = ax.patches[i]
    #         female_rect = ax.patches[i + 1]
    #         male_width = int(male_rect.get_x())
    #         female_width = int(female_rect.get_width()) - male_width
    #         total = male_width + female_width
    #         ax.annotate(f'{total}', 
    #                     xy=(male_rect.get_x() + male_rect.get_width(), 
    #                         male_rect.get_y() + male_rect.get_height() / 2),
    #                     xytext=(offset, 0),  # 3 points horizontal offset
    #                     textcoords="offset points",
    #                     ha='left', va='center', fontsize=10)
    
    # annotate_bars(ax1_twin, -12)
    # annotate_bars(ax2_twin, 3)
    
    # add legend
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=10, linestyle='') for color, label in zip(["black", "dimgray"], ["Male", "Female"])]
    ax2_twin.legend(handles=legend_handles, loc="lower right", frameon=False)
    
    gs.tight_layout(fig)
    fig.savefig(os.path.join(figure_dir, filename))




def plot_gene_counts_percent_miscalled(df_somatic, df_ch, subset_to_1_perc, figure_dir, filename):
    """
    Each gene is a percentage - percentage of somatic vs CH mutations.
    """
    # Calculate percentage that is CH vs ctDNA for every gene.
    # df_ch = baseline.copy()
    if subset_to_1_perc:
        df_ch = df_ch[df_ch["VAF_t"] > 1].reset_index(drop = True)
    # df_somatic = baseline_somatic.copy()
    df_somatic = df_somatic[["Diagnosis", "Gene"]].sort_values("Gene").value_counts().reset_index()
    df_somatic.columns = ["Diagnosis", "Gene", "Counts_ctDNA"]
    df_ch = df_ch[["Diagnosis", "Gene"]].sort_values("Gene").value_counts().reset_index()
    df_ch.columns = ["Diagnosis", "Gene", "Counts_ch"]
    
    # Get percentage values
    combined = df_ch.merge(df_somatic, on = ["Diagnosis", "Gene"], how = "outer").fillna(0)
    combined["summed"] = combined["Counts_ch"] + combined["Counts_ctDNA"]
    combined["perc_ch"] = (combined["Counts_ch"]/combined["summed"])*100
    combined["perc_ctDNA"] = (combined["Counts_ctDNA"]/combined["summed"])*100
    # Determine gene order
    bladder_df = combined[combined["Diagnosis"] == "Bladder"].sort_values(by = "perc_ch", ascending = False)
    kidney_df = combined[combined["Diagnosis"] == "Kidney"]
    kidney_df = bladder_df["Gene"].reset_index().merge(kidney_df, on = "Gene") # same order as bladder
    kidney_df["perc_ch"] = kidney_df["perc_ch"]*-1
    kidney_df["perc_ctDNA"] = kidney_df["perc_ctDNA"]*-1
    # plotting
    fig = plt.figure(figsize=(7, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], wspace = 0)
    ax1 = plt.subplot(gs[0]) # kidney
    ax2 = plt.subplot(gs[1], sharey = ax1) # bladder
    ax1_twin = ax1.twiny()
    ax2_twin = ax2.twiny()
    
    ax1_twin.barh(kidney_df["Gene"], kidney_df["perc_ch"], color='black', edgecolor = "None")
    ax1_twin.barh(kidney_df["Gene"], kidney_df["perc_ctDNA"], color='dimgray', left = kidney_df["perc_ch"].fillna(0), edgecolor = "None")
    ax2_twin.barh(bladder_df["Gene"], bladder_df["perc_ch"], color='black', edgecolor = "None")
    ax2_twin.barh(bladder_df["Gene"], bladder_df["perc_ctDNA"], color='dimgray', left = bladder_df["perc_ch"], edgecolor = "None")
    
    ax1.spines[["left", "bottom"]].set_visible(False)
    ax2.spines[["bottom", "right"]].set_visible(False)
    ax1_twin.spines[["left", "bottom"]].set_visible(False)
    ax2_twin.spines[["bottom", "right"]].set_visible(False)
    
    ax1_twin.set_ylabel("")
    ax2_twin.set_ylabel("")
    ax1_twin.set_xlabel("Kidney")
    ax2_twin.set_xlabel("Bladder")
    ax1_twin.set_xticks([0, -20, -40, -60, -80, -100])
    ax1_twin.set_xticklabels(["0", "20", "40", "60", "80", "100"])
    ax2_twin.set_xticks([0, 20, 40, 60, 80, 100])
    ax2_twin.set_xticklabels(["0", "20", "40", "60", "80", "100"])  
    
    ax1_twin.set_ylim((-1, 41.5))
    ax2_twin.set_ylim((-1, 41.5))
    ax1.set_ylim((-1, 41.5))
    ax2.set_ylim((-1, 41.5))
    
    # set the y tick labels
    gene_names = bladder_df["Gene"].tolist()
    ax2.set_yticklabels(gene_names, ha='center')
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax2.yaxis.set_tick_params(pad=30)   
    ax1.set_xticks([])
    ax2.set_xticks([])
    
    # add legend
    # legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=10, linestyle='') for color, label in zip(["black", "dimgray"], ["CH", "ctDNA"])]
    # ax2_twin.legend(handles=legend_handles, loc="lower right", frameon=False)
    
    gs.tight_layout(fig)
    fig.savefig(os.path.join(figure_dir, filename))

def plot_gene_counts_percent_miscalled_together(df_somatic, df_ch, subset_to_1_perc, figure_dir, filename):
    """
    Each gene is a percentage - percentage of somatic vs CH mutations.
    Together just means we don't separate bladder and kidney.
    """
    # Calculate percentage that is CH vs ctDNA for every gene.
    df_ch = baseline.copy()
    if subset_to_1_perc:
        df_ch = df_ch[df_ch["VAF_t"] > 1].reset_index(drop = True)
    df_somatic = baseline_somatic.copy()
    df_somatic = df_somatic[["Diagnosis", "Gene"]].sort_values("Gene").value_counts().reset_index()
    df_somatic.columns = ["Diagnosis", "Gene", "Counts_ctDNA"]
    df_ch = df_ch[["Diagnosis", "Gene"]].sort_values("Gene").value_counts().reset_index()
    df_ch.columns = ["Diagnosis", "Gene", "Counts_ch"]
    
    # Get percentage values
    combined = df_ch.merge(df_somatic, on = ["Diagnosis", "Gene"], how = "outer").fillna(0)
    combined["summed"] = combined["Counts_ch"] + combined["Counts_ctDNA"]
    combined["perc_ch"] = (combined["Counts_ch"]/combined["summed"])*100
    combined["perc_ctDNA"] = (combined["Counts_ctDNA"]/combined["summed"])*100
    # Determine gene order
    bladder_df = combined[combined["Diagnosis"] == "Bladder"].sort_values(by = "Counts_ch", ascending = False)
    kidney_df = combined[combined["Diagnosis"] == "Kidney"]
    kidney_df = bladder_df["Gene"].reset_index().merge(kidney_df, on = "Gene") # same order as bladder
    kidney_df["perc_ch"] = kidney_df["perc_ch"]*-1
    kidney_df["perc_ctDNA"] = kidney_df["perc_ctDNA"]*-1
    # plotting
    fig = plt.figure(figsize=(7, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], wspace = 0)
    ax1 = plt.subplot(gs[0]) # kidney
    ax2 = plt.subplot(gs[1], sharey = ax1) # bladder
    ax1_twin = ax1.twiny()
    ax2_twin = ax2.twiny()
    
    ax1_twin.barh(kidney_df["Gene"], kidney_df["perc_ch"], color='black', edgecolor = "None")
    ax1_twin.barh(kidney_df["Gene"], kidney_df["perc_ctDNA"], color='dimgray', left = kidney_df["perc_ch"].fillna(0), edgecolor = "None")
    ax2_twin.barh(bladder_df["Gene"], bladder_df["perc_ch"], color='black', edgecolor = "None")
    ax2_twin.barh(bladder_df["Gene"], bladder_df["perc_ctDNA"], color='dimgray', left = bladder_df["perc_ch"], edgecolor = "None")
    
    ax1.spines[["left", "bottom"]].set_visible(False)
    ax2.spines[["bottom", "right"]].set_visible(False)
    ax1_twin.spines[["left", "bottom"]].set_visible(False)
    ax2_twin.spines[["bottom", "right"]].set_visible(False)
    
    ax1_twin.set_ylabel("")
    ax2_twin.set_ylabel("")
    ax1_twin.set_xlabel("Kidney")
    ax2_twin.set_xlabel("Bladder")
    ax1_twin.set_xticks([0, -20, -40, -60, -80, -100])
    ax1_twin.set_xticklabels(["0", "20", "40", "60", "80", "100"])
    ax2_twin.set_xticks([0, 20, 40, 60, 80, 100])
    ax2_twin.set_xticklabels(["0", "20", "40", "60", "80", "100"])  
    
    ax1_twin.set_ylim((-1, 41.5))
    ax2_twin.set_ylim((-1, 41.5))
    ax1.set_ylim((-1, 41.5))
    ax2.set_ylim((-1, 41.5))
    
    # set the y tick labels
    gene_names = bladder_df["Gene"].tolist()
    ax2.set_yticklabels(gene_names, ha='center')
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax2.yaxis.set_tick_params(pad=30)   
    ax1.set_xticks([])
    ax2.set_xticks([])
    
    # add legend
    # legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=10, linestyle='') for color, label in zip(["black", "dimgray"], ["CH", "ctDNA"])]
    # ax2_twin.legend(handles=legend_handles, loc="lower right", frameon=False)
    
    gs.tight_layout(fig)
    fig.savefig(os.path.join(figure_dir, filename))






    
    
    

def plot_gene_counts_baseline_vs_prog(baseline_df, prog_df, figure_dir, figure_name):
    """
    Plots the number of mutations each individual gene has in the cohort as a stacked BAR CHART.
    Bars are colored to reflect the portion of mutations called in just baseline vs baseline and progression. Only considering pts that have both baseline and prog available.
    No option to do unique here. I mean we could but do we really need it? 
    --merged: baseline and prog merged.
    """
    baseline_df = prog_df[["Patient_id", "Diagnosis"]].drop_duplicates().merge(baseline_df, how = "left")
    combined = baseline_df.merge(prog_df, how = "outer", on = ["Patient_id", "Diagnosis", "Gene", "Chrom", "Position", "Protein_annotation", "Consequence"], indicator = True)
    # Create a pivot table for plotting
    pivot_table = combined.pivot_table(index='Gene', columns='_merge', aggfunc='size', fill_value=0)
    pivot_table.columns = ["Baseline only", "Progression only", "Both"]
    # if pivot_table['right_only'].eq(0).all():
    #     pivot_table = pivot_table.drop(columns='right_only')
    # if pivot_table['left_only'].eq(0).all():
    #     pivot_table = pivot_table.drop(columns='left_only')
    pivot_table = pivot_table.sort_values("Both", ascending = False)
    fig, ax = plt.subplots(figsize = (8, 4))
    pivot_table.plot(kind='bar', stacked=True, colormap='copper', ax=ax)
    ax.set_xlabel('Genes')
    ax.set_ylabel('Number of CH events')
    ax.set_title('CH events in baseline and progression samples\nSome samples may have more than one progression sample.\nOnly showing patients where at least one progression sample was available.')
    ax.legend(title='')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # ax.set_yticks([0, 5, 10, 15, 20, 25])
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
    fig.tight_layout()
    fig.savefig(os.path.join(figure_dir, figure_name))

def make_pie_chart(df, PATH_sample_info, figure_dir, filename):
    """
    Given a df with mutations (baseline or progression) and a list of patients in the cohort, generates a pie chart showing the percentage of the cohort CH positive for each given gene.
    """
    baseline_info = pd.read_csv(PATH_sample_info, "\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"]).query('Timepoint == "Baseline"').reset_index(drop = True)
    all_pts = baseline_info.assign(Type="all_pts")[["Patient_id", "Type"]]
    df = df.assign(Type="CH_positive_pts")[["Patient_id", "Type", "Gene"]]
    combined = all_pts.merge(df, how = "left", on = "Patient_id", indicator = True)
    combined = combined[combined["_merge"] == "both"]
    gene_counts = combined["Gene"].value_counts()
    # PLOTTING
    fig, ax = plt.subplots()
    ax.pie(gene_counts, labels=gene_counts.index, autopct='%1.1f%%', startangle=140)
    ax.set_title('Distribution of Genes in the DataFrame')
    fig.savefig(os.path.join(figure_dir, filename))

def plot_pie_for_ch_presence_in_df(df, PATH_sample_info, figure_dir, figure_title, figure_name, color_dict = None, fontsize = 10):
    """
    Plots a simple pie chart to show the percentage of the cohort CH+. Can work with both baseline and progression data frames.
    """
    timepoint = df["Timepoint"].unique()[0]
    diagnosis = df["Diagnosis"].unique()[0]
    sample_df = pd.read_csv(PATH_sample_info, sep="\t", header=None, names = ["Patient_id", "Date", "Diagnosis", "Timepoint"])
    sample_df = sample_df[(sample_df["Timepoint"] == timepoint) & (sample_df["Diagnosis"] == diagnosis)].reset_index(drop = True)
    
    # Get unique patients from sample information and patients with CH+ from the given DataFrame
    all_pts = sample_df["Patient_id"].unique()
    ch_pts = df["Patient_id"].unique()
    
    # Count CH+ and CH- patients
    ch_positive_count = len(set(ch_pts).intersection(all_pts))
    ch_negative_count = len(all_pts) - ch_positive_count
    
    # Create a pie chart using ax
    fig, ax = plt.subplots(figsize=(1, 1))
    labels = ['CH+', 'CH-']
    sizes = [ch_positive_count, ch_negative_count]
    if color_dict is None: 
        colors = ['lightcoral', 'lightblue']
    else:
        colors = color_dict
    ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=colors)
    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    ax.set_title(figure_title)
    
    # Save the figure
    figure_path = os.path.join(figure_dir, figure_name)
    fig.tight_layout()
    plt.savefig(figure_path)

def plot_per_patient_counts(df, figure_dir, figure_name, colorby=None, ax = None, bar_color = None):
    """
    Plots a simple bar chart indicating the number of patients having n number of mutations.
    y axis is the number of patients. 
    x axis is the number of mutations.
    If colorby is given as "diagnosis", plot a grouped bar chart separated by diagnosis.
    """
    timepoint = df["Timepoint"].unique()[0]
    diagnosis = df["Diagnosis"].unique()[0]
    if "Status_n" in df.columns: 
        keyword = "CH"
    elif "Status" in df.columns:
        keyword = "ctDNA"
    
    # Plotting stacked bar chart if colorby is "diagnosis"
    if colorby == "Diagnosis":
        # Bladder
        count_table_bladder = df[df["Diagnosis"] == "Bladder"].groupby("Patient_id").size().value_counts().reset_index()
        count_table_bladder.columns = ["Number of mutations", "Number of patients bladder"]
        count_table_bladder = count_table_bladder.sort_values(by="Number of mutations")
        # Kidney
        count_table_kidney = df[df["Diagnosis"] == "Kidney"].groupby("Patient_id").size().value_counts().reset_index()
        count_table_kidney.columns = ["Number of mutations", "Number of patients kidney"]
        count_table_kidney = count_table_kidney.sort_values(by="Number of mutations")
        merged = count_table_bladder.merge(count_table_kidney, how = "outer")
        merged.set_index('Number of mutations', inplace=True)
        merged["Number of patients kidney"] = merged["Number of patients kidney"].fillna(0)
        merged["Number of patients bladder"] = merged["Number of patients bladder"].fillna(0)
        merged = merged.rename(columns = {"Number of patients kidney": "Kidney", "Number of patients bladder": "Bladder"})
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(3.5, 3))
        else: 
            merged.plot(kind='bar', stacked=True, color=['deepskyblue', 'orangered'], ax=ax, edgecolor = None, linewidth=0.0)
            # bars1 = ax.bar(count_table_bladder["Number of mutations"], count_table_bladder["Number of patients bladder"], color='deepskyblue', label='Bladder', edgecolor = None, linewidth=0.0)
            # bars2 = ax.bar(count_table_kidney["Number of mutations"], count_table_kidney["Number of patients kidney"], bottom=merged["Bladder"], color='orangered', label='Kidney', edgecolor = None, linewidth=0.0)
            ax.legend(frameon = False)
            # ax.set_xticks(range(0, max(merged.index.tolist())))
            ax.set_xticklabels(merged.index.tolist(), rotation = 0, fontsize = 8)
    else:
        count_table = df.groupby("Patient_id").size().value_counts().reset_index()
        count_table.columns = ["Number of mutations", "Number of patients"]
        count_table = count_table.sort_values(by="Number of mutations")
        if ax is None:
            fig, ax = plt.subplots()
        if bar_color is None:
            ax.bar(count_table["Number of mutations"], count_table["Number of patients"], color="black")
        else: 
            ax.bar(count_table["Number of mutations"], count_table["Number of patients"], color=bar_color)
        ax.set_title(diagnosis, loc = "left")
    # Aesthetics
    ax.set_xticks(range(1, max(count_table["Number of mutations"])+1))
    ax.set_xticklabels(range(1, max(count_table["Number of mutations"])+1), fontsize = 6, rotation = 0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Number of CH mutations")
    ax.set_ylabel("Number of patients")
    ax.tick_params(axis="both", direction="out", which="both", left=True, bottom=True , colors='k')
    
    # y tick labels
    # max_y = max(ax.get_yticks())
    # new_max_y = round_to_nearest_10(max_y)
    # ax.set_ylim((0, new_max_y))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))  # Adjust the number of ticks as needed
    
    # Save the figure
    if ax is not None:
        return(ax)
    else:
        figure_path = os.path.join(figure_dir, figure_name)
        fig.tight_layout()
        fig.savefig(figure_path)

        

def plot_ctDNA_bins_and_VAF_correlations(df, ctDNA_fractions_path, dir_figures, ext = "png"):
    ct_frac = pd.read_csv(PATH_mutation_ctfractions)[["Patient_id", "Mutation_ctDNA_fraction"]]
    df = df.merge(ct_frac, how = "left")
    df["Mutation_ctDNA_fraction"] = df["Mutation_ctDNA_fraction"]*100
    
    bins = [0, 2, 30, 100]
    labels = ['0-2', '2-30', '30-100']
    df['Mutation_ctDNA_bin'] = pd.cut(df['Mutation_ctDNA_fraction'], bins=bins, labels=labels, right=False)
    
    # Correlate WBC VAf with cfDNA VAF in each bin
    correlations = {}
    for label in labels:
        bin_data = df[df['Mutation_ctDNA_bin'] == label] # Filter the DataFrame for the current bin
        correlation = bin_data['VAF_t'].corr(bin_data['VAF_n']) # Calculate Pearson correlation between 'VAF_t' and 'VAF_n' in the current bin
        correlations[label] = correlation # Store the correlation in the dictionary
    
    # Set up the plot
    fig = plt.figure(figsize=(3, 3))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
    
    ax1 = plt.subplot(gs[0]) # dotplot
    ax2 = plt.subplot(gs[1]) # boxplot
    
    # Plot the connected line plot for the pearson
    ax1.scatter(labels, correlations.values(), s=70, color = "black")
    ax1.plot(labels, correlations.values(), color='black', linewidth=1)
    ax1.set_ylabel("Pearson R")
    ax1.set_xlabel("")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_yticks([1, 0.98, 0.96, 0.94, 0.92])
    ax1.set_yticklabels(["1", "0.98", "0.96", "0.94", "0.92"])
    ax1.set_xticks([])
    ax1.tick_params(axis="both", direction="out", which="both", left=True, bottom=False, colors='k')
    
    # Plot the boxplots with scatters
    sns.boxplot(x='Mutation_ctDNA_bin', y="Mutation_ctDNA_fraction", data=df, order=['0-2', '2-30', '30-100'], ax=ax2, color='black', medianprops={'color': 'black'}, fliersize = 0, boxprops={"alpha": .2, "edgecolor": 'black'})
    sns.stripplot(x='Mutation_ctDNA_bin', y="Mutation_ctDNA_fraction", data=df, order=['0-2', '2-30', '30-100'], ax=ax2, jitter=True, size = 3, color='black')
    ax2.set_ylabel("ctDNA%")
    ax2.set_xlabel("Binned ctDNA fractions")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.set_yticks([0, 25, 50, 75, 100])
    ax2.set_yticklabels(["0", "25", "50", "75", "100"])
    ax2.tick_params(axis="both", direction="out", which="both", left=True, bottom=True, colors='k')
    fig.tight_layout()
    fig.savefig(os.path.join(dir_figures, f"binned_VAF_correlation_plot.{ext}"))
    
    

def plot_vaf_scatter(df, figure_dir, figure_name, color_dict = None, colorby = None, fontsize = 10):
    """
    Plots a scatter plot comparing the cfDNA VAF to WBC VAF.
    """  
    # helps remove dependent calls with vaf 0 for both wbc and cfdna
    df = df[df["Dependent"] == False]
    df = df.sort_values(by="VAF_n", ascending=False)
    df2 = df[(df["VAF_n"] <= 2) & (df["VAF_t"] <= 2)]
        
    # Plotting
    fig = plt.figure(figsize=(6, 2.8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])  # 1 row, 2 columns, with same widths
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    # Scatter plots
    if colorby is None:
        sns.regplot(x="VAF_n", y="VAF_t", data=df, ax=ax1, scatter_kws={'s': 9, 'color': 'black'}, line_kws={'color': 'orange', 'linewidth': 1, "linestyle": "--"})
        sns.scatterplot(x="VAF_n", y="VAF_t", data=df2, ax=ax2, s=9)
    elif colorby == "Diagnosis":
        sns.regplot(x="VAF_n", y="VAF_t", data=df, ax=ax1, scatter = False, line_kws={'color': 'black', 'linewidth': 1, "linestyle": "--"})
        sns.scatterplot(x="VAF_n", y="VAF_t", data=df, ax=ax1, s=9, hue="Diagnosis", edgecolor = None, palette = color_dict)
        sns.scatterplot(x="VAF_n", y="VAF_t", data=df2, ax=ax2, s=9, hue="Diagnosis", edgecolor = None, palette = color_dict)
    # Adding y = x line
    # ax1.plot([0, 100], [0, 100], linestyle='--', color='blue', label='y = x')  # y = x line
    # ax1.annotate('x = y', xy=(65, 90), xytext=(65, 90), color='blue', fontsize=10)
    
    # ax2.plot([0, 2.2], [0, 2.2], linestyle='--', color='blue', label='y = x')  # y = x line
    # ax2.annotate('y = x', xy=(60, 80), xytext=(60, 80), color='blue', fontsize=10)
    
    spearman_corr = df[["VAF_n", "VAF_t"]].corr(method="spearman").iloc[0, 1]
    ax1.text(0.05, 0.9, f'Spearman\'s\ncorrelation: {spearman_corr:.2f}', transform=ax1.transAxes, fontsize=fontsize, color='black')
    ax2.set_xlim((0, 2.1))
    ax2.set_ylim((0, 2.1))
    ax2.set_xticks([0, 2])
    ax2.set_yticks([0, 2])
    ax2.set_xticklabels(["0", "2"], fontsize=fontsize)
    ax2.set_yticklabels(["0", "2"], fontsize=fontsize)
    for ax in [ax1, ax2]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    ax1.set_xlabel("WBC VAF%", fontsize=fontsize)
    ax1.set_ylabel("Plasma cfDNA VAF%", fontsize=fontsize)
    ax2.set_xlabel("")
    ax2.set_ylabel("")
    # Save the figure
    figure_path = os.path.join(figure_dir, figure_name)
    fig.suptitle(f"")
    ax1.tick_params(axis="both", direction="out", which="both", left=True, bottom=True, colors='k')
    ax2.tick_params(axis="both", direction="out", which="both", left=True, bottom=True, colors='k')
    ax1.legend().set_visible(False)
    ax2.legend().set_visible(False)
    
    # Draw the highlighted bit at the bottom left
    # ax1.plot([2, 2], [0, 2], color='r', linestyle='--', linewidth = 0.5)
    # ax1.plot([0, 2], [2, 2], color='r', linestyle='--', linewidth = 0.5)
    
    ax1.set_yticks([0, 20, 40, 60, 80])
    ax1.set_yticklabels(["0", "20", "40", "60", "80"], fontsize = fontsize)
    ax1.set_xticks([0, 20, 40, 60, 80])
    ax1.set_xticklabels(["0", "20", "40", "60", "80"], fontsize = fontsize)
    fig.tight_layout()
    fig.savefig(figure_path)
    
def plot_boxplot(df, numerical_col_to_plot, categorical_col_to_plot, color_dict_dots, color_dict_boxes, figure_dir, figure_title):
    """
    color_dict_dots: Dict. Keys should be values in the categorical_col_to_plot.
    color_dict_boxes: Dict. Keys should be values in the categorical_col_to_plot.
    """
    fig, ax = plt.subplots(figsize=(3, 2.8))
    sns.boxplot(x=categorical_col_to_plot, y=numerical_col_to_plot, data=df, ax=ax, medianprops={'color': 'black'}, fliersize = 0, boxprops={"alpha": 1, "edgecolor": 'black'})
    sns.stripplot(x=categorical_col_to_plot, y=numerical_col_to_plot, data=df, ax=ax, jitter=True, size = 4)
    # Customize the plot
    ax.set_xlabel(categorical_col_to_plot)
    ax.set_ylabel(numerical_col_to_plot)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(figure_title)
    figure_path = os.path.join(figure_dir, figure_title + ".pdf")
    # Statistical test
    value1 = df[categorical_col_to_plot].unique()[0]
    value2 = df[categorical_col_to_plot].unique()[1]
    g1 = df[df[categorical_col_to_plot] == value1][numerical_col_to_plot].tolist()
    g2 = df[df[categorical_col_to_plot] == value2][numerical_col_to_plot].tolist()
    stat, p_value = mannwhitneyu(g1, g2)
    formatted_p_value = "{:.2e}".format(p_value)
    # Display p-value on the plot
    # ax.text(0.5, max(combined[numerical_col_to_plot]) + 1, f'P-Value: {formatted_p_value}', ha='center', va='center', backgroundcolor='white')
    
    # display the numerical values too
    g1_median = df[df[categorical_col_to_plot] == value1][numerical_col_to_plot].median()
    g2_median = df[df[categorical_col_to_plot] == value2][numerical_col_to_plot].median()
    ax.text(0.25, g1_median, f'{round(g1_median)}', ha='center', va='bottom', color='black', fontsize=8)
    ax.text(1.25, g2_median, f'{round(g2_median)}', ha='center', va='bottom', color='black', fontsize=8)
    ax.tick_params(axis="both", direction="out", which="both", left=True, bottom=True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.suptitle(f'P-Value: {formatted_p_value}', fontsize=8, y=0.89)
    fig.savefig(figure_path)
    print(f"Test Statistic: {stat}, P-Value: {formatted_p_value}")


def age_vs_CH_presence(df_CH, numerical_col_to_plot, PATH_sample_information, clin, figure_title, figure_dir, gene = None, min_WBC_VAF_threshold = None, test_to_use = "MWU", fontsize = 10):
    """
    Boxplot. Plotting the ages on the y-axis of patients, grouped into CH+ and CH-. Give a df consisting of only one timepoint, baseline or OT. 
    test_to_use: use either 'MWU' or 'T test'
    """
    # Data manipulation
    diagnosis = df_CH["Diagnosis"].unique()[0]
    timepoint = df_CH["Timepoint"].unique()[0]
    all_pts = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
    all_pts = all_pts[(all_pts["Diagnosis"] == diagnosis) & (all_pts["Timepoint"] == timepoint)].reset_index(drop = True)
    # filtering
    if gene is None:
        plotting_df = df_CH.copy()
    else:
        plotting_df = df_CH[df_CH["Gene"] == gene]
    plotting_df["CH_status"] = "CH+"
    combined = all_pts.merge(plotting_df[["Patient_id", "CH_status", "VAF_n"]], how="left")
    combined["CH_status"] = combined["CH_status"].fillna("CH-")
    combined = combined.merge(clin[["Patient_id", numerical_col_to_plot]], how="left", on="Patient_id")
    combined = combined[~pd.isnull(combined[numerical_col_to_plot])]
    if min_WBC_VAF_threshold is not None:
        ch_pos_df = combined[(combined["CH_status"] == "CH+") & (combined["VAF_n"] > min_WBC_VAF_threshold)]
        ch_neg_df = combined[combined["CH_status"] == "CH-"]
        combined = pd.concat([ch_pos_df, ch_neg_df]).reset_index(drop = True)
    del combined["VAF_n"]    
    combined = combined.drop_duplicates().reset_index(drop = True)
    # Plotting
    fig, ax = plt.subplots(figsize=(3, 2.8))
    sns.boxplot(x='CH_status', y=numerical_col_to_plot, data=combined, order=['CH-', 'CH+'], ax=ax, color='black', medianprops={'color': 'black'}, palette={'CH-': '#4D93C3', 'CH+': '#FF7F0E'}, fliersize = 0, boxprops={"alpha": .5, "edgecolor": 'black'})
    sns.stripplot(x='CH_status', y=numerical_col_to_plot, data=combined, order=['CH-', 'CH+'], ax=ax, jitter=True, size = 4, color='black', palette={'CH-': '#4D93C3', 'CH+': '#FF7F0E'})
    # Customize the plot
    ax.set_xlabel('CH status')
    ax.set_ylabel(numerical_col_to_plot)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(figure_title)
    figure_path = os.path.join(figure_dir, figure_title + "_boxplot_scatter.pdf")
    # Statistical test
    ch_pos_values = combined[combined['CH_status'] == 'CH+'][numerical_col_to_plot].tolist()
    ch_neg_values = combined[combined['CH_status'] == 'CH-'][numerical_col_to_plot].tolist()
    if test_to_use == "MWU":
        stat, p_value = mannwhitneyu(ch_pos_values, ch_neg_values)
    else:
        stat, p_value = ttest_ind(ch_pos_values, ch_neg_values)
    formatted_p_value = "{:.2e}".format(p_value)
    # Display p-value on the plot
    # ax.text(0.5, max(combined[numerical_col_to_plot]) + 1, f'P-Value: {formatted_p_value}', ha='center', va='center', backgroundcolor='white')
    
    # display the median ages too
    median_ch_pos = combined[combined['CH_status'] == 'CH+'][numerical_col_to_plot].median()
    median_ch_neg = combined[combined['CH_status'] == 'CH-'][numerical_col_to_plot].median()
    ax.text(0.25, median_ch_neg, f'{round(median_ch_neg)}', ha='center', va='bottom', color='black', fontsize=fontsize)
    ax.text(1.25, median_ch_pos, f'{round(median_ch_pos)}', ha='center', va='bottom', color='black', fontsize=fontsize)
    ax.tick_params(axis="both", direction="out", which="both", left=True, bottom=True)
    fig.suptitle(f'P-Value: {formatted_p_value}', fontsize=fontsize, y=0.89)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(figure_path)
    print(f"Test Statistic: {stat}, P-Value: {formatted_p_value}")

def age_vs_CH_presence_grouped(df_CH, age_df, numerical_col_to_plot, PATH_sample_information, figure_title, figure_dir, gene = None, min_WBC_VAF_threshold = None, test_to_use = "MWU", fontsize = 10):
    """
    VERY SIMILAR TO age_vs_CH_presence. ONLY DIFFERENCE IS PLOTS BLADDER AND KIDNEY TOGETHER.
    Boxplot. Plotting the ages on the y-axis of patients, grouped into CH+ and CH-. Give a df consisting of only one timepoint, baseline or OT. 
    test_to_use: use either 'MWU' or 'T test'
    """
    # Data manipulation
    timepoint = df_CH["Timepoint"].unique()[0]
    all_pts = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
    all_pts = all_pts[(all_pts["Timepoint"] == timepoint)].reset_index(drop = True)
    # filtering
    if gene is None:
        plotting_df = df_CH.copy()
    else:
        plotting_df = df_CH[df_CH["Gene"] == gene]
    plotting_df["CH_status"] = "CH+"
    combined = all_pts.merge(plotting_df[["Patient_id", "CH_status", "VAF_n"]], how="left")
    combined["CH_status"] = combined["CH_status"].fillna("CH-")
    combined = combined.merge(age_df, how = "left")
    combined = combined[~pd.isnull(combined[numerical_col_to_plot])]
    if min_WBC_VAF_threshold is not None:
        ch_pos_df = combined[(combined["CH_status"] == "CH+") & (combined["VAF_n"] > min_WBC_VAF_threshold)]
        ch_neg_df = combined[combined["CH_status"] == "CH-"]
        combined = pd.concat([ch_pos_df, ch_neg_df]).reset_index(drop = True)
    del combined["VAF_n"]    
    combined = combined.drop_duplicates().reset_index(drop = True)
    combined["group"] = "dummy_value"
    combined.loc[(combined["Diagnosis"] == "Bladder") & (combined["CH_status"] == "CH+"), "group"] = "Bladder CH+"
    combined.loc[(combined["Diagnosis"] == "Bladder") & (combined["CH_status"] == "CH-"), "group"] = "Bladder CH-"
    combined.loc[(combined["Diagnosis"] == "Kidney") & (combined["CH_status"] == "CH+"), "group"] = "Kidney CH+"
    combined.loc[(combined["Diagnosis"] == "Kidney") & (combined["CH_status"] == "CH-"), "group"] = "Kidney CH-"    
    
    fig, ax = plt.subplots(figsize=(3, 3))
    sns.boxplot(x="group", y=numerical_col_to_plot, order = ["Bladder CH+", "Bladder CH-", "Kidney CH+", "Kidney CH-"], data=combined, ax=ax, medianprops={'color': 'black'}, fliersize = 0, whiskerprops={'color':'black'}, capprops={'color':'black'}, boxprops={"edgecolor": 'black'}, hue = "Diagnosis", palette = {"Bladder": (143/255, 215/255, 239/255, 255/255), "Kidney": (239/255, 169/255, 143/255, 255/255)}, linewidth=0.5, width=0.7, dodge = False)
    sns.stripplot(x="group", y=numerical_col_to_plot, order = ["Bladder CH+", "Bladder CH-", "Kidney CH+", "Kidney CH-"], data=combined, ax=ax, jitter=True, size = 3, hue = "Diagnosis", palette = color_dict)
    # Customize the plot
    ax.set_xlabel('CH status')
    ax.set_ylabel("Age at baseline blood draw")
    ax.set_xticklabels(["Bladder\nCH+", "Bladder\nCH-", "Kidney\nCH+", "Kidney\nCH-"], fontsize = fontsize -2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("")
    ax.legend().set_visible(False)
    figure_path = os.path.join(figure_dir, figure_title + "_boxplot_scatter.pdf")
    # Statistical test
    ch_pos_values_bladder = combined[(combined['CH_status'] == 'CH+') & (combined["Diagnosis"] == "Bladder")][numerical_col_to_plot].tolist()
    ch_neg_values_bladder = combined[(combined['CH_status'] == 'CH-') & (combined["Diagnosis"] == "Bladder")][numerical_col_to_plot].tolist()
    stat_bl, p_value_bl = mannwhitneyu(ch_pos_values_bladder, ch_neg_values_bladder)
    formatted_p_value_bladder = "{:.2e}".format(p_value_bl)
    
    ch_pos_values_kidney = combined[(combined['CH_status'] == 'CH+') & (combined["Diagnosis"] == "Kidney")][numerical_col_to_plot].tolist()
    ch_neg_values_kidney = combined[(combined['CH_status'] == 'CH-') & (combined["Diagnosis"] == "Kidney")][numerical_col_to_plot].tolist()
    stat_ki, p_value_ki = mannwhitneyu(ch_pos_values_kidney, ch_neg_values_kidney)
    formatted_p_value_kidney = "{:.2e}".format(p_value_ki)
    
    # Display p-value on the plot
    ax.text(0.5, max(combined[numerical_col_to_plot]) + 3, f'p={formatted_p_value_bladder}', ha='center', va='center', fontsize = fontsize - 2)
    ax.text(2.5, max(combined[numerical_col_to_plot]) + 3, f'p={formatted_p_value_kidney}', ha='center', va='center', fontsize = fontsize - 2)
    
    # display the median ages too
    # BLADDER
    median_ch_pos_bladder = combined[(combined['CH_status'] == 'CH+') & (combined["Diagnosis"] == "Bladder")][numerical_col_to_plot].median()
    median_ch_neg_bladder = combined[(combined['CH_status'] == 'CH-') & (combined["Diagnosis"] == "Bladder")][numerical_col_to_plot].median()
    ax.text(0.25, median_ch_pos_bladder, f'{round(median_ch_pos_bladder)}', ha='center', va='bottom', color='black', fontsize=fontsize-3)
    ax.text(1.25, median_ch_neg_bladder, f'{round(median_ch_neg_bladder)}', ha='center', va='bottom', color='black', fontsize=fontsize-3)
    # KIDNEY
    median_ch_pos_kidney = combined[(combined['CH_status'] == 'CH+') & (combined["Diagnosis"] == "Kidney")][numerical_col_to_plot].median()
    median_ch_neg_kidney = combined[(combined['CH_status'] == 'CH-') & (combined["Diagnosis"] == "Kidney")][numerical_col_to_plot].median()
    ax.text(2.25, median_ch_pos_kidney, f'{round(median_ch_pos_kidney)}', ha='center', va='bottom', color='black', fontsize=fontsize-3)
    ax.text(3.25, median_ch_neg_kidney, f'{round(median_ch_neg_kidney)}', ha='center', va='bottom', color='black', fontsize=fontsize-3)
    ax.tick_params(axis="both", direction="out", which="both", left=True, bottom=True)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    # fig.suptitle(f'P-Value: {formatted_p_value}', fontsize=fontsize, y=0.89)
    fig.savefig(figure_path)
    # print(f"Test Statistic: {stat}, P-Value: {formatted_p_value}")

def categorical_vs_CH_vaf(df_CH, categorical_col_to_plot, clin, figure_title, figure_dir):
    """
    Boxplot. Separates the patients by sex and plots the highest VAF per patient on the y axis.
    """
    # Data manipulation
    diagnosis = df_CH["Diagnosis"].unique()[0]
    all_pts = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
    all_pts = all_pts[all_pts["Diagnosis"] == diagnosis].reset_index(drop = True)
    combined = all_pts.merge(df_CH[["Patient_id", "VAF_n"]], how="left")
    combined = combined.merge(clin[["Patient_id", categorical_col_to_plot]], how="left", on="Patient_id")
    combined = combined[~pd.isnull(combined[categorical_col_to_plot])]
    combined["VAF_n"] = combined["VAF_n"].fillna(0)
    combined = combined.loc[combined.groupby("Patient_id")["VAF_n"].idxmax()].reset_index(drop = True)
    combined = combined.drop_duplicates().reset_index(drop = True)
    # figure out what cetagories we are plotting
    category1 = combined[categorical_col_to_plot].unique()[0]
    category2 = combined[categorical_col_to_plot].unique()[1]   
    # Plotting
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.boxplot(x=categorical_col_to_plot, y="VAF_n", data=combined, order = [category1, category2], ax=ax, color='lightgray', medianprops={'color': 'black'})
    sns.stripplot(x=categorical_col_to_plot, y="VAF_n", data=combined, order = [category1, category2], ax=ax, jitter=True, color='black', alpha=0.7)
    # Customize the plot
    ax.set_xlabel(categorical_col_to_plot)
    ax.set_ylabel("WBC VAF %")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(figure_title)
    figure_path = os.path.join(figure_dir, figure_title + "_boxplot_scatter.png")
    # Statistical test
    cat1 = combined[combined[categorical_col_to_plot] == category1]["VAF_n"].tolist()
    cat2 = combined[combined[categorical_col_to_plot] == category2]["VAF_n"].tolist()
    t_stat, p_value = ttest_ind(cat1, cat2)
    formatted_p_value = "{:.2e}".format(p_value)
    # Display p-value on the plot
    ax.text(0.5, max(combined["VAF_n"]) + 1, f'P-Value: {formatted_p_value}', ha='center', va='center', backgroundcolor='white')
    
    # display the median ages too
    median_cat1 = combined[combined[categorical_col_to_plot] == category1]["VAF_n"].median()
    median_cat2 = combined[combined[categorical_col_to_plot] == category2]["VAF_n"].median()
    ax.text(0.25, median_cat1, f'{round(median_cat1, 3)}', ha='center', va='bottom', color='red', fontsize=12)
    ax.text(1.25, median_cat2, f'{round(median_cat2, 3)}', ha='center', va='bottom', color='red', fontsize=12)
    fig.tight_layout()
    plt.savefig(figure_path)
    print(f"T-Statistic: {t_stat}, P-Value: {formatted_p_value}")
    # ax.set_title(f"{figure_title}\nP-Value: {p_value}")


def bar_chart_effects(PATH_baseline_chip): 
    """

    """
    # Count the occurrences of each combination of 'Gene' and 'Consequence'
    counts = df.groupby(['Gene', 'Consequence']).size().reset_index(name='Count')

    # Pivot the DataFrame for plotting
    pivot_df = counts.pivot(index='Gene', columns='Consequence', values='Count').fillna(0)

    # Order rows (genes) by the sum of counts in descending order
    row_order = pivot_df.sum(axis=1).sort_values(ascending=False).index
    pivot_df = pivot_df.loc[row_order]

    sns.set(style="white")
    ax = pivot_df.plot(kind='bar', stacked=True, colormap='viridis', figsize=(10, 6))

    # Adding labels and title
    ax.set_xlabel("Gene")
    ax.set_ylabel("Count")
    ax.set_title("Mutations effects")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 0)

    plt.tight_layout()
    plt.savefig('/groups/wyattgrp/users/amunzur/therapy_chip/results/figures/mutation_effects_and_genes.png')



# def mut_counts_per_patient(PATH_baseline_chip): 
#     """
#     """
#     df = pd.read_csv(PATH_baseline_chip).rename(columns={"Sample_name_t": "Sample_name", "VAF_t": "VAF"})["Patient_id"]
#     df["Patient_id"].table()
# # Assuming df is your DataFrame
# df = pd.read_csv(PATH_baseline_chip).rename(columns={"Sample_name_t": "Sample_name", "VAF_t": "VAF"})["Patient_id"]

# # Count the number of mutations per patient
# mutation_counts = df.value_counts()

# # Plotting using pandas
# fig, ax = plt.subplots(figsize=(6, 6))
# # Plotting using pandas
# mutation_counts.value_counts().sort_index().plot(kind='bar', ax=ax, color='skyblue')

# # Adding labels and title
# ax.set_xlabel("Number of Mutations per Patient")
# ax.set_ylabel("Number of Patients")
# ax.set_title("Distribution of Mutations per Patient")
# ax.set_xticklabels(ax.get_xticklabels(), rotation = 0)
# ax.spines["top"].set_visible(False)
# ax.spines["right"].set_visible(False)

# plt.tight_layout()
# plt.savefig("/groups/wyattgrp/users/amunzur/therapy_chip/results/figures/number of muts per pt.png")
    
def do_survival_analysis(PATH_clinical, PATH_sample_information, mutations_df, annotate_what, figure_name, plot_title, annotate_gene = False, diagnosis = "Bladder", figure_dir = "/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/survival_analysis"):
    """
    PATH_clinical: path to clinical data. Should have these columns: Patient_id, Date of last follow-up or death and Death at last follow up.
    """
    surv_df = pd.read_csv(PATH_clinical)[["Patient_id", 'Date of last follow-up or death', 'Death at last follow up', 'Age at last follow-up or death']]
    pts = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    if diagnosis == "Both":
        pts = pts[pts["Timepoint"] == "Baseline"].reset_index(drop = True)
    else: 
        pts = pts[(pts["Diagnosis"] == diagnosis) & (pts["Timepoint"] == "Baseline")].reset_index(drop = True)
    surv_df = surv_df.merge(pts, how = "inner")
    surv_df["Date of last follow-up or death"] = pd.to_datetime(surv_df["Date of last follow-up or death"])
    surv_df["Date_collected"] = pd.to_datetime(surv_df["Date_collected"], format='%Y%b%d', errors='coerce')
    surv_df["Overall survival"] = (surv_df["Date of last follow-up or death"] - surv_df["Date_collected"]).astype('timedelta64[M]')
    surv_df["Overall survival"] = surv_df["Overall survival"].round().astype(int)
    
    base_mut_status = annotate_mutation_status(mutations_df, diagnosis, PATH_sample_information, annotate_what = annotate_what, annotate_gene = annotate_gene)
    base_mut_status = base_mut_status[base_mut_status["Timepoint"] == "Baseline"]
    surv_df_merged = base_mut_status.merge(surv_df[["Patient_id", "Overall survival", "Death at last follow up", "Age at last follow-up or death"]].drop_duplicates())
    # if diagnosis == "Bladder":
    surv_df_merged["Death at last follow up"] = surv_df_merged["Death at last follow up"].map({True: True, False: False, 'True': True, 'False': False, 'Lost to follow-up': False}) # Map values to boolean, treating 'Lost to follow-up' as False
    
    # exclude if all neg or all pos
    if (surv_df_merged[annotate_what + " status"].str.contains('Negative').any() and surv_df_merged[annotate_what + " status"].str.contains('Positive').any()):
        if diagnosis == "Bladder":
            output_path = os.path.join(figure_dir, figure_name)
        elif diagnosis == "Kidney":
            output_path = os.path.join(figure_dir, figure_name)
        else: 
            output_path = os.path.join(figure_dir, figure_name)
        make_survival_curve(surv_df_merged, annotate_what + " status", output_path, plot_title)
        return(output_path)

def make_survival_curve(surv_df_merged, stratify_by, output_path, plot_title):
    """
    stratify_by = Name of the column in df to stratify the curves by.
    """
    label_keyword = re.search("(ctDNA|CHIP)", stratify_by, flags=re.IGNORECASE).group()
    kmf_positive = KaplanMeierFitter()
    kmf_negative = KaplanMeierFitter()
    groups = surv_df_merged[stratify_by]
    ix_positive = (groups == 'Positive')
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Negative group
    kmf_negative.fit(durations=surv_df_merged["Overall survival"][~ix_positive], event_observed=surv_df_merged["Death at last follow up"][~ix_positive], label=f"{label_keyword}(-)")
    kmf_negative.plot_survival_function(ax=ax, ci_show=False, show_censors=True)
    
    # Calculate median survival for negative group and show on the plot
    median_survival_negative = kmf_negative.median_survival_time_
    # ax.annotate(f'Median survival {label_keyword}(-): {median_survival_negative} months', xy=(28, 0.5), xytext=(100, 60), textcoords='offset points', fontsize=10)
    ax.text(0.8, 0.7, f'Median survival {label_keyword}(-): {median_survival_negative} months', transform=ax.transAxes, fontsize=10, ha='center', va='center')
    
    # Positive group
    kmf_positive.fit(durations=surv_df_merged["Overall survival"][ix_positive], event_observed=surv_df_merged["Death at last follow up"][ix_positive], label=f"{label_keyword}(+)")
    kmf_positive.plot_survival_function(ax=ax, ci_show=False, show_censors=True)
    
    # add risk counts
    add_at_risk_counts(kmf_positive, kmf_negative, ax=ax)
    
    # Calculate median survival for positive group
    median_survival_positive = kmf_positive.median_survival_time_
    ax.text(0.8, 0.75, f'Median survival {label_keyword}(+): {median_survival_positive} months', transform=ax.transAxes, fontsize=10, ha='center', va='center')
    
    # check if the median survivals are significantly different using log rank test
    results = logrank_test(surv_df_merged["Overall survival"][~ix_positive], surv_df_merged["Overall survival"][ix_positive], event_observed_A=surv_df_merged["Death at last follow up"][~ix_positive], event_observed_B=surv_df_merged["Death at last follow up"][ix_positive])
    p_value_logrank = round(results.p_value, 4)
    ax.text(0.8, 0.64, f'Log rank p = {p_value_logrank}', transform=ax.transAxes, fontsize=10, ha='center', va='center')
    # ax.annotate(f'p value: {p_value}', xy=(1, 0.5), xytext=(100, 80), textcoords='offset points', fontsize=12)
    ax.set_xlabel("Months")
    ax.set_ylabel("Survival Probability")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc='best', frameon=False)
    ax.set_title(plot_title)
    
    # Run CPH here and print the hazard ratios to the plot.
    try:
        run_cox_proportional_hazards(surv_df_merged, "Overall survival", "Death at last follow up", stratify_by, ax)
    except Exception as e:
        print(f"An error occurred while running Cox Proportional Hazards.")
        
    fig.tight_layout()
    fig.savefig(output_path)
    return median_survival_negative, median_survival_positive

def run_cox_proportional_hazards(surv_df_merged, duration_col, event_col, stratify_by, ax):
    """
    Fits the CPH model to survival data.
    duration_col: Column that has OS data.
    event_col: Column that has indicates if the event of interest has happened.
    stratify_by: The categorical variable that we use to split the curves into two.
    ax: Axis to plot hazard ratios if provided.
    """
    cph = CoxPHFitter()
    surv_df_merged = pd.get_dummies(surv_df_merged, columns=[stratify_by], drop_first=True)
    cph.fit(surv_df_merged[[duration_col, event_col, stratify_by+"_Positive", "Age at last follow-up or death"]], duration_col=duration_col, event_col=event_col)
    
    # Get hazard ratios, confidence intervals, and p-values
    hazard_ratio = cph.summary['exp(coef)'][stratify_by+"_Positive"]
    ci_lower = cph.summary['exp(coef) lower 95%'][stratify_by+"_Positive"]
    ci_upper = cph.summary['exp(coef) upper 95%'][stratify_by+"_Positive"]
    p_value = cph.summary['p'][stratify_by+"_Positive"]
    
    # Print hazard ratios, confidence intervals, and p-values on the plot
    ax.text(0.9, 0.5, f"{stratify_by}\nHR={hazard_ratio:.2f} [{ci_lower:.2f}, {ci_upper:.2f}]\np={p_value:.4f}", transform=ax.transAxes, fontsize=10, ha='center', va='center')

def annotate_mutation_status(mutations_df, diagnosis, PATH_sample_information, annotate_what, annotate_gene = False, drop_dependent = True): 
    """
    Given a mutation mutations_df and a list of all pts in the cohort, annotate the status of patients in the mutation list.
    Intended for survival analysis.
    annotate_what: Either provide ctDNA or CHIP.
    annotate_gene: If a gene name is provided it will annotate the status of that gene. A list can also be given.
    """
    if isinstance(annotate_gene, str):
        mutations_df = mutations_df[mutations_df["Gene"] == annotate_gene].reset_index(drop = True)
    elif isinstance(annotate_gene, list):
        mutations_df = mutations_df[mutations_df["Gene"].isin(annotate_gene)].reset_index(drop=True)
    if drop_dependent:
        mutations_df = mutations_df[mutations_df["Dependent"] == False].reset_index(drop = True)[["Patient_id", "Timepoint", "Sample_name_t", "Gene", "VAF_n", "Protein_annotation"]]
    else: 
        mutations_df = mutations_df[["Patient_id", "Timepoint", "Sample_name_t", "Gene", "VAF_n", "Protein_annotation"]]
    mutations_df[annotate_what + " status"] = "Positive"
    mutations_df = mutations_df[["Patient_id", annotate_what + " status"]]
    all_pts = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
    mutations_df = all_pts[["Patient_id", "Diagnosis", "Timepoint"]].drop_duplicates().merge(mutations_df, how = "left")
    if diagnosis != "Both":
        mutations_df = mutations_df[mutations_df["Diagnosis"] == diagnosis]
    mutations_df[annotate_what + " status"] = mutations_df[annotate_what + " status"].fillna("Negative")
    mutations_df = mutations_df.drop_duplicates().reset_index(drop = True)[["Patient_id", "Diagnosis", "Timepoint", annotate_what+" status"]]
    return(mutations_df)


def add_gene_category_to_df(df, PATH_gene_categories):
    """
    Given a df containing the column GENES, add another col to indicate which category the gene belongs to.
    """
    genes_df = pd.read_csv(PATH_gene_categories, sep = "\t").rename(columns = {"Panel genes": "Gene"})
    df = df.merge(PATH_gene_categories, how = "left", on = "Gene")
    return(df)

def make_mutation_count_histograms(mutations_df, figure_dir):
    """
    This is mostly for QC purposes. Helps identify patients with unexpectedly high number of mutations, likely indicative of oxidative damage. 
    Note that the mutation numbers are taken from samples, not patients.
    """
    diagnosis = mutations_df["Diagnosis"].unique()[0]
    mut_type = mutations_df["Status_n"].unique()[0]
    # Get mutation counts per sample
    mutation_counts = mutations_df["Sample_name_t"].value_counts()
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot histogram
    ax.hist(mutation_counts, bins=20, color='black', edgecolor='black')
    ax.set_title('Distribution of Mutation Counts per Sample')
    ax.set_xlabel('Mutation Count')
    ax.set_ylabel('Frequency')
    # Save the plot
    fig.savefig(os.path.join(figure_dir, f'{diagnosis}_mutation_count_histogram.png'))


def make_violin_gene_cn(coverage_df, ax, ax_title):
    """
    Given the output from cnvkit and an ax to plot, plots the CN violins based on the median log2 values for each gene.
    """
    sns.set_style("whitegrid")
    
    # Calculate median log2 value for each gene
    median_log2 = coverage_df.groupby('gene')['log2'].median().reset_index(name='median_log2')
    
    # Apply get_cna and get_color functions to median log2 values
    median_log2['cna'] = median_log2['median_log2'].apply(get_cna)
    median_log2['color'] = median_log2['cna'].apply(get_color)
    
    # Merge median log2 values and color information back to coverage_df
    coverage_df = coverage_df.merge(median_log2[['gene', 'color']], on='gene')
    colors = coverage_df[["gene", "color"]].drop_duplicates()["color"]
    
    sns.violinplot(ax=ax, x='gene', y='log2', data=coverage_df, palette=colors, inner=None, linewidth=0)
    ax.set_xlabel('')
    ax.set_ylabel('Log2 coverage ratio')
    ax.set_title(ax_title)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ymax = max(ax.get_ylim())
    ax.set_ylim((-2.5, ymax))
    return ax

def get_cna(lr):
    """
    Based on a given log2 ratio, return the most likely CNA event in terms of numbers. Then these numbers can be mapped to colors and used in violin plots.
    """
    if lr <= -1:
        return -2
    if lr <= -0.3 and lr > -1:
        return -1
    if lr <= 0.3 and lr > -0.3:
        return 0
    if lr <= 1 and lr > 0.3:
        return 1
    if lr > 1:
        return 2

def get_color(cna):
    """
    Assigns a color to each CNA event value.
    """
    if cna == 2:
        return 'darkred'
    if cna == 1:
        return 'red'
    if cna == 0:
        return 'black'
    if cna == -1:
        return 'blue'
    if cna == -2:
        return 'darkblue'

def generate_upset_plot(df, path_figure, color_dict, genes=None):
    """
    Given a variant calls df, make an upset plot with colors based on patient diagnosis.
    """
    # If genes or diagnosis colors are not provided, use defaults
    if genes is None:
        genes = df["Gene"].unique().tolist()
    # Create a dictionary to store patients per gene
    my_dict = {}
    for gene in genes:
        subset_df = df[df["Gene"] == gene]
        pts = subset_df["Patient_id"].unique()
        my_dict[gene] = pts
    # Create an upset plot
    upset_df = upsetplot.from_contents(my_dict).rename(columns = {"id": "Patient_id"})
    idx = upset_df.index
    upset_df = upset_df.merge(df[["Patient_id", "Diagnosis"]].drop_duplicates(), how = "left")
    upset_df.index = idx
    upset = upsetplot.UpSet(upset_df, intersection_plot_elements=0, shading_color="white", show_percentages=True)  # disable the default bar chart
    upset.add_stacked_bars(by="Diagnosis", colors=color_dict, title="Count by diagnosis", elements=6)
    upset.style_subsets(linewidth = 1)
    ax = upset.plot()["matrix"]    
    ax.legend().set_visible(False)
    plt.savefig(path_figure)  


# Example usage:
# generate_upset_plot(df, "upset_plot.png", genes=["Gene1", "Gene2", "Gene3"], diagnosis_colors={"Kidney": "blue", "Bladder": "green"})

def make_density_plot(df, mut_effect, ax, fontsize = 8):
    df = df[df["Consequence"] == mut_effect]
    sns.kdeplot(data=df["VAF_n"], fill=True, color='black', alpha=0.6, ax = ax, bw_adjust= 0.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.set_title(mut_effect, fontsize = fontsize, loc = "left")
    return(ax)

def plot_detection_limit(df, PATH_sample_information, path_figure, ax, fontsize = 10, legendon = True):
    """
    Plots a simple regression, what percentage of our cohort would be CH+ if our LOD was 0.25, 0.50, 0.75 etc.
    """
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    sample_info = sample_info[sample_info["Diagnosis"] != "Healthy"].reset_index(drop = True)
    n_total = sample_info.shape[0]
    all_genes = df["Gene"].unique()
    gene_groups = {"DTA genes": ["DNMT3A", "TET2", "ASXL1"], "DDR genes": ["ATM", "BRCA2", "PPM1D", "CHEK2", "TP53"], "All genes": all_genes}
    results = {}
    for name, genes in gene_groups.items():
        lod = 0.25
        subset = df[df["Gene"].isin(genes)]
        perc_pos_dict = {}
        while lod <= 2:
            n_pos = subset[subset["VAF_n"] >= lod]["Patient_id"].drop_duplicates().shape[0]
            perc_pos = n_pos/n_total*100
            perc_pos_dict[lod] = perc_pos
            lod += 0.25
        results[name] = perc_pos_dict
    # Now we plot
    color_dict = {"DTA genes": "hotpink", "DDR genes":"limegreen", "All genes":"black"}
    for name in gene_groups.keys():
        color = color_dict[name]
        df = pd.DataFrame.from_dict(results[name], orient='index').reset_index().rename(columns = {"index": "Limit of detection (%)", 0: "CH+ %"})
        # sns.regplot(x="Limit of detection (%)", y="CH+ %", data=df, ax=ax, scatter_kws={'s': 5, 'color': color}, line_kws={'color': color, 'linewidth': 1})
        ax.plot(df["Limit of detection (%)"], df["CH+ %"], marker='o', color=color, markersize=5, linestyle='-', linewidth=1)
    ax.legend(frameon = False, loc = "right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2])
    ax.set_xticklabels(["0.25", "0.5", "0.75", "1", "1.25", "1.5", "1.75", "2"], fontsize = fontsize-2)
    ax.set_ylabel("CH+ %", fontsize = fontsize)
    ax.set_xlabel("Limit of detection (%)", fontsize = fontsize)
    ax.set_ylim(0, 70)
    ax.set_yticks([0, 35, 70])
    ax.set_yticklabels(["0", "35", "70"])
    # legend
    if legendon:
        legend_colors = color_dict.values()
        legend_labels = color_dict.keys()
        legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
        ax.legend(handles=legend_handles, loc="upper right", frameon=False)
    return(ax)

def plot_detection_limit_per_gene(df, genes_list, ax, fontsize = 10):
    """
    Given a list of genes and variant calls, make a barchart showing the percentage of patients we would misrepresent as mutated.
    """
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    sample_info = sample_info[sample_info["Diagnosis"] != "Healthy"].reset_index(drop = True)
    n_total = sample_info.shape[0]
    df = df[df["Gene"].isin(genes_list)][["Patient_id", "Gene", "VAF_t"]]
    df.loc[df["Gene"] == "BRCA1", "Gene"] = "BRCA1/2" # We will represent them as BRCA1/2
    df.loc[df["Gene"] == "BRCA2", "Gene"] = "BRCA1/2" # We will represent them as BRCA1/2
    df = df[df["VAF_t"] >= 1][["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts()
    df = df.reset_index()
    df["Gene perc"] = df["Gene"]/n_total*100
    df["xticklabels"] = df["index"].astype(str) + "\n" + "n=" + df["Gene"].astype(str)
    results = {}
    # Plotting the bar chart
    bars = ax.bar(df['index'], df['Gene perc'], color='hotpink')
    ax.set_xlabel('')
    ax.set_ylabel('Mischaracterized\npatients %', fontsize = fontsize)
    ax.set_title('')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim((0, 4))
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_yticklabels(["0", "1", "2", "3", "4"])
    ax.set_xlabel("DDR genes on the panel")
    ax.set_xticklabels(df["xticklabels"], fontsize = fontsize-3)
    ax.tick_params(axis='x', rotation=0)  # Rotate x-axis labels for better visibility
    return(ax)

def check_mutation_prevalence_and_sex(muts_df, sex_df, diagnosis, PATH_sample_information, annotate_what, annotate_gene = False):
    """
    Checks whether the prevalence of a gene changes with sex.
    """
    mut_status = annotate_mutation_status(muts_df, diagnosis, PATH_sample_information, annotate_what = annotate_what, annotate_gene = annotate_gene)
    df = mut_status.merge(sex_df, how = "left").dropna(subset=['Sex'])
    
    contingency_table = pd.crosstab(df[f'{annotate_what} status'], df['Sex'])
    odds_ratio, p_value = fisher_exact(contingency_table)
    if p_value < 0.05:
        print(annotate_what)
        print(f'Odds Ratio: {odds_ratio}')
        print(f'P-value: {p_value}')
        return(contingency_table)

def calculate_OR(muts_df, clinical_df, genes_list, variable, vaf_cutoff = None):
    """
    Given a continuous variable calculate the odds ratio of having a mutation in these list of genes. Logistic.
    """
    if vaf_cutoff is not None:
        df = muts_df[muts_df["VAF_n"] > vaf_cutoff]
    df = df[(df["Gene"].isin(genes_list)) & (df["Dependent"] == False)].reset_index(drop = True)[["Patient_id", "Gene"]].drop_duplicates().reset_index(drop = True)
    clinical_subset = clinical_df[["Patient_id", variable]].dropna(subset=[variable])
    
    # Merge mutation data with clinical data
    merged_df = df.merge(clinical_subset, on='Patient_id', how='inner')
    results = []
    unique_genes = merged_df['Gene'].unique()   
    # Loop through each gene to fit a logistic regression model
    for gene in genes_list:
        merged_df[gene] = (merged_df['Gene'] == gene).astype(int) # Create a binary column indicating the presence of the gene mutation
        # Fit the logistic regression model
        X = merged_df[[variable]]
        X = sm.add_constant(X)  # Add a constant term for the intercept
        y = merged_df[gene]
        model = sm.Logit(y, X)
        result = model.fit(disp=False)
        # Get the coefficient for age
        coef = result.params[variable]
        odds_ratio = np.exp(coef)
        p_value = result.pvalues[variable]
        # Apped the results
        results.append({
            'Gene': gene,
            'Coefficient': coef,
            'Odds_ratio': odds_ratio,
            'P_value': p_value
        })  
    # Convert the results to a DataFrame
    results_df = pd.DataFrame(results)
    # FDR correction
    reject, corrected_p_values, _, _ = multipletests(results_df['P_value'], method='fdr_bh')
    results_df['FDR_corrected_p_value'] = corrected_p_values
    # Display the results
    return(results_df)   

def plot_OR_forest_plots(results_df, dir_figures, filename, plot_title):
    # Calculate the standard errors and confidence intervals
    results_df = results_df.reset_index().sort_values(by = "index", ascending = False)
    results_df['SE'] = np.abs(results_df['Coefficient'] / np.sqrt(results_df['P_value']))
    results_df['CI_lower'] = np.exp(results_df['Coefficient'] - 1.96 * results_df['SE'])
    results_df['CI_upper'] = np.exp(results_df['Coefficient'] + 1.96 * results_df['SE'])
    # Plotting the forest plot
    fig, ax = plt.subplots(figsize=(2.8, 3))
    ax.errorbar(results_df['Odds_ratio'], results_df['Gene'], 
                xerr=[results_df['Odds_ratio'] - results_df['CI_lower'], results_df['CI_upper'] - results_df['Odds_ratio']], 
                fmt='o', color='black', ecolor='black', elinewidth=0.5, capsize=0, markersize = 3)
    ax.axvline(x=1, linestyle='--', color='black', lw = 0.5)
    ax.set_xlabel('OR')
    ax.set_title(plot_title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    # Annotate significant points with asterisks
    for i, (x, y, ci_lower, ci_upper) in enumerate(zip(results_df['Odds_ratio'], results_df['Gene'], results_df['CI_lower'], results_df['CI_upper'])):
        if ci_lower <= 1 <= ci_upper:  # Check if 1 is within the confidence interval
            ax.annotate('*', xy=(x, y), xytext=(3, 0), textcoords='offset points', ha='left', va='center', fontsize=8)
    
    fig.tight_layout()
    fig.savefig(os.path.join(dir_figures, filename))
    return(ax)

def calculate_OR_two_categorical(independent_variable_df, independent_variable_name, dependent_variable_df, independent_variable_df):
    """
    Runs logistic regression between two variables.
    """


        cols = 



genes_list = ["DNMT3A", "PPM1D", "TET2", "ASXL1", "CHEK2", "ATM", "TP53", "GNAS", "STAG2", "RAD21", "KMT2D", "SF3B1", "CBL", "JAK2", "BRCC3"]
calculate_OR(muts_df = all_vars_chip, kidney_clinical = kidney_clinical, bladder_clinical = bladder_clinical, diagnosis = "Kidney", genes_list = genes_list, variable = "Age")


    kidney_clinical = pd.read_csv(PATH_kidney_clinical)
bladder_clinical = pd.read_csv(PATH_clinical_bladder)

base_kidney_chip
base_bladder_chip
