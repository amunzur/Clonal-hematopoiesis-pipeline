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

def CH_presence_line_thresholds(df_CH, age_df, PATH_sample_information, ax, axforest, fontsize = 10):
    """
    Plos the fraction of the cohort that is CH positive in each age.
    """
    import scipy.stats as stats
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    sample_info = sample_info[(sample_info["Timepoint"] == "Baseline")].reset_index(drop = True)
    
    # Generate dfs of varying thresholds and annotate their CHIP mutation status
    ch1 = annotate_mutation_status(df_CH[(df_CH["VAF_n"] <= 2)], "Both", PATH_sample_information, annotate_what = "CHIP").merge(age_df, how = "left").dropna()
    ch2 = annotate_mutation_status(df_CH[(df_CH["VAF_n"] > 2) & (df_CH["VAF_n"] <= 10)], "Both", PATH_sample_information, annotate_what = "CHIP").merge(age_df, how = "left").dropna()
    ch3 = annotate_mutation_status(df_CH[df_CH["VAF_n"] > 10], "Both", PATH_sample_information, annotate_what = "CHIP").merge(age_df, how = "left").dropna()
    
    ch1 = ch1[ch1["Timepoint"] == "Baseline"].reset_index(drop = True)
    ch2 = ch2[ch2["Timepoint"] == "Baseline"].reset_index(drop = True)
    ch3 = ch3[ch3["Timepoint"] == "Baseline"].reset_index(drop = True)
    
    colors = {"0.25%-2%": "limegreen", "2%-10%": "#709bd0", ">10%": "#224193"}
    
    bins = [20, 40, 50, 60, 70, 80, 90]
    labels = ['20-39', '40-49', '50-59', '60-69', '70-79', '80-89']
    
    for j, (df, df_annot) in enumerate(zip([ch1, ch2, ch3], ["0.25%-2%", "2%-10%", ">10%"])):
        # Bin the ages
        df['Age_bin'] = pd.cut(df['Age'], bins=bins, labels=labels, right=False)
        df['Age_bin'] = pd.Categorical(df['Age_bin'], categories=labels, ordered=True)
        
        # Calculate fraction that is CH+
        fractions = []
        for i, group in df.groupby("Age_bin"):
            n_pos = group[group["CHIP status"] == "Positive"].shape[0]
            n_neg = group[group["CHIP status"] == "Negative"].shape[0]
            perc_pos = n_pos/(n_pos+n_neg)
            fractions.append(perc_pos)
        
        # LOGISTIC REGRESSOIN
        df['CHIP_binary'] = df['CHIP status'].apply(lambda x: 1 if x == 'Positive' else 0)
        df['Age_bin'] = pd.Categorical(df['Age_bin'], ordered=True) # If 'Age_Bin' is already categorical, convert it to numeric codes (ordered)
        df['Age_bin_numeric'] = df['Age_bin'].cat.codes  # Convert Age Bin to numeric
        
        # Step 1: Prepare X (predictor) and y (outcome)
        X = df['Age_bin_numeric']  # Age bins as predictor
        y = df['CHIP_binary']      # Binary CHIP status as outcome
        
        # Step 2: Add constant to X for the intercept
        X = sm.add_constant(X)
        
        # Step 3: Fit logistic regression model
        logit_model = sm.Logit(y, X)
        result = logit_model.fit()
        
        odds_ratios = pd.DataFrame({"Odds Ratio": np.exp(result.params),"p-value": result.pvalues})
        print(f"CH category: {df_annot}")
        print(odds_ratios)
        
        # Step 4: Calculate Odds Ratios and Confidence Intervals
        odds_ratios = pd.DataFrame({
            "Odds Ratio": np.exp(result.params),
            "p-value": result.pvalues,
            "Lower 95% CI": np.exp(result.params - 1.96 * result.bse),  # Lower bound
            "Upper 95% CI": np.exp(result.params + 1.96 * result.bse)   # Upper bound
        })
        OR = odds_ratios[odds_ratios.index == "Age_bin_numeric"]["Odds Ratio"][0]
        p_value = odds_ratios[odds_ratios.index == "Age_bin_numeric"]["p-value"][0]
        lower_ci = odds_ratios[odds_ratios.index == "Age_bin_numeric"]["Lower 95% CI"][0]
        upper_ci = odds_ratios[odds_ratios.index == "Age_bin_numeric"]["Upper 95% CI"][0]
        
        # Step 4. Plotting LR results as an inset plot
        age_bins_numeric_range = np.linspace(df['Age_bin_numeric'].min(), df['Age_bin_numeric'].max(), 100)
        X_plot = sm.add_constant(age_bins_numeric_range)
        y_pred_probs = result.predict(X_plot)
        ax.scatter(sorted(df['Age_bin_numeric'].unique()), fractions, color=colors[df_annot], label='Actual Fraction Positive', s = 10) #Plot the actual fraction of CH-positive individuals by age bin
        ax.plot(age_bins_numeric_range, y_pred_probs, color=colors[df_annot], label='Logistic Regression Fit', linewidth = 0.8) # Step 4: Plot the logistic regression predicted probabilities
        
        # Calculating error bars for each dot
        age_bin_stats = df.groupby('Age_bin').agg(count_positive=('CHIP_binary', 'sum'),count_total=('CHIP_binary', 'size')).reset_index()
        age_bin_stats['proportion_positive'] = age_bin_stats['count_positive'] / age_bin_stats['count_total']
        conf_intervals = []
        for index, row in age_bin_stats.iterrows():
            n = row['count_total']
            p = row['proportion_positive']
            
            if n > 0:
                ci_low, ci_high = stats.binom.interval(0.95, n=n, p=p)
                # Calculate the error margins
                conf_intervals.append((p - ci_low/n, ci_high/n - p))
            else:
                conf_intervals.append((0, 0))  # No error if no total
        conf_intervals = np.array(conf_intervals).T
        
        # Plotting the error bar
        df_plotting = pd.DataFrame({'Age_bin': age_bin_stats['Age_bin'],'Positive': age_bin_stats['proportion_positive']})
        ax.errorbar(df_plotting.index, df_plotting['Positive'], yerr=conf_intervals, fmt='none', color=colors[df_annot], capsize=2, elinewidth=0.5, capthick = 0.5)
        
        # Plotting the forest
        if axforest is not None:
            axforest.scatter(OR, j, color=colors[df_annot], s = 15)
            axforest.text(3, j, p_value, ha='center', va='center', fontsize=7, color = "black")
            axforest.errorbar(OR, j, xerr = [[OR - lower_ci], [upper_ci - OR]], fmt='none', color=colors[df_annot], capsize=2, elinewidth=1, capthick = 1)
    
    # Annotate the number of people in each age group below the x-tick labels
    patient_ages = sample_info.merge(age_df)
    patient_ages["Age_bin"] = pd.cut(patient_ages['Age'], bins=bins, labels=labels, right=False)
    patient_ages['Age_bin'] = pd.Categorical(patient_ages['Age_bin'], categories=labels, ordered=True)
    patient_ages = patient_ages["Age_bin"].value_counts().reset_index()
    patient_ages['index'] = pd.Categorical(patient_ages['index'], categories=labels, ordered=True)
    patient_ages = patient_ages.sort_values('index')
    x_ticklabellist = patient_ages["index"].astype(str) + "\nn=" + patient_ages["Age_bin"].astype(str).tolist()    
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_xticklabels(x_ticklabellist)
    # AES
    ax.set_ylim((-0.1, 1.1))
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=10)
    ax.set_xlabel('Age at baseline draw')
    ax.set_ylabel('% of patients with CH')
    ax.spines[["top", "right"]].set_visible(False)
    
    # add legend
    legend_colors = colors.values()
    legend_labels = colors.keys()
    legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="best", frameon=False, fontsize = 8, handletextpad=0.1, title = "CH VAF interval")
    
    if axforest is not None:
        axforest.axvline(1, color='black', linestyle='--', linewidth = 0.5)
        axforest.set_ylim((-0.1, 2.1))
        axforest.set_yticks([0, 1, 2])
        axforest.tick_params(axis='y', left=False, labelleft=False)
        axforest.spines[["top", "right", "left"]].set_visible(False)
        axforest.set_yticklabels(["0.25%-2%", "2%-10%", ">10%"])
        axforest.set_xlabel("OR")
        return(ax, axforest)
    else:
        return(ax)

def plot_age_plots_stacked(df, age_df, PATH_sample_information, figure_dir, fontsize, color_dict, ext = "png"):
    """
    SIMILAR TO THE plot_age_plots FUNCTION - ONLY DIFFERENCE IS IT PLOTS BLADDER AND KIDNEY TOGETHER.
    Makes a figure with two subplots stacked on top of each other.
    Top one is a line plot showing the percentage of the cohort that is CH positive for each age bin.
    Bottom is a histogram showing the number of patients we have in each bin.
    """
    color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}
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
    
    # Get the WBC vafs for each age group 
    max_vaf_per_patient = df.groupby("Patient_id")["VAF_n"].max().reset_index().merge(age_df, how = "left").dropna()
           
    # Create figure and gridspec
    fig = plt.figure(figsize=(3, 5))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])
    
    # Plot boxplot in the first subplot - boxplot of vafs per age group
    ax0 = plt.subplot(gs[0])
    for diagnosis, color in color_dict.items():
        subset = max_vaf_per_patient[max_vaf_per_patient["Diagnosis"] == diagnosis]
        ax0.scatter(subset["Age"], subset["VAF_n"], color=color, label=diagnosis, s = 3)
        # sns.regplot(x="Age", y="VAF_n", data=subset, ax=ax0, scatter_kws={'s': 2, 'color': color}, line_kws={'color': color, 'linewidth': 2, "linestyle": "-"})
    ax0.set_xlim((20, 100))
    ax0.set_xticks([20, 40, 60, 80, 100])
    ax0.set_xticklabels(["20", "40", "60", "80", "100"], fontsize=fontsize-2)
    ax0.set_ylim((-5, 100))
    ax0.set_yticks([0, 25, 50, 75, 100])
    ax0.set_yticklabels(["0", "25", "50", "75", "100"], fontsize=fontsize-2)
    ax0.set_xlabel('Age at baseline blood draw', fontsize=fontsize)
    ax0.set_ylabel('Highest\nWBC VAF (%)', fontsize=fontsize)
    ax0.spines["top"].set_visible(False)
    ax0.spines["right"].set_visible(False)
    
    # Plot dot plot in the first subplot
    ax1 = plt.subplot(gs[1])
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
    ax1.set_ylim((-10, 110))
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    
    # Plot histogram in the second subplot
    ax2 = plt.subplot(gs[2], sharex=ax1)
    patient_counts = all_pts['Age_bin'].value_counts().reindex(labels).fillna(0)
    ax2.bar(labels, patient_counts_bladder, color=color_dict["Bladder"], label='Bladder', edgecolor = None, linewidth=0.0)
    ax2.bar(labels, patient_counts_kidney, bottom=patient_counts_bladder, color=color_dict["Kidney"], label='Kidney', edgecolor = None, linewidth=0.0)
    ax2.set_xlabel('Age at baseline blood draw', fontsize=fontsize)
    ax2.set_ylabel('Number of\npatients', fontsize=fontsize)
    ax2.set_xticklabels(labels, fontsize=fontsize-2)
    ax2.set_yticks([0, 20, 40, 60, 80, 100])
    ax2.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=fontsize)
    # ax2.set_ylim((-10, 110))
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

def generate_patient_profiles_dot_plot(patient_df, figure_dir, PATH_sample_information, show_dates_on_x_axis, ax_main, color_map_name, patient_id, diagnosis, variant_type, title, show_legend = True, add_figure_title = True, spine_to_hide = ["top", "right"], add_somatic_to_legend = False, gene=None, markersize=8, gene_dict = None):    
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
    # if gene_dict is not None:
        
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
                return f"{months_from_baseline}m"
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
                ax_main.plot(baseline_date, [y_values_list[0][j]][0], marker='o', color=mutation_color, markersize=markersize, linestyle='-', linewidth=1, zorder = zorder)
            elif len(y_values_list) > 1:
                x_values = [vals[j] for vals in x_values_list]
                y_values = [vals[j] for vals in y_values_list]
                ax_main.plot(x_values, y_values, marker='o', color=mutation_color, markersize=markersize, linestyle='-', linewidth=1, zorder = zorder)
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

def plot_treatment_landscape(ax_treatment, ax_chip, patient_id, PATH_treatment_landscape, sample_info, show_nac = True):
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
    
    if not show_nac:
        pt_treatment = pt_treatment[pt_treatment["Treatment line"] != "NAC"]
    
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

def plot_gene_counts_grouped(all_vars_chip, figure_dir, filename, PATH_sample_information, ax, gene_list = None):
    """
    Plots the number of mutations each individual gene has in the cohort as a BAR CHART. If unique, non unique mutations in the same gene are counted only once.
    Plots kidney and bladder separatelyin grouped bar charts. 
    """
    sample_info = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
    n_bladder = sample_info[sample_info["Diagnosis"] == "Bladder"].drop_duplicates("Patient_id").shape[0]
    n_kidney = sample_info[sample_info["Diagnosis"] == "Kidney"].drop_duplicates("Patient_id").shape[0]
    
    df = all_vars_chip[all_vars_chip["Gene"].isin(gene_list)].reset_index(drop = True)
    # subset
    bladder_df = df[df["Diagnosis"] == "Bladder"].reset_index(drop = True)[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "Counts_bladder"}).reset_index()
    kidney_df = df[df["Diagnosis"] == "Kidney"].reset_index(drop = True)[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "Counts_kidney"})
    # Go from counts to fractions
    bladder_df["frac_bladder"] = (bladder_df["Counts_bladder"]/n_bladder)*100
    kidney_df["frac_kidney"] = (kidney_df["Counts_kidney"]/n_kidney)*100
    
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
    # fig, ax = plt.subplots(figsize = (5.6, 2.8))
    ax2 = ax.twinx()
    for i, row in combined.iterrows():
        ax.bar(row["index"]+0.2, row["frac_kidney"], color=kidney_color, width = 0.4, edgecolor = "None")
        ax.bar(row["index"]-0.2, row["frac_bladder"], color=bladder_color, width = 0.4, edgecolor = "None")
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
    ax.set_xticklabels(tick_df["Gene"].tolist(), rotation = 90, fontstyle = "italic")
    ax.tick_params(axis = "x", which='both', length=0)
    ax.set_ylabel("% of the cohort mutated")
    ax.set_ylim((0, 75))
    ax.set_yticks([0, 25, 50, 75])
    ax.set_yticklabels(["0", "25", "50", "75"])
    ax.set_xlim((-0.6, 14.6))
    ax2.set_ylabel("WBC VAF%")
    ax2.set_ylim((0, 50))
    ax2.set_yticks([0, 10, 20, 30, 40, 50])
    ax2.set_yticklabels(["0", "10", "20", "30", "40", "50"])
    ax2.tick_params(axis='x', pad=1)
    
    # Calculate proportions for each gene in bladder and kidney cancers
    combined['Proportion_bladder'] = combined['Counts_bladder'] / combined['Counts_bladder'].sum()
    combined['Proportion_kidney'] = combined['Counts_kidney'] / combined['Counts_kidney'].sum()
    # Perform the Z-test for each gene
    p_values_ztest = []
    # Perform the Z-test for each gene
    for i, (bladder_count, kidney_count) in enumerate(zip(combined["Counts_bladder"], combined["Counts_kidney"])):
        counts = np.array([bladder_count, kidney_count])
        nobs = np.array([combined["Counts_bladder"].sum(), combined["Counts_kidney"].sum()])
        _, p_val = proportions_ztest(counts, nobs, alternative="larger")  # alternative="larger" for one-sided test
        p_values_ztest.append(p_val)
    
    # Plotting with significance marks
    max_val_ax = ax.get_ylim()[1]
    for i, p_val in enumerate(p_values_ztest):
        if p_val < 0.001:
            ax.text(i + 0.2, max_val_ax-1, "***", ha='center', va='bottom', color='black')
        elif p_val < 0.01:
            ax.text(i + 0.2, max_val_ax-1, "**", ha='center', va='bottom', color='black')
        elif p_val < 0.05:
            ax.text(i + 0.2, max_val_ax-1, "*", ha='center', va='bottom', color='black')
    
    print(p_values_ztest)
    return([ax, ax2])

def plot_gene_counts_grouped_single_ax(chip_df, ax, gene_list = None, horizontal = False):
    """
    Plots the number of mutations each individual gene has in the cohort as a BAR CHART. If unique, non unique mutations in the same gene are counted only once.
    Plots kidney and bladder separatelyin grouped bar charts. 
    """
    # Generate our denominators, which will be the number of CH+ people.
    n_kidney = chip_df[chip_df["Diagnosis"] == "Kidney"]["Patient_id"].unique().shape[0]
    n_bladder = chip_df[chip_df["Diagnosis"] == "Bladder"]["Patient_id"].unique().shape[0]
    
    if gene_list is not None: 
        df = chip_df[chip_df["Gene"].isin(gene_list)].reset_index(drop = True)
    
    # subset
    bladder_df = df[df["Diagnosis"] == "Bladder"].reset_index(drop = True)[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "Counts_bladder"})
    kidney_df = df[df["Diagnosis"] == "Kidney"].reset_index(drop = True)[["Patient_id", "Gene"]].drop_duplicates()["Gene"].value_counts().reset_index().rename(columns = {"index": "Gene", "Gene": "Counts_kidney"})
    
    # Go from counts to fractions
    bladder_df["frac_bladder"] = (bladder_df["Counts_bladder"]/n_bladder)*100
    kidney_df["frac_kidney"] = (kidney_df["Counts_kidney"]/n_kidney)*100
    
    combined = bladder_df.merge(kidney_df, on = "Gene", how = "inner").reset_index() # gene order will be determined based on bladder df
    gene_rank = {gene: rank for rank, gene in enumerate(gene_list)}
    combined['index'] = combined['Gene'].map(gene_rank)
    combined = combined.sort_values(by = "index")
    
    kidney_color = "orangered"
    bladder_color = "deepskyblue"
    
    # plotting
    for i, row in combined.iterrows():
        # Plot bars
        if horizontal:
            ax.barh(row["index"]+0.2, row["frac_kidney"], color=kidney_color, height=0.4, edgecolor="None")
            ax.barh(row["index"]-0.2, row["frac_bladder"], color=bladder_color, height=0.4, edgecolor="None")
            
            ax.text(row["frac_kidney"]+4, row["index"]+0.2, row["Counts_kidney"], va='center', color='black', fontsize=5)
            ax.text(row["frac_bladder"]+4, row["index"]-0.2, row["Counts_bladder"], va='center', color='black', fontsize=5)
        else:
            ax.bar(row["index"]+0.2, row["frac_kidney"], color=kidney_color, width = 0.4, edgecolor = "None")
            ax.bar(row["index"]-0.2, row["frac_bladder"], color=bladder_color, width = 0.4, edgecolor = "None")
            
            ax.text(row["index"]+0.2, row["frac_kidney"]+2, row["Counts_kidney"], ha='center', color='black', fontsize = 5)
            ax.text(row["index"]-0.2, row["frac_bladder"]+2, row["Counts_bladder"], ha='center', color='black', fontsize = 5)    
    
    # aesthetics
    if horizontal:
        # ax.tick_params(axis='x', bottom=False, labelbottom=False)
        ax.tick_params(axis='y', left=False, labelleft=False)
        ax.spines[["top", "left"]].set_visible(False)
        tick_df = combined[["index", "Gene"]].sort_values(by = "index")
        ax.set_yticks(tick_df["index"].tolist())
        ax.set_yticklabels(tick_df["Gene"].tolist(), rotation = 0, fontstyle = "italic", va = "center", fontsize = 8, ha = "center")
        ax.tick_params(axis='y', which='both', pad=20)  # Adjust padding between tick labels and axis
        ax.tick_params(axis = "y", which='both', length=0)
        ax.set_xlabel("% mutated in CH+ pts")
        ax.set_xlim((0, 60))
        ax.set_xticks([0, 15, 30, 45, 60])
        ax.set_xticklabels(["0", "15", "30", "45", "60"])
        ax.set_ylim((-0.6, 14.6))
    else:
        ax.spines[["top", "right"]].set_visible(False)
        tick_df = combined[["index", "Gene"]].sort_values(by = "index")
        ax.set_xticks(tick_df["index"].tolist())
        ax.set_xticklabels(tick_df["Gene"].tolist(), rotation = 90, fontstyle = "italic", va = "center", fontsize = 8)
        ax.tick_params(axis='x', which='both', pad=20)  # Adjust padding between tick labels and axis
        ax.tick_params(axis = "x", which='both', length=0)
        ax.set_ylabel("% mutated in CH+ pts")
        ax.set_ylim((0, 60))
        ax.set_yticks([0, 15, 30, 45, 60])
        ax.set_yticklabels(["0", "15", "30", "45", "60"])
        ax.set_xlim((-0.6, 14.6))
    
    # Calculate proportions for each gene in bladder and kidney cancers
    combined['Proportion_bladder'] = combined['Counts_bladder'] / combined['Counts_bladder'].sum()
    combined['Proportion_kidney'] = combined['Counts_kidney'] / combined['Counts_kidney'].sum()
    
    # Perform the Z-test for each gene
    p_values_ztest = []
    nobs = np.array([n_bladder, n_kidney])
    for i, (bladder_count, kidney_count) in enumerate(zip(combined["Counts_bladder"], combined["Counts_kidney"])):
        counts = np.array([bladder_count, kidney_count])
        _, p_val = proportions_ztest(counts, nobs, alternative="larger")  # alternative="larger" for one-sided test
        p_values_ztest.append(p_val)
    
    adjusted_p_values_ztest = multipletests(p_values_ztest, method="fdr_bh")[1]
    adjusted_p_values_ztest = [float(p) for p in adjusted_p_values_ztest]
    3.064581261712751e-05
    # Plotting with significance marks
    # max_val_ax = ax.get_ylim()[1]
    for i, p_val in enumerate(adjusted_p_values_ztest):
        # Get the bladder count, the p values will be annotated right on top of it.
        bladder_count = combined[combined["index"] == i]["frac_bladder"].iloc[0]
        if horizontal:
            if p_val < 0.05:
                formatted_p_value = round(float(p_val), 3)
                formatted_text = f"p={formatted_p_value}"
                ax.text(bladder_count+4, i-0.2, formatted_text, ha='right', va='bottom', color='black', fontsize = 6, fontdict={'family': 'Arial'})
            elif p_val < 0.001: 
                coefficient = f"{value:.2e}".split('e')[0]  # e.g., 2.59
                exponent = int(f"{value:.2e}".split('e')[1])  # e.g., -5
                formatted_p_value = rf"${coefficient} \times 10^{{{exponent}}}$"
                formatted_text = f"p={formatted_p_value}"
                ax.text(bladder_count+4, i-0.2, formatted_text, ha='right', va='bottom', color='black', fontsize = 6, fontdict={'family': 'Arial'})
        else:
            if p_val < 0.001:
                ax.text(i + 0.2, bladder_count+5, "***", ha='center', va='bottom', color='black')
            elif p_val < 0.01:
                ax.text(i + 0.2, bladder_count+5, "**", ha='center', va='bottom', color='black')
            elif p_val < 0.05:
                ax.text(i + 0.2, bladder_count+5, "*", ha='center', va='bottom', color='black')
    
    print(p_values_ztest)
    
    # Generate the gene_order dict for the subsequent plot
    gene_order = combined[["Gene", "index"]]
    
    # Add legend
    legend_colors = ["deepskyblue", "orangered"]
    legend_labels = ["mUC", "mRCC"]
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="lower left", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.5)
    
    if horizontal:
        ax.invert_xaxis()
        ax.set_xlabel("% mutated in CH+ pts")
    return(ax, gene_order)

def plot_vafs_swarm(chip_df, ax_swarm, gene_list = None, revert_vafs = False, gene_order = None, add_violins = True, horizontal = False):
    """
    Plots the VAFs as swarm.
    """
    # Subset to genes of interest
    if gene_list is not None: 
        df = chip_df[chip_df["Gene"].isin(gene_list)].reset_index(drop = True)
    
    df = df[["Gene", "VAF_n", "Diagnosis"]]
    
    if gene_order is not None:
        df = df.merge(gene_order, how = "left")
    
    # if revert_vafs:
    #     df["VAF_n"] = df["VAF_n"]*-1
    # Convert VAFs to a log scale, keeping them negative
    df["VAF_n_log"] = np.log10(df["VAF_n"].replace(0, np.nan))  # Replace 0 with NaN to avoid log(0)
    
    kidney_vafs = df[df["Diagnosis"] == "Kidney"]
    bladder_vafs = df[df["Diagnosis"] == "Bladder"]
        
    for i, row in kidney_vafs.iterrows():
        jitter = np.random.uniform(-0.05, 0.05, 1)
        ax_swarm.scatter(row["VAF_n_log"], row["index"]+0.2+jitter, color="black", alpha = 0.7, s = 0.15, zorder = 1) # kidney
    
    for i, row in bladder_vafs.iterrows():
        jitter = np.random.normal(-0.05, 0.05, 1)
        ax_swarm.scatter(row["VAF_n_log"], row["index"]-0.2+jitter, color="black", alpha = 0.7, s = 0.15, zorder = 1) # Bladder
    
    # ax_swarm.set_yscale('log')
    ax_swarm.spines[["right", "top"]].set_visible(False)
    ax_swarm.set_xlim((np.log10(0.25), np.log10(50)))
    ax_swarm.set_xticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(50)])
    ax_swarm.set_xticklabels([".25", "1", "2", "10", "50"])
    ax_swarm.tick_params(axis='y', labelleft=False, left = False)
    ax_swarm.set_xlabel("WBC VAF%")
    
    # Shades behind the swarmplot
    if add_violins:
        positions_kidney = np.arange(len(gene_list)) + 0.2
        positions_bladder = np.arange(len(gene_list)) - 0.2
        width = 0.3 
        x_min, x_max = ax_swarm.get_xlim() # Get the current y limits (assuming your axis has been set up already
        bar_height = x_max - x_min #  Calculate the height for the bars
        
        # Add rectangles
        for pos in positions_kidney:
            ax_swarm.barh(pos, bar_height, height = width, left = x_min, color='orangered', alpha=0.1, edgecolor = None, zorder = 0)
        
        for pos in positions_bladder:
            ax_swarm.barh(pos, bar_height, height = width, left = x_min, color='deepskyblue', alpha=0.1, edgecolor = None, zorder = 0)
    
    return(ax_swarm)

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


def plot_fraction_of_CH_calls(df_somatic_main, df_ch_main, ax0, ax1):
    """
    For each gene show the fraction of calls that are CHIP. 
    Kidney and bladder shown separately.
    Top plot is a line plot showing the percentage of calls that are CH, bottom plot is a barchart showing the count for total number of mutations.
    """
    # Calculate percentage that is CH vs ctDNA for every gene.
    # df_ch_main = baseline.copy()        
    # df_somatic_main = baseline_somatic.copy()
    df_somatic = df_somatic_main[["Gene"]].sort_values("Gene").value_counts().reset_index()
    df_somatic.columns = ["Gene", "Counts_ctDNA"]
    df_ch = df_ch_main[["Gene"]].sort_values("Gene").value_counts().reset_index()
    df_ch.columns = ["Gene", "Counts_ch"]
    
    # Get percentage values
    combined = df_ch.merge(df_somatic, on = ["Gene"], how = "outer").fillna(0)
    combined["summed"] = combined["Counts_ch"] + combined["Counts_ctDNA"]
    combined["perc_ch"] = (combined["Counts_ch"]/combined["summed"])*100
    combined["perc_ctDNA"] = (combined["Counts_ctDNA"]/combined["summed"])*100
    # combined = combined[combined["summed"] >= 5]
    combined = combined.sort_values(by = "perc_ch", ascending = False).reset_index(drop = True)
    myd88_row = combined[combined['Gene'] == 'MYD88']
    remaining_rows = combined[combined['Gene'] != 'MYD88']
    combined = pd.concat([myd88_row, remaining_rows]) 
    # get the number of patients with CH mutations in each gene
    patient_counts = df_ch_main[["Patient_id", "Gene"]].drop_duplicates().reset_index()["Gene"].value_counts().reset_index()
    patient_counts.columns = ["Gene", "n_patients"]
    
    combined = combined.merge(patient_counts, how = "left", on = "Gene").fillna(0)
    
    # Some genes require additional margin
    additional_margin = {
        "KMT2D": -1,
        "SETD2": 15,
        "PBRM1": 10, 
        "KRAS": 2, 
        "CEBPA": 3,
        "MPL": 5,
        "ETV6": 6,
        "LPAR6": 7, 
        "SPOP": 7, 
        "PTPN11": 7, 
        "NRAS": 7}
    combined["Text margin"] = combined["Gene"].map(additional_margin).fillna(0).astype(int)
    
    for i, row in combined.iterrows():
        if row["Gene"] in ["CHEK2", "ATM", "TERT", "TP53", "BRCA1", "BRCA2", "ARID1A", "ERBB2", "AR", "KMT2D", "PBRM1", "BAP1", "SETD2"]:
            row_color = "#eb4034"
        else:
            row_color = "dimgray"
        ax1.bar(i, row["summed"], color = row_color, edgecolor = "None")
        if row["summed"] > 0:
            ax1.text(i, row["summed"] + 0.3 + row["Text margin"], str(int(row["summed"])), ha='center', va='bottom', fontsize = 7)
    
    ax0.scatter(combined.index, combined["perc_ch"], color = "black", edgecolor = "None", zorder = 999999, s = 15)
    ax0.plot(combined.index, combined["perc_ch"], color = "black", linewidth = 0.5, zorder = 999999)
    
    for ax in [ax0, ax1]:
        ax.spines[["top"]].set_visible(False)
        ax.set_xlabel("")
    
    ax0.set_ylim((-3, 102))
    ax0.set_ylabel("% of mutations with CH origin")
    ax0.set_yticks([0, 20, 40, 60, 80, 100])
    ax0.set_yticklabels(["0", "20", "40", "60", "80", "100"])
    ax0.set_xticks(combined.index.tolist())
    ax0.set_xlim((combined.index.min()-1, combined.index.max()+1))
    ax1.set_xticklabels(combined["Gene"], rotation=90, fontsize = 8)
    ax1.tick_params(axis='x', which='both', bottom=False, top=False)
    ax1.xaxis.set_tick_params(pad=0)  # Adjust padding to bring labels closer to the axis
    
    ax1.set_ylabel("Number of somatic mutations")
    ax1.set_yticks([0, 50, 100, 150, 200])
    ax1.set_yticklabels(["0", "50", "100", "150", "200"])
    ax1.set_ylim((0, 210))
    
    # ax0.figure.canvas.draw()  # Ensure the figure is drawn
    # ax0.set_zorder(ax1.get_zorder() + 1)  # Put ax0 on top of ax1
    # fig = ax0.figure
    # fig.canvas.draw()  # Ensure the figure is drawn
    # fig.subplots_adjust(bottom=0.2)  # Adjust bottom margin if necessary
    # fig.axes[0].set_zorder(1)  # Ensure ax0 is on top of ax1
    return(ax0, ax1)


def plot_variance_per_pt(baseline, ax, scatter_color, bar_color):
    # baseline.loc[baseline["Gene"] == "ZRSR2", "VAF_n"] = baseline.loc[baseline["Gene"] == "ZRSR2", "VAF_n"]/2  
    
    # Calculate coefficient of variance per pt
    grouped = baseline.groupby("Patient_id")["VAF_n"]
    cv_per_pt = (grouped.std() / grouped.mean()).reset_index(name='CV')
    
    cv_per_pt.columns = ["Patient_id", "CV"]
    cv_per_pt.sort_values(by = "CV", ascending = False, inplace = True)
    cv_per_pt = cv_per_pt.dropna(subset = "CV").reset_index(drop = True)
    cv_per_pt = cv_per_pt.reset_index()
    
    # Pull the WBC VAF of all mutations in these patients
    baseline["VAF_n"] = np.log10(baseline["VAF_n"])
    vaf_df = baseline[["Patient_id", "VAF_n"]]
    vaf_df = vaf_df[vaf_df["Patient_id"].isin(cv_per_pt["Patient_id"])].reset_index(drop = True)
    
    # merge
    cv_per_pt = cv_per_pt.merge(vaf_df)
    
    # plotting
    ax_twin = ax.twinx()
    for i, group in cv_per_pt.groupby("Patient_id"):
        var = group["CV"].unique()[0]
        x_pos = group["index"].unique()[0]
        ax.bar(x_pos, var, color=bar_color, edgecolor = "None")
        ax_twin.scatter(np.repeat(x_pos, len(group)), group["VAF_n"], s = 2, edgecolor = None, color = scatter_color, alpha = 0.7, zorder = 100)
    
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_zorder(1000)
    ax.set_ylabel("Coefficient of variance")
    
    ax_twin.set_ylim((np.log10(0.25), np.log10(60)))
    ax_twin.set_yticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(50)])
    ax_twin.set_yticklabels([".25", "1", "2", "10", "50"])
    
    ax_twin.spines["top"].set_visible(False)
    ax_twin.set_ylabel("WBC VAF%")
    ax_twin.spines["bottom"].set_zorder(1000)
    
    # Set x limits
    xmin = cv_per_pt["index"].min() - 1
    xmax = cv_per_pt["index"].max() + 1
    ax.set_xlim((xmax, xmin))
    ax_twin.set_xlim((xmax, xmin))
    ax_twin.set_xlabel("Baseline samples with >1 CH mutation")
    
    return(ax, ax_twin)

def plot_fraction_of_CH_calls_with_bars(df_somatic_main, df_ch_main, ax0, ax1, ax2, stacked = True, yval_max = -200):
    """
    For each gene show the fraction of calls that are CHIP. 
    Kidney and bladder shown separately.
    Top plot is a line plot showing the percentage of calls that are CH, bottom plot is a barchart showing the count for total number of mutations.
    """
    # Calculate percentage that is CH vs ctDNA for every gene.
    # df_ch_main = baseline.copy()        
    # df_somatic_main = baseline_somatic.copy()
    df_somatic = df_somatic_main[["Gene"]].sort_values("Gene").value_counts().reset_index()
    df_somatic.columns = ["Gene", "Counts_ctDNA"]
    df_ch = df_ch_main[["Gene"]].sort_values("Gene").value_counts().reset_index()
    df_ch.columns = ["Gene", "Counts_ch"]   
    
    # Get percentage values
    combined = df_ch.merge(df_somatic, on = ["Gene"], how = "outer").fillna(0)
    combined["summed"] = combined["Counts_ch"] + combined["Counts_ctDNA"]
    combined["perc_ch"] = (combined["Counts_ch"]/combined["summed"])*100
    combined["perc_ctDNA"] = (combined["Counts_ctDNA"]/combined["summed"])*100
    # combined = combined[combined["summed"] >= 5]
    combined = combined.sort_values(by = "perc_ch", ascending = False).reset_index(drop = True)
    myd88_row = combined[combined['Gene'] == 'MYD88']
    remaining_rows = combined[combined['Gene'] != 'MYD88']
    combined = pd.concat([myd88_row, remaining_rows]) 
    # get the number of patients with CH mutations in each gene
    patient_counts = df_ch_main[["Patient_id", "Gene"]].drop_duplicates().reset_index()["Gene"].value_counts().reset_index()
    patient_counts.columns = ["Gene", "n_patients"] 
    
    combined = combined.merge(patient_counts, how = "left", on = "Gene").fillna(0)
    
    for i, row in combined.iterrows():
        if stacked:
            ax0.bar(i, row["perc_ch"], color = "black", edgecolor = "None")
            ax0.bar(i, row["perc_ctDNA"], bottom = row["perc_ch"], color = "gray", edgecolor = "None")
        else:
            ax0.bar(i, row["perc_ch"], color = "black", edgecolor = "None")
    
    # Do the broken ax for the bar chart
    ax1.bar(combined.index, -combined["summed"], color = "gray", edgecolor = "None")    
    ax2.bar(combined.index, -combined["summed"], color = "gray", edgecolor = "None")  
    ax1.set_ylim((-40, 0))  
    ax2.set_ylim((yval_max, -60))  
    
    ax0.set_ylim((0, 100))
    if stacked: 
        ax0.set_ylabel("% Origin")
    else:
        ax0.set_ylabel("% CH Origin")
    ax0.set_yticks([0, 50, 100])
    ax0.set_yticklabels(["0", "50", "100"])
    ax0.set_xticks(combined.index.tolist())
    ax0.set_xlim((combined.index.min()-1, combined.index.max()+1))
    ax0.spines[["top", "right"]].set_visible(False)  
    ax0.tick_params(axis='x', bottom=False)
    ax0.set_xticklabels(combined["Gene"], rotation=90, fontsize = 8, va = "center")
    ax0.xaxis.set_tick_params(pad=15)
    ax0.axhline(50, linestyle='--', color='gray', lw = 0.5, zorder = 9999999)
    
    # Bold some of the x tick labels
    bolded_genes = ["CHEK2", "ATM", "BRCA1", "BRCA2", "ERBB2", "TP53", "BAP1", "SETD2", "MTOR", "KMT2D", "PBRM1"]
    # Modify labels: make specific ones bold
    for tick in ax0.get_xticklabels():
        if tick.get_text() in bolded_genes:
            tick.set_fontweight('bold')
    
    ax1.tick_params(axis='x', which='both', bottom=False, top=False)
    ax1.xaxis.set_tick_params(pad=0)  # Adjust padding to bring labels closer to the axis   
    ax1.set_ylabel("Count")
    ax1.set_yticks([-40, -20, 0])
    ax1.set_yticklabels(["40", "20", "0"])
    ax1.spines[["bottom", "right"]].set_visible(False)
    ax1.tick_params(axis='x', bottom=False, labelbottom=False)
    ax1.axhline(-40, linestyle='--', color='gray', lw = 0.5)       
    
    ax2.spines[["bottom", "right", "top"]].set_visible(False) 
    ax2.set_yticks([-60, yval_max]) 
    ax2.set_yticklabels(["60", str(abs(yval_max))]) 
    ax2.tick_params(axis='x', bottom=False, labelbottom=False)
    ax2.axhline(-60, linestyle='--', color='gray', lw = 0.5)       
    
    return(ax0, ax1, ax2)

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

def plot_per_patient_counts_grouped_bar_chart(muts_df, PATH_sample_information, ax):
    """
    In a given ax plots the fraction of the bladder and kidney cohorts separately that have CH mutations.
    """
    kidney_color = "orangered"
    bladder_color = "deepskyblue"
    
    sample_info = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
    n_bladder = sample_info[sample_info["Diagnosis"] == "Bladder"].drop_duplicates("Patient_id").shape[0]
    n_kidney = sample_info[sample_info["Diagnosis"] == "Kidney"].drop_duplicates("Patient_id").shape[0]
    
    # get mut counts
    mutation_counts = muts_df.groupby(['Diagnosis', 'Patient_id']).size().reset_index(name='Mutation count')
    count_table = mutation_counts.groupby(['Diagnosis', 'Mutation count']).size().reset_index(name='Number of patients')
    pivot_table = count_table.pivot(index='Mutation count', columns='Diagnosis', values='Number of patients').fillna(0).astype(int)
    pivot_table.columns = [f'Number of patients {col.lower()}' for col in pivot_table.columns]
    pivot_table.reset_index(inplace=True)
    pivot_table["Bladder_fraction"] = pivot_table["Number of patients bladder"] / n_bladder
    pivot_table["Kidney_fraction"] = pivot_table["Number of patients kidney"] / n_kidney
    
    # Get the number of patients with 0 mutations
    n_ch_positive_bladder = muts_df[muts_df["Diagnosis"] == "Bladder"]["Patient_id"].unique().shape[0]
    n_ch_positive_kidney = muts_df[muts_df["Diagnosis"] == "Kidney"]["Patient_id"].unique().shape[0]
    n_0_bladder = n_bladder - n_ch_positive_bladder
    n_0_kidney = n_kidney - n_ch_positive_kidney
    
    pts_with_0_muts_dict = {"Mutation count": 0, "Number of patients bladder": n_0_bladder, "Number of patients kidney": n_0_kidney, "Bladder_fraction": n_0_bladder/n_bladder, "Kidney_fraction": n_0_kidney/n_kidney}
    pts_with_0_muts_df = pd.DataFrame.from_dict(pts_with_0_muts_dict, orient='index').T
    pivot_table = pd.concat([pivot_table, pts_with_0_muts_df]).sort_values(by = "Mutation count").reset_index(drop = True)
    pivot_table["xpos"] = pivot_table.index
    pivot_table.loc[pivot_table["xpos"] == 12, "xpos"] = 13
    
    # For the swarm plot get the vafs of all muts
    muts_df = muts_df.merge(mutation_counts, how = "left")
    muts_df.loc[muts_df["Mutation count"] == 17, "Mutation count"] = 13
    kidney_muts = muts_df[muts_df["Diagnosis"] == "Kidney"].reset_index(drop = True)
    bladder_muts = muts_df[muts_df["Diagnosis"] == "Bladder"].reset_index(drop = True)
    
    # Plot them
    ax2 = ax.twinx()
    for i, row in pivot_table.iterrows():
        ax.bar(row["xpos"]+0.2, row["Kidney_fraction"], color=kidney_color, width = 0.4, edgecolor = "None")
        ax.text(row["xpos"] + 0.2, 0.3, int(row["Number of patients kidney"]), ha='center', va='center', fontsize=5, color='black')
        
        ax.bar(row["xpos"]-0.2, row["Bladder_fraction"], color=bladder_color, width = 0.4, edgecolor = "None")
        ax.text(row["xpos"] - 0.2, 0.3, int(row["Number of patients bladder"]), ha='center', va='center', fontsize=5, color='black')
    
    # plot the vafs in logscale
    muts_df["VAF_n_log"] = np.log10(muts_df["VAF_n"].replace(0, np.nan))
    for i, row in muts_df.iterrows():
        if row["Diagnosis"] == "Kidney": 
            offset = 0.2
        else:
            offset = -0.2
        jitter = np.random.uniform(-0.08, 0.08, 1)
        ax2.scatter(row["Mutation count"]+offset+jitter[0], row["VAF_n_log"], color="black", s = 0.15, alpha = 0.7)
    
    # plot the vafs for bladder
    # for i, row in bladder_muts.iterrows():
    #     jitter = np.random.normal(-0.05, 0.05, 1)
    #     print(row["Mutation count"]-0.2+jitter, row["VAF_n"])
    #     ax2.scatter(row["Mutation count"]-0.2+jitter, row["VAF_n"], color="black", s = 2)
    
    # Aes
    for a in [ax, ax2]:
        a.spines["top"].set_visible(False)
    # x ticks
    labels = [f"{i}" for i in range(12)]  # Create labels for 0 to 11
    labels.append("17")
    ax.set_xticks(np.append(np.arange(0, 12), 13))
    ax.set_xticklabels(labels)
    ax.set_xlabel("Number of CH mutations")
    ax.set_ylabel("% of patients with CH")
    ax.tick_params(axis='x', bottom=False)
    ax2.tick_params(axis='x', pad=-5)
    # ax.tick_params(axis="both", direction="out", which="both", left=True, bottom=True , colors='k')    
    ax.set_ylim((0, 0.3))
    ax.set_yticks([0, 0.1, 0.2, 0.3])
    ax.set_yticklabels(["0", "10", "20", "30"])
    ax2.set_ylim((np.log10(0.25), np.log10(60)))
    ax2.set_yticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(50)])
    ax2.set_yticklabels([".25", "1", "2", "10", "50"])
    ax2.set_ylabel("WBC VAF %")

    # Add legend
    legend_colors = ["deepskyblue", "orangered"]
    legend_labels = ["mUC", "mRCC"]
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.5)
    return([ax, ax2])

def plot_ch_presence_absence_bars(muts_df, PATH_sample_information, ax, annotate_what):
    """
    Plots the presence / absence of CH in the cohort as bar charts.
    """
    all_pts = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
    all_pts = all_pts[all_pts["Timepoint"] == "Baseline"]
    
    kidney_pts_n = all_pts[all_pts["Diagnosis"] == "Kidney"].shape[0]
    bladder_pts_n = all_pts[all_pts["Diagnosis"] == "Bladder"].shape[0]
    color_dict = {"Bladder": "deepskyblue", "Kidney": "orangered"}
    
    if annotate_what.lower() == "chip":
        col_to_filter = "VAF_n"
    else:
        col_to_filter = "VAF_t"
    
    results_dict = {}
    for min_vaf in [0, 2, 10]:
        muts_df_filtered  = muts_df[muts_df[col_to_filter] >= min_vaf]
        # percentage of the bladder patients that are CH+
        bladder_pts = muts_df_filtered[muts_df_filtered["Diagnosis"] == "Bladder"]["Patient_id"].unique().shape[0]
        bladder_perc = round((bladder_pts/bladder_pts_n)*100)
        # percentage of the kidney patients that are CH+
        kidney_pts = muts_df_filtered[muts_df_filtered["Diagnosis"] == "Kidney"]["Patient_id"].unique().shape[0]
        kidney_perc = round((kidney_pts/kidney_pts_n)*100)
        # save results
        results_dict[min_vaf] = {"Kidney": kidney_perc, "Bladder": bladder_perc, "Kidney npts": kidney_pts, "Bladder npts": bladder_pts}
        sum_perc = (bladder_pts + kidney_pts)/301
        print(f"Min vaf={min_vaf}, {sum_perc}")
    
    # Now plotting the bar charts
    df = pd.DataFrame.from_dict(results_dict, orient='index').reset_index().rename(columns = {"index": "min_vaf"})
    for i, row in df.iterrows():
        ax.bar(i-0.2, row["Bladder"], color = color_dict["Bladder"], edgecolor = None, width = 0.4)
        ax.bar(i+0.2, row["Kidney"], color = color_dict["Kidney"], edgecolor = None, width = 0.4)
        
        # Annotate the number of patients on top of the bars
        ax.text(i - 0.2, row["Bladder"]+1, str(row["Bladder npts"]), ha='center', va='bottom', fontsize=7, color='black')
        ax.text(i + 0.2, row["Kidney"]+1, str(row["Kidney npts"]), ha='center', va='bottom', fontsize=7, color='black')
    
    ax.set_xlim((-0.7, 2.5))
    ax.set_xticks([0, 1, 2])
    ax.tick_params(axis='x', bottom=False)
    ax2.tick_params(axis='x', pad=2)
    ax.set_xticklabels(["≥0.25%", "≥2%", "≥10%"], fontsize = 8)
    ax.set_ylabel("% of patients with CH")
    ax.set_xlabel("Minimum CH VAF%")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    # ADD LEGEND
    legend_colors = ["deepskyblue", "orangered"]
    legend_labels = ["mUC", "mRCC"]
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.5)
    return(ax)

def plot_per_patient_counts(df, figure_dir = None, figure_name = None, colorby=None, ax = None, bar_color = None, superimpose_swarm = False):
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
    
    if superimpose_swarm:
        ax2 = ax.twinx()
        mutcounts = df["Patient_id"].value_counts().reset_index()
        mutcounts.columns = ["Patient_id", "count"]
        mutcounts = df[["Patient_id", "VAF_t"]].merge(mutcounts)
        for i, group in mutcounts.groupby("count"):
            xpos = group["count"].unique()[0]
            jitter = np.random.uniform(0.05, -0.05, group.shape[0])
            ax2.scatter(np.repeat(xpos, group.shape[0])+jitter, group["VAF_t"], s = 0.5, alpha = 0.6, color = "black")
        # Aes
        ax2.spines["top"].set_visible(False)
        ax2.set_ylabel("Plasma VAF%")
    
    # y tick labels
    # max_y = max(ax.get_yticks())
    # new_max_y = round_to_nearest_10(max_y)
    # ax.set_ylim((0, new_max_y))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))  # Adjust the number of ticks as needed
    
    # Save the figure
    if ax is not None and superimpose_swarm: 
        return(ax, ax2)
    elif ax is not None and not superimpose_swarm: 
        return(ax)
    if figure_dir is not None:
        figure_path = os.path.join(figure_dir, figure_name)
        fig.tight_layout()
        fig.savefig(figure_path)
    

def plot_ctDNA_bins_and_VAF_correlations(PATH_sample_information, ch_df, ctdna_df, ax1, ax2):
    """
    Shows how the correlation between the WBC VAF and cfDNA VAF changes with ctDNA fraction. 
    """
    # group patients based on the vaf of their ctDNA mutation
    ct_fr = ctdna_df.groupby("Patient_id")["VAF_t"].max().reset_index().rename(columns = {"VAF_t": "Max VAF"})
    ct_fr["ct bin"] = ["≥10%" if vaf >= 10 else "<10%" for vaf in ct_fr["Max VAF"]]
    
    # Add the ctDNA negative patients
    ctdna_status = annotate_mutation_status(ctdna_df, "Both", PATH_sample_information, annotate_what = "ctDNA")
    ctdna_status = ctdna_status[ctdna_status["Timepoint"] == "Baseline"]
    ctdna_neg = ctdna_status[ctdna_status["ctDNA status"] == "Negative"]
    ctdna_neg["Max VAF"] = 0
    ctdna_neg["ct bin"] = 0
    ct_fr = pd.concat([ct_fr, ctdna_neg[["Patient_id", "Max VAF", "ct bin"]]])
    
    # Counts, to visualize in the plots
    bin_high_count = ct_fr[ct_fr["ct bin"] == "≥10%"].shape[0]
    bin_low_count = ct_fr[ct_fr["ct bin"] == "<10%"].shape[0]
    bin_zero_count = ct_fr[ct_fr["ct bin"] == 0].shape[0]
    
    import scipy.stats as stats
    ct_frac = pd.read_csv(PATH_mutation_ctfractions)
    ct_frac = ct_frac[ct_frac["Sample_name_t"].str.contains("Baseline")].reset_index(drop = True)
    ct_frac = ct_frac[["Patient_id", "Mutation_ctDNA_fraction"]]
    ct_frac["Mutation_ctDNA_fraction"] = ct_frac["Mutation_ctDNA_fraction"]*100
    
    ch_df = ch_df.merge(ct_fr, how = "left")
    
    # Correlate WBC VAF with cfDNA VAF in each bin
    slopes = {}
    errors = {}
    labels = [0, "<10%", "≥10%"]
    for label in labels:
        bin_data = ch_df[ch_df['ct bin'] == label] # Filter the DataFrame for the current bin
        slope, intercept, r_value, p_value, std_err = stats.linregress(bin_data['VAF_t'], bin_data['VAF_n']) # Perform linear regression
        intercept
        slopes[label] = slope # Store the correlation in the dictionary
        errors[label] = std_err  # Store the standard error in the dictionary
    
    # Convert dictionaries to lists for plotting
    slope_values = [slopes.get(label, 0) for label in labels]
    error_values = [errors.get(label, 0) for label in labels]
    
    # Plot the connected line plot for the pearson
    ax1.scatter([0, 1, 2], slope_values, s=10, color = "black")
    # ax1.plot(labels, slopes.values(), color='black', linewidth=1)
    ax1.errorbar([0, 1, 2], slope_values, yerr=error_values, fmt='o', color='black', capsize=2, linewidth=0.6, capthick=0.6, markersize=1)
    ax1.set_ylabel("Slope")
    ax1.set_xlabel("")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_ylim((0.84, 1.1))
    ax1.set_yticks([1.02, 0.96, 0.9, 0.84])
    ax1.set_yticklabels(["1.02", "0.96", "0.9", "0.84"])
    ax1.set_xticks([])
    ax1.tick_params(axis="both", direction="out", which="both", left=True, bottom=False, labelbottom = False, colors='k')
    
    # Plot the boxplots with scatters
    # Plot the violin plots
    sns.violinplot(x='ct bin', y="Max VAF", data=ch_df, order=[0, '<10%', '≥10%'], ax=ax2, color='dimgray', inner=None, linewidth = 0.3)  # This removes the internal boxplot    
    sns.boxplot(x='ct bin', y="Max VAF", data=ch_df, order=[0, '<10%', '≥10%'], ax=ax2, fliersize = 0, width = 0.025, color='black', 
                medianprops={'color': 'black', "lw": 0.3},
                boxprops={"edgecolor": 'black', 'zorder': 2}, 
                capprops={'linewidth': 0.3, 'zorder': 2},
                whiskerprops={'linewidth': 0.3, 'zorder': 2})
    # sns.stripplot(x='Mutation_ctDNA_bin', y="Mutation_ctDNA_fraction", data=df, order=['0-29', '30-100'], ax=ax2, jitter=True, size = 2, color='black')
    ax2.set_ylabel("ctDNA%")
    ax2.set_xlabel("ctDNA% bins")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.set_yticks([0, 25, 50, 75, 100])
    ax2.set_yticklabels(["0", "25", "50", "75", "100"])
    ax2.tick_params(axis="both", direction="out", which="both", left=True, bottom=True, colors='k')
    
    # Set x ticklabels
    ax2.set_xticks([0, 1, 2])
    x_tick_labels = [f"0\nn={bin_zero_count}", f"<10%\nn={bin_low_count}", f"≥10%\nn={bin_high_count}"]    
    ax2.set_xticklabels(x_tick_labels)
    return(ax1, ax2)

def plot_vaf_scatter(df_chip, ax1, annotate_genes = True, fontsize = 10):
    """
    Plots a scatter plot comparing the cfDNA VAF to WBC VAF.
    """  
    # helps remove dependent calls with vaf 0 for both wbc and cfdna
    if "Dependent" in df_chip.columns:
        df_chip = df_chip[df_chip["Dependent"] == False]
    
    df_chip = df_chip.sort_values(by="VAF_n", ascending=False)
    df_chip["DTA genes"] = df_chip["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])
    df_chip["zorder"] = np.random.randint(1, 3, df_chip.shape[0])
    df_chip["color"] = df_chip["DTA genes"].map({True: "blue", False: "black"})
    
    # Plotting
    sns.regplot(x="VAF_n", y="VAF_t", data=df_chip, ax=ax1, scatter = False, line_kws={'color': 'red', 'linewidth': 0.5, "linestyle": "--"})
    ax1.scatter(df_chip["VAF_n"], df_chip["VAF_t"], s = 1, c = df_chip["color"].tolist(), alpha = 0.6)
    
    # Run spearman's corr
    import scipy.stats as stats
    spearman_corr, p_value = stats.spearmanr(df_chip["VAF_n"], df_chip["VAF_t"])
    spearman_corr = round(spearman_corr, 2)
    p_value_str = f"{p_value:.2e}"
    ########################################################################################################
    # Add an inset ax to show the zoomed in values between 0-2.
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    # inset_ax = inset_axes(ax1, width="35%", height="35%", loc=2, bbox_to_anchor=(0.15, 0.9, 0.35, 0.35), bbox_transform=ax1.transAxes)
    inset_ax = inset_axes(ax1, width="55%", height="55%", loc='upper left', bbox_to_anchor=(0.1, 0.5, 0.5, 0.5), bbox_transform=ax1.transAxes)
    
    # Move the inset plot to the right by adjusting the bbox_to_anchor
    inset_ax.set_position([inset_ax.get_position().x0 + 0.25,  # Shift right by increasing x0
                           inset_ax.get_position().y0,         # Keep the same y0
                           inset_ax.get_position().width, 
                           inset_ax.get_position().height])
    
    inset_ax.tick_params(axis='both', which='major', labelsize=5, width = 0.25, length = 1)  # Set tick size
    inset_ax.spines[['top', 'right']].set_visible(False)
    inset_ax.spines[['bottom', 'left']].set_linewidth(0.25)
    
    df_chip_2_perc = df_chip[(df_chip["VAF_t"] < 2) & (df_chip["VAF_n"] < 2)]
    sns.regplot(x="VAF_n", y="VAF_t", data=df_chip_2_perc, ax=inset_ax, scatter = False, line_kws={'color': 'red', 'linewidth': 0.5, "linestyle": "--"})
    inset_ax.scatter(df_chip_2_perc["VAF_n"], df_chip_2_perc["VAF_t"], s = 0.5, c = df_chip_2_perc["color"].tolist(), alpha = 0.6)
    inset_ax.set_ylim((0, 2))
    inset_ax.set_xlim((0, 2))
    inset_ax.set_yticks([0, 1, 2])
    inset_ax.set_xticks([0, 1, 2])
    inset_ax.set_yticklabels(["0", "1", "2"], fontsize = 6)
    inset_ax.set_xticklabels(["0", "1", "2"], fontsize = 6)
    inset_ax.set_ylabel("")
    inset_ax.set_xlabel("")
    inset_ax.xaxis.set_tick_params(pad=2)  # Reduces padding for x-ticks
    inset_ax.yaxis.set_tick_params(pad=2)  # Reduces padding for y-ticks
    ########################################################################################################
    
    # ax1.text(0.45, 0.85, f'Spearman\'s\ncorr: {spearman_corr:.2f}\np={p_value_str}', transform=ax1.transAxes, fontsize=8, color='red')
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.set_xlabel("WBC VAF%", fontsize=fontsize)
    ax1.set_ylabel("Plasma cfDNA VAF%", fontsize=fontsize)
    
    if annotate_genes:
        top_points = df_chip.nlargest(3, 'VAF_t').reset_index(drop = True)
        texts = []
        line_length = 10
        for i, row in top_points.iterrows():
            y_pos_text = row["VAF_t"]-3
            if i < 2:
                # Top 2 points will have a line extending to the left
                x_right = row["VAF_n"]
                x_left = x_right-line_length
                ax1.text(x_left, y_pos_text, f'$\it{{{row["Gene"]}}}$', fontsize=8, color='black', ha = "right")
            else:
                x_left = row["VAF_n"]
                x_right = x_right+line_length
                ax1.text(x_right, y_pos_text, f'$\it{{{row["Gene"]}}}$', fontsize=8, color='black', ha = "left")
            ax1.plot([x_left, x_right], [row["VAF_t"], row["VAF_t"]], color='black', linewidth=0.5)
    
        # Couple of misc points to label, do it manually
        ax1.text(44.8+10, 34.3-3, f'$\it{"DNMT3A"}$', fontsize=8, color='blue', ha = "left")
        ax1.plot([44.8, 44.8+10], [34.3, 34.3], color='blue', linewidth=0.5)
        
        ax1.text(38.55-10, 41.48-3, f'$\it{"ASXL1"}$', fontsize=8, color='blue', ha = "right")
        ax1.plot([38.55, 38.55-10], [41.48, 41.48], color='blue', linewidth=0.5)
    
    ax1.set_ylim([0, 60])
    ax1.set_xlim([0, 60])
    ax1.set_yticks([0, 20, 40, 60])
    ax1.set_yticklabels(["0", "20", "40", "60"], fontsize = fontsize)
    ax1.set_xticks([0, 20, 40, 60])
    ax1.set_xticklabels(["0", "20", "40", "60"], fontsize = fontsize)
    
    # add legend for the colors
    legend_colors = ["blue", "black"]
    legend_labels = ["DTA", "non-DTA"]
    legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=2, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax1.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize = 8, handlelength=2, handletextpad=0.1)
    return(ax1)


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

def age_vs_CH_presence_grouped_ch_thresholds(df_CH, age_df, PATH_sample_information, ax, fontsize = 10):
    """
    VERY SIMILAR TO age_vs_CH_presence. ONLY DIFFERENCE IS PLOTS BLADDER AND KIDNEY TOGETHER.
    Boxplot. Plotting the ages on the y-axis of patients, grouped into CH+ and CH-. Give a df consisting of only one timepoint, baseline or OT. 
    test_to_use: use either 'MWU' or 'T test'
    """
    dict_box_colors = {"Bladder": (143/255, 215/255, 239/255, 255/255), "Kidney": (239/255, 169/255, 143/255, 255/255)}
    dict_scatter_colors = {"Bladder": "deepskyblue", "Kidney": "orangered"}
    dict_offset = {"Bladder": -0.2, "Kidney": 0.2}
    
    # Data manipulation
    timepoint = df_CH["Timepoint"].unique()[0]
    all_pts = pd.read_csv(PATH_sample_information, sep="\t", names=["Patient_id", "Date", "Diagnosis", "Timepoint"])
    all_pts = all_pts[(all_pts["Timepoint"] == timepoint)].reset_index(drop = True)
    
    all_ch = annotate_mutation_status(df_CH, "Both", PATH_sample_information, annotate_what = "CHIP")
    all_ch = all_ch[all_ch["Timepoint"] == "Baseline"].reset_index(drop = True).merge(age_df, how = "left").dropna()
    ch_neg = all_ch[(all_ch["CHIP status"] == "Negative")]
    
    max_ch_df = df_CH[["Patient_id", "VAF_n"]].groupby("Patient_id")["VAF_n"].max().reset_index()
    
    ch1 = maxvaf_df[maxvaf_df["VAF_n"] < 2].merge(age_df).assign(**{"CHIP status": "Positive"})
    ch2 = maxvaf_df[(maxvaf_df["VAF_n"] >= 2) & (maxvaf_df["VAF_n"] < 10)].merge(age_df).assign(**{"CHIP status": "Positive"})
    ch3 = maxvaf_df[maxvaf_df["VAF_n"] >= 10].merge(age_df).assign(**{"CHIP status": "Positive"})
    
    flierprops_dict = dict(marker='o', markersize=5, markeredgecolor='black', linestyle='None')
    whiskerprops_dict =dict(color='black')
    medianprops_dict = dict(color='black')
    capprops_dict = dict(color='black')
    
    # Plot the boxes
    x_tick_labels_list = []
    for i, df in enumerate([ch_neg, ch1, ch2, ch3]):
        df["CHIP status"] = df["CHIP status"].replace({"Negative": "CH-", "Positive": "CH+"})
        if i != 0:
            df = df[df["CHIP status"] == "CH+"]
        for diagnosis in ["Bladder", "Kidney"]:
            df_diagnosis = df[df["Diagnosis"] == diagnosis]
            boxprops_dict = dict(facecolor=dict_box_colors[diagnosis], edgecolor='black', linewidth = 0.7)  
            boxplot = ax.boxplot(df_diagnosis["Age"], positions = [i+dict_offset[diagnosis]], flierprops = flierprops_dict, boxprops = boxprops_dict, medianprops = medianprops_dict, capprops = capprops_dict, widths = 0.3, showfliers = False, patch_artist = True)
            ax.scatter(np.random.uniform(i+dict_offset[diagnosis]-0.08, i+dict_offset[diagnosis]+0.08, len(df_diagnosis["Age"])), df_diagnosis["Age"], s = 1, color = dict_scatter_colors[diagnosis], alpha = 1, zorder = 100)
            n_patients = len(df_diagnosis["Age"])
            x_tick_labels_list.append(f"n={n_patients}")
    # Aes
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Age at baseline draw")
    ax.set_ylim((0, 120))
    ax.set_yticks([20, 40, 60, 80, 100])
    ax.set_yticklabels(["20", "40", "60", "80", "100"])
    ax.set_xticklabels(x_tick_labels_list, fontsize = 7)
    
    # add text
    # text_pos_dict = {0: "CH-", 1: "0.25%-2%", 2: "2%-10%", 3: ">10%"}
    # for index, (key, value) in enumerate(text_pos_dict.items()):
    #     ax.text(key, 0, value, ha='center', va='top', fontsize=9)
    
    # do MWU between various groups
    for diagnosis in ["Kidney", "Bladder"]:
        age_ch_neg = ch_neg[ch_neg["Diagnosis"] == diagnosis]["Age"]
        age_ch1 = ch1[ch1["Diagnosis"] == diagnosis]["Age"]
        age_ch2 = ch2[ch2["Diagnosis"] == diagnosis]["Age"]
        age_ch3 = ch3[ch3["Diagnosis"] == diagnosis]["Age"]
        
        mwu1 = "{:.2e}".format(mannwhitneyu(age_ch_neg, age_ch1)[1])
        mwu2 = "{:.2e}".format(mannwhitneyu(age_ch_neg, age_ch2)[1])
        mwu3 = "{:.2e}".format(mannwhitneyu(age_ch_neg, age_ch3)[1])     
        
        # np.median(age_ch_neg)
        # np.median(age_ch1)
        # np.median(age_ch2)
        # np.median(age_ch3)
        
        print(diagnosis)
        print(f"CH neg vs 0.25%-2% p={mwu1}")
        print(f"CH neg vs 2%-10% p={mwu2}")
        print(f"CH neg vs >10% p={mwu3}")
        # "{:.2e}".format(mannwhitneyu(age_ch1, age_ch3)[1])
        # "{:.2e}".format(mannwhitneyu(age_ch_neg, age_ch3)[1])
    
    # Add legend
    legend_colors = ["deepskyblue", "orangered"]
    legend_labels = ["mUC", "mRCC"]
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize = 8, handlelength=2, handletextpad = 0.5)
    
    return(ax)


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

def return_median_survival(surv_df_merged, stratify_by, event_col, duration_col):
    """
    Returns median survival values.
    """
    import math
    kmf_positive = KaplanMeierFitter()
    kmf_negative = KaplanMeierFitter()
    groups = surv_df_merged[stratify_by]
    ix_positive = (groups == 'Positive') | (groups == 1)
    
    # Negative group
    kmf_negative.fit(durations=surv_df_merged[duration_col][~ix_positive], event_observed=surv_df_merged[event_col][~ix_positive], label=f"{stratify_by}(-)")
    kmf_positive.fit(durations=surv_df_merged[duration_col][ix_positive], event_observed=surv_df_merged[event_col][ix_positive], label=f"{stratify_by}(+)")
    
    if kmf_positive.median_survival_time_ == math.inf: 
        median_survival_positive = "Not reached"
    else:
        median_survival_positive = round(kmf_positive.median_survival_time_)
    if kmf_negative.median_survival_time_ == math.inf:
        median_survival_negative = "Not reached"
    else:
        median_survival_negative = round(kmf_negative.median_survival_time_)
   
    # check if the median survivals are significantly different using log rank test
    if median_survival_negative != "Not reached" and median_survival_positive != "Not reached":
        results = logrank_test(surv_df_merged[duration_col][~ix_positive], surv_df_merged[duration_col][ix_positive], event_observed_A=surv_df_merged[event_col][~ix_positive], event_observed_B=surv_df_merged[event_col][ix_positive])
        p_value_logrank = f"{results.p_value:.2e}"
        dict_result = {"Positive": median_survival_positive, "Negative": median_survival_negative, "logrank p": p_value_logrank}
    else:
        dict_result = {"Positive": median_survival_positive, "Negative": median_survival_negative}
    return dict_result

def make_survival_curve(surv_df_merged, stratify_by, event_col, duration_col, output_path, plot_title, print_HRs = True, pos_curve_label_positions = None, neg_curve_label_positions = None, xmax = None, add_legend = True, print_survival = True, xlabel = None, ax = None):
    """
    stratify_by = Name of the column in df to stratify the curves by.
    """
    import math
    label_keyword = re.search("(ctDNA|CHIP|CH)", stratify_by, flags=re.IGNORECASE).group()
    kmf_positive = KaplanMeierFitter()
    kmf_negative = KaplanMeierFitter()
    groups = surv_df_merged[stratify_by]
    ix_positive = (groups == 'Positive') | (groups == 1)
    
    neg_curve_color = "black"
    if label_keyword == "ctDNA": 
        pos_curve_color = "forestgreen"
    else:
        pos_curve_color = "red"
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    
    # Negative group
    kmf_negative.fit(durations=surv_df_merged[duration_col][~ix_positive], event_observed=surv_df_merged[event_col][~ix_positive], label=f"{label_keyword}(-)")
    kmf_negative.plot_survival_function(ax=ax, ci_show=False, show_censors=True, linewidth=0.8, censor_styles={'marker': '|', 'ms': 5, 'markerfacecolor': neg_curve_color, 'mew': 0.8}, color = neg_curve_color)
    
    # Positive group
    kmf_positive.fit(durations=surv_df_merged[duration_col][ix_positive], event_observed=surv_df_merged[event_col][ix_positive], label=f"{label_keyword}(+)")
    kmf_positive.plot_survival_function(ax=ax, ci_show=False, show_censors=True, linewidth=0.8, censor_styles={'marker': '|', 'ms': 5, 'markerfacecolor': pos_curve_color, 'mew': 0.8}, color = pos_curve_color)
        
    if print_survival is not None:
        if isinstance(print_survival, bool) and print_survival:
            ax_to_plot = ax
        elif isinstance(print_survival, plt.Axes):
            ax_to_plot = print_survival
            if kmf_positive.median_survival_time_ == math.inf: 
                median_survival_positive = "Not reached"
            else:
                median_survival_positive = round(kmf_positive.median_survival_time_)
            if kmf_negative.median_survival_time_ == math.inf:
                median_survival_negative = "Not reached"
            else:
                median_survival_negative = round(kmf_negative.median_survival_time_)
                # check if the median survivals are significantly different using log rank test
                results = logrank_test(surv_df_merged[duration_col][~ix_positive], surv_df_merged[duration_col][ix_positive], event_observed_A=surv_df_merged[event_col][~ix_positive], event_observed_B=surv_df_merged[event_col][ix_positive])
                p_value_logrank = f"{results.p_value:.2e}"
                ax_to_plot.text(0.3, 0.37, f"Log rank p={p_value_logrank}", transform=ax_to_plot.transAxes, fontsize=10, ha='left', va='center')
            
            ax_to_plot.text(1.1, 0.8, "Median", transform=ax_to_plot.transAxes, fontsize=10, ha='center', va='center')
            ax_to_plot.text(1.1, 0.67, "mo", transform=ax_to_plot.transAxes, fontsize=10, ha='center', va='center', fontstyle = "italic")
            
            ax_to_plot.text(0.3, 0.57, f'{label_keyword}+', transform=ax_to_plot.transAxes, fontsize=10, ha='left', va='center')
            ax_to_plot.text(0.3, 0.47, f'{label_keyword}-', transform=ax_to_plot.transAxes, fontsize=10, ha='left', va='center')
            
            ax_to_plot.text(1.1, 0.57, f'{median_survival_positive}', transform=ax_to_plot.transAxes, fontsize=10, ha='center', va='center')
            ax_to_plot.text(1.1, 0.47, f'{median_survival_negative}', transform=ax_to_plot.transAxes, fontsize=10, ha='center', va='center')
            
            ax_to_plot.text(0.3, 0.17, "Hazard ratio,", transform=ax_to_plot.transAxes, fontsize=10, ha='left', va='center')
    
    # ax.annotate(f'p value: {p_value}', xy=(1, 0.5), xytext=(100, 80), textcoords='offset points', fontsize=12)
    if xlabel is None:
        ax.set_xlabel("Months since cfDNA collection")
    else: 
        ax.set_xlabel(xlabel)
    ax.set_ylabel("Survival fraction")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc='best', frameon=False)
    ax.set_ylim((0, 1))
    ax.set_yticks((0, 0.2, 0.4, 0.6, 0.8, 1))
    ax.set_yticklabels(("0", "0.2", "0.4", "0.6", "0.8", "1"))
    if xmax is None:
        xmax = ax.get_xlim()[1]
    ax.set_xlim((0, xmax))
    ax.set_xticks(np.arange(0, xmax + 10, 10))
    if xmax >= 100:
        ax.set_xticks(np.arange(0, xmax + 10, 20))
    if xmax < 40:
        ax.set_xticks(np.arange(0, xmax + 10, 5))
    
    # add risk counts
    add_at_risk_counts(kmf_positive, kmf_negative, ax=ax, rows_to_show = ["At risk"])
    
    if not add_legend:
        ax.legend_.remove()
    ax.set_title(plot_title)
    
    # label the curves
    if neg_curve_label_positions is not None and pos_curve_label_positions is not None:
        ax.text(neg_curve_label_positions[0], 
                   neg_curve_label_positions[1], 
                   f'{label_keyword}-', transform=ax.transAxes, fontsize=10, ha='center', va='center', color = neg_curve_color)
        ax.text(pos_curve_label_positions[0], 
                   pos_curve_label_positions[1], 
                   f'{label_keyword}+', transform=ax.transAxes, fontsize=10, ha='center', va='center', color = pos_curve_color)
    
    # Run CPH here and print the hazard ratios to the plot.
    if print_HRs:
        try:
            ax = run_cox_proportional_hazards(surv_df_merged, duration_col, event_col, stratify_by, ax_to_plot)
        except Exception as e:
            print(f"An error occurred while running Cox Proportional Hazards.")
    
    if output_path is not None:
        fig.tight_layout()
        fig.savefig(output_path)
    return ax

def run_cox_proportional_hazards(surv_df_merged, duration_col, event_col, stratify_by, ax = None):
    """
    Fits the CPH model to survival data.
    duration_col: Column that has OS data.
    event_col: Column that has indicates if the event of interest has happened.
    stratify_by: The categorical variable that we use to split the curves into two.
    ax: Axis to plot hazard ratios if provided.
    """
    cph = CoxPHFitter()
    surv_df_merged= surv_df_merged[[duration_col, event_col, stratify_by]]
    # surv_df_merged[stratify_by]=surv_df_merged[stratify_by].astype(bool)
    # surv_df_merged = pd.get_dummies(surv_df_merged, columns=[stratify_by], drop_first=True)
    cph.fit(surv_df_merged[[duration_col, event_col, stratify_by]], duration_col=duration_col, event_col=event_col)
    
    # Get hazard ratios, confidence intervals, and p-values
    hazard_ratio = cph.summary['exp(coef)'][stratify_by]
    ci_lower = cph.summary['exp(coef) lower 95%'][stratify_by]
    ci_upper = cph.summary['exp(coef) upper 95%'][stratify_by]
    p_value = cph.summary['p'][stratify_by]
    p_value = f"{p_value:.2e}"
    
    # Print hazard ratios, confidence intervals, and p-values on the plot
    if ax is None:
        result_dict = {"HR": hazard_ratio, "CI_upper": ci_upper, "CI_lower": ci_lower, "p": p_value}
        return(result_dict)
    else:
        return(ax)

def make_survival_curve_with_forest(merged_df, stratify_by, ax_km, ax_cox, pos_curve_label_positions = None, neg_curve_label_positions = None, add_legend = False):
    """
    stratify_by = Name of the column in df to stratify the curves by.
    pos_curve_label_positions: x and y positions to label the curve
    neg_curve_label_positions: x and y positions to label the curve
    """
    label_keyword = re.search("(ctDNA+|CH+)", stratify_by, flags=re.IGNORECASE).group()
    kmf_positive = KaplanMeierFitter()
    kmf_negative = KaplanMeierFitter()
    groups = merged_df[stratify_by]
    ix_positive = (groups == True)
    
    neg_curve_color = "black"
    if label_keyword == "ctDNA": 
        pos_curve_color = "forestgreen"
    else:
        pos_curve_color = "red"
    
    # fig = plt.figure(figsize=(8, 4))
    # gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.7], hspace=0.4)
    # ax_km = fig.add_subplot(gs[0])
    # ax_cox = fig.add_subplot(gs[1])
    
    # Negative group
    kmf_negative.fit(durations=merged_df["Overall survival"][~ix_positive], event_observed=merged_df["Death at last follow up"][~ix_positive], label=f"{label_keyword}(-)")
    kmf_negative.plot_survival_function(ax=ax_km, ci_show=False, show_censors=True, censor_styles={'marker': '|', 'ms': 7, 'markerfacecolor': neg_curve_color}, color = neg_curve_color)
    
    # Positive group
    kmf_positive.fit(durations=merged_df["Overall survival"][ix_positive], event_observed=merged_df["Death at last follow up"][ix_positive], label=f"{label_keyword}(+)")
    kmf_positive.plot_survival_function(ax=ax_km, ci_show=False, show_censors=True, censor_styles={'marker': '|', 'ms': 7, 'markerfacecolor': pos_curve_color}, color = pos_curve_color)
    ax_km.set_xticks([0, 10, 20, 30, 40, 50, 60])
    ax_km.set_xticklabels(["0", "10", "20", "30", "40", "50", "60"])
    # remove legend
    if not add_legend:
        ax_km.legend_.remove()
    
    # Calculate median survival for negative group and show on the plot
    median_survival_negative = kmf_negative.median_survival_time_
    median_survival_positive = kmf_positive.median_survival_time_
    
    # check if the median survivals are significantly different using log rank test
    results = logrank_test(merged_df["Overall survival"][~ix_positive], merged_df["Overall survival"][ix_positive], event_observed_A=merged_df["Death at last follow up"][~ix_positive], event_observed_B=merged_df["Death at last follow up"][ix_positive])
    p_value_logrank = round(results.p_value, 4)
    
    # Show on the plot
    ax_km.text(0.9, 0.95, "Median survival (months)", transform=ax_km.transAxes, fontsize=8, ha='right', va='center')
    ax_km.text(0.9, 0.9, f'{label_keyword}(-): {median_survival_negative}', transform=ax_km.transAxes, fontsize=8, ha='right', va='center')
    ax_km.text(0.9, 0.85, f'{label_keyword}(+): {median_survival_positive}', transform=ax_km.transAxes, fontsize=8, ha='right', va='center')
    ax_km.text(0.9, 0.8, f'Log rank p = {p_value_logrank}', transform=ax_km.transAxes, fontsize=8, ha='right', va='center')
    
    # add risk counts
    add_at_risk_counts(kmf_positive, kmf_negative, ax=ax_km, rows_to_show = ["At risk"])    
    
    # ax.annotate(f'p value: {p_value}', xy=(1, 0.5), xytext=(100, 80), textcoords='offset points', fontsize=12)
    ax_km.set_xlabel("OS since cfDNA collection (mo.)")
    ax_km.set_ylabel("Survival Probability")
    ax_km.spines[["right", "top"]].set_visible(False)
    # ax_km.legend(loc='best', frameon=False)
    # ax_km.set_title("Overall survival")
    
    # label the curves
    if neg_curve_label_positions is not None and pos_curve_label_positions is not None:
        ax_km.text(neg_curve_label_positions[0], 
                   neg_curve_label_positions[1], 
                   f'{label_keyword}(-)', transform=ax_km.transAxes, fontsize=10, ha='center', va='center', color = neg_curve_color)
        ax_km.text(pos_curve_label_positions[0], 
                   pos_curve_label_positions[1], 
                   f'{label_keyword}(+)', transform=ax_km.transAxes, fontsize=10, ha='center', va='center', color = pos_curve_color)
    
    # Run a univariate cox analysis to print the hazard ratios to the KM plot
    cph = CoxPHFitter()
    cph.fit(merged_df[["Overall survival", "Death at last follow up", stratify_by]], duration_col="Overall survival", event_col="Death at last follow up")
    
    # Get hazard ratios, confidence intervals, and p-values
    hazard_ratio = cph.summary['exp(coef)'][stratify_by]
    ci_lower = cph.summary['exp(coef) lower 95%'][stratify_by]
    ci_upper = cph.summary['exp(coef) upper 95%'][stratify_by]
    p_value = cph.summary['p'][stratify_by]
    
    # Print hazard ratios, confidence intervals, and p-values on the plot
    ax_km.text(0.9, 0.7, f"HR={hazard_ratio:.2f} [{ci_lower:.2f}, {ci_upper:.2f}]\np={p_value:.3f}", transform=ax_km.transAxes, fontsize=8, ha='right', va='center')
    
    # Run CPH here and print the hazard ratios to the second ax.
    # cph = CoxPHFitter()
    # cph.fit(merged_df, duration_col='Overall survival', event_col='Death at last follow up')
    # cph.plot(hazard_ratios=True, ax=ax_cox)
    # ax_cox.spines[["right", "top"]].set_visible(False)
    # ax_cox.set_xlabel("HR (95% CI)")
    # ax_cox.yaxis.set_tick_params(labelsize=6)
    return ax_km

def make_survival_curve_with_forest_3_groups(merged_df, ax_km, event_col, duration_col, diagnosis = None, pos_curve_label_positions = None, neg_curve_label_positions = None, add_legend = False):
    """
    stratify_by = Name of the column in df to stratify the curves by.
    pos_curve_label_positions: x and y positions to label the curve
    neg_curve_label_positions: x and y positions to label the curve
    """
    
    kmf_ctdna_neg = KaplanMeierFitter()
    kmf_ctdna_pos_ch_neg = KaplanMeierFitter()
    kmf_ctdna_pos_ch_pos = KaplanMeierFitter()
    
    ix_ctdna_neg = merged_df["ctDNA-"] == True
    ix_ctdna_pos_ch_neg = merged_df["ctDNA+ CH-"] == True
    ix_ctdna_pos_ch_pos = merged_df["ctDNA+ CH+"] == True
    
    curve_color1 = "black"
    curve_color2 = "darkorange"
    curve_color3 = "blue"
    
    # Negative group
    kmf_ctdna_neg.fit(durations=merged_df[duration_col][ix_ctdna_neg], event_observed=merged_df[event_col][ix_ctdna_neg], label="ctDNA-")
    kmf_ctdna_neg.plot_survival_function(ax=ax_km, ci_show=False, show_censors=True, linewidth=0.8, censor_styles={'marker': '|', 'ms': 5, 'markerfacecolor': curve_color1, 'mew': 0.8}, color = curve_color1)
    
    kmf_ctdna_pos_ch_neg.fit(durations=merged_df[duration_col][ix_ctdna_pos_ch_neg], event_observed=merged_df[event_col][ix_ctdna_pos_ch_neg], label="ctDNA+ CH-")
    kmf_ctdna_pos_ch_neg.plot_survival_function(ax=ax_km, ci_show=False, show_censors=True, linewidth=0.8, censor_styles={'marker': '|', 'ms': 5, 'markerfacecolor': curve_color2, 'mew': 0.8}, color = curve_color2)
    
    kmf_ctdna_pos_ch_pos.fit(durations=merged_df[duration_col][ix_ctdna_pos_ch_pos], event_observed=merged_df[event_col][ix_ctdna_pos_ch_pos], label="ctDNA+ CH+")
    kmf_ctdna_pos_ch_pos.plot_survival_function(ax=ax_km, ci_show=False, show_censors=True, linewidth=0.8, censor_styles={'marker': '|', 'ms': 5, 'markerfacecolor': curve_color3, 'mew': 0.8}, color = curve_color3)
    
    # remove legend
    if not add_legend:
        ax_km.legend_.remove()
    
    # Calculate median survival for negative group and show on the plot
    median_survival_ch = kmf_ctdna_neg.median_survival_time_
    median_survival_ctdna = kmf_ctdna_pos_ch_neg.median_survival_time_
    median_survival_both = kmf_ctdna_pos_ch_pos.median_survival_time_
    
    # check if the median survivals are significantly different using log rank test
    results = logrank_test(merged_df[duration_col][ix_ctdna_pos_ch_neg], merged_df[duration_col][ix_ctdna_pos_ch_pos], event_observed_A=merged_df[event_col][ix_ctdna_pos_ch_neg], event_observed_B=merged_df[event_col][ix_ctdna_pos_ch_pos])
    p_value_logrank = round(results.p_value, 4)
    
    # Show on the plot
    # ax_km.text(0.9, 0.95, "Median survival (months)", transform=ax_km.transAxes, fontsize=8, ha='right', va='center')
    # ax_km.text(0.9, 0.9, f'{label_keyword}(-): {median_survival_negative}', transform=ax_km.transAxes, fontsize=8, ha='right', va='center')
    # ax_km.text(0.9, 0.85, f'{label_keyword}(+): {median_survival_positive}', transform=ax_km.transAxes, fontsize=8, ha='right', va='center')
    # ax_km.text(0.9, 0.8, f'Log rank p = {p_value_logrank}', transform=ax_km.transAxes, fontsize=8, ha='right', va='center')
    
    # xticks
    if diagnosis == "Bladder":
        ax_km.set_xticks([0, 10, 20, 30, 40, 50])
        ax_km.set_xticklabels(["0", "10", "20", "30", "40", "50"])
    elif diagnosis == "Kidney":
        ax_km.set_xticks([0, 10, 20, 30, 40, 50, 60, 70])
        ax_km.set_xticklabels(["0", "10", "20", "30", "40", "50", "60", "70"])
    
    # add risk counts
    add_at_risk_counts(kmf_ctdna_neg, kmf_ctdna_pos_ch_neg, kmf_ctdna_pos_ch_pos, ax=ax_km, rows_to_show = ["At risk"])    
    
    # ax.annotate(f'p value: {p_value}', xy=(1, 0.5), xytext=(100, 80), textcoords='offset points', fontsize=12)
    ax_km.set_xlabel("OS since cfDNA collection (mo.)")
    ax_km.set_ylabel("Survival Probability")
    ax_km.spines[["right", "top"]].set_visible(False)
    # ax_km.legend(loc='best', frameon=False)
    # ax_km.set_title("Overall survival", ha = "left")
    ax_km.text(0.0, 1.02, "ctDNA-", color=curve_color1, transform=ax_km.transAxes, fontsize=10, ha='left')
    ax_km.text(0.0, 0.97, "ctDNA+ CH-", color=curve_color2, transform=ax_km.transAxes, fontsize=10, ha='left')
    ax_km.text(0.0, 0.92, "ctDNA+ CH+", color=curve_color3, transform=ax_km.transAxes, fontsize=10, ha='left')
    
    # label the curves
    if neg_curve_label_positions is not None and pos_curve_label_positions is not None:
        ax_km.text(neg_curve_label_positions[0], 
                   neg_curve_label_positions[1], 
                   f'{label_keyword}(-)', transform=ax_km.transAxes, fontsize=10, ha='center', va='center', color = neg_curve_color)
        ax_km.text(pos_curve_label_positions[0], 
                   pos_curve_label_positions[1], 
                   f'{label_keyword}(+)', transform=ax_km.transAxes, fontsize=10, ha='center', va='center', color = pos_curve_color)
    
    return(ax_km)

def plot_cph_forest(merged_df, ax_cox, event_col, duration_col):
    """
    Runs CPH and plots the forest.
    """
    cph = CoxPHFitter()
    cph.fit(merged_df, duration_col=duration_col, event_col=event_col)
    
    # Extract hazard ratios and confidence intervals
    hazard_ratios = cph.hazard_ratios_
    conf_int = cph.confidence_intervals_
    
    # # Manually plot hazard ratios and confidence intervals
    # for i, variable in enumerate(hazard_ratios.index):
    #     hr = hazard_ratios[variable]
    #     ci_low = conf_int.loc[variable, '95% lower-bound']
    #     ci_high = conf_int.loc[variable, '95% upper-bound']
        
    #     # Plot the hazard ratio as a point
    #     ax_cox.plot([hr, hr], [i - 0.4, i + 0.4], 'ro-', markersize=8, label=variable if i == 0 else "")
    #     ax_cox.plot([ci_low, ci_high], [i, i], 'k-', lw=2)
    cph.plot(hazard_ratios=True, ax=ax_cox)
    
    # Customize the plot
    # ax_cox.set_yticks(range(len(hazard_ratios)))
    # ax_cox.set_yticklabels(hazard_ratios.index, fontsize = 7)
    ax_cox.set_xlabel("HR (95% CI)")
    # ax_cox.set_title("")
    ax_cox.yaxis.set_tick_params(labelsize=6)
    
    # Remove right and top spines
    ax_cox.spines[["right", "top"]].set_visible(False)  
    return ax_cox, cph.summary


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
    mutations_df = mutations_df[["Patient_id", annotate_what + " status"]].drop_duplicates().reset_index(drop = True)
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
    mut_status = mut_status[mut_status["Timepoint"] == "Baseline"]
    
    df = mut_status.merge(sex_df, how = "left").dropna(subset=['Sex'])
    
    contingency_table = pd.crosstab(df[f'{annotate_what} status'], df['Sex'])
    if contingency_table.shape == (2, 2):
        odds_ratio, p_value = fisher_exact(contingency_table)
        if p_value < 0.05:
            print(annotate_what)
            print(f'Odds Ratio: {odds_ratio}')
            print(f'P-value: {p_value}')
        return(p_value)

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

def ch_presence_absence_bars(baseline_chip, PATH_sample_information, sex_df, smoking_df, initial_diagnosis_df, bladder_staging, dir_figures, fig_name, gs, title, show_legend = True, show_ax_titles= True, show_x_ticks = True):
    """
    Simple stacked bar charts showing the fraction of patients at baseline that are CH+ and CH-. 
    Bars are colored according to sex or diagnosis.
    sex_df = two cols, patient id and sex.
    """
    # fig = plt.figure(figsize=(6, 2.8))
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
    ax4 = plt.subplot(gs[3])
    ax5 = plt.subplot(gs[4])
    
    base_mut_status = annotate_mutation_status(baseline_chip, "Both", PATH_sample_information, annotate_what="CHIP", annotate_gene=False)
    base_mut_status = base_mut_status[base_mut_status["Timepoint"] == "Baseline"]
    base_mut_status["Diagnosis"] = base_mut_status["Diagnosis"].replace({"Bladder": "mUC", "Kidney": "mRCC"})
    diagnosis_counts = base_mut_status[["Diagnosis", "CHIP status"]].value_counts().reset_index().rename(columns={0: "Count"})
    sex_counts = base_mut_status.merge(sex_df)[["CHIP status", "Sex"]].value_counts().reset_index().rename(columns={0: "Count"})
    smoking_counts = base_mut_status.merge(smoking_df)[["CHIP status", "Previous smoking history"]].value_counts().reset_index().rename(columns={0: "Count"})
    initial_diagnosis_counts = base_mut_status.merge(initial_diagnosis_df)[["CHIP status", "Disease stage at initial diagnosis"]].value_counts().reset_index().rename(columns={0: "Count"})
    staging_counts = base_mut_status.merge(bladder_staging)[["CHIP status", "Disease stage at initial diagnosis"]].value_counts().reset_index().rename(columns={0: "Count"})
    
    # Pivot for plotting
    diagnosis_counts_pivot = diagnosis_counts.pivot(index='Diagnosis', columns='CHIP status', values='Count').fillna(0)
    sex_counts_pivot = sex_counts.pivot(index='Sex', columns='CHIP status', values='Count').fillna(0)
    smoking_counts_pivot = smoking_counts.pivot(index='Previous smoking history', columns='CHIP status', values='Count').fillna(0)
    initial_diagnosis_counts_pivot = initial_diagnosis_counts.pivot(index='Disease stage at initial diagnosis', columns='CHIP status', values='Count').fillna(0)
    staging_counts_pivot = staging_counts.pivot(index='Disease stage at initial diagnosis', columns='CHIP status', values='Count').fillna(0)
    
    # Fisher's exact
    oddsratio, p_value_sex = fisher_exact(sex_counts_pivot)
    oddsratio, p_value_diagnosis = fisher_exact(diagnosis_counts_pivot)
    oddsratio, p_value_smoking = fisher_exact(smoking_counts_pivot)
    oddsratio, p_value_initial_diagnosis = fisher_exact(initial_diagnosis_counts_pivot)
    oddsratio, p_value_staging = fisher_exact(staging_counts_pivot)
    
    # Colors for CHIP status
    colors = {"Positive": "black", "Negative": "silver"}
    
    # Stacked bar chart for diagnosis
    bottom = np.zeros(len(diagnosis_counts_pivot))
    for status in diagnosis_counts_pivot.columns:
        bars = ax1.bar(diagnosis_counts_pivot.index, diagnosis_counts_pivot[status], bottom=bottom, label=status, color=colors[status], width=0.7)
        for bar, count, color in zip(bars, diagnosis_counts_pivot[status], [colors[status] for _ in range(len(bars))]):
            text_color = 'white' if color == 'black' else 'black'
            # ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height()-10, f'{count}', ha='center', va='center', color=text_color, fontsize=6)
        bottom += diagnosis_counts_pivot[status]
    
    # add p values for Fisher's exact
    ax1.text(0.5, 1, f'p = {p_value_sex:.2f}', ha='center', va='center', color="black", fontsize=7, transform=ax1.transAxes)
    ax2.text(0.5, 1, f'p = {p_value_diagnosis:.2f}', ha='center', va='center', color="black", fontsize=7, transform=ax2.transAxes)
    ax3.text(0.5, 1, f'p = {p_value_initial_diagnosis:.2f}', ha='center', va='center', color="black", fontsize=7, transform=ax3.transAxes)
    ax4.text(0.5, 1, f'p = {p_value_smoking:.2f}', ha='center', va='center', color="black", fontsize=7, transform=ax4.transAxes)
    ax5.text(0.5, 1, f'p = {p_value_staging:.2f}', ha='center', va='center', color="black", fontsize=7, transform=ax5.transAxes)
    
    if show_ax_titles:
        ax1.set_title(f'Cancer\n{title}', loc = "left", fontsize = 8)
    else:
        ax1.set_title(f'{title}', loc = "left", fontsize = 8)
    ax1.set_xlabel('')
    ax1.set_ylabel('Number of patients')
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.set_xticklabels(["mRCC", "mUC"], fontsize = 8)
    ax1.set_yticks([0, 100, 200])
    ax1.set_yticklabels(["0", "100", "200"], fontsize = 8)
    
    # Stacked bar chart for sex
    bottom = np.zeros(len(sex_counts_pivot))
    for status in sex_counts_pivot.columns:
        bars = ax2.bar(sex_counts_pivot.index, sex_counts_pivot[status], bottom=bottom, label=status, color=colors[status], width = 0.7)
        for bar, count, color in zip(bars, sex_counts_pivot[status], [colors[status] for _ in range(len(bars))]):
            text_color = 'white' if color == 'black' else 'black'
            # ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height()-10, f'{count}', ha='center', va='center', color=text_color, fontsize=6)
        bottom += sex_counts_pivot[status]
    
    if show_ax_titles:
        ax2.set_title('Sex',  loc = "left", fontsize = 8)
    ax2.set_xlabel('')
    # ax2.set_ylabel('Fraction')
    ax2.set_xticklabels(["Female", "Male"], fontsize = 8)
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.set_yticks([0, 100, 200])
    ax2.set_yticklabels(["0", "100", "200"], fontsize = 8)
    
    # Stacked bar chart for initial diagnosis stage
    bottom = np.zeros(len(initial_diagnosis_counts_pivot))
    for status in initial_diagnosis_counts_pivot.columns:
        bars = ax3.bar(initial_diagnosis_counts_pivot.index, initial_diagnosis_counts_pivot[status], bottom=bottom, label=status, color=colors[status], width=0.7)
        for bar, count, color in zip(bars, initial_diagnosis_counts_pivot[status], [colors[status] for _ in range(len(bars))]):
            text_color = 'white' if color == 'black' else 'black'
            # ax3.text(bar.get_x() + bar.get_width() / 2, bar.get_height() -10, f'{count}', ha='center', va='center', color=text_color, fontsize=6)
        bottom += initial_diagnosis_counts_pivot[status]
    
    if show_ax_titles:
        ax3.set_title('Stage at initial diagnosis', loc = "left", fontsize = 8)
    
    ax3.set_xlabel('')
    ax3.spines[["top", "right"]].set_visible(False)
    ax3.set_xticklabels(["Localized", "Metastatic"], fontsize = 8)
    ax3.set_yticks([0, 100, 200])
    ax3.set_yticklabels(["0", "100", "200"], fontsize = 8)
       
    # Stacked bar chart for smoking
    bottom = np.zeros(len(smoking_counts_pivot))
    for status in smoking_counts_pivot.columns:
        bars = ax4.bar(smoking_counts_pivot.index, smoking_counts_pivot[status], bottom=bottom, label=status, color=colors[status], width = 0.7)
        for bar, count, color in zip(bars, smoking_counts_pivot[status], [colors[status] for _ in range(len(bars))]):
            text_color = 'white' if color == 'black' else 'black'
            # ax4.text(bar.get_x() + bar.get_width() / 2, bar.get_height()-10, f'{count}', ha='center', va='center', color=text_color, fontsize=6)
        bottom += smoking_counts_pivot[status]
    
    if show_ax_titles:
        ax4.set_title('Smoking status', loc = "left", fontsize = 8)
    
    ax4.set_xlabel('')
    ax4.spines[["top", "right"]].set_visible(False)
    ax4.set_xticklabels(["Never\nsmoker", "Prev./current\nsmoker"], fontsize = 8)
    ax4.set_yticks([0, 35, 70])
    ax4.set_yticklabels(["0", "35", "70"], fontsize = 8)
    
    # Stacked bar chart for bladder staging
    bottom = np.zeros(len(staging_counts_pivot))
    for status in staging_counts_pivot.columns:
        bars = ax5.bar(staging_counts_pivot.index, staging_counts_pivot[status], bottom=bottom, label=status, color=colors[status], width = 0.7)
        for bar, count, color in zip(bars, staging_counts_pivot[status], [colors[status] for _ in range(len(bars))]):
            text_color = 'white' if color == 'black' else 'black'
            # ax5.text(bar.get_x() + bar.get_width() / 2, bar.get_height()-10, f'{count}', ha='center', va='center', color=text_color, fontsize=6)
        bottom += staging_counts_pivot[status]
    
    if show_ax_titles:
        ax5.set_title('MIBC/NMIBC status', loc = "left", fontsize = 8)
    
    ax5.set_xlabel('')
    ax5.spines[["top", "right"]].set_visible(False)
    ax5.set_xticklabels(["MIBC", "NMIBC"], fontsize = 8)
    ax5.set_yticks([0, 25, 50])
    ax5.set_yticklabels(["0", "25", "50"], fontsize = 8)
        
    if show_legend:
        ax6 = plt.subplot(gs[5])
        handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in colors.values()]
        legend_labels = ["CH+", "CH-"]
        ax6.legend(handles=handles, labels=legend_labels, loc='best', frameon=False)
        ax6.spines[["top", "right", "bottom", "left"]].set_visible(False)
        ax6.set_yticks([])
        ax6.set_xticks([])
    
    if not show_x_ticks:
        ax1.set_xticks([])
        ax1.set_xticklabels([])
        ax2.set_xticks([])
        ax2.set_xticklabels([])
        ax3.set_xticks([])
        ax3.set_xticklabels([])
        ax4.set_xticks([])
        ax4.set_xticklabels([])
        ax5.set_xticks([])
        ax5.set_xticklabels([])
    return(gs)

def CH_gene_and_age_plot_per_category(baseline_chip, age_df, PATH_sample_information, figure_dir, figure_name):
    """
    Scatter plot showing the percentage of the cohort that carries at least one mutation in a given set of genes.
    Plots multiple panels for multiple genes of interest.
    age_df is a df with 2 cols: patient id and age. Bladder and kidney combined.
    """
    fig = plt.figure(figsize=(5, 7))  # Adjusted figure size for better layout
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1], hspace=0.4)
    gene_lists = {
        "DTA genes": ["DNMT3A", "TET2", "ASXL1"], 
        "DDR genes": ["TP53", "PPM1D", "CHEK2", "ATM"], 
        "Splicing genes": ['SF3B1', 'SRSF2', 'U2AF1', 'ZRSR2']}
    
    # Define different color palettes for each group of genes
    color_palettes = {
        "DTA genes": plt.cm.tab10,
        "DDR genes": plt.cm.Accent,
        "Splicing genes": plt.cm.Dark2
    }
    
    for index, (key, value) in enumerate(gene_lists.items()):
        ax = plt.subplot(gs[index])
        color_palette = color_palettes[key]  # Get the corresponding color palette
        max_y_value = 0  # Initialize max y value
        
        for i, gene in enumerate(value):
            df = annotate_mutation_status(baseline_chip, "Both", PATH_sample_information, annotate_what="CHIP", annotate_gene=gene, drop_dependent=True)
            df = df.merge(age_df, how="left").dropna()
            
            # Bin into age groups and return the percentage of the cohort positive for mutation at that bin
            bins = [20, 40, 50, 60, 70, 80, 90]
            labels = ['20-40', '40-50', '50-60', '60-70', '70-80', '80-90']
            df['Age_bin'] = pd.cut(df['Age'], bins=bins, labels=labels, right=False)
            age_groups = df.groupby('Age_bin')['CHIP status'].value_counts(normalize=True).unstack().fillna(0)
            age_groups['Positive_percentage'] = age_groups['Positive'] * 100
            
            # Find the highest ymax value
            max_y_value = max(max_y_value, age_groups['Positive_percentage'].max())
            
            
            # Plotting with the specific color palette for the group
            ax.plot(age_groups.index, age_groups['Positive_percentage'], marker='o', color=color_palette(i / len(value)), label=gene)
        
        # Round up the highest ymax to the nearest number ending with 0
        ymax_rounded = np.ceil(max_y_value / 10) * 10
        ax.set_ylim(-2, ymax_rounded+2)  # Set the y-axis limit to the rounded ymax
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))  # Adjust y-axis ticks
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(key)
        ax.set_xlabel("Age Group")
        ax.set_ylabel("CH+ (%)")
        ax.legend(title="Genes", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
        
    gs.tight_layout(fig)
    fig.savefig(os.path.join(figure_dir, figure_name + ".pdf"))

def make_ppm1d_plot():
    """
    Makes a plot showing the number of PPM1D mutations present in patients with and without chemo exposure.
    """
    chemo_yes_mutation_yes = 13
    chemo_yes_mutation_no = 41
    chemo_no_mutation_yes = 5
    chemo_no_mutation_no = 49
    
    # Stacked bar data
    categories = ['Chemo\nExposed', 'Not Chemo\nExposed']
    mutation_yes = [chemo_yes_mutation_yes, chemo_no_mutation_yes]
    mutation_no = [chemo_yes_mutation_no, chemo_no_mutation_no]
    
    # Bar chart
    fig, ax = plt.subplots(figsize = (3, 3))
    bar_width = 0.5
    x = np.arange(len(categories))
    
    ax.bar(x, mutation_no, bar_width, label='Mutation No', color='skyblue')
    ax.bar(x, mutation_yes, bar_width, bottom=mutation_no, label='Mutation Yes', color='salmon')
    
    # Labels and title
    ax.set_xlabel('Chemo Exposure')
    ax.set_ylabel('Number of Patients')
    ax.set_title('PPM1D Mutations in Patients with and without Chemo Exposure')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend(frameon = False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig("/groups/wyattgrp/users/amunzur/pipeline/results/figures/amazing_figures/bladder/ppm1d_chemo_plot.pdf")

def make_tmb_plot(df_ch, df_ctdna, ax):
    """
    Scatter plot. 
    """
    ch_count = df_ch[~(df_ch["Gene"].isin(["DNMT3A", "TET2"]) & (df_ch["VAF_t"] >= 0.005))]["Patient_id"].value_counts().reset_index().rename(columns = {"index": "Patient_id", "Patient_id":"CH_count"})
    ctdna_count = df_ctdna[~df_ctdna["Gene"].isin(["DNMT3A", "TET2"])]["Patient_id"].value_counts().reset_index().rename(columns = {"index": "Patient_id", "Patient_id":"ctDNA_count"})
    merged = ch_count.merge(ctdna_count, how = "outer").fillna(0)
    merged["Summed count"] = merged["CH_count"] + merged["ctDNA_count"]
    merged = merged.sort_values(by = "Summed count", ascending = False)
    
    ax.scatter(merged["Summed count"] + np.random.normal(0, 0.15, size=len(merged)), merged["ctDNA_count"] + np.random.normal(0, 0.15, size=len(merged)), s = 4, edgecolor = None, color = "black", alpha = 0.7, zorder = 100)
    # ax.axhline(10, linestyle='--', color='red')       
    # ax.axvline(10, linestyle='--', color='red')
    
    # Plot the shaded region for x > 10 and y > 10
    ax.set_xlim((-2, 25))
    ax.set_ylim((-2, 25))
    ax.set_xticks([0, 5, 10, 15, 20, 25])
    ax.set_yticks([0, 5, 10, 15, 20, 25])
    ax.set_xticklabels(["0", "5", "10", "15", "20", "25"])
    ax.set_yticklabels(["0", "5", "10", "15", "20", "25"])
    
    ax.axvspan(10, max(ax.get_ylim()), ymin=0, ymax=1, alpha=0.2, color='red')
    ax.axhspan(10, max(ax.get_xlim()), xmin=0, xmax=1, alpha=0.2, color='blue')
    
    ax.set_ylabel("With WBC control")      
    ax.set_xlabel("Without WBC control")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # label the coordinates
    ax.text(11, 23.5, 'I', fontsize=9, ha='center', va='center')
    ax.text(-1, 23.5, 'II', fontsize=9, ha='center', va='center')
    ax.text(-1, 8.5, 'III', fontsize=9, ha='center', va='center')
    ax.text(11, 8.5, 'IV', fontsize=9, ha='center', va='center')
    
    return(ax)

def make_gene_cooccurrence_plot(mutations_df, genes_list, ax):
    """
    Given a mutations_df and a list of genes, make a gene co-occurrence plot.
    """
    # Filter the dataframe to include only the specified genes
    mutations_df = mutations_df[mutations_df["Gene"].isin(genes_list)][["Patient_id", "Gene"]]
    # Create a binary matrix of patients (rows) and genes (columns)
    binary_matrix = mutations_df.pivot_table(index='Patient_id', columns='Gene', aggfunc='size', fill_value=0)
    binary_matrix = (binary_matrix > 0).astype(int)
    # Initialize a DataFrame to hold the p-values from Fisher's Exact Test
    p_values = pd.DataFrame(index=genes_list, columns=genes_list, dtype=float)
    # Perform pairwise Fisher's Exact Test
    for gene1 in genes_list:
        for gene2 in genes_list:
            if gene1 != gene2:
                # Create contingency table
                contingency_table = pd.crosstab(binary_matrix[gene1], binary_matrix[gene2])
                # Perform Fisher's Exact Test
                _, p_value = fisher_exact(contingency_table)
                # Store the p-value
                p_values.loc[gene1, gene2] = p_value
            else:
                p_values.loc[gene1, gene2] = np.nan  # Avoid self-comparison
    # Flatten the lower triangle of the p-values DataFrame to apply multiple test correction
    lower_triangle_p_values = p_values.where(np.tril(np.ones(p_values.shape), -1).astype(bool)).stack().values
    # Perform multiple test correction (Benjamini-Hochberg)
    _, corrected_p_values, _, _ = multipletests(lower_triangle_p_values, method='fdr_bh')
    
    # Assign corrected p-values back to the DataFrame
    corrected_p_values_df = p_values.copy()
    corrected_p_values_df.values[np.tril_indices_from(corrected_p_values_df, -1)] = corrected_p_values
    # Define significance threshold
    significance_threshold = 0.05
    # Apply mask for the upper triangle
    mask = np.triu(np.ones_like(corrected_p_values_df, dtype=bool))
    # Create a mask for significant p-values
    significant_mask = corrected_p_values_df < significance_threshold
    # Plotting the heatmap with significant p-values only
    sns.heatmap(corrected_p_values_df, ax=ax, annot=True, mask=mask | ~significant_mask, cmap='viridis', cbar=True, square=True, linewidths=0.5)
    ax.set_title('Gene Co-occurrence Significant P-values (Fisher Exact Test)')
    return ax

# Example usage:
# fig, ax = plt.subplots(figsize=(10, 8))
# make_gene_cooccurrence_plot(baseline, ["DNMT3A", "TET2", "PPM1D", "ASXL1", "ATM"], ax)
# fig.savefig(os.path.join(figure_dir, "test.png"))

def age_and_n_muts(baseline, PATH_kidney_clinical, PATH_clinical_bladder, ax):
    """
    Correlation between age and the number of mutations per patient.
    """
    # clean up clinical data
    kidney_clin = pd.read_csv(path_kidney_clin, sep = "\t")[["GUBB ID", "Age at GUBB draw"]].rename(columns= {"GUBB ID": "Patient_id", "Age at GUBB draw": "Age"})
    kidney_clin["Diagnosis"] = "mRCC"
    kidney_clin["Color"] = "orangered"
    
    bladder_clin = pd.read_csv(path_bladder_clin)
    bladder_clin = bladder_clin[bladder_clin["First sample?"] == True]
    bladder_clin = bladder_clin[["Patient_id", "Age at blood draw"]].rename(columns = {"Age at blood draw": "Age"})
    bladder_clin["Diagnosis"] = "mUC"
    bladder_clin["Color"] = "deepskyblue"
    
    # Prepare age data
    age_df = pd.concat([kidney_clin, bladder_clin]).reset_index(drop = True)
    bins = [20, 40, 50, 60, 70, 80, 90]
    labels = [f'{bins[i]}-{bins[i+1]}' for i in range(len(bins) - 1)]
    age_df['Age_Bin'] = pd.cut(age_df['Age'], bins=bins, labels=labels, right=False)
    age_df['Age_Bin'] = pd.Categorical(age_df['Age_Bin'], categories=labels, ordered=True)
    
    n_muts = baseline["Patient_id"].value_counts().reset_index().rename(columns = {"index": "Patient_id", "Patient_id": "n_muts"})
    merged = n_muts.merge(age_df, how = "inner")
    merged = merged.sort_values(by='Age_Bin') # Sort by Age_Bin to ensure proper order
    merged['Age_Bin_Num'] = merged['Age_Bin'].cat.codes
    
    # Fit and plot second-degree polynomial for mRCC (kidney)
    kidney_data = merged[merged["Diagnosis"] == "mRCC"]
    if not kidney_data.empty:
        kidney_fit = np.polyfit(kidney_data["Age_Bin_Num"], kidney_data["n_muts"], 2)
        kidney_poly = np.poly1d(kidney_fit)
        ax.plot(kidney_data["Age_Bin"], kidney_poly(kidney_data["Age_Bin_Num"]), color="orangered", label="Kidney (mRCC)")
    
    # Fit and plot second-degree polynomial for mUC (bladder)
    bladder_data = merged[merged["Diagnosis"] == "mUC"]
    if not bladder_data.empty:
        bladder_fit = np.polyfit(bladder_data["Age_Bin_Num"], bladder_data["n_muts"], 2)
        bladder_poly = np.poly1d(bladder_fit)
        ax.plot(bladder_data["Age_Bin"], bladder_poly(bladder_data["Age_Bin_Num"]), color="deepskyblue", label="Bladder (mUC)")
    
    ax.scatter(merged["Age_Bin_Num"].astype(int) + np.random.uniform(-0.25, 0.25, size=len(merged)), merged["n_muts"], color=merged["Color"], s=2)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
    ax.set_yticks([0, 4, 8, 12, 16, 20])
    ax.set_yticklabels(["0", "4", "8", "12", "16", "20"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Number of mutations\nper patient")
    return(ax)

def plot_ctDNA_and_CH_VAF_box(ch_df, ctdna_df, PATH_sample_information, ax, ax_title, hide_xticks = False, show_legend = True, annotate_median = True, run_mwu = True, all_ctdna_muts_df = None, scatter_zorder = 0, scatter_alpha = 0.5, scatter_size = 2):
    """
    Boxp. Shows overlap of ctDNA and CH vafs.
    """
    import math
    x_ax_pos = {"<10%": 1, '≥10%': 2}
    
    # Log the VAFs
    ch_df["VAF_t_log"] = np.log10(ch_df["VAF_t"].replace(0, np.nan))
    ctdna_df["VAF_t_log"] = np.log10(ctdna_df["VAF_t"].replace(0, np.nan))
    
    # group patients based on the vaf of their ctDNA mutation
    ct_fr = ctdna_df.groupby("Patient_id")["VAF_t"].max().reset_index().rename(columns = {"VAF_t": "Max VAF"})
    ct_fr["ct bin"] = ["≥10%" if vaf >= 10 else "<10%" for vaf in ct_fr["Max VAF"]]
    
    # Add the ctDNA negative patients
    ctdna_status = annotate_mutation_status(ctdna_df, "Both", PATH_sample_information, annotate_what = "ctDNA")
    ctdna_status = ctdna_status[ctdna_status["Timepoint"] == "Baseline"]
    ctdna_neg = ctdna_status[ctdna_status["ctDNA status"] == "Negative"]
    ctdna_neg["Max VAF"] = 0
    ctdna_neg["ct bin"] = "<10%"
    ct_fr = pd.concat([ct_fr, ctdna_neg[["Patient_id", "Max VAF", "ct bin"]]])
    
    # Counts, to visualize in the plots
    bin_high_count = ct_fr[ct_fr["ct bin"] == "≥10%"].shape[0]
    bin_low_count = ct_fr[ct_fr["ct bin"] == "<10%"].shape[0]
    
    for pt, group in ct_fr.groupby("Patient_id"):
        
        # Get relevant CH and ctDNA mutations
        pt_ch_mutations = ch_df[ch_df["Patient_id"] == pt]
        pt_ctdna_mutations = ctdna_df[ctdna_df["Patient_id"] == pt]
        
        jitter = np.random.normal(-0.05, 0.05, 1)
        # Plot the CH mutations
        for i, row in pt_ch_mutations.iterrows():
            x_pos = x_ax_pos[group["ct bin"].unique()[0]] # position along the x axis
            ax.scatter(x_pos+0.2+jitter, row["VAF_t_log"], color="red", s = scatter_size, zorder = scatter_zorder, alpha = scatter_alpha)
        
        # Plot the ctDNA mutations
        for i, row in pt_ctdna_mutations.iterrows():
            jitter = np.random.normal(-0.05, 0.05, 1)
            x_pos = x_ax_pos[group["ct bin"].unique()[0]] # position along the x axis
            ax.scatter(x_pos-0.2-jitter, row["VAF_t_log"], color="forestgreen", s = scatter_size, zorder = scatter_zorder, alpha = scatter_alpha)
    
    # Now add the boxplots
    ch_df["Date_collected"] = pd.to_datetime(ch_df["Date_collected"], format = '%Y%b%d')
    ctdna_df["Date_collected"] = pd.to_datetime(ctdna_df["Date_collected"], format = '%Y%b%d')
    ch_df = ch_df.merge(ct_fr, how = "left")
    ctdna_df = ctdna_df.merge(ct_fr, how = "left")
    
    # Plotting    
    for mut_df, x_pos_offset in zip([ch_df, ctdna_df], [0.15, -0.15]):
        for i, group in mut_df.groupby("ct bin"):
            if not group.empty:
                x_pos = x_ax_pos[group["ct bin"].unique()[0]] + x_pos_offset
                boxplot = ax.boxplot(group["VAF_t_log"], 
                positions = [x_pos], widths = 0.3, zorder = 1, 
                showfliers = False, patch_artist = True, 
                boxprops=dict(facecolor='none'), 
                medianprops = dict(color='black'), 
                capprops = dict(color='black'), 
                whiskerprops = dict(color='black'))
    if annotate_median:
        for mut_df, x_pos_offset_median, ha in zip([ch_df, ctdna_df], [0.32, -0.32], ["left", "right"]):
            for i, group in mut_df.groupby("ct bin"):
                median_value = np.median(group["VAF_t"])
                median_position = np.median(group["VAF_t_log"])
                median_xpos = x_ax_pos[group["ct bin"].unique()[0]] + x_pos_offset_median
                ax.text(median_xpos, median_position, f'{median_value:.2f}', ha=ha, va='center', fontsize=5, color='black', zorder=2)
    
    # Running the significance test
    if run_mwu:
        for [ct_bin, x_pos_p_value] in zip(["<10%", "≥10%"], [1, 2]): 
            ch_values = ch_df[ch_df["ct bin"] == ct_bin]["VAF_t"]
            ctdna_values = ctdna_df[ctdna_df["ct bin"] == ct_bin]["VAF_t"]
            if not ch_values.empty and not ctdna_values.empty: 
                stat, p = mannwhitneyu(ch_values, ctdna_values)
            
                # format the p values for plotting
                p_base, p_exp = np.format_float_scientific(p, precision=1).split('e')
                p_exp = int(p_exp)  # Convert exponent to an integer
                ax.text(x_pos_p_value, np.log10(100), f'p={p_base}×10^{{{p_exp}}}', ha='center', va='bottom', fontsize=5, color='black', zorder=2)
    
    # Now to annotate the number of dots in each boxplot
    for [ct_bin, x_pos_value] in zip(["<10%", "≥10%"], [1, 2]): 
        n_ch = ch_df[ch_df["ct bin"] == ct_bin].shape[0]
        n_ctdna = ctdna_df[ctdna_df["ct bin"] == ct_bin].shape[0]
        ax.text(x_pos_value+0.2, np.log10(0.12), f'n={n_ch}', ha='center', va='bottom', fontsize=5, color='black', zorder=2)
        ax.text(x_pos_value-0.2, np.log10(0.12), f'n={n_ctdna}', ha='center', va='bottom', fontsize=5, color='black', zorder=2)
    
    # AES
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["<10%", "≥10%"])
    ax.set_xlabel("ctDNA bin")
    ax.set_ylabel("cfDNA VAF%")
    ax.set_ylim((np.log10(0.10), np.log10(90)))
    ax.set_yticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(50), np.log10(90)])
    ax.set_yticklabels([".25", "1", "2", "10", "50", "90"])
    ax.set_title(ax_title, fontstyle='italic')
    
    if hide_xticks:
        ax.tick_params(axis='x', bottom=False, labelbottom=False)
    
    # add legend
    if show_legend:
        legend_colors = ["red", "forestgreen"]
        legend_labels = ["CH", "ctDNA"]
        legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=3, linestyle='') for color, label in zip(legend_colors, legend_labels)]
        ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize = 8, handletextpad=0.3)
    
    return ax

def plot_ctDNA_and_CH_VAF_split_violin(PATH_sample_information, ch_df, ctdna_df, ax, ax_title, hide_xticks = False, show_legend = True):
    """
    Boxp. Shows overlap of ctDNA and CH vafs.
    """
    import math
    
    # group patients based on the vaf of their ctDNA mutation
    ct_fr = ctdna_df.groupby("Patient_id")["VAF_t"].max().reset_index().rename(columns = {"VAF_t": "Max VAF"})
    ct_fr["ct bin"] = ["≥10%" if vaf >= 10 else "<10%" for vaf in ct_fr["Max VAF"]]
    
    # Add the ctDNA negative patients
    ctdna_status = annotate_mutation_status(ctdna_df, "Both", PATH_sample_information, annotate_what = "ctDNA")
    ctdna_status = ctdna_status[ctdna_status["Timepoint"] == "Baseline"]
    ctdna_neg = ctdna_status[ctdna_status["ctDNA status"] == "Negative"]
    ctdna_neg["Max VAF"] = 0
    ctdna_neg["ct bin"] = 0
    ct_fr = pd.concat([ct_fr, ctdna_neg[["Patient_id", "Max VAF", "ct bin"]]])
    
    # Counts, to visualize in the plots
    bin_high_count = ct_fr[ct_fr["ct bin"] == "≥10%"].shape[0]
    bin_low_count = ct_fr[ct_fr["ct bin"] == "<10%"].shape[0]
    bin_zero_count = ct_fr[ct_fr["ct bin"] == 0].shape[0]
    
    # Add ctDNA fraction to the CH df
    ch_df["Mut_type"] = "CH"
    ch_df["Date_collected"] = pd.to_datetime(ch_df["Date_collected"], format = '%Y%b%d')
    ch_df = ch_df[["Patient_id", "Date_collected", "VAF_t", "Mut_type"]].merge(ct_fr, how = "left")
    
    # Add ctDNA fraction to the ctDNA df
    ctdna_df["Mut_type"] = "ctDNA"
    ctdna_df["Date_collected"] = pd.to_datetime(ctdna_df["Date_collected"], format = '%Y%b%d')
    ctdna_df = ctdna_df[["Patient_id", "Date_collected", "VAF_t", "Mut_type"]].merge(ct_fr, how = "left")
    
    combined = pd.concat([ch_df, ctdna_df]).reset_index(drop = True)
    
    sns.violinplot(data=combined, x="ct bin", y="VAF_t", hue="Mut_type", split=True, order = [0, "<10%", "≥10%"],
                   inner=None, fill=False, palette={"CH": "red", "ctDNA": "forestgreen"}, 
                   ax = ax, legend = False, gap=1, cut=0)
    
    # The line around the violins will be removed:
    from matplotlib.collections import PolyCollection
    for art in ax.get_children():
         if isinstance(art, PolyCollection):
             art.set_edgecolor((0.3, 0.3, 0.3, 0))
    
    # AES
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("Highest ctDNA mutation VAF%")
    ax.set_ylabel("cfDNA VAF%")
    ax.set_title(ax_title)
    
    # add legend
    legend_colors = ["red", "forestgreen"]
    legend_labels = ["CH", "ctDNA"]
    legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=3, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize = 8, handletextpad=0.3)
    
    # Set x tick labels
    x_tick_labels = [f"0\nn={bin_zero_count}", f"<10%\nn={bin_low_count}", f"≥10%\nn={bin_high_count}"]
    ax.set_xticklabels(x_tick_labels)
    # ax.text(0, 100, bin_0_29_count, ha='center', va='bottom', color='black', fontsize = 5)
    # ax.text(1, 100, bin_30_100_count, ha='center', va='bottom', color='black', fontsize = 5)
    return(ax)

def clean_up_df(muts_df):
    """
    Subsetting the df to the relevant cols before analysis.
    """
    muts_df = muts_df[(muts_df["Dependent"] == False) & (muts_df["Timepoint"] == "Baseline")].reset_index(drop = True)
    
    if "Status" in muts_df.columns:
        variant_type = "ctDNA"
    else:
        variant_type = "CH"
    
    muts_df_subsetted = muts_df[["Patient_id", "Gene", "VAF_t", "Consequence"]]
    muts_df_subsetted["Variant_type"] = variant_type
    
    return(muts_df_subsetted)

def calculate_retention_fraction(muts_df_filtered, muts_df_total, gene, mut_type):
    """
    Given a gene name calculates the fraction of ctDNA mutations we correctly identified.
    Provide CH or ctDNA for mut_type.
    """
    muts_df_filtered_subset = muts_df_filtered[muts_df_filtered["Variant_type"] == mut_type]
    muts_df_total_subset = muts_df_total[muts_df_total["Variant_type"] == mut_type]
    
    n_total_gene = muts_df_total_subset[muts_df_total_subset["Gene"] == gene].shape[0]
    n_identified = muts_df_filtered_subset[muts_df_filtered_subset["Gene"] == gene].shape[0]
    
    if n_total_gene == 0:
        identification_fraction = np.nan
    else:
        identification_fraction = round(n_identified/n_total_gene, 3)
    
    result_dict = {"Gene": gene, "n total": n_total_gene, "n retained": n_identified, "retention fraction": identification_fraction}
    
    return(result_dict)

def make_variant_retention_barchart(results, ax, color_dict):
    """
    For each gene shows the fraction and number of variants retained and missed, for both CH and ctDNA variants
    """
    # Plot ctDNA stats
    bar_height = 0.8
    ax.bar(results.index, results["n retained_ctDNA"], width=bar_height, color=color_dict["True positive"], alpha=0.8)
    ax.bar(results.index, results["n missed_ctDNA"], bottom = results["n retained_ctDNA"], width=bar_height, color=color_dict["False negative"], alpha=0.3)
    
    # Plot CH stats
    ax.bar(results.index, results["n retained_CH"], bottom = 80, width=bar_height, color=color_dict["False positive"], alpha=0.8)
    ax.bar(results.index, results["n missed_ctDNA"], bottom = 80+results["n retained_CH"], width=bar_height, color=color_dict["True negative"], alpha=0.3)
    
    # Aes
    ax.axhline(80, color='black', linestyle='--', linewidth = 0.25)
    ax.set_xticks(results.index)
    ax.set_xticklabels(results["Gene"], fontsize = 8, rotation = 90)
    ax.tick_params(axis='x', bottom=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim((results.index.min() - bar_height*0.75, results.index.max() + bar_height*0.75))    
    ax.set_yticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110])
    ax.set_yticklabels(["0", "10", "20", "30", "40", "50", "60", "70", "0", "10", "20", "30"])
    ax.text(-0.5, 110, "CH variants", ha='left', va='bottom', color='black', fontsize = 10)
    ax.text(-0.5, 70, "ctDNA variants", ha='left', va='bottom', color='black', fontsize = 10)
    
    # add legend
    legend_handles = [plt.Line2D([0], [0], marker='s', color=color_dict[label], label=label, markersize=4, linestyle='') for label in color_dict.keys()]
    ax.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize=8, handletextpad=0.1, bbox_to_anchor=(1, 0.3))
    return(ax)

def filter_like_industry(all_vars_chip, all_vars_somatic):
    """
    Filters like indutry, reports findings as a df.
    """
    ctdna_vars = clean_up_df(all_vars_somatic)
    chip_vars = clean_up_df(all_vars_chip)
    combined = pd.concat([ctdna_vars, chip_vars]).reset_index(drop = True)
    
    # Now filtering based on a few criteria.
    # 1. Excluding DTA genes
    combined_filtered = combined[~combined["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])]
    
    # 2. Exclude all variants less than 1% in VAF
    combined_filtered = combined_filtered[combined_filtered["VAF_t"] > 1]
    
    # 3. Excluding variants that are <25% of the sample's maximum VAF (i.e. a surrogate for tumor fraction), i.e. highly "subclonal" variants.
    # Compute sample maximum VAF
    subclonal_df = combined.groupby("Patient_id")["VAF_t"].max().reset_index()
    subclonal_df["Subclonal VAF threshold"] = subclonal_df["VAF_t"]*0.25
    del subclonal_df["VAF_t"]
    combined_filtered = combined_filtered.merge(subclonal_df, how = "left")
    combined_filtered = combined_filtered[combined_filtered["VAF_t"] > combined_filtered["Subclonal VAF threshold"]]    
    
    # Combining results
    ctDNA_results_list = []
    CH_results_list = []
    for gene in ctdna_vars["Gene"].unique(): 
        ctDNA_retention = calculate_retention_fraction(combined_filtered, combined, gene = gene, mut_type = "ctDNA")
        CH_retention = calculate_retention_fraction(combined_filtered, combined, gene = gene, mut_type = "CH")
        
        # Append to list
        ctDNA_results_list.append(ctDNA_retention)
        CH_results_list.append(CH_retention)    
    
    # Make a ctDNA df
    ctDNA_df = pd.DataFrame(ctDNA_results_list)
    ctDNA_df = ctDNA_df.sort_values(by = "n total", ascending = False).reset_index(drop = True)
    ctDNA_df["missed fraction"] = 1 - ctDNA_df["retention fraction"]
    ctDNA_df["n missed"] = ctDNA_df["n total"] - ctDNA_df["n retained"] 
    
    # Make a CH df
    CH_df = pd.DataFrame(CH_results_list)
    CH_df = CH_df.sort_values(by = "n total").reset_index(drop = True)
    CH_df["missed fraction"] = 1 - CH_df["retention fraction"]
    CH_df["n missed"] = CH_df["n total"] - CH_df["n retained"]  
    
    filter_like_industry_results = ctDNA_df.merge(CH_df, on = "Gene", suffixes = ["_ctDNA", "_CH"])
    filter_like_industry_results = filter_like_industry_results[~filter_like_industry_results["Gene"].isin(["DNMT3A", "TET2", "ASXL1"])]
    
    return(filter_like_industry_results)

def line_plot_with_gene_categories(muts_df, ax, PATH_sample_information):
    """
    CH prevalence in the cohort in certain genes with certain CH positivity thresholds.
    """
    sample_info = pd.read_csv(PATH_sample_information, sep = "\t", names = ["Patient_id", "Date_collected", "Diagnosis", "Timepoint"])
    sample_info = sample_info[(sample_info["Timepoint"] == "Baseline")].reset_index(drop = True)
    n_pts_total = sample_info.shape[0]
    
    thresholds = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
    gene_groups = [muts_df["Gene"].unique(), ["DNMT3A", "TET2", "ASXL1"], ["TP53", "ATM", "CHEK2", "BRCA1", "BRCA2", "ARID1A"]]
    colors = ["black", "hotpink", "limegreen"]
    
    for gene_group, color in zip(gene_groups, colors):
        values = []
        for threshold in thresholds:
                muts = muts_df[(muts_df["VAF_n"] >=threshold) & (muts_df["Gene"].isin(gene_group))]
                n_pts = muts["Patient_id"].unique().size
                fr = n_pts/n_pts_total
                values.append(fr)
        # Plotting
        ax.plot(thresholds, values, marker='o', color=color, markersize = 4)
    
    # Add legend
    legend_colors = colors
    legend_labels = ["All genes", "DTA genes", "DDR genes"]
    legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize = 8, handletextpad=0.3, title = "")
    
    # Aes
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_yticklabels(["0", "20", "40", "60", "80"])
    ax.set_ylabel("% of patients with CH")
    ax.set_xlabel("Minimum CH VAF%")
    return(ax)

def CH_plot_VAFs_dot(df_CH, age_df, ax):
    """
    Divides the CH into 3 categories based on VAF, and plots the VAFs of mutations for each category with a different color.
    """    
    # Generate dfs of varying thresholds and annotate their CHIP mutation status   
    colors = {"0.25%-2%": "limegreen", "2%-10%": "#709bd0", ">10%": "#224193"}
    
    bins = [20, 40, 50, 60, 70, 80, 90]
    labels = ['20-39', '40-49', '50-59', '60-69', '70-79', '80-89']
    label_pos = {label: i for i, label in enumerate(labels)}
    
    # Add the age bin info
    df = df_CH[["Patient_id", "VAF_n"]].merge(age_df)
    df['Age_bin'] = pd.cut(df['Age'], bins=bins, labels=labels, right=False)
    df["label x pos"] = df["Age_bin"].map(label_pos)
    df["VAF_n_log"] = np.log10(df["VAF_n"].replace(0, np.nan))  # Replace 0 with NaN to avoid log(0)
    
    for i, row in df.iterrows():
        if row["VAF_n"] >= 0.25 and row["VAF_n"] < 2:
            color = "limegreen"
        elif row["VAF_n"] >= 2 and row["VAF_n"] < 10:
            color = "#709bd0"
        elif row["VAF_n"] >= 10: 
            color = "#224193"
        # Plotting
        jitter = np.random.uniform(-0.2, 0.2, 1)[0]
        ax.scatter(row["label x pos"]+jitter, row["VAF_n_log"], marker='o', color=color, s = 2)
    
    # AES
    ax.set_xlim((np.log10(0.25), np.log10(50)))
    ax.set_xticks([np.log10(0.25), np.log10(1), np.log10(2), np.log10(10), np.log10(50)])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel("WBC VAF%")
    return(ax)

def make_irAE_association_plots(all_vars_chip, PATH_kidney_clinical, PATH_sample_information, ax1, ax2, ax3, irae_legend_ax, gene = None):
    """
    Checks for association between irAEs and chip status
    """
    clin_df = pd.read_csv(PATH_kidney_clinical)[["Patient_id", "irAE", "irAE AI", "irAE colitis"]]
    clin_df = clin_df[~clin_df["irAE"].isna()].reset_index(drop = True) # drop patients that didn't receive systemic treatment
    clin_df.columns = ["Patient_id", "Any irAE", "Adrenal insufficiency", "Colitis"]
    
    # Annotate CH status
    ch1 = annotate_mutation_status(all_vars_chip, "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = gene).merge(clin_df, how = "inner").dropna()
    ch2 = annotate_mutation_status(all_vars_chip[all_vars_chip["VAF_n"] >= 2], "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = gene).merge(clin_df, how = "inner").dropna()
    ch3 = annotate_mutation_status(all_vars_chip[all_vars_chip["VAF_n"] >= 10], "Kidney", PATH_sample_information, annotate_what = "CHIP", annotate_gene = gene).merge(clin_df, how = "inner").dropna()
    
    # Subset to baseline
    ch1 = ch1[ch1["Timepoint"] == "Baseline"].reset_index(drop = True)
    ch2 = ch2[ch2["Timepoint"] == "Baseline"].reset_index(drop = True)
    ch3 = ch3[ch3["Timepoint"] == "Baseline"].reset_index(drop = True)
    
    # Plotting
    n_total = clin_df.shape[0] # will be the denom
    for irAE_name, ax in zip(["Any irAE", "Adrenal insufficiency", "Colitis"], [ax1, ax2, ax3]):
        for i, (label, df) in enumerate(zip(["Any CH", "CH≥2%"], [ch1, ch2])):
            n_irae_neg_ch_pos = df[(df["CHIP status"] == "Positive") & (df[irAE_name] == 0)].shape[0]
            n_irae_neg_ch_neg = df[(df["CHIP status"] == "Negative") & (df[irAE_name] == 0)].shape[0]
            n_irae_pos_ch_pos = df[(df["CHIP status"] == "Positive") & (df[irAE_name] == 1)].shape[0]
            n_irae_pos_ch_neg = df[(df["CHIP status"] == "Negative") & (df[irAE_name] == 1)].shape[0]
            
            perc_irae_neg_ch_pos = n_irae_neg_ch_pos/n_total*100
            perc_irae_neg_ch_neg = n_irae_neg_ch_neg/n_total*100
            perc_irae_pos_ch_pos = n_irae_pos_ch_pos/n_total*100
            perc_irae_pos_ch_neg = n_irae_pos_ch_neg/n_total*100
            
            # Bars for irAE absent
            ax.bar(i-0.1, perc_irae_neg_ch_pos, color="black", width = 0.2, edgecolor = "None", zorder = 10)
            ax.bar(i-0.1, perc_irae_neg_ch_neg, bottom = perc_irae_neg_ch_pos, color="gray", width = 0.2, edgecolor = "None", zorder = 10)
            
            # Bars for irAE present
            ax.bar(i+0.1, perc_irae_pos_ch_pos, color="black", width = 0.2, edgecolor = "None", zorder = 10)
            ax.bar(i+0.1, perc_irae_pos_ch_neg, bottom = perc_irae_pos_ch_pos, color="gray", width = 0.2, edgecolor = "None", zorder = 10)
            
            # Run fishers and annotate on the plot
            import scipy.stats as stats
            oddsratio, p_value = stats.fisher_exact([[n_irae_neg_ch_pos, n_irae_pos_ch_pos], [n_irae_neg_ch_neg, n_irae_pos_ch_neg]])
            p_val = round(p_value, 3)
            ax.text(i, 100, f"P={p_val}", ha='center', va='center', color='black', fontsize = 8, zorder = 15)
            
            # Annotate numbers on plot
            ax.text(i-0.1, perc_irae_neg_ch_pos-1, n_irae_neg_ch_pos, ha='center', va='top', color='white', fontsize = 7, zorder = 15)
            ax.text(i-0.1, perc_irae_neg_ch_pos+perc_irae_neg_ch_neg-1, n_irae_neg_ch_neg, ha='center', va='top', color='black', fontsize = 7, zorder = 15)
            ax.text(i+0.1, perc_irae_pos_ch_pos-1, n_irae_pos_ch_pos, ha='center', va='top', color='white', fontsize = 7, zorder = 15)
            ax.text(i+0.1, perc_irae_pos_ch_pos+perc_irae_pos_ch_neg-1, n_irae_pos_ch_neg, ha='center', va='top', color='black', fontsize = 7, zorder = 15)
            
            # Annotate CH status on the plot
            ax.bar(i, 115, width = 0.8, color = "whitesmoke", edgecolor = "None")
            ax.text(i, 110, label, ha='center', va='center', color='black', fontsize = 8, zorder = 15)
        
        # AES
        ax.set_title(irAE_name)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xticks([-0.1, 0.1, 0.9, 1.1])
        ax.set_xticklabels(["Absent", "Present", "Absent", "Present"], rotation = 90)
        ax.tick_params(axis='x', bottom=False, pad = 0)        
        ax.set_xlabel("irAE status")
        ax.set_ylim((0, 120))
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_yticklabels(["0", "25", "50", "75", "100"])
        ax.set_xlim((-0.7, 1.7))
        if irAE_name == "Any irAE":
            ax.set_ylabel("% of patients")
        else:
            ax.set_ylabel("")
        
        # Add legend
        legend_colors = ["black", "gray"]
        legend_labels = ["CH+", "CH-"]
        legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=5, linestyle='') for color, label in zip(legend_colors, legend_labels)]
        irae_legend_ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize = 8, handletextpad=0.3, title = "")
        irae_legend_ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
        irae_legend_ax.set_yticks([])
        irae_legend_ax.set_yticklabels([])
        irae_legend_ax.set_xticks([])
        irae_legend_ax.set_xticklabels([])
        
    return(ax1, ax2, ax3)

def run_correlation_age_and_nCH(age_df, all_vars_chip, ax):
    """
    Runs Spearman's correlation between age and the number of CH mutations
    """
    df = all_vars_chip["Patient_id"].value_counts().reset_index().rename(columns = {"index": "Patient_id", "Patient_id": "n_muts"})
    df = df.merge(age_df, how = "inner")
    correlation, p_value = spearmanr(df['Age'], df['n_muts'])
    
    # Plot the scatter
    ax.scatter(df["Age"], df["n_muts"], s = 2)
    
    return(ax)


