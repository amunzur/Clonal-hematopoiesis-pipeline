import numpy as np
import pandas as pd
from scipy import stats

def run_bootstrap(distribution_1, distribution_2, n_bootstraps):
	"""
	Calculates confidence intervals using bootsrapping across two distributions.
	"""
	# Function to calculate the difference in means
	def bootstrap_mean_difference(data1, data2):
		return np.mean(np.random.choice(data1, size=len(data1), replace=True)) - \
		np.mean(np.random.choice(data2, size=len(data2), replace=True))
	
	# Generate bootstrap samples
	bootstrap_differences = [bootstrap_mean_difference(distribution_1, distribution_2) for _ in range(n_bootstraps)]
	observed_difference = np.mean(distribution_1) - np.mean(distribution_2)
	
	# Calculate the 95% confidence interval
	confidence_interval = np.percentile(bootstrap_differences, [2.5, 97.5])
	p_value = np.sum(np.abs(bootstrap_differences) >= np.abs(observed_difference)) / n_bootstraps
	
	return(bootstrap_differences, observed_difference, confidence_interval, p_value)

CH_mutant_fragments = np.array(result_dict["CH_mutant_fragments"])
CH_WT_fragments = np.array(result_dict["CH_WT_fragments"])
ctDNA_mutant_fragments = np.array(result_dict["ctDNA_mutant_fragments"])
ctDNA_WT_fragments = np.array(result_dict["ctDNA_WT_fragments"])

# Define the range
lower_bound = 100
upper_bound = 500

# Subset the arrays
CH_mutant_fragments_subset = CH_mutant_fragments[(CH_mutant_fragments >= lower_bound) & (CH_mutant_fragments <= upper_bound)]
CH_WT_fragments_subset = CH_WT_fragments[(CH_WT_fragments >= lower_bound) & (CH_WT_fragments <= upper_bound)]
ctDNA_mutant_fragments_subset = ctDNA_mutant_fragments[(ctDNA_mutant_fragments >= lower_bound) & (ctDNA_mutant_fragments <= upper_bound)]
ctDNA_WT_fragments_subset = ctDNA_WT_fragments[(ctDNA_WT_fragments >= lower_bound) & (ctDNA_WT_fragments <= upper_bound)]

# Generate KDE
x_range = np.linspace(100, 200, 1000)
ch_mutant_kde = gaussian_kde(CH_mutant_fragments_subset, bw_method=0.01)(x_range)
ch_wt_kde = gaussian_kde(CH_WT_fragments_subset, bw_method=0.01)(x_range)
ctdna_mutant_kde = gaussian_kde(ctDNA_mutant_fragments_subset, bw_method=0.01)(x_range)
ctdna_wt_kde = gaussian_kde(ctDNA_WT_fragments_subset, bw_method=0.01)(x_range)

kde_CH_bootstrap_differences, kde_CH_observed_difference, kde_CH_confidence_interval, kde_CH_p_value = run_bootstrap(ch_mutant_kde, ch_wt_kde, n_bootstraps = 10000)
kde_ctDNA_bootstrap_differences, kde_ctDNA_observed_difference, kde_ctDNA_confidence_interval, kde_ctDNA_p_value = run_bootstrap(ctdna_mutant_kde, ctdna_wt_kde, n_bootstraps = 10000)

# Using raw data:
# CH_bootstrap_differences, CH_observed_difference, CH_confidence_interval, CH_p_value = run_bootstrap(CH_mutant_fragments_subset, CH_WT_fragments_subset, n_bootstraps = 10000)
# ctDNA_bootstrap_differences, ctDNA_observed_difference, ctDNA_confidence_interval, ctDNA_p_value = run_bootstrap(ctDNA_mutant_fragments_subset, ctDNA_WT_fragments_subset, n_bootstraps = 10000)

# Conduct the Kruskal-Wallis Test - this is what I ended up using
stats.kruskal(ctdna_mutant_kde, ctdna_wt_kde)
stats.kruskal(ch_mutant_kde, ch_wt_kde)

