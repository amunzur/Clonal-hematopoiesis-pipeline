# -*- coding: utf-8 -*-
"""
@author: Elie
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
import math
from scipy.stats import chi2
import gc
import datetime
pd.options.mode.chained_assignment = None

dir = os.path.dirname(__file__)
parser = argparse.ArgumentParser()
parser.add_argument("-in", help="input file" ,dest="input", type=str, required=True)
parser.add_argument("-out", help="output filename" ,dest="output", type=str, required=True)

args = parser.parse_args()
args.input = os.path.expanduser(args.input)
args.output = os.path.expanduser(args.output)
basename = os.path.splitext(args.input)[0]

def main():
	def load_vstat(path):
		#tnvstat columns needed
		# change to one dictionary, typing gets dict and cols gets keys.
		print(f"loading file {path} at {str(datetime.datetime.now())}")
		col_type = {"chrom": "category", "pos": np.int64, "ref": "category", "reads_all": np.int32, "matches": np.int32, "mismatches": np.int32, "deletions": np.int32, "insertions": np.int32, "A": np.int32, "C": np.int32, "T": np.int32, "G": np.int32, "N": np.int32, "sample": "category"}
		df = pd.read_csv(path, sep='\t', names=col_type.keys(), dtype=col_type, low_memory=True, engine='c', header=None)
		gc.collect()
		return df
		
	def compute_means(df):
		print(f"finding alternates for ref=A at {str(datetime.datetime.now())}")
		#find the number of each non reference base
		df['alt_AtoC'] = df['C'].where(df['ref'] == "A", 0)
		df['alt_AtoT'] = df['T'].where(df['ref'] == "A", 0)
		df['alt_AtoG'] = df['G'].where(df['ref'] == "A", 0)
		print(f"finding alternates for ref=T at {str(datetime.datetime.now())}")
		df['alt_TtoC'] = df['C'].where(df['ref'] == "T", 0)
		df['alt_TtoG'] = df['G'].where(df['ref'] == "T", 0)
		df['alt_TtoA'] = df['A'].where(df['ref'] == "T", 0)
		print(f"finding alternates for ref=G at {str(datetime.datetime.now())}")
		df['alt_GtoC'] = df['C'].where(df['ref'] == "G", 0)
		df['alt_GtoA'] = df['A'].where(df['ref'] == "G", 0)
		df['alt_GtoT'] = df['T'].where(df['ref'] == "G", 0)
		print(f"finding alternates for ref=C at {str(datetime.datetime.now())}")
		df['alt_CtoG'] = df['G'].where(df['ref'] == "C", 0)
		df['alt_CtoA'] = df['A'].where(df['ref'] == "C", 0)
		df['alt_CtoT'] = df['T'].where(df['ref'] == "C", 0)
		print(f"totals per position at {str(datetime.datetime.now())}")
		df['reads_total_for_pos'] = df['reads_all'].groupby(df['pos']).transform('sum')
		gc.collect()
			
		print(f"sum alt A {str(datetime.datetime.now())}")
		df['sum_alt_AtoC'] = df['alt_AtoC'].groupby(df['pos']).transform('sum')
		df['sum_alt_AtoT'] = df['alt_AtoT'].groupby(df['pos']).transform('sum')
		df['sum_alt_AtoG'] = df['alt_AtoG'].groupby(df['pos']).transform('sum')
		print(f"sum alt T {str(datetime.datetime.now())}")
		df['sum_alt_TtoC'] = df['alt_TtoC'].groupby(df['pos']).transform('sum')
		df['sum_alt_TtoG'] = df['alt_TtoG'].groupby(df['pos']).transform('sum')
		df['sum_alt_TtoA'] = df['alt_TtoA'].groupby(df['pos']).transform('sum')
		print(f"sum alt G {str(datetime.datetime.now())}")
		df['sum_alt_GtoC'] = df['alt_GtoC'].groupby(df['pos']).transform('sum')
		df['sum_alt_GtoA'] = df['alt_GtoA'].groupby(df['pos']).transform('sum')
		df['sum_alt_GtoT'] = df['alt_GtoT'].groupby(df['pos']).transform('sum')
		print(f"sum alt C {str(datetime.datetime.now())}")
		df['sum_alt_CtoG'] = df['alt_CtoG'].groupby(df['pos']).transform('sum')
		df['sum_alt_CtoA'] = df['alt_CtoA'].groupby(df['pos']).transform('sum')
		df['sum_alt_CtoT'] = df['alt_CtoT'].groupby(df['pos']).transform('sum')
		print(f"sum alt indels {str(datetime.datetime.now())}")
		df['sum_alt_del'] = df['deletions'].groupby(df['pos']).transform('sum')
		df['sum_alt_ins'] = df['insertions'].groupby(df['pos']).transform('sum')
		gc.collect()

		print(f"compute means for A>N {str(datetime.datetime.now())}")
		df['mean_errorAtoC'] = df['sum_alt_AtoC'] / df['reads_total_for_pos']
		df['mean_errorAtoT'] = df['sum_alt_AtoT'] / df['reads_total_for_pos']
		df['mean_errorAtoG'] = df['sum_alt_AtoG'] / df['reads_total_for_pos']
		print(f"compute means for T>N {str(datetime.datetime.now())}")
		df['mean_errorTtoC'] = df['sum_alt_TtoC'] / df['reads_total_for_pos']
		df['mean_errorTtoG'] = df['sum_alt_TtoG'] / df['reads_total_for_pos']
		df['mean_errorTtoA'] = df['sum_alt_TtoA'] / df['reads_total_for_pos']
		print(f"compute means for G>N {str(datetime.datetime.now())}")
		df['mean_errorGtoC'] = df['sum_alt_GtoC'] / df['reads_total_for_pos']
		df['mean_errorGtoA'] = df['sum_alt_GtoA'] / df['reads_total_for_pos']
		df['mean_errorGtoT'] = df['sum_alt_GtoT'] / df['reads_total_for_pos']
		print(f"compute means for C>N {str(datetime.datetime.now())}")
		df['mean_errorCtoG'] = df['sum_alt_CtoG'] / df['reads_total_for_pos']
		df['mean_errorCtoA'] = df['sum_alt_CtoA'] / df['reads_total_for_pos']
		df['mean_errorCtoT'] = df['sum_alt_CtoT'] / df['reads_total_for_pos']
		print(f"compute means for indel {str(datetime.datetime.now())}")
		df['mean_errordel'] = df['sum_alt_del'] / df['reads_total_for_pos']
		df['mean_errorins'] = df['sum_alt_ins'] / df['reads_total_for_pos']
		gc.collect()
		return df

	def output_file(df, path):
		df = df[["chrom", "pos", "ref", "mean_errorAtoC", "mean_errorAtoT", "mean_errorAtoG", "mean_errorTtoC", "mean_errorTtoG", "mean_errorTtoA", "mean_errorGtoC", "mean_errorGtoT", "mean_errorGtoA", "mean_errorCtoG", "mean_errorCtoT", "mean_errorCtoA", "mean_errordel", "mean_errorins"]]
		error_cols = ["mean_errorAtoC", "mean_errorAtoT", "mean_errorAtoG", "mean_errorTtoC", "mean_errorTtoG", "mean_errorTtoA", "mean_errorGtoC", "mean_errorGtoT", "mean_errorGtoA", "mean_errorCtoG", "mean_errorCtoT", "mean_errorCtoA", "mean_errordel", "mean_errorins"]
		df[error_cols] = df[error_cols].round(decimals = 4)
		df = df.drop_duplicates()
		print(f"writing file {path} at {str(datetime.datetime.now())}")
		df.to_csv(f"{path}", sep='\t', index=False)
	
	input_table = load_vstat(args.input)
	table = compute_means(input_table)
	output_file(table, args.output)
	
if __name__ == "__main__":
	main()
