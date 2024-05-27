#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 15:17:05 2024

@author: jbacon
"""

import pandas as pd

error=pd.read_csv("/groups/wyattgrp/users/jbacon/err_rates/prostate_exons/nextseq/error_rate/error_rates.tsv", sep='\t')

error['mean_errorindel'] = error[['mean_errordel', 'mean_errorins']].max(axis=1)
#error.insert(2,'end',error['pos'])
#error=error.astype({"pos":"Int64"})
#error=error.astype({"end":"Int64"})

error.rename(columns={'chrom': 'CHROM'}, inplace=True)
error.rename(columns={'pos': 'POS'}, inplace=True)
error.rename(columns={'ref': 'REF'}, inplace=True)

error.to_csv("/groups/wyattgrp/users/jbacon/err_rates/prostate_exons/nextseq/error_rate/error_rates_fixed.tsv", sep='\t')