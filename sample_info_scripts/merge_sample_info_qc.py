#!/usr/bin/env python3
# usage: merge_sample_info_qc.py <sample_info_csv> <qc_tsv> <output_tsv>
# example dir: results/total_qc
# example: merge_sample_info_qc.py sample_info/sample_info_csv/220809_total_sample_info.csv 220809_total_qc.tsv 220809_total_qc_sample_info.tsv


import os
import sys
import pandas as pd


sample_info_table = pd.read_csv(sys.argv[1], usecols=['Sample_ID','Name','type','tag'])
qc_table = pd.read_csv(sys.argv[2], sep='\t')
output_file = sys.argv[3]

sample_info_table.drop_duplicates(subset=['Sample_ID'], inplace=True)
output_table = pd.merge(sample_info_table, qc_table,right_on='sample_id', left_on='Sample_ID')
output_table.drop(columns=['sample_id'], inplace=True)
output_table.to_csv(output_file, sep='\t', index=False)
