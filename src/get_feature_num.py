#!/usr/bin/env python3
# get the qc info for WGS data
# usage: stat_scripts/get_feature_num.py <input_dir> <output_tsv>
# example dir: results/end_motif_220809/all_feature_selection
# example: stat_scripts/get_feature_num.py . feature_num.tsv


import os
import sys
import pandas as pd
from glob import glob

input_dir = os.path.abspath(sys.argv[1])
output_file = os.path.abspath(sys.argv[2])
FILE_ENDS = 'feature_selection.csv'
tables = []

for _file in glob(os.path.join(input_dir, "**", FILE_ENDS), recursive=True):
    print(f'Parsing {_file}')
    feature_table = pd.read_csv(_file, usecols=['method','feature_number']).set_index('method').T
    tag = os.path.dirname(_file).split('/')[-1]
    feature_table['tag'] = tag
    tables.append(feature_table)

output_table = pd.concat(tables)
cols = output_table.columns.tolist()
output_table = output_table[cols[-1:] + cols[:-1]]
output_table.to_csv(output_file, sep='\t', index=False)
