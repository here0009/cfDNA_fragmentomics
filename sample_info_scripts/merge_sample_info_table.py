#!/usr/bin/env python3
# usage: get_frag_features_from_bed.py <input_bed> <output_dir>
# example dir: sample_info
# example: merge_sample_info_table.py 220628_cancer_sample_info.csv 220708_control_sample_info.csv 220708_total_sample_info.csv


import os
import sys
import pandas as pd

cancer_table = pd.read_csv(sys.argv[1])
control_table = pd.read_csv(sys.argv[2])
output_file = sys.argv[3]

cancer_table['type'] = 'cancer'
control_table['type'] = 'control'

output_table = pd.concat([cancer_table, control_table])
output_table.to_csv(output_file, index=False)