#!/usr/bin/env python3
# get the qc info for WGS data
# usage: total_qc_info.py <sample_info_table> <working_dir> <output_tsv>
# example dir: Project/cfDNA_fragmentomics/220615_control
# example: batch_qc_info.py sample_info/220615_sample_info.csv . total_qc.tsv


import os
import sys
import pandas as pd
import utils

input_file = os.path.abspath(sys.argv[1])
working_dir = os.path.abspath(sys.argv[2])
output_file = os.path.abspath(sys.argv[3])

input_table = utils.read_csv(input_file)

input_sample_lst = input_table['Sample_ID'].tolist()
qc_table_lst = [pd.read_csv(os.path.join(working_dir, _sample_id, "qc_summary.tsv"), sep="\t") for _sample_id in input_sample_lst]
total_table = pd.concat(qc_table_lst)
print(total_table.head())
total_table.to_csv(output_file, sep='\t', index=False)