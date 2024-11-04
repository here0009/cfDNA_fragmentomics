#!/usr/bin/env python3
# convert excel table to csv file
# usage: convert_sample_info_table.py <input_table> <output_table>
# example: convert_sample_info_table.py cfDNA_fragmentomics/sample_info/raw/CRC_frag_screen_cancer_sample_info_20220426.xlsx cfDNA_fragmentomics/sample_info/cancer_sample_info_220610.csv


import os
import sys
import pandas as pd


input_file = os.path.abspath(sys.argv[1])
output_file = os.path.abspath(sys.argv[2])

selected_cols = ['sed_id','sample_id','lib_number','panel','name','data_id']
converted_cols = ['Sed_ID','Sample_ID','Lib_number','Name','Panel','Data_ID']
input_table = pd.read_excel(input_file, sheet_name=0, usecols=selected_cols)
input_table.columns = converted_cols
input_table.to_csv(output_file, sep=',', index=False)

