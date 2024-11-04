#!/usr/bin/env python3
# usage: sample_info/rename_samples.py <sample_info_excel> <path> <rename_script>
# example working dir: sample_info/TCR/raw
# example: sample_info/rename_samples.py sequence_220706-E100047228_2.xls 220706-E100047228_local_path.tsv 220706-E100047228_rename.sh


import pandas as pd
import os
import sys
import re


sample_info_file = os.path.abspath(sys.argv[1])
local_path_file = os.path.abspath(sys.argv[2])
rename_sh_file = os.path.abspath(sys.argv[3])
rename_sh_fhand = open(rename_sh_file, 'w')
oss_dir = "oss://sz-hapseq/rawfq/MGI_merge/220706-E100047228/"

sample_info_table = pd.read_excel(sample_info_file, sheet_name=0)
print(sample_info_table.head())
local_path_table = pd.read_csv(local_path_file, sep='\t')
local_path_table['Sample_ID'] = local_path_table['Sample_ID'].apply(lambda x: x.replace('tmp_', ''))
print(local_path_table.head())
local_R1_dict = dict(zip(local_path_table['Sample_ID'], local_path_table['fq_R1']))
local_R2_dict = dict(zip(local_path_table['Sample_ID'], local_path_table['fq_R2']))


rename_dict = dict(zip(sample_info_table['Adapter'], sample_info_table['编号 ']))
print(rename_dict)
for key, val in rename_dict.items():
    if key not in local_R1_dict or key not in local_R2_dict:
        print('{} not in local path'.format(key))
        continue
    R1 = local_R1_dict[key]
    R2 = local_R2_dict[key]
    S_tag = re.search(r'S\d+', val).group(0)
    R1_path = os.path.dirname(R1)
    R2_path = os.path.dirname(R2)
    renamed_R1 = os.path.join(R1_path, f'{val}_{S_tag}_R1_001.fastq.gz')
    renamed_R2 = os.path.join(R2_path, f'{val}_{S_tag}_R2_001.fastq.gz') 
    rename_sh_fhand.write('mv {} {}\n'.format(R1, renamed_R1))
    rename_sh_fhand.write('mv {} {}\n'.format(R2, renamed_R2))
    rename_sh_fhand.write('ossutil cp {} {}\n'.format(renamed_R1, oss_dir))
    rename_sh_fhand.write('ossutil cp {} {}\n'.format(renamed_R2, oss_dir))
rename_sh_fhand.close()