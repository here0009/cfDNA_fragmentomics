#!/usr/bin/env python3
# usage: merge_samples.py <input_sample_info> <output_sample_info>
# example dir: Project/cfDNA_fragmentomics/220718_merged_cancer
# example: merge_samples.py fq_matched.tsv fq_matched_merged.tsv rawfq/cfDNA_fragmentomics/220718_merged_cancer


import os
import sys
import pandas as pd

input_table = pd.read_csv(sys.argv[1], sep='\t')
output_file = os.path.abspath(sys.argv[2])
working_dir = os.path.abspath(sys.argv[3])
merge_script = os.path.join(working_dir, 'merge_samples.sh')
merge_script_fhand = open(merge_script, 'w')

print(input_table.head())
grouped_table = input_table.groupby('Lib_number').agg({
    'Sed_ID':'first',
    'Sample_ID':'first',
    'Lib_number':'first',
    'Name':'first',
    'Panel':'first',
    'Data_ID':'first',
    'fq_R1':lambda x: ';'.join(x),
    'fq_R2':lambda x: ';'.join(x)
})
print(grouped_table.head())
grouped_table['fq_R1_sep'] = grouped_table['fq_R1']
grouped_table['fq_R2_sep'] = grouped_table['fq_R2']
grouped_table['fq_R1'] = grouped_table['Lib_number'].apply(lambda x: os.path.join(working_dir, f'total_{x}_R1_001.fastq.gz'))
grouped_table['fq_R2'] = grouped_table['Lib_number'].apply(lambda x: os.path.join(working_dir, f'total_{x}_R2_001.fastq.gz'))

for _, row in grouped_table.iterrows():
    merge_script_fhand.write(f'cat {row["fq_R1_sep"].replace(";"," ")} > {row["fq_R1"]}\n')
    merge_script_fhand.write(f'cat {row["fq_R2_sep"].replace(";"," ")} > {row["fq_R2"]}\n')

merge_script_fhand.close()
grouped_table.to_csv(output_file, index=False, sep='\t')
