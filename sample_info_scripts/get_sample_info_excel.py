#!/usr/bin/env python3
# usage: get_sample_info_excel.py <template> <sample_lst> <input_sample_info> <output_sample_info>
# example working dir: sample_info/TCR/raw/220718
# example: sample_info/get_sample_info_excel.py 220718_cfdna_fragmentomics_tmplate.xlsx 220718_total_sample.xls 220718_cfDNA_fragmentomics_sample_info.csv
import pandas as pd
import os
import sys


template_file = os.path.abspath(sys.argv[1])
# sample_lst_file = os.path.abspath(sys.argv[2])
input_sample_info_file = os.path.abspath(sys.argv[2])
output_file = os.path.abspath(sys.argv[3])


template_table = pd.read_excel(template_file, sheet_name=0)
sample_names = template_table['样本姓名'].tolist()

# sample_names = pd.read_csv(sample_lst_file, header=None)[0].tolist()
print(f'There are {len(sample_names)} samples in the sample list, there are {len(set(sample_names))} unique sample names')
input_sample_table = pd.read_excel(input_sample_info_file, sheet_name=0)
print(template_table.head())
print(input_sample_table.head())
input_names = input_sample_table['受检者姓名'].tolist()
print(f'There are {len(input_names)} samples in the input sample info table') 
print(f'there are {len(set(input_names))} unique sample names')
age_dict = dict(zip(input_sample_table['受检者姓名'], input_sample_table['年龄']))
sample_id_dict = dict(zip(input_sample_table['受检者姓名'], input_sample_table['送检单号']))
gender_dict = dict(zip(input_sample_table['受检者姓名'], input_sample_table['性别']))
telephone_dict = dict(zip(input_sample_table['受检者姓名'], input_sample_table['联系电话']))
health_status_dict = dict(zip(input_sample_table['受检者姓名'], input_sample_table['情况']))

# template_table['送检单号'] = template_table['样本姓名'].apply(lambda x: sample_id_dict[x])
template_table['年龄'] = template_table['样本姓名'].apply(lambda x: age_dict[x])
template_table['性别'] = template_table['样本姓名'].apply(lambda x: gender_dict[x])
template_table['电话'] = template_table['样本姓名'].apply(lambda x: telephone_dict[x])
template_table['癌种'] = template_table['样本姓名'].apply(lambda x: health_status_dict[x])

template_table.to_csv(output_file, index=False, encoding='gbk') # use gbk encoding to avoid as required by hapyun https://mall.hapyun.com/mall/order/manage
