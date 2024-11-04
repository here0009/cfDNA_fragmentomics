#!/usr/bin/env python3
# usage: aggregate_features.py <input_dir> <output_dir>
# example dir: Project/cfDNA_fragmentomics/220615_control/
# example: aggregate_features.py . total_features/


import os
import sys
import pandas as pd
import utils
import glob


def df_lst_aggregate(df_lst, id_col):
    total_df = pd.concat(df_lst, axis=0)
    _cols_set = set(total_df.columns)
    assert id_col in _cols_set
    _cols_set.remove(id_col)
    # print(f'columns of the output file is:')
    # print(sorted(list(_cols_set)))
    total_df = total_df.reindex([id_col] + sorted(list(_cols_set)), axis=1)
    
    return total_df

input_dir = os.path.abspath(sys.argv[1])
output_dir = os.path.abspath(sys.argv[2])
os.makedirs(output_dir, exist_ok=True)

READ_LEN_MIN = 1
READ_LEN_MAX = 500

read_length_csv = "read_length_summary.csv"
end_motif_4_csv = "end_motif_4_summary.csv"
end_motif_6_csv = "end_motif_6_summary.csv"
breakpoint_motif_csv = "break_point_motif_summary.csv"

# read length
count = 0
read_length_df_lst = []
for _file in glob.glob(os.path.join(input_dir, '**', read_length_csv), recursive=True):
    print(f'Parsing {_file}')
    # print(os.path.split(_file))
    sample_id = _file.split(os.path.sep)[-3]
    df = pd.read_csv(_file).set_index("read_length").transpose()
    df['sample_id'] = sample_id
    read_length_df_lst.append(df)
    count += 1
read_length_total_df = df_lst_aggregate(read_length_df_lst,  'sample_id')
for _col in read_length_total_df.columns:
    if _col != 'sample_id' and (_col < READ_LEN_MIN or _col > READ_LEN_MAX):
        read_length_total_df.drop(_col, axis=1, inplace=True)
read_length_total_df.to_csv(os.path.join(output_dir, f'total_{read_length_csv}'), index=False)
print(f"There are {count} samples of read length df in total.")
# end motif 4
end_motif_4_df_lst = []
count = 0
for _file in glob.glob(os.path.join(input_dir, '**', end_motif_4_csv), recursive=True):
    print(f'Parsing {_file}')
    sample_id = _file.split(os.path.sep)[-3]
    df = pd.read_csv(_file)
    total_counts = df['count'].sum()
    # df['count'] = df['count'] / total_counts * 4**4
    df2=df.set_index("end_motif_4").transpose()
    df2['sample_id'] = sample_id
    df2['total_count'] = total_counts
    end_motif_4_df_lst.append(df2)
    count += 1
end_motif_4_total_df = df_lst_aggregate(end_motif_4_df_lst, 'sample_id')
end_motif_4_total_df.to_csv(os.path.join(output_dir, f'total_{end_motif_4_csv}'), index=False)
print(f"There are {count} samples of end motif 4 df in total.")
# end motif 6
end_motif_6_df_lst = []
count = 0
for _file in glob.glob(os.path.join(input_dir, '**', end_motif_6_csv), recursive=True):
    print(f'Parsing {_file}')
    sample_id = _file.split(os.path.sep)[-3]
    df = pd.read_csv(_file)
    total_counts = df['count'].sum()
    # df['count'] = df['count'] / total_counts * 4**6
    df2=df.set_index("end_motif_6").transpose()
    df2['sample_id'] = sample_id
    df2['total_count'] = total_counts
    end_motif_6_df_lst.append(df2)
    count += 1
end_motif_6_total_df = df_lst_aggregate(end_motif_6_df_lst, 'sample_id')
end_motif_6_total_df.to_csv(os.path.join(output_dir, f'total_{end_motif_6_csv}'), index=False)
print(f"There are {count} samples of end motif 6 df in total.")
# beakpoint motif
breakpoint_motif_lst = []
count = 0
for _file in glob.glob(os.path.join(input_dir, '**', breakpoint_motif_csv), recursive=True):
    print(f'Parsing {_file}')
    sample_id = _file.split(os.path.sep)[-3]
    df = pd.read_csv(_file)
    total_counts = df['count'].sum()
    # df['count'] = df['count'] / total_counts * 4**6
    df2=df.set_index("break_point").transpose()
    df2['sample_id'] = sample_id
    df2['total_count'] = total_counts
    breakpoint_motif_lst.append(df2)
    count += 1
break_point_total_df = df_lst_aggregate(breakpoint_motif_lst, 'sample_id')
break_point_total_df.to_csv(os.path.join(output_dir, f'total_{breakpoint_motif_csv}'), index=False)
print(f"There are {count} samples of breakpoint motif df in total.")