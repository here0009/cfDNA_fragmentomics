#!/usr/bin/env python3
# usage: dask_get_frag_features_from_bed.py <input_bed> <output_dir> <threads>
# example dir: Project/cfDNA_fragmentomics/220629_control/S093_SZ20220426013WHB-d_cfdna_genome_75028/results/test
# example: dask_get_frag_features_from_bed.py test.bed . 6


import itertools
import os
import sys
import dask.dataframe as dd
import pandas as pd
import utils
from collections import Counter, defaultdict
import itertools
import dask
from multiprocessing.pool import ThreadPool
import time


def get_break_point(break_point, end_motif):
    # concate the seq before and after the break point
    return break_point[-BREAKPOINT_LEN:] + end_motif[:BREAKPOINT_LEN]

def add_dict(dict_a, dict_b):
    for key,value in dict_b.items():
        if key in dict_a:
            dict_a[key] += value
    return dict_a

# os['environ'] += os.path.sep + '/miniconda3/envs/bioinfo/bin'
start_time = time.time()
bedfile = os.path.abspath(sys.argv[1])
output_dir = os.path.abspath(sys.argv[2])
threads = int(sys.argv[3])
dask.config.set(pool=ThreadPool(threads))
os.makedirs(output_dir, exist_ok=True)


bed_table = dd.read_csv(bedfile, sep='\t', header=None)
bed_table.columns = ["chr", "start", "end", "name", "score", "strand", "breakpoint_F", "end_motif_F", "breakpoint_R", "end_motif_R", "gc"]
bed_table = bed_table.drop(["name","score","strand","gc"], axis=1) # drop unnecessary columns
os.makedirs(output_dir, exist_ok=True)
# hg19_fa.fetch("chr1", 750101, 750107)
bed_table = bed_table.repartition(npartitions=threads)
bed_table = bed_table.persist()

END_MOTIF_4 = 4
END_MOTIF_6 = 6
READ_LEN_MIN = 1
READ_LEN_MAX = 500
MAP_QUAL_MIN = 30
FRAG_SIZE_BIN = 5*10**6
FRAG_SIZE_CUT_OFF = 150
BREAKPOINT_LEN = 3 # the length of breakpoint on each side, total length is 2*BREAKPOINT_LEN
CHROMS = set([f'chr{i}' for i in range(1, 23)] + ['chrX'])

read_length_counter = Counter()
combs_4 = itertools.product('ACGT', repeat=END_MOTIF_4)
combs_6 = itertools.product('ACGT', repeat=END_MOTIF_6)
comb_bp = itertools.product('ACGT', repeat=BREAKPOINT_LEN * 2)
end_motif_6_counter = {''.join(i):0 for i in combs_6}
end_motif_4_counter = {''.join(i):0 for i in combs_4}
break_point_counter = {''.join(i):0 for i in comb_bp}


# end motif 6
end_motif_6_counter = add_dict(end_motif_6_counter, dict(bed_table['end_motif_F'].value_counts().compute()))
end_motif_6_counter = add_dict(end_motif_6_counter, dict(bed_table['end_motif_R'].value_counts().compute()))

end_motif_6_counter_df = pd.DataFrame.from_dict(end_motif_6_counter, orient='index')
end_motif_6_counter_df.index.name = "end_motif_6"
end_motif_6_counter_df.columns = ["count"]
end_motif_6_counter_df.sort_index(inplace=True)
end_motif_6_counter_df.to_csv(os.path.join(output_dir,"end_motif_6_summary.csv"))
print(f'There are {len(end_motif_6_counter)} for {END_MOTIF_6} end motifs.')

# breakpoint
bed_table['breakpoint_F'] = bed_table.apply(lambda x: get_break_point(x['breakpoint_F'], x['end_motif_F']), axis=1, meta=('breakpoint_F', 'object')) # left breakpoint
bed_table['breakpoint_R'] = bed_table.apply(lambda x: get_break_point(x['breakpoint_R'], x['end_motif_R']), axis=1, meta=('breakpoint_R', 'object')) # right breakpoint
break_point_counter = add_dict(break_point_counter, dict(bed_table['breakpoint_F'].value_counts().compute()))
break_point_counter = add_dict(break_point_counter, dict(bed_table['breakpoint_R'].value_counts().compute()))
break_point_counter_df = pd.DataFrame.from_dict(break_point_counter, orient='index')
break_point_counter_df.index.name = "break_point"
break_point_counter_df.columns = ["count"]
break_point_counter_df.sort_index(inplace=True)
break_point_counter_df.to_csv(os.path.join(output_dir,"break_point_motif_summary.csv"))
print(f'There are {len(break_point_counter)} for {BREAKPOINT_LEN * 2} breakpoint motifs.')

# end motif 4
bed_table['end_motif_F'] = bed_table['end_motif_F'].apply(lambda x: x[:END_MOTIF_4], meta=('end_motif_F', 'object'))
bed_table['end_motif_R'] = bed_table['end_motif_R'].apply(lambda x: x[:END_MOTIF_4], meta=('end_motif_R', 'object'))
end_motif_4_counter = add_dict(end_motif_4_counter, dict(bed_table['end_motif_F'].value_counts().compute()))
end_motif_4_counter = add_dict(end_motif_4_counter, dict(bed_table['end_motif_R'].value_counts().compute()))
end_motif_4_counter_df = pd.DataFrame.from_dict(end_motif_4_counter, orient='index')
end_motif_4_counter_df.index.name = "end_motif_4"
end_motif_4_counter_df.columns = ["count"]
end_motif_4_counter_df.sort_index(inplace=True)
end_motif_4_counter_df.to_csv(os.path.join(output_dir,"end_motif_4_summary.csv"))
print(f'There are {len(end_motif_4_counter)} for {END_MOTIF_4} end motifs.')


# read length
bed_table['read_len'] = bed_table['end'] - bed_table['start']
read_length_counter = dict(bed_table['read_len'].value_counts().compute())
read_length_counter_df = pd.DataFrame.from_dict(read_length_counter, orient="index")
read_length_counter_df.columns = ["count"]
read_length_counter_df.index.name = "read_length"
read_length_counter_df.sort_index(inplace=True)
read_length_counter_df.to_csv(os.path.join(output_dir,"read_length_summary.csv"))


end_time = time.time()
minutes, seconds = divmod(end_time - start_time, 60)
print(f'Time used: {minutes} minutes {seconds:.2f} seconds.')