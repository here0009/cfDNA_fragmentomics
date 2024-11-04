#!/usr/bin/env python3
# usage: get_frag_features_from_bed.py <input_bed> <output_dir>
# example dir: Project/cfDNA_fragmentomics/220629_control/S091_SZ20220426016WHB-1_cfdna_genome_75026/results
# example: get_frag_features_from_bed.py results/S097_SZ20220425072WHB-d_cfdna_genome_71498.dedup.bam fragmentomics_features/


import itertools
import os
import sys
import pandas as pd
import utils
import pysam
from collections import Counter, defaultdict
import itertools
import time


def get_break_point(break_point, end_motif):
    # concate the seq before and after the break point
    return break_point[-BREAKPOINT_LEN:] + end_motif[:BREAKPOINT_LEN]

def add_dict(dict_a, dict_b):
    for key,value in dict_b.items():
        if key in dict_a:
            dict_a[key] += value
    return dict_a


start_time = time.time()
# os['environ'] += os.path.sep + '/miniconda3/envs/bioinfo/bin'
bedfile = os.path.abspath(sys.argv[1])
output_dir = os.path.abspath(sys.argv[2])
hg19_fa_file = os.path.abspath(sys.argv[3])
# hg19_fa_file = "/data/hg19_fragmentomics/hg19.fa"

bed_table = pd.read_csv(bedfile, sep='\t', header=None)
bed_table.columns = ["chr", "start", "end", "name", "score", "strand", "breakpoint_F", "end_motif_F", "breakpoint_R", "end_motif_R", "gc"]
bed_table.drop(["name","score","strand","gc"], axis=1, inplace=True) # drop unnecessary columns
os.makedirs(output_dir, exist_ok=True)
# hg19_fa.fetch("chr1", 750101, 750107)

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
end_motif_6_counter = add_dict(end_motif_6_counter, dict(bed_table['end_motif_F'].value_counts()))
end_motif_6_counter = add_dict(end_motif_6_counter, dict(bed_table['end_motif_R'].value_counts()))

end_motif_6_counter_df = pd.DataFrame.from_dict(end_motif_6_counter, orient='index')
end_motif_6_counter_df.index.name = "end_motif_6"
end_motif_6_counter_df.columns = ["count"]
end_motif_6_counter_df.sort_index(inplace=True)
end_motif_6_counter_df.to_csv(os.path.join(output_dir,"end_motif_6_summary.csv"))
print(f'There are {len(end_motif_6_counter)} for {END_MOTIF_6} end motifs.')

# end motif 4
bed_table['end_motif_F'] = bed_table['end_motif_F'].apply(lambda x: x[:END_MOTIF_4])
bed_table['end_motif_R'] = bed_table['end_motif_R'].apply(lambda x: x[:END_MOTIF_4])
end_motif_4_counter = add_dict(end_motif_4_counter, dict(bed_table['end_motif_F'].value_counts()))
end_motif_4_counter = add_dict(end_motif_4_counter, dict(bed_table['end_motif_R'].value_counts()))
end_motif_4_counter_df = pd.DataFrame.from_dict(end_motif_4_counter, orient='index')
end_motif_4_counter_df.index.name = "end_motif_4"
end_motif_4_counter_df.columns = ["count"]
end_motif_4_counter_df.sort_index(inplace=True)
end_motif_4_counter_df.to_csv(os.path.join(output_dir,"end_motif_4_summary.csv"))
print(f'There are {len(end_motif_4_counter)} for {END_MOTIF_4} end motifs.')


# read length
bed_table['read_len'] = bed_table['end'] - bed_table['start']
read_length_counter = dict(bed_table['read_len'].value_counts())
read_length_counter_df = pd.DataFrame.from_dict(read_length_counter, orient="index")
read_length_counter_df.columns = ["count"]
read_length_counter_df.index.name = "read_length"
read_length_counter_df.sort_index(inplace=True)
read_length_counter_df.to_csv(os.path.join(output_dir,"read_length_summary.csv"))


# breakpoint
bed_table['breakpoint_F'] = bed_table.apply(lambda x: get_break_point(x['breakpoint_F'], x['end_motif_F']), axis=1) # left breakpoint
bed_table['breakpoint_R'] = bed_table.apply(lambda x: get_break_point(x['breakpoint_R'], x['end_motif_R']), axis=1) # right breakpoint
break_point_counter = add_dict(break_point_counter, dict(bed_table['breakpoint_F'].value_counts()))
break_point_counter = add_dict(break_point_counter, dict(bed_table['breakpoint_R'].value_counts()))
break_point_counter_df = pd.DataFrame.from_dict(break_point_counter, orient='index')
break_point_counter_df.index.name = "break_point"
break_point_counter_df.columns = ["count"]
break_point_counter_df.sort_index(inplace=True)
break_point_counter_df.to_csv(os.path.join(output_dir,"break_point_motif_summary.csv"))
print(f'There are {len(break_point_counter)} for {BREAKPOINT_LEN * 2} breakpoint motifs.')

end_time = time.time()
minutes, seconds = divmod(end_time - start_time, 60)
print(f'Time used: {minutes} minutes {seconds:.2f} seconds.')