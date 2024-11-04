#!/usr/bin/env python3
# usage: fragmentomics_features.py <input_bam> <output_dir>
# example dir: Project/cfDNA_fragmentomics/test/S097_SZ20220425072WHB-d_cfdna_genome_71498
# example: fragmentomics_features.py results/S097_SZ20220425072WHB-d_cfdna_genome_71498.dedup.bam fragmentomics_features/


import itertools
import os
import sys
import pandas as pd
import utils
import pysam
from collections import Counter, defaultdict
import itertools

# os['environ'] += os.path.sep + '/miniconda3/envs/bioinfo/bin'
samfile = pysam.AlignmentFile(sys.argv[1], "rb")
output_dir = os.path.abspath(sys.argv[2])
os.makedirs(output_dir, exist_ok=True)

END_MOTIF_LEN = 5
READ_LEN_MIN = 1
READ_LEN_MAX = 500
MAP_QUAL_MIN = 30
FRAG_SIZE_BIN = 5*10**6
FRAG_SIZE_CUT_OFF = 150
CHROMS = set([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM'])

read_length_counter = Counter()
combs = itertools.product('ACGT', repeat=END_MOTIF_LEN)
end_motif_counter = {''.join(i):0 for i in combs}
frag_size_counter = dict()

print(f'There are {len(end_motif_counter)} for {END_MOTIF_LEN} end motifs.')
for read in samfile:
    if (not read.is_proper_pair) or (read.mapping_quality < MAP_QUAL_MIN) or (read.is_unmapped) or (read.reference_name not in CHROMS):
        continue
    ends = read.query_sequence[:END_MOTIF_LEN]
    if ends in end_motif_counter:
        end_motif_counter[ends] += 1
    if read.is_read1:
        read_len = abs(read.template_length)
        # read_len = read.inferred_length
        # print(read_len)
        if read_len > READ_LEN_MIN and read_len < READ_LEN_MAX:
            # total fragment length distribution
            read_length_counter[read_len] += 1
            chrom, pos = read.reference_name, read.reference_start
            # short to long reads ratio in each fragment size bin
            frag_bin = read.reference_start // FRAG_SIZE_BIN
            if (chrom, frag_bin) not in frag_size_counter:
                frag_size_counter[(chrom, frag_bin)] = [0, 0]
            if read_len > FRAG_SIZE_CUT_OFF:
                frag_size_counter[(chrom, frag_bin)][1] += 1
            else:
                frag_size_counter[(chrom, frag_bin)][0] += 1

print(end_motif_counter)
end_motif_counter_df = pd.DataFrame.from_dict(end_motif_counter, orient='index')
end_motif_counter_df.index.name = "end_motif"
end_motif_counter_df.columns = ["count"]
end_motif_counter_df.sort_index(inplace=True)
end_motif_counter_df.to_csv(os.path.join(output_dir,"end_motif_summary.csv"))

print(read_length_counter)
read_length_counter_df = pd.DataFrame.from_dict(read_length_counter, orient="index")
read_length_counter_df.columns = ["count"]
read_length_counter_df.index.name = "read_length"
read_length_counter_df.sort_index(inplace=True)
read_length_counter_df.to_csv(os.path.join(output_dir,"read_length_summary.csv"))

print(frag_size_counter)
frag_size_fhand = open(os.path.join(output_dir,"frag_size_summary.csv"), "w")
# frag_size_fhand.write("chrom,frag_bin,short_reads,long_reads\n")
for key in frag_size_counter.keys():
    frag_size_fhand.write('\t'.join([str(key[0]), str(key[1]), str(frag_size_counter[key][0]), str(frag_size_counter[key][1])]) + '\n')
    
frag_size_fhand.close()
samfile.close()
