#!/usr/bin/env python3
# usage: get_insert_size.py <sample_info_tsv> <input_dir> <output_shell_file>
# example dir: Project/cfDNA_fragmentomics/220805
# example: get_insert_size.py fq_matched.tsv . get_insert_size.sh

import os
import sys
import pandas as pd
import utils


sample_info_table = pd.read_csv(sys.argv[1], usecols=['Sample_ID'], sep='\t')
print(sample_info_table)
input_dir = os.path.abspath(sys.argv[2])
output_shell = sys.argv[3]
output_shell_fhand = open(output_shell, 'w')

sample_ids = sample_info_table['Sample_ID'].unique().tolist()
for sample_id in sample_ids:
    sample_dir = os.path.join(input_dir, sample_id, 'results')
    dedup_bam = os.path.join(sample_dir, f'{sample_id}.dedup.bam')
    if not utils.check_file(dedup_bam):
        continue
    insert_size_file = os.path.join(sample_dir, f'{sample_id}.insert_size.tsv')
    output_shell_fhand.write("samtools view -f66 {}  | cut -f9 | awk '{{print sqrt($0^2)}}' | sort | uniq -c  > {} &\n".format(dedup_bam, insert_size_file))
output_shell_fhand.write('wait\n')
output_shell_fhand.close()