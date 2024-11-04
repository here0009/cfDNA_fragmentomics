#!/usr/bin/env python3
# usage: filter_bam2bed_v2.py <input_bam> <filtered_bam_bed> <breakpoint_bed>
# example dir: Project/cfDNA_fragmentomics/test/S097_SZ20220425072WHB-d_cfdna_genome_71498
# example: filter_bam2bed_v2.py results/S097_SZ20220425072WHB-d_cfdna_genome_71498.dedup.bam filtered_bam.bed
# add break point and gc content information to bed file


from builtins import breakpoint
import os
import sys
import pysam
from collections import defaultdict
import utils


samfile = pysam.AlignmentFile(sys.argv[1], "rb")
filtered_bam_bed = os.path.abspath(sys.argv[2])
hg19_fa_file = os.path.abspath(sys.argv[3])
filtered_bam_bed_fhand = open(filtered_bam_bed, 'w')
# hg19_fa_file = "/data/hg19_fragmentomics/hg19.fa"
hg19_fa = pysam.FastaFile(hg19_fa_file)

READ_LEN_MIN = 1
READ_LEN_MAX = 500
MAP_QUAL_MIN = 20
END_MOTIF_LEN = 6
CHROMS = set([f'chr{i}' for i in range(1, 23)] + ['chrX'])
BREAKPOINT_LEN = 3 # the length of breakpoint on each side, total length is 2*BREAKPOINT_LEN
# bed header: chr, start, end, name, score, strand, breakpoint_F, end_motif_F, breakpoint_R, end_motif_R, gc
# print(CHROMS)
filtered_bam_bed_fhand.write('\t'.join(['#chr', 'start', 'end', 'name', 'score', 'strand', 'breakpoint_F', 'end_motif_F', 'breakpoint_R', 'end_motif_R', 'gc']) + '\n')

def get_seq(chromsome, start, end):
    try:
        seq = hg19_fa.fetch(chromsome, start, end).upper()
    except:
        print(chromsome, start, end, "got no feteched seq")
        return None
    return seq

def isSoftClipped(cigar):
    """
    see here for more information about this function
    references:
        https://pysam.readthedocs.io/en/latest/api.html
        https://davetang.org/wiki/tiki-index.php?page=SAM
    """
    for (op, _) in cigar:
        if op in [4, 5, 6]:
            return True
    return False

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    reference:
        https://www.biostars.org/p/306041/
    with small modifications:
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        # filter reads
        if read.is_unmapped or read.mate_is_unmapped or read.is_qcfail:
            continue
        if not read.is_paired:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if read.next_reference_name != read.reference_name or read.reference_name not in CHROMS : # not same chromosome
            # print(read.reference_id)
            continue
        if abs(read.template_length) < READ_LEN_MIN or abs(read.template_length) > READ_LEN_MAX:
            continue
        if isSoftClipped(read.cigar):
            continue
        qname = read.query_name
        # print(qname, read.is_read1, read.is_read2)
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


for read1, read2 in read_pair_generator(samfile):
    r1_start, r1_end = read1.reference_start, read1.reference_end
    r2_start, r2_end = read2.reference_start, read2.reference_end
    r1_mapq, r2_mapq = read1.mapping_quality, read2.mapping_quality
    if r1_mapq < MAP_QUAL_MIN or r2_mapq < MAP_QUAL_MIN: # map quality filter
        continue
    if read1.is_reverse:
        start, end = r2_start, r1_end
        strand = '-'
        end_motif_F = read2.query_sequence[:END_MOTIF_LEN]
        end_motif_R = utils.reverse_complement(read1.query_sequence[-END_MOTIF_LEN:])
    else:
        start, end = r1_start, r2_end
        strand = '+'
        end_motif_F = read1.query_sequence[:END_MOTIF_LEN]
        end_motif_R = utils.reverse_complement(read2.query_sequence[-END_MOTIF_LEN:])
    if start < 0 or end < 0 or start >= end: # check if this happens
        continue
    mapq = min(r1_mapq, r2_mapq)
    bp_template = get_seq(read1.reference_name, start - BREAKPOINT_LEN, end + BREAKPOINT_LEN) # breakpoint + template + breakpoint
    if bp_template is None:
        continue
    
    breakpoint_F = bp_template[:BREAKPOINT_LEN]
    breakpoint_R = utils.reverse_complement((bp_template[-BREAKPOINT_LEN:]))
    tmplate = bp_template[BREAKPOINT_LEN:-BREAKPOINT_LEN]
    # end_motif_F = tmplate[:END_MOTIF_LEN]
    # end_motif_R = utils.reverse_complement(tmplate[-END_MOTIF_LEN:])
    gc = utils.gc_content(tmplate)
    filtered_bam_bed_fhand.write('\t'.join([read1.reference_name, str(start), str(end), read1.query_name, str(mapq), strand, breakpoint_F, end_motif_F, breakpoint_R, end_motif_R, f'{gc:.3f}']) + '\n')

filtered_bam_bed_fhand.close()
samfile.close()
