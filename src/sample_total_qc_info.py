#!/usr/bin/env python3
# get the qc info for WGS data
# usage: total_qc_info.py <sample_id> <fastp_qc_json> <bwa_dedup_metrics> <mosdepth_summary> < mosdepth_threshold> <output_file>
# example dir: Project/cfDNA_fragmentomics/220615_control/S002_SZ20220425057WHB-5_cfdna_genome_72630
# example: total_qc_info.py S002_SZ20220425057WHB-5_cfdna_genome_72630 results/S002_SZ20220425057WHB-5_cfdna_genome_72630_fastp.json results/S002_SZ20220425057WHB-5_cfdna_genome_72630.dedup_metrics.txt QC/S002_SZ20220425057WHB-5_cfdna_genome_72630.mosdepth.summary.txt QC/S002_SZ20220425057WHB-5_cfdna_genome_72630.thresholds.bed.gz qc_summary.tsv


import os
import sys
from collections import OrderedDict
import pandas as pd
import json

class Fastp(object):

        def __init__(self, json_f):
                data = json.load(json_f)
                self.raw_reads_num = data["summary"]["before_filtering"]["total_reads"]
                self.raw_bases_num = data["summary"]["before_filtering"]["total_bases"]
                self.raw_q20_rate = data["summary"]["before_filtering"]["q20_rate"]
                self.raw_q30_rate = data["summary"]["before_filtering"]["q30_rate"]
                self.raw_gc_content = data["summary"]["before_filtering"]["gc_content"]
                self.clean_reads_num = data["summary"]["after_filtering"]["total_reads"]
                self.clean_bases_num = data["summary"]["after_filtering"]["total_bases"]
                self.clean_q20_rate = data["summary"]["after_filtering"]["q20_rate"]
                self.clean_q30_rate = data["summary"]["after_filtering"]["q30_rate"]
                self.clean_gc_content = data["summary"]["after_filtering"]["gc_content"]
                self.command = data["command"]
                try:
                        self.dup_rate = data["duplication"]["rate"]
                except:
                        self.dup_rate = data["summary"]["duplication"]["rate"]
                self.insert_size_peak = data["insert_size"]["peak"]

                self.r1_before_quality = data["read1_before_filtering"]["quality_curves"]
                self.r1_before_mean_quality = data["read1_before_filtering"]["quality_curves"]["mean"]
                self.r1_after_quality = data["read1_after_filtering"]["quality_curves"]
                self.r1_after_mean_quality = data["read1_after_filtering"]["quality_curves"]["mean"]
                self.r2_before_quality = data["read2_before_filtering"]["quality_curves"]
                self.r2_before_mean_quality = data["read2_before_filtering"]["quality_curves"]["mean"]
                self.r2_after_quality = data["read2_after_filtering"]["quality_curves"]
                self.r2_after_mean_quality = data["read2_after_filtering"]["quality_curves"]["mean"]

                self.r1_before_content = data["read1_before_filtering"]["content_curves"]
                self.r1_before_gc_content = data["read1_before_filtering"]["content_curves"]["GC"]
                self.r1_after_content = data["read1_after_filtering"]["content_curves"]
                self.r1_after_gc_content = data["read1_after_filtering"]["content_curves"]["GC"]
                self.r2_before_content = data["read2_before_filtering"]["content_curves"]
                self.r2_before_gc_content = data["read2_before_filtering"]["content_curves"]["GC"]
                self.r2_after_content = data["read2_after_filtering"]["content_curves"]
                self.r2_after_gc_content = data["read2_after_filtering"]["content_curves"]["GC"]

        def percent(self, num):
                return float(num) * 100

        def filter_rate(self, raw, clean):
                return (float(clean) / float(raw)) * 100


sample_id = sys.argv[1]
fastp_json = os.path.abspath(sys.argv[2])
bwa_dedup_metrics = os.path.abspath(sys.argv[3])
mosdepth_summary = os.path.abspath(sys.argv[4])
mosdepth_threshold = os.path.abspath(sys.argv[5])
output_file = os.path.abspath(sys.argv[6])

qc_info = OrderedDict()
qc_info['sample_id'] = sample_id

# fastp
fastp = Fastp(open(fastp_json, "r"))
qc_info['Raw_Yield(G)'] = fastp.raw_bases_num / 10**9
qc_info['Raw_Reads_Num(M)'] = fastp.raw_reads_num / 10**6
qc_info['Raw_Q30(%)'] = fastp.raw_q30_rate * 100
qc_info['Raw_Q20(%)'] = fastp.raw_q20_rate * 100
qc_info['Raw_GC(%)'] = fastp.raw_gc_content * 100
qc_info['Effective(%)'] = fastp.percent(fastp.clean_bases_num / fastp.raw_bases_num)
qc_info['Duplication_Rate_Raw(%)'] = fastp.dup_rate * 100
qc_info['Insert_Size_Peak(bp)'] = fastp.insert_size_peak

# bwa dedup metrics
bwa_dedup_metrics_table = pd.read_csv(bwa_dedup_metrics, sep="\t", header=1, nrows=2)
# print(bwa_dedup_metrics_table)
qc_info['Toatl_clean_read(M)'] = (bwa_dedup_metrics_table['UNPAIRED_READS_EXAMINED'] + 2*bwa_dedup_metrics_table['READ_PAIRS_EXAMINED'] + bwa_dedup_metrics_table['SECONDARY_OR_SUPPLEMENTARY_RDS']) / 10**6
qc_info['Mapping_rate(%)'] = (qc_info['Toatl_clean_read(M)'] - bwa_dedup_metrics_table['SECONDARY_OR_SUPPLEMENTARY_RDS'] / 10**6) / qc_info['Toatl_clean_read(M)'] * 100
qc_info["Duplication_Rate(%)"] = bwa_dedup_metrics_table['PERCENT_DUPLICATION'] * 100 

# mosdepth_summary
mosdepth_summary_table = pd.read_csv(mosdepth_summary, sep="\t")
# print(mosdepth_summary_table.head())
mosdepth_summary_table.set_index("chrom", inplace=True)
qc_info["Average_depth(X)"] = mosdepth_summary_table.loc['total', 'mean']

# mosdepth threshold
mosdepth_threshold_table = pd.read_csv(mosdepth_threshold, sep="\t")
# print(mosdepth_threshold_table.head())
# print(mosdepth_threshold_table.columns)
mosdepth_total_length = (mosdepth_threshold_table['end'] - mosdepth_threshold_table['start']).sum()
qc_info['Covreage_1X(%)'] = mosdepth_threshold_table['1X'].sum() / mosdepth_total_length * 100
qc_info['Covreage_5X(%)'] = mosdepth_threshold_table['5X'].sum() / mosdepth_total_length * 100
qc_info['Covreage_10X(%)'] = mosdepth_threshold_table['10X'].sum() / mosdepth_total_length * 100


output_table = pd.DataFrame(qc_info, index=[0])
output_table.to_csv(output_file, sep='\t', index=False, float_format='%.2f')

