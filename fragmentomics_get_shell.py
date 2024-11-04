#!/usr/bin/env python3
# Version 1.0
# usage: fragmentomics_get_shell.py <sample_info_file> <config_file> <output_dir>
# sample_info_file : tsv file of sample info, columns: Sample_ID, fq_R1, fq_R2
# config_file: configuration file for the software
# output_dir: output dir




import os
import sys
import pandas as pd
import utils
import json


sample_file = sys.argv[1]
config_file = sys.argv[2]
output_dir = os.path.abspath(sys.argv[3])

sample_table = utils.read_csv(sample_file)
config_fhand = open(config_file, "r")
config_dict = json.load(config_fhand)
# print(config_dict)
all_run_script = os.path.join(output_dir, "all_run_script.sh")
all_run_fhand = open(all_run_script, "w")
batch_size = config_dict.get("batch_size", 3) # default batch_size is 3
sample_counts = 0
print('\n'.join([f'#{_k}\t{_v}'  for _k, _v in config_dict.items()]) + '\n')

for _, row in sample_table.iterrows():
    sample_id = row['Sample_ID']
    fq1 = row.get('fq_R1', None)
    fq2 = row.get('fq_R2', None)
    fq1_base = os.path.basename(fq1) if fq1 else None
    fq2_base = os.path.basename(fq2) if fq2 else None
    utils.check_file(fq1)
    utils.check_file(fq2)
    sample_dir = os.path.join(output_dir, sample_id)
    sub_dir_lst = ['results', 'CNA', 'QC', 'frag_features', 'length_ratio']
    sub_dir_dict = {_dir: os.path.join(sample_dir, _dir) for _dir in sub_dir_lst}
    for sub_dir in sub_dir_dict.values():
        os.makedirs(sub_dir, exist_ok=True)

    sample_counts += 1
    sample_shell_script = os.path.join(sample_dir, "run.sh")
    all_run_fhand.write(f'sh {sample_shell_script} > {sample_shell_script}.o 2> {sample_shell_script}.e &\n')
    if sample_counts % batch_size == 0:
        all_run_fhand.write("wait\n")
    sample_shell_fhand = open(sample_shell_script, "w")
    sample_shell_fhand.write(f'# Generating shell script of cfDNA fragmentomics for: {sample_id}\n')
    
    # fastp
    sample_shell_fhand.write("# Step 1 Remove adapter\n")
    sample_shell_fhand.write(f"{config_dict['fastp']} -i {fq1} -I {fq2} -o {sub_dir_dict['results']}/{fq1_base} -O {sub_dir_dict['results']}/{fq2_base} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -h {sub_dir_dict['results']}/{sample_id}_fastp.html -j {sub_dir_dict['results']}/{sample_id}_fastp.json\n")
    
    # # bwa
    sample_shell_fhand.write("# Step 2 Mapping\n")
    sample_shell_fhand.write(f"{config_dict['export_sentieon_license']}\n")
    bwa_header_info = f"@RG\\tID:{sample_id}\\tLB:{sample_id}\\tPL:ILLUMINA\\tSM:{sample_id}"
    sample_shell_fhand.write(f"{config_dict['sentieon']} bwa mem -R \"{bwa_header_info}\" -t {config_dict['threads']} -k 32 -M {config_dict['hg19']} {sub_dir_dict['results']}/{fq1_base} {sub_dir_dict['results']}/{fq2_base} | {config_dict['sentieon']} util sort -o {sub_dir_dict['results']}/{sample_id}.sort.bam -t {config_dict['threads']} --sam2bam -i - \n")
    
    # # dedup
    sample_shell_fhand.write("# Step 3 Deduplication\n")
    sample_shell_fhand.write(f"{config_dict['sentieon']} driver -r {config_dict['hg19']} -t {config_dict['threads']} -i {sub_dir_dict['results']}/{sample_id}.sort.bam  --algo LocusCollector --fun score_info {sub_dir_dict['results']}/{sample_id}.score.txt\n")
    sample_shell_fhand.write(f"{config_dict['sentieon']} driver -t {config_dict['threads']} -i {sub_dir_dict['results']}/{sample_id}.sort.bam --algo Dedup --rmdup --score_info {sub_dir_dict['results']}/{sample_id}.score.txt --metrics {sub_dir_dict['results']}/{sample_id}.dedup_metrics.txt {sub_dir_dict['results']}/{sample_id}.dedup.bam 2>{sub_dir_dict['results']}/{sample_id}.bam.dup.log\n")
    
    # QC
    sample_shell_fhand.write("# Step 3 Quality control\n")
    sample_shell_fhand.write(f"cd {sub_dir_dict['QC']} && {config_dict['mosdepth']} -n -t {min(4, config_dict['threads'])} -b {config_dict['hg19_bed']} --thresholds 1,5,10 {sample_id} {sub_dir_dict['results']}/{sample_id}.dedup.bam\n") # more than 4 threads is useless according to the documentation of mosdepth
    # # total QC
    sample_shell_fhand.write(f'{config_dict["total_qc"]} {sample_id} {sub_dir_dict["results"]}/{sample_id}_fastp.json {sub_dir_dict["results"]}/{sample_id}.dedup_metrics.txt {sub_dir_dict["QC"]}/{sample_id}.mosdepth.summary.txt {sub_dir_dict["QC"]}/{sample_id}.thresholds.bed.gz {sample_dir}/qc_summary.tsv\n')


    # CNV
    sample_shell_fhand.write("# Step 5 CNV analysis\n")
    # wig and ichorCNA
    sample_shell_fhand.write(f"export PATH={config_dict['r_env']}:$PATH\n")
    sample_shell_fhand.write(f'export ichorcna_exdata={config_dict['ichorcna_exdata']}\n')
    sample_shell_fhand.write(f'readCounter --window 1000000 --quality 20 --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"  {sub_dir_dict["results"]}/{sample_id}.dedup.bam > {sub_dir_dict["CNA"]}/{sample_id}.wig\n')
    sample_shell_fhand.write(f'runIchorCNA.R --id {sample_id} --WIG {sub_dir_dict["CNA"]}/{sample_id}.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 3 --gcWig $ichorcna_exdata/gc_hg19_1000kb.wig --mapWig $ichorcna_exdata/map_hg19_1000kb.wig --centromere $ichorcna_exdata/GRCh37.p13_centromere_UCSC-gapTable.txt --normalPanel $ichorcna_exdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" --genomeStyle UCSC --scStates "c()" --txnE 0.9999 --txnStrength 10000 --outDir {sub_dir_dict["CNA"]} &/\n')

    # bam2bed
    sample_shell_fhand.write("# Step 5 Length ratio analysis and get features\n")
    # sample_shell_fhand.write(f'rm {sub_dir_dict["results"]}/{sample_id}.dedup.filtered.bed.gz\n')
    sample_shell_fhand.write(f'{config_dict["bam2bed"]} {sub_dir_dict["results"]}/{sample_id}.dedup.bam {sub_dir_dict["results"]}/{sample_id}.dedup.bed {config_dict["hg19"]}\n')
    sample_shell_fhand.write(f'{config_dict["bedtools"]} intersect -v -a {sub_dir_dict["results"]}/{sample_id}.dedup.bed -b {config_dict["hg19_blacklist_bed"]} > {sub_dir_dict["results"]}/{sample_id}.dedup.filtered.bed\n')
    #for gz bed file
    # Lenght ratio
    sample_shell_fhand.write(f'{config_dict["length_ratio_R"]} --bedfile {sub_dir_dict["results"]}/{sample_id}.dedup.filtered.bed --outdir {sub_dir_dict["length_ratio"]} --id {sample_id}\n')
    # sample_shell_fhand.write(f'{config_dict["length_ratio_R"]} --bedfile {sub_dir_dict["results"]}/{sample_id}.dedup.filtered.bed.gz --outdir {sub_dir_dict["length_ratio"]} --id {sample_id}\n')
    # get feautres from bed file
    # sample_shell_fhand.write(f'{config_dict["features_from_bed"]} {sub_dir_dict["results"]}/{sample_id}.dedup.filtered.bed.gz {sub_dir_dict["frag_features"]}\n')
    # sample_shell_fhand.write(f'{config_dict["features_from_bed"]} {sub_dir_dict["results"]}/{sample_id}.dedup.filtered.bed {sub_dir_dict["frag_features"]} {config_dict["threads"]}\n')
    sample_shell_fhand.write(f'{config_dict["features_from_bed"]} {sub_dir_dict["results"]}/{sample_id}.dedup.filtered.bed {sub_dir_dict["frag_features"]} {config_dict["hg19"]}\n')
    # clean up
    sample_shell_fhand.write("# Step 6 Clean Up\n")
    sample_shell_fhand.write("wait\n")
    sample_shell_fhand.write(f'gzip {sub_dir_dict["results"]}/{sample_id}.dedup.filtered.bed\n')
    sample_shell_fhand.write(f"rm -f {sub_dir_dict['results']}/{sample_id}.sort.bam {sub_dir_dict['results']}/{fq1_base} {sub_dir_dict['results']}/{fq2_base} {sub_dir_dict['results']}/{sample_id}.dedup.bed\n")


print(f'Total samples: {sample_counts}')
all_run_fhand.write("wait\n") # wait for all samples to finish
config_fhand.close()
all_run_fhand.close()