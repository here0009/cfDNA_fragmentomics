#!/usr/bin/bash
input_dir=$1
batch_size=$2
counts=0
for folder in `ls $input_dir`; do
    echo $folder
    if [ -f $input_dir/$folder/results/$folder.dedup.bam ]; then
        echo "input bam file is : $folder.dedup.bam"
        mkdir -p $input_dir/$folder/length_ratio
        echo "output dir is : $input_dir/$folder/length_ratio"
        let "counts+=1"
        src/filter_bam2bed.py $input_dir/$folder/results/$folder.dedup.bam $input_dir/$folder/results/filtered_bam.bed 
        bedtools intersect -v -a $input_dir/$folder/results/filtered_bam.bed -b data/hg19_fragmentomics/hg19-blacklist.v2.bed > $input_dir/$folder/results/no_blacklist.bed
        rm -f $input_dir/$folder/results/filtered_bam.bed
        src/length_ratio.R --bedfile $input_dir/$folder/results/no_blacklist.bed --outdir $input_dir/$folder/length_ratio --id $folder
        if [ $counts -eq $batch_size ]; then
            echo "batch is done"
            let "counts=0"
            wait
        fi
    fi
done
