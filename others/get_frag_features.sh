#!/usr/bin/bash
input_dir=$1
batch_size=$2
counts=0
for folder in `ls $input_dir`; do
    echo $folder
    if [ -f $input_dir/$folder/results/$folder.dedup.bam ]; then
        echo "input bam file is : $folder.dedup.bam"
        mkdir -p $input_dir/$folder/frag_features
        echo "output dir is : $input_dir/$folder/frag_features"
        let "counts+=1"
        src/fragmentomics_features.py $input_dir/$folder/results/$folder.dedup.bam $input_dir/$folder/frag_features/ &
        if [ $counts -eq $batch_size ]; then
            echo "batch is done"
            let "counts=0"
            wait
        fi
    fi
done
