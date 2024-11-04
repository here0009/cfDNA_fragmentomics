# download rawfq
tag=$1
sample_info=$2
rawfq_dir=$3
working_dir=$4

mkdir -p $rawfq_dir/$tag
mkdir -p $working_dir/$tag

/common_tools/download_fq.py $sample_info $rawfq_dir/$tag
cd $rawfq_dir/$tag
nohup bash download_fq.sh
# match local sample with sample info

/haima_pipeline_tools/match_sample.py $sample_info $rawfq_dir/$tag $rawfq_dir/$tag/fq_matched.tsv
cp download_fq.sh $working_dir/$tag
cp fq_matched.tsv $working_dir/$tag
cp $sample_info $working_dir/$tag

/haima_pipeline_tools/merge_samples.py $rawfq_dir/$tag/fq_matched.tsv $rawfq_dir/$tag/fq_matched_merged.tsv $rawfq_dir/$tag
cp $rawfq_dir/$tag/fq_matched_merged.tsv $working_dir/$tag
bash $working_dir/$tag/merge_samples.sh
# get running shell
fragmentomics_get_shell.py $rawfq_dir/$tag/fq_matched_merged.tsv config.json $working_dir/$tag
# run shell
cd $working_dir/$tag
bash all_run_script.sh
# batch QC
wait
batch_qc_info.py $sample_info $working_dir/$tag $working_dir/$tag/total_qc.tsv