#!/usr/bin/env python3
# usage: get_prediction_shell.py output_dir
# example dir: results/selected_features/glmnet
# example: get_prediction_shell.py r_scripts/caret_gbm_predict.R . . 10

import os
import sys
import pandas as pd
import utils


predict_method = os.path.abspath(sys.argv[1])
input_dir = os.path.abspath(sys.argv[2])
output_dir = os.path.abspath(sys.argv[3])
threads = int(sys.argv[4])
output_shell = os.path.join(output_dir, 'prediction.sh')
output_shell_fh = open(output_shell, 'w')
# predict_method = 'r_scripts/caret_feature_selection.R'
sample_info_file = 'sample_info/sample_info_csv/220809_total_sample_info.csv'
threads = 10

input_files = ['end_motif_6.csv', 'end_motif_4.csv', 'breakpoint_motif.csv', 'read_length.csv']
# input_dir = 'results/end_motif_2'
for _file in input_files:
    _input_file = os.path.join(input_dir, _file)
    if not utils.check_file(_input_file):
        continue
    _tag = _file.split('.')[0]
    _output_dir = os.path.join(output_dir, _tag)
    _log_file = os.path.join(_output_dir, f'{_tag}_predict.log')
    os.makedirs(_output_dir, exist_ok=True)
    output_shell_fh.write(f'{predict_method} {_input_file} {_output_dir} {sample_info_file} {threads} > {_log_file} 2>&1\n')
output_shell_fh.close()