#!/usr/bin/env python3
# usage: get_filtered_data.py output_dir
# example dir: sample_info
# example: get_filtered_data.py 

import os
import sys
import pandas as pd
import utils


data_dir = os.path.abspath(sys.argv[1])
filter_features_dir = os.path.abspath(sys.argv[2])
output_dir = os.path.abspath(sys.argv[3])
threads = 10
SELECTED_METHODS = ['glmnet', 'rfe_fifty_predictors']


input_feautre_pairs ={'end_motif_6.csv':'end_motif_6/feature_selection.csv', 'end_motif_4.csv':'end_motif_4/feature_selection.csv', 'breakpoint_motif.csv':'breakpoint_motif/feature_selection.csv', 'read_length.csv':'read_length/feature_selection.csv'}

for method in SELECTED_METHODS:
    _output_dir = os.path.join(output_dir, method)
    os.makedirs(_output_dir, exist_ok=True)
    for _data_file, _feature_file in input_feautre_pairs.items():
        _file = os.path.join(data_dir, _data_file)
        _feature = os.path.join(filter_features_dir, _feature_file)
        if not utils.check_file(_file) or not utils.check_file(_feature):
            continue
        feature_df = pd.read_csv(_feature, sep=',')
        feature_df.set_index('method', inplace=True)
        selected_features = feature_df.loc[method, 'features'].split(';')
        _data = pd.read_csv(_file, sep=',')
        _out_file = os.path.join(_output_dir, _data_file)
        _data[['sample_id'] + selected_features  + ['type']].to_csv(_out_file, sep=',', index=False)
        print('{} is done'.format(_out_file))