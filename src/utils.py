#!/usr/bin/env python3
import os
import sys
import pandas as pd
from collections import Counter

complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

def read_csv(table_file):
    _sep = ',' if table_file.endswith(".csv") else '\t'
    table = pd.read_csv(table_file, sep=_sep)
    return table

def check_file(_file):
    if _file is None or not os.path.isfile(_file):
        print(f"[ERROR] File {_file} does not exist!")
        return False
    return True

def check_dir(_dir):
    if not os.path.isdir(_dir):
        print(f"[ERROR] Directory {_dir} does not exist!")
        return False
    return True

def concate_table(dir_lst, table_name):
    _sep = ',' if table_name.endswith(".csv") else '\t'
    table_lst = [pd.read_csv(f"{dir}/{table_name}", sep=_sep) for dir in dir_lst]
    total_table = pd.concat(table_lst, ignore_index=True)
    return total_table

def reverse_seq(seq):
    return seq[::-1].upper()

def reverse_complement(seq):
    return ''.join([complement_dict.get(base, 'N') for base in reverse_seq(seq)])

def gc_content(seq):
    counts = Counter(seq.upper())
    return (counts['G'] + counts['C']) / (counts['G'] + counts['C'] + counts['A'] + counts['T'])