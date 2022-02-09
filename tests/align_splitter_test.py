#! /usr/bin/python3

import os
import argparse
import subprocess
import time
# import pytest
import pandas as pd
# from ...modules.alignment_splitter import split_alignment
import sys
sys.path.append('../modules')
from alignment_splitter import split_alignment

def parse_args():
    parser = argparse.ArgumentParser(prog='align_splitter test', \
        description='This program tests the alignment_splitter module of intensiphy. \
        Test currently checks that individual files were output and structured \
        as expected \
        EXAMPLE COMMAND: align_splitter_test.py \
        EXAMPLE PYTEST COMMAND: pytest align_splitter_test.py')
    # parser.add_argument('--ep_path', help='Absolute path to your Extensiphy directory.')
    return parser.parse_args()

def main():
    args = parse_args()
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]
    path = '/' + absolute_path.strip('/tests')

    split_alignment(path + '/combo.fas', path)

def test_output():
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]
    path = '/' + absolute_path.strip('/tests')
    for file in path:
        if file.endswith('_.fas'):
            print(file)
