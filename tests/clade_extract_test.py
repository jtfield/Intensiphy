#! /usr/bin/python3

import os
import argparse
import subprocess
import time
# import pytest
import pandas as pd
import sys
import filecmp

# sys.path.append('../modules')
# from alignment_splitter import split_alignment

def parse_args():
    parser = argparse.ArgumentParser(prog='clade extract test', \
        description='This program tests the clade extraction software packaged with Intensiphy. \
        EXAMPLE COMMAND: clade_extract_test.py \
        EXAMPLE PYTEST COMMAND: pytest clade_extract_test.py')
    # parser.add_argument('--ep_path', help='Absolute path to your Extensiphy directory.')
    return parser.parse_args()

def main():
    args = parse_args()
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]
    path = '/' + absolute_path.strip('/tests')
    print(path)
    print(absolute_path)

    tree_file = absolute_path + '/combo.tre'
    metadata_csv = absolute_path + '/combo.csv'
    output_list = absolute_path + '/test_clade.txt'
    lower_bound = '4'
    upper_bound = '7'

    subprocess.run([path + '/clade_filter.py', '--tree_file', tree_file, '--metadata_csv', metadata_csv, '--output_file', output_list, '--lower_bound', lower_bound, '--upper_bound', upper_bound])
    # print(absolute_path + '/clade_filter.py', '--tree_file', tree_file, '--metadata_csv', metadata_csv, '--output_file', output, '--lower_bound', lower_bound, '--upper_bound', upper_bound)

    subprocess.run([path + '/tree_trim.py', '--tree_file', tree_file, '--taxa_list', output_list])

def test_clade_select():
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]
    path = '/' + absolute_path.strip('/tests')
    output_tree = absolute_path + '/trimmed_tree.tre'
    output_clade = absolute_path + '/test_clade.txt'
    example_clade = absolute_path + '/example_clade.txt'
    example_tree = absolute_path + '/example_trimmed_tree.tre'

    # if 







# def test_output():
#     split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
#     absolute_path = split_path_and_name[0]
#     path = '/' + absolute_path.strip('/tests')
#     for file in path:
#         if file.endswith('_.fas'):
#             print(file)

if __name__ == '__main__':
    main()
