#!/usr/bin/env python3

import os
import argparse
import pathlib
import shutil
from numpy.lib.shape_base import split
import pandas as pd
import subprocess
import datetime
import dateutil
# import sys
# sys.path.append('./modules')
# import modules.seq_similarity_assessment
# import modules.alignment_splitter
from modules.seq_similarity_assessment import *
from modules.alignment_splitter import split_alignment
from modules.fetch_and_align import *
from modules.tree_assess import *
# import tests.accessions_tests.accession_download_test
# import tests.assembly_tests.gon_phy_test

def parse_args():
    parser = argparse.ArgumentParser(prog='Table Filter', \
        description='Filter a table by only keeping .')
    parser.add_argument('--csv_input', default=False, help='CSV file with table you want to filter.')
    parser.add_argument('--filter_column_label', help='Name of the column you wish to filter on')
    parser.add_argument('--output_file', default='filtered_table.csv', help='Name of file with new, filtered table.')
    return parser.parse_args()

def main():

    args = parse_args()

    # Read csv input file
    input = pd.read_csv(args.csv_input)

    filtered_table = filter_table(input, args.filter_column_label)

    filtered_table.to_csv(args.output_file)




def filter_table(table_, filter_column_):
    """
    Filter the CSV file by excluding all rows that have and empty value for the chosen column
    """

    new_df = table_.dropna(subset=[filter_column_])

    # print(new_df)

    return new_df






if __name__ == '__main__':
    main()
