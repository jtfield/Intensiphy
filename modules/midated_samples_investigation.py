#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
# import subprocess
# import dendropy
# import pathlib
# path_root = pathlib.Path(__file__).parents[0]
# sys.path.append(str(path_root) + '/modules')
# from modules.analysis_functions import *


def parse_args():
    parser = argparse.ArgumentParser(prog='misdated tree investigation', \
        description='After running EP on multiple reference sequences, analyze the results.')
    parser.add_argument('--tree_file', default=False, help='input phylogeny option.')
    parser.add_argument('--csv_file', help='table of isolate metadata.')
    parser.add_argument('--output_dir', default='sorted_snp_clusters.csv', help='table of sorted by SNP cluster values.')
    # parser.add_argument('--taxa_list', default=False, help='File that contains the names of each taxon to keep.')
    # parser.add_argument('--output_tree_file', default='trimmed_tree.tre', help='output tree file. (DEFAULT: trimmed_tree.tre)')

    return parser.parse_args()

def main():
    args = parse_args()

    df = pd.read_csv(args.csv_file)

    multi_tables = sep_snp_tables(df, args.output_dir)

    # df = df.sort_values(by=['SNP cluster'])
    #
    # # print(df['SNP cluster'])
    #
    # df.to_csv(args.output_csv)


def sep_snp_tables(csv_, out_dir):
    """
    Takes in a table of isolate metadata.
    Separates the table based on SNP cluster values.
    Returns those tables.
    """

    uniq_snp_clusters = csv_['SNP cluster'].dropna().unique()

    print(uniq_snp_clusters)

    num_tables = len(uniq_snp_clusters)

    df_columns = csv_.columns
    print(df_columns)

    df_container = []

    for cluster in uniq_snp_clusters:
        # df = pd.DataFrame(columns=df_columns)

        new_df = csv_[csv_['SNP cluster'] == cluster]

        # print(new_df)

        file_name = out_dir + '/' + cluster + '.csv'

        new_df.to_csv(file_name)








if __name__ == '__main__':
    main()
