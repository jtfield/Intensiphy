#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file', help='alignment file with samples you want to pull from a csv metadata file.')
    parser.add_argument('--csv_file', help='large metadata file that you want data for a specific set of taxa from.')
    parser.add_argument('--output_csv', default='filtered_dataset.csv', help='output csv file name.')
    return parser.parse_args()


def main():
    args = parse_args()

    read_align = open(args.align_file, 'r').read()

    read_csv = pd.read_csv(args.csv_file)

    seq_names = get_names(read_align)

    output_table = filter_df(read_csv, seq_names)
    # print(output_table)

    output_table.to_csv(args.output_csv)


def get_names(align):
    """
    Pull the sequence names from the alignment file.
    """
    names_list = []

    split_align = align.split('>')

    for chunk in split_align:
        if len(chunk) > 1:
            name_and_seq = chunk.split('\n', 1)

            names_list.append(name_and_seq[0])

    print(len(names_list))

    return names_list


def filter_df(csv, sample_names):
    """
    Take in a list of names and pull the samples from the metadata file, making a new csv.
    """
    #pull column names of big metadata df
    df_columns = csv.columns

    #initialize dataframe with columns
    clade_df = pd.DataFrame(columns=df_columns)

    #Find sample ids in the big dataframe
    for id in sample_names:
        try:
            # print(id)
            # print(big_df_drop_na_runs['Run'])

            found_id = csv.loc[csv['Run'].str.contains(id, na=False)]
            # print(found_id)

            clade_df = clade_df.append(found_id)

        except IndexError:
            print("Missing tip in chosen clade: ", id)

    return clade_df









if __name__ == '__main__':
    main()
