#!/usr/bin/env python3

import os
import argparse
import numpy
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(prog='Table Filter', \
        description='Filter a table by only keeping .')
    parser.add_argument('--samples_csv', default=False, help='CSV file with the samples you want to search for.')
    parser.add_argument('--data_csv', help='CSV file with the rows you want for the selected samples')
    parser.add_argument('--output_csv', default='filtered_table.csv', help='Name of file with new, filtered table.')
    return parser.parse_args()

def main():

    args = parse_args()

    # Read csv input file
    samples_input = pd.read_csv(args.samples_csv)

    big_csv = pd.read_csv(args.data_csv)

    sample_names = samples_input['Run'].tolist()

    big_csv_columns = big_csv.columns

    new_df = pd.DataFrame(columns=big_csv_columns)
    # print(new_df)

    for name in sample_names:
        # print(name)
        # print(type(name))
        try:
            found_row = big_csv[big_csv["Run"] == name]
            # print(found_row)
            new_df = new_df.append(found_row)
        except ValueError:
            print("not found ", name)

    new_df.to_csv(args.output_csv)







if __name__ == '__main__':
    main()
