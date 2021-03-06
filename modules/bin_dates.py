#!/usr/bin/env python3

import os
import argparse
import pathlib
import shutil
import subprocess
import dendropy
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(prog='Bin Dates', \
        description='Parse dates in a csv file from PathoDB and break up the samples by time periods of 6 months.')
    parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    parser.add_argument('--output_dir', default='dir_of_binned_dates', help='output directory containing output separate csv files binned by their collection dates.')
    return parser.parse_args()

def main():
    args = parse_args()

    #Read in metadata file
    metadata_csv = pd.read_csv(args.metadata_csv)

    # print(metadata_csv['Collection date'])

    col_dates = metadata_csv['Collection date']

    # print(col_dates)

    str_col_dates = col_dates.astype('str')

    # str_col_dates.replace({r'(\d\d\d\d)(-\d\d)(-\d\d)' : '\\1\\2'}, regex=True)

    metadata_csv['Collection date'].replace({r'(\d\d\d\d)(-\d\d)(-\d\d)' : '\\1\\2'}, regex=True, inplace=True)

    # print(metadata_csv['Collection date'])

    year_only_df = metadata_csv.copy()

    # year_only_df_dates = year_only_df['Collection date'].astype(str)

    yead_only_df = year_only_df['Collection date'].astype(str)

    # year_only_df_dates.replace({r'(\d\d\d\d)(-\d\d)' : '\\1'}, regex=True, inplace=True)

    year_only_df.replace({r'(\d\d\d\d)(-\d\d)' : '\\1'}, regex=True, inplace=True)

    int_series = pd.to_numeric(year_only_df['Collection date'])

    print(int_series.value_counts())

    max_year = int_series.max()
    min_year = int_series.min()

    print(min_year)
    print(max_year)

    # year_max = year_only_df['Collection date'].loc[year_only_df['Collection date'].idmax()]
    #
    # year_min = year_only_df['Collection date'].loc[year_only_df['Collection date'].idmin()]
    #
    # print(year_min)
    # print(year_max)

    # for row in int_series:
    #     print(row)
    #     print(type(row))

    # year_only_df['Collection date'].replace({r'(\d\d\d\d)(-\d\d)(-\d\d)' : '\\1'}, regex=True, inplace=True)
    #
    # print(year_only_df['Collection date'])



if __name__ == '__main__':
    main()
