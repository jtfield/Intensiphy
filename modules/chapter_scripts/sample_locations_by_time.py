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
import re


def parse_args():
    parser = argparse.ArgumentParser(prog='sample locations by time', \
        description='Run this program to reformat a csv and build a table tracking samples collected in a collection over time.')
    parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    parser.add_argument('--output_table_file', default='dated_location_data.csv', help='output csv file with location data split by dates.')
    return parser.parse_args()

def main():
    args = parse_args()

    metadata_csv = pd.read_csv(args.metadata_csv)

    locations = metadata_csv['Location'].unique()

    years = metadata_csv['Collection date'].unique()

    for i in range(0, len(years)):
        years[i] = int(years[i])

    years.sort()

    counts_df = pd.DataFrame(columns=locations, index=years).fillna(0)

    for idx, row in metadata_csv.iterrows():
        if row['Collection date'] != 0 or row['Collection date'] != '0':
            year = row['Collection date']
            location = row['Location']
            # print(year)
            # print(location)
            updated_value = counts_df.at[year, location] = counts_df.at[year, location] + 1


    counts_df.to_csv(args.output_table_file)












if __name__ == '__main__':
    main()
