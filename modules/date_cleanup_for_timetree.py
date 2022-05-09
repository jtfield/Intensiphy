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
    parser = argparse.ArgumentParser(prog='Reformat csv', \
        description='Run this program to reformat a csv downloaded from PathoDB to contain info in the format required for TimeTree.')
    parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    parser.add_argument('--output_file', default='timetree_data.csv', help='Alignment file.')
    return parser.parse_args()

def main():
    args = parse_args()

    input_csv = pd.read_csv(args.metadata_csv).fillna('0')

    columns = input_csv.columns

    new_df = pd.DataFrame(columns=columns)

    year_regex = '\d\d\d\d'

    compile_year_regex = re.compile(year_regex)

    for idx, row in input_csv.iterrows():
        if row['Collection date'] != '0':

            pull_found_year = re.match(compile_year_regex, row['Collection date'])

            if pull_found_year:
                input_csv.loc[idx, 'Collection date'] = pull_found_year[0]

                new_df = new_df.append(row)
            # print(row['Run'])

        elif row['Collection date'] == '0':
            new_date = row['SRA release date']

            pull_year = re.match(compile_year_regex, new_date)

            if pull_year:

                input_csv.loc[idx, 'Collection date'] = pull_year[0]
                print(row['Collection date'])
                new_df = new_df.append(row)

    #             new_df.append(replacement_date_row)
    #
    # print(new_df)

    new_df.to_csv(args.output_file)


    # timetree_df = input_csv[['Run', 'Collection date']].copy()
    #
    # timetree_df.rename(columns= {'Run':'accession'}, inplace=True)
    # timetree_df.rename(columns= {'Collection date':'date'}, inplace=True)
    #
    # timetree_df.to_csv(args.output_file)








if __name__ == '__main__':
    main()
