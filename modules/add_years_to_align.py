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
    parser = argparse.ArgumentParser(prog='add dates to align', \
        description='Run this program to reformat a sequence names in an alignment file and add the relevant date found in a csv.')
    parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    parser.add_argument('--align_file', default=False, help='Multiple sequence alignment file.')
    parser.add_argument('--date_column_name', default='Collection date', help='name of date column to use.')
    parser.add_argument('--output_file', default='date_added.aln', help='Alignment file with dates added.')
    return parser.parse_args()

def main():
    args = parse_args()

    output = []

    input_csv = pd.read_csv(args.metadata_csv)

    read_align = open(args.align_file, 'r').read()

    split_align = read_align.split('>')

    year_regex = '\d\d\d\d'

    compile_year_regex = re.compile(year_regex)

    for chunk in split_align:

        if len(chunk) > 0:

            split_name_and_seq = chunk.split('\n', 1)

            name = split_name_and_seq[0]
            seq = split_name_and_seq[1]

            try:
                found_row = input_csv.loc[input_csv['Run'] == name].fillna('0')
                if name == 'ERR2525602':
                    print(found_row['Collection date'])
                    print(found_row['SRA release date'])
                date = found_row[args.date_column_name].values[0]
                backup_date = found_row['SRA release date'].values[0]
                # print(date)

                find_date = re.match(compile_year_regex, date)
                name_and_date = ''
                if find_date:

                    find_date = find_date[0]

                    if find_date != '0':
                        name_and_date = name + '_' + find_date
                        # print(name_and_date)
                        # output.append('>' + name_and_date + '\n' + seq)

                else:
                    find_backup_date = re.match(compile_year_regex, backup_date)

                    if find_backup_date:

                        find_backup_date = find_backup_date[0]
                        name_and_date = name + '_' + find_backup_date

                    # name_and_date = name

                if len(name_and_date) == 0:
                    print(name)
                    print(name_and_date)

                output.append('>' + name_and_date + '\n' + seq)

            except ValueError:
                print("not found ", name)


    output_file = open(args.output_file, 'w')
    for item in output:
        output_file.write(item)

    output_file.close()







if __name__ == '__main__':
    main()
