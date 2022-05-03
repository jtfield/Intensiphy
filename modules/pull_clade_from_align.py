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


def parse_args():
    parser = argparse.ArgumentParser(prog='Build Alignment', \
        description='Run this program and specify the output directory of an Intensiphy run to build an alignment from the sequence database.')
    # parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    parser.add_argument('--align_file', help='Alignment file.')
    parser.add_argument('--csv_file', help='CSV file.')
    parser.add_argument('--output_file', default='filtered_alignment.fas', help='Alignment file.')
    return parser.parse_args()

def main():
    args = parse_args()

    output = []

    align_file = open(args.align_file, 'r').read()

    csv_file = pd.read_csv(args.csv_file)

    split_align = align_file.split('>')

    clade_sras = list(csv_file['Run'])

    for chunk in split_align:

        if len(chunk) > 0:

            split_name_and_seq = chunk.split('\n', 1)

            name = split_name_and_seq[0]
            seq = split_name_and_seq[1]

            # print(name)

            if name in clade_sras:
                output.append(chunk)

    output_file = open(args.output_file, 'w')

    for seq in output:
        output_file.write('>' + seq)

    output_file.close()





if __name__ == '__main__':
    main()
