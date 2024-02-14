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

    timetree_df = input_csv[['Run', 'Collection date']].copy()

    timetree_df.rename(columns= {'Run':'accession'}, inplace=True)
    timetree_df.rename(columns= {'Collection date':'date'}, inplace=True)

    timetree_df.to_csv(args.output_file)








if __name__ == '__main__':
    main()
