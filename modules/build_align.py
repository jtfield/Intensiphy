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
from modules.tree_assess import *

def parse_args():
    parser = argparse.ArgumentParser(prog='Build Alignment', \
        description='Run this program and specify the output directory of an Intensiphy run to build an alignment from the sequence database.')
    # parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    parser.add_argument('--ip_dir', default='ip_output', help='Directory containing the outputs of an Intensiphy run.')
    return parser.parse_args()

def main():
    args = parse_args()


    # This whole progam is jut to give the user a "make an alignment whenever i want" button
    # So if the run stops for any reason, they can use what they have in the sequence database
    construct_align_and_place(args.ip_dir)




if __name__ == '__main__':
    main()
