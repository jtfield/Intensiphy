#! /usr/bin/python3

import os
import argparse
import pathlib
import shutil
# from numpy.lib.shape_base import split
# import pandas as pd
import subprocess
# import datetime
# import dateutil
# import tests.accessions_tests.accession_download_test
# import tests.assembly_tests.gon_phy_test
# import sys
# sys.path.append('./modules')
# import modules.seq_similarity_assessment.py
# import modules.alignment_splitter.py
from modules.seq_similarity_assessment import *
from modules.alignment_splitter import split_alignment

def parse_args():
    parser = argparse.ArgumentParser(prog='Intensiphy', \
        description='Automate downloading of high-throughput sequence data and updating of alignments using Extensiphy.')
    parser.add_argument('--align_file', default=False, help='input alignment option. \
    #     If an alignment is added using this option, samples will be added to this alignment.')
    # parser.add_argument('--read_dir', help='Directory of reads')
    parser.add_argument('--ep_out_dir', help='Absolute path and folder name to create for outputs')
    # parser.add_argument('--ep_path', help='Path to extensiphy folder.')
    # parser.add_argument('--seq_1_file', help='Path to file one.')
    # parser.add_argument('--seq_2_file', help='Path to file two.')
    return parser.parse_args()

def main():
    args = parse_args()

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    # ep_path = args.ep_path
    alignment = args.align_file
    # read_dir = args.read_dir
    ep_outdir = args.ep_out_dir
    #####################################################################
    # file_1 = args.seq_1_file
    # file_2 = args.seq_2_file
    #
    # read_file_1 = open(file_1, 'r').read()
    # file_1_split = read_file_1.split('\n', 1)
    # seq_1 = file_1_split[1].replace('\n','')
    #
    # read_file_2 = open(file_2, 'r').read()
    # file_2_split = read_file_2.split('\n', 1)
    # seq_2 = file_2_split[1].replace('\n','')
    #
    #
    # check_sequence_similarities(seq_1, seq_2)
    ##########################################################################

    # split_alignment(alignment, ep_outdir)

    build_or_update_df(ep_outdir, False, '_.fas')

    check_sims_and_remove(ep_outdir, .99)

    # seq_compare(ep_outdir, '_.fas', None)

if __name__ == '__main__':
    main()
