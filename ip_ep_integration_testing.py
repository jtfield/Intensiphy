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

def parse_args():
    parser = argparse.ArgumentParser(prog='Intensiphy', \
        description='Automate downloading of high-throughput sequence data and updating of alignments using Extensiphy.')
    parser.add_argument('--align_file', default=False, help='input alignment option. \
        If an alignment is added using this option, samples will be added to this alignment.')
    parser.add_argument('--read_dir', help='Directory of reads')
    parser.add_argument('--ep_out_dir', help='Absolute path and folder name to create for outputs')
    parser.add_argument('--ep_path', help='Path to extensiphy folder.')
    return parser.parse_args()

def main():
    args = parse_args()

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    ep_path = args.ep_path
    alignment = args.align_file
    read_dir = args.read_dir
    ep_outdir = args.ep_out_dir

    print(ep_path + "/extensiphy.sh", "-a", alignment, "-d", read_dir, "-1", "_1.fastq", "-2", "_2.fastq", "-o", ep_outdir)
    ep_process = subprocess.Popen([ep_path + "/extensiphy.sh", "-a", alignment, "-d", read_dir, "-1", "_1.fastq", "-2", "_2.fastq", "-o", ep_outdir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(ep_process.communicate()[0].decode('utf-8'))


if __name__ == '__main__':
    main()
