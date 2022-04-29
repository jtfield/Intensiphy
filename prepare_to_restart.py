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
from modules.fetch_and_align import *

def parse_args():
    parser = argparse.ArgumentParser(prog='Prepare to Restart', \
        description='Run this program and specify the output directory of an Intensiphy run to clean up files so you can continue youre Intensiphy run.')
    # parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    parser.add_argument('--ip_dir', default='ip_output', help='Directory containing the outputs of an Intensiphy run.')
    return parser.parse_args()

def main():
    args = parse_args()

    absolute_dir_path = os.path.abspath(args.ip_dir)

    clean_downloads(absolute_dir_path)

    clean_ep_run(absolute_dir_path)


def clean_downloads(outdir):
    """Delete files and folders left behind by fasterq-dump if a download is interrupted"""
    read_dir = outdir + '/read_files'

    # List files in read dir
    files_list = os.listdir(read_dir)

    # Check if file_dir is empty
    # if empty, no need to remove anything
    if len(files_list) == 0:
        print("Read directory is empty and ready for new reads.")

    else:

        # Loop over files and check if any are directories
        for file in files_list:
            file_path = read_dir + '/' + file
            check_dir = os.path.isdir(file_path)
            check_file = os.path.isfile(file_path)

            # If file is a directory, try to remove
            # throw an exception if this doesnt work
            if check_dir:
                try:
                    print("Removing incomplete data download.")
                    shutil.rmtree(file_path)
                except:
                    print("Could not delete dir :", file_path)


            if check_file:
                try:
                    print("Removing data files.")
                    os.remove(file_path)
                except:
                    print("Could not delete file :", file_path)


def clean_ep_run(outdir):

    ep_dir = outdir + "/intermediate_files/ep_output"

    check_dir = os.path.isdir(ep_dir)

    if check_dir:
        try:
            print("Removing old Extensiphy run output folder.")
            shutil.rmtree(ep_dir)
        except:
            print("Could not delete dir :", ep_dir)

    else:
        print("No leftover Extensiphy directory was found.")




if __name__ == '__main__':
    main()
