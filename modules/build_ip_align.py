#!/usr/bin/env python3

import os
import argparse
import pathlib
import numpy as np
import pandas as pd
import subprocess
import datetime

def parse_args():
    parser = argparse.ArgumentParser(prog='Tree to clades', \
        description='Split a large tree into more manageable monophyletic clades.')
    parser.add_argument('--ip_dir', default=False, help='Directory output by an IP run.')
    return parser.parse_args()

def main():
    args = parse_args()

    construct_align_(args.ip_dir)


def construct_align_(outdir):
    """Builds an alignment using all present sequences. \
    Uses starting tree to place new sequences in the tree (RAxML EPA)"""

    # Find current date and make folder with the date in the name
    # for record keeping purposes
    now = datetime.datetime.now()
    seq_storage_path = outdir + '/sequence_storage'
    # starting_tree_path = outdir + '/intermediate_files/RAxML_bestTree.starting_tree.tre'
    taxa_list = os.listdir(seq_storage_path)
    suffix = '_.fas'

    align_dir = 'ip_alignment_' + now.strftime('%Y-%m-%d-%H-%M-%S')
    align_dir_full_path = outdir + '/' + align_dir

    output_alignment_path = align_dir_full_path + '/extended.aln'

    os.mkdir(phylo_dir_full_path)

    # cat_command_start = ['cat']
    output = open(output_alignment_path, 'a')

    for dir in taxa_list:
        taxon_dir_path = seq_storage_path + '/' + dir
        taxon_seq_file_path = taxon_dir_path + '/' + dir + suffix
        # print(taxon_seq_file_path)

        # Add file path to cat command list
        # cat_command_start.append(taxon_seq_file_path)
        seq_file = open(taxon_seq_file_path, 'r').read()
        # print(type(seq_file))
        output.write(seq_file)
        output.write('\n')


    output.close()

    print("Updated alignment construction complete.")


if __name__ == '__main__':
    main()
