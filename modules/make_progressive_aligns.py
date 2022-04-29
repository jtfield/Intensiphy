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
    parser.add_argument('--ip_dir', help='ip_dir.')
    parser.add_argument('--split_number', default=4, help='number of alignments you want to split the new sequences into.')
    parser.add_argument('--output_dir', default='split_placement_alignments_dir', help='output location to make output folders.')
    return parser.parse_args()

def main():
    args = parse_args()

    align_file = open(args.align_file, 'r').read()


    starting_align_names_list = read_original_align_names(align_file)

    # print(starting_align_names_list)

    splits = make_splits(args.ip_dir, starting_align_names_list, args.split_number)

    # print(splits)

    make_aligns_and_dirs(splits, args.output_dir, args.ip_dir)


def read_original_align_names(alignment):
    names_list = []

    split_align_file = alignment.split('>')

    for chunk in split_align_file:

        if len(chunk) > 0:

            split_name_and_seq = chunk.split('\n', 1)
            name = split_name_and_seq[0]
            seq = split_name_and_seq[1]

            names_list.append(name)

    return names_list


def make_splits(ip_dir_, list_of_names_, splits):


    list_of_seqs = os.listdir(ip_dir_ + '/sequence_storage')
    # print(list_of_seqs)

    for name in list_of_names_:
        if name in list_of_seqs:
            list_of_seqs.remove(name)


    seq_per_split = int(len(list_of_seqs) / int(splits))

    batches = [list_of_seqs[x:x+seq_per_split] for x in range(0, len(list_of_seqs), seq_per_split)]

    # print(batches)

    for chunk in batches:
        for name in list_of_names_:
            chunk.append(name)
        # print(chunk)

    return batches



def make_aligns_and_dirs(batch_list, outdir, ip_dir):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    num_of_aligns = len(batch_list)

    dir_name_template = 'alignment_dir_'
    dir_count = 0
    list_of_dirs = []


    # print(num_of_alins)

    os.chdir(outdir)

    for num in range(0, num_of_aligns):
        dir_name = dir_name_template + str(num)
        list_of_dirs.append(dir_name)

        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

    for num, seq_names in enumerate(batch_list):

        print(num)
        print(seq_names)

        make_align(seq_names, ip_dir + '/sequence_storage', outdir + '/' + list_of_dirs[num])


def make_align(list_of_seqs, seq_storage, outdir):

    alignment = []

    for seq in list_of_seqs:

        print(seq)

        read_seq_file = open(seq_storage + '/' + seq + '/' + seq + '_.fas', 'r').read()

        split_align_file = read_seq_file.split('>')

        for chunk in split_align_file:

            if len(chunk) > 0:

                split_name_and_seq = chunk.split('\n', 1)
                name = split_name_and_seq[0]
                seq = split_name_and_seq[1]

                alignment.append(chunk)

    # print(chunk)








if __name__ == '__main__':
    main()
