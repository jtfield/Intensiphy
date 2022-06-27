#!/usr/bin/env python3
# program produces the best reference file from the original alignment
# currently, this program just isolates the first sequence in the file
# this handles if there are new line breaks and the sequence isn't read as a single line (some genome files)
import os
import argparse
import re
from .fasta_manipulation import *

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file')
    parser.add_argument('--out_file')
    # parser.add_argument('--ref_select')
    # parser.add_argument('-m', action='store_true', help='Turn on multi-fasta output option')
    # parser.add_argument('-r', action='store_true', help='Turn on single taxon RANDOM reference fasta output option')
    # parser.add_argument('-s', action='store_true', help='Turn on single taxon SEARCH reference fasta output option')
    return parser.parse_args()


def main():
    args = parse_args()

    split_fasta_into_seqs(args.align_file, args.out_file)




if __name__ == '__main__':
    main()
