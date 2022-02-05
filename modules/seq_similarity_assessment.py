#! /usr/bin/python3

import os
import argparse
import pathlib

def check_sequence_similarities(_seq_1, _seq_2):
    """Checks similarity of two sequences"""
    print("waffle")

    # standard_nucleotides = ['A', 'C', 'G', 'T']
    # gaps = ['-']
    # degen_nucleotides = ['W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N']

    ident_nucs = 0
    ident_gaps = 0
    ident_degens = 0
    non_ident_bases = 0

    list_seq_1 = list(_seq_1)
    list_seq_2 = list(_seq_2)

    zipped_seqs = list(zip(list_seq_1, list_seq_2))

    results = map(check_nucs, zipped_seqs)

    for i in list(results):
        if i == 1:
            ident_nucs+=1
        elif i == 2:
            ident_gaps+=1
        elif i == 3:
            ident_degens+=1
        elif i == 4:
            non_ident_bases+=1

    summed_ident = ident_nucs + ident_gaps + ident_degens
    print(summed_ident)
    print(summed_ident / len(list_seq_1))
    print(summed_ident / len(list_seq_2))


def check_nucs(nuc_tuple):
    standard_nucleotides = ['A', 'C', 'G', 'T']
    gaps = ['-']
    degen_nucleotides = ['W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N']
    nuc_1 = nuc_tuple[0].upper()
    nuc_2 = nuc_tuple[1].upper()

    if nuc_1 == nuc_2:
        if nuc_1 in standard_nucleotides:
            return 1
        elif nuc_1 in gaps:
            return 2
        elif nuc_1 in degen_nucleotides:
            return 3
    else:
        return 4
