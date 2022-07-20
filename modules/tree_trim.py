#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import dendropy
import pathlib
import pandas as pd
# path_root = pathlib.Path(__file__).parents[0]
# sys.path.append(str(path_root) + '/modules')
# from modules.analysis_functions import *


def parse_args():
    parser = argparse.ArgumentParser(prog='trim around subtree', \
        description='Provide a list of taxa to keep. Prune the rest of the taxa.')
    parser.add_argument('--tree_file', default=False, help='input phylogeny option.')
    parser.add_argument('--taxa_set', default=False, help='File that contains the names of each taxon to keep.')
    parser.add_argument('--taxa_file_type', default='csv', help='File type of list of taxa. Either csv or a staight list .txt File.')
    parser.add_argument('--output_tree_file', default='trimmed_tree.tre', help='output tree file. (DEFAULT: trimmed_tree.tre)')

    return parser.parse_args()


def main():
    args = parse_args()

    if args.taxa_file_type != 'csv':

        if_list(args.tree_file, args.taxa_set, args.output_tree_file)

    elif args.taxa_file_type == 'csv':

        if_csv(args.tree_file, args.taxa_set, args.output_tree_file)


def if_csv(tree_file, taxa_list, outfile):

    names_to_keep = []

    # # establish taxon namespace
    tns = dendropy.TaxonNamespace()
    #
    in_tree = dendropy.Tree.get(path=tree_file, schema='newick', taxon_namespace=tns, preserve_underscores=True)
    #
    # print(in_tree)

    # input_taxa_list = open(taxa_list, 'r')
    input_taxa_list = pd.read_csv(taxa_list)

    input_taxa_list = input_taxa_list['Run']
    #
    # split_names = args.taxa_list.split(',')

    for name in input_taxa_list:
        # print(name)
        names_to_keep.append(name.strip())
    # print(names_to_keep)

    taxa_to_retain = set([taxon for taxon in in_tree.taxon_namespace if taxon.label in names_to_keep])

    # print(taxa_to_retain)


    filtered_tree = in_tree.extract_tree_with_taxa(taxa=taxa_to_retain)
    #
    #
    #
    print(filtered_tree.as_ascii_plot())

    filtered_tree.write(path=outfile, schema="newick")


def if_list(tree_file, taxa_list, outfile):
    names_to_keep = []

    # # establish taxon namespace
    tns = dendropy.TaxonNamespace()
    #
    in_tree = dendropy.Tree.get(path=tree_file, schema='newick', taxon_namespace=tns, preserve_underscores=True)
    #
    # print(in_tree)

    input_taxa_list = open(taxa_list, 'r')
    #
    # split_names = args.taxa_list.split(',')

    for name in input_taxa_list:
        # print(name)
        names_to_keep.append(name.strip())
    # print(names_to_keep)

    taxa_to_retain = set([taxon for taxon in in_tree.taxon_namespace if taxon.label in names_to_keep])

    # print(taxa_to_retain)


    filtered_tree = in_tree.extract_tree_with_taxa(taxa=taxa_to_retain)
    #
    #
    #
    print(filtered_tree.as_ascii_plot())

    filtered_tree.write(path=outfile, schema="newick")

if __name__ == '__main__':
    main()
