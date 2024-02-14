#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import dendropy
import pathlib
# path_root = pathlib.Path(__file__).parents[0]
# sys.path.append(str(path_root) + '/modules')
# from modules.analysis_functions import *


def parse_args():
    parser = argparse.ArgumentParser(prog='remove taxa from tree', \
        description='Provide a list of taxa to keep. Prune the rest of the taxa.')
    parser.add_argument('--tree_file', default=False, help='input phylogeny option.')
    parser.add_argument('--tree_format', default='newick', help='input phylogeny file type.')
    parser.add_argument('--taxa_list', default=False, help='File that contains the names of each taxon to remove.')
    parser.add_argument('--output_tree_file', default='trimmed_tree.tre', help='output tree file. (DEFAULT: trimmed_tree.tre)')

    return parser.parse_args()

def main():
    args = parse_args()

    names_to_remove = []

    # # establish taxon namespace
    tns = dendropy.TaxonNamespace()
    #
    in_tree = dendropy.Tree.get(path=args.tree_file, schema=args.tree_format, taxon_namespace=tns, preserve_underscores=True)
    #
    # print(in_tree)

    input_taxa_list = open(args.taxa_list, 'r')
    #
    # split_names = args.taxa_list.split(',')

    for name in input_taxa_list:
        # print(name)
        names_to_remove.append(name.strip())
    # print(names_to_keep)

    taxa_to_keep = set([taxon for taxon in in_tree.taxon_namespace if taxon.label not in names_to_remove])

    # print(taxa_to_retain)


    filtered_tree = in_tree.extract_tree_with_taxa(taxa=taxa_to_keep)
    #
    #
    #
    print(filtered_tree.as_ascii_plot())

    filtered_tree.write(path=args.output_tree_file, schema="newick")





if __name__ == '__main__':
    main()
