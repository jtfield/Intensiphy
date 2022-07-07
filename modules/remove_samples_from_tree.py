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
        description='After running EP on multiple reference sequences, analyze the results.')
    parser.add_argument('--tree_file', default=False, help='input phylogeny option.')
    # parser.add_argument('--taxa_list', default=False, help='File that contains the names of each taxon to keep.')
    parser.add_argument('--output_tree_file', default='trimmed_tree.tre', help='output tree file. (DEFAULT: trimmed_tree.tre)')

    return parser.parse_args()

def main():
    args = parse_args()

    names_to_keep = []

    names_to_remove = ['SRR8833557', 'SRR17007432', 'SRR16996016', 'SRR17541083', 'SRR9877804', 'SRR9985644', 'SRR11904725', 'SRR11724687', 'SRR11698932', 'SRR11724686']

    # # establish taxon namespace
    tns = dendropy.TaxonNamespace()
    #
    in_tree = dendropy.Tree.get(path=args.tree_file, schema='newick', taxon_namespace=tns, preserve_underscores=True)

    in_tree.prune_taxa_with_labels(names_to_remove)


    print(in_tree.as_ascii_plot())

    in_tree.write(path=args.output_tree_file, schema="newick")





if __name__ == '__main__':
    main()
