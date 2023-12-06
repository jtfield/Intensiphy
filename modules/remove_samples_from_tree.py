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
    parser.add_argument('--tree_dir', default=False, help='input phylogeny option.')
    parser.add_argument('--rm_taxa', default=False, help='taxa you want to remove. if multiple taxa, separate with comma (,).')
    parser.add_argument('--output_tree_file', default='trimmed_tree.tre', help='output tree file. (DEFAULT: trimmed_tree.tre)')

    return parser.parse_args()

def main():
    args = parse_args()

    names_to_keep = []

    # names_to_remove = args.rm_taxa.split(',')

    tree_file_list = os.listdir(args.tree_dir)

    for file in tree_file_list:

        ref_name = file.replace('RAxML_bestTree_', '')
        ref_name = ref_name.replace('_tree', '')

        names_to_remove = [ref_name]

        # # establish taxon namespace
        tns = dendropy.TaxonNamespace()
        #
        in_tree = dendropy.Tree.get(path=args.tree_dir + '/' + file, schema='newick', taxon_namespace=tns, preserve_underscores=True)

        in_tree.prune_taxa_with_labels(names_to_remove)


        print(in_tree.as_ascii_plot())

        in_tree.write(path='ref_name_removed_RAxML_bestTree_' + ref_name + '_ref_tre', schema="newick")





if __name__ == '__main__':
    main()
