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
import dendropy


def parse_args():
    parser = argparse.ArgumentParser(prog='monophyly checker', \
        description='Reads a tree file and finds the clade-of-interest. \
        Then checks if there are taxa not in the clade table that are present in this version of the tree. \
        Also checks if there are taxa from the clade table missing in this version of the tree.')
    # parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    parser.add_argument('--clade_csv_file', help='CSV file containing only the clade of interest metadata.')
    parser.add_argument('--tree', help='Phylogenetic tree in newick format.')
    # parser.add_argument('--output_file', default='updated_query_clade.csv', help='CSV file.')
    return parser.parse_args()

def main():
    args = parse_args()

    output = []

    clade_csv = pd.read_csv(args.clade_csv_file)

    clade_names = list(clade_csv['Run'])

    tree = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True)

    mrca = tree.mrca(taxon_labels=clade_names)

    for leaf_node in mrca.leaf_iter():
        # print(dir(leaf_node))
        tip_name = str(leaf_node.taxon)

        output.append(tip_name.strip("'"))

    # print(output)
    # print(clade_names)

    names_in_one_not_in_another = 0

    for name in clade_names:
        if name not in output:
            print("Name in input clade NOT found in tree: ", name)
            names_in_one_not_in_another+=1

    for name2 in output:
        if name2 not in clade_names:
            print("Name in tree NOT found in provided clade: ", name2)
            names_in_one_not_in_another+=1

    #Just another check to see if the lists are equivalent
    if names_in_one_not_in_another == 0:
        assert clade_names.sort() == output.sort()




if __name__ == '__main__':
    main()
