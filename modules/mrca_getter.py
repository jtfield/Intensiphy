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
    parser = argparse.ArgumentParser(prog='MRCA getter', \
        description='Takes in a tree and a csv of taxa representing a clade you want to collect all the members of. \
        Outputs a list of all taxa in the chosen clade. Specifically for finding new taxa added during placement.')
    # parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    parser.add_argument('--clade_csv_file', help='CSV file containing only the clade of interest metadata.')
    parser.add_argument('--big_csv_file', help='CSV file containing all samples being analyzed.')
    parser.add_argument('--tree_dir', help='directory of tree files.')
    parser.add_argument('--tree_suffix', default='.tre', help='Newick phylogeny file.')
    parser.add_argument('--output_file', default='updated_query_clade.csv', help='CSV file.')
    return parser.parse_args()

def main():
    args = parse_args()

    # tree = dendropy.Tree.get(path=args.tree_file, schema="newick", preserve_underscores=True)

    output = []

    clade_csv = pd.read_csv(args.clade_csv_file)

    big_csv = pd.read_csv(args.big_csv_file)

    clade_names = list(clade_csv['Run'])

    big_csv_columns = big_csv.columns

    # mrca = tree.mrca(taxon_labels=clade_names)

    list_of_files = os.listdir(args.tree_dir)

    for file_name in list_of_files:
        if file_name.endswith(args.tree_suffix):
            tree = dendropy.Tree.get(path=args.tree_dir + '/' + file_name, schema="newick", preserve_underscores=True)

            mrca = tree.mrca(taxon_labels=clade_names)

            for leaf_node in mrca.leaf_iter():
                # print(dir(leaf_node))
                tip_name = str(leaf_node.taxon)

                output.append(tip_name.strip("'"))


    output = set(output)
    print(len(output))

    updated_clade = pd.DataFrame(columns=big_csv_columns)

    for name in output:
        name = name.replace('QUERY___','')

        try:
            found_id = big_csv.loc[big_csv['Run'].str.contains(name, na=False)]

            updated_clade = updated_clade.append(found_id)

        except IndexError:
            print("Missing tip in chosen clade: ", name)

    # print(updated_clade)

    updated_clade.to_csv(args.output_file)





if __name__ == '__main__':
    main()
