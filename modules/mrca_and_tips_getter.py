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
    parser = argparse.ArgumentParser(prog='MRCA and tips getter', \
        description='Takes a list of taxa you wish to get the mcra for and all tips connected to that MRCA.')
    # parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    parser.add_argument('--big_csv', help='CSV file containing only the clade of interest metadata.')
    parser.add_argument('--taxa_set', default=False, help='File with the taxa of interest within the clade of interest.')
    parser.add_argument('--tree_file', default=False, help='input phylogeny option.')
    parser.add_argument('--taxa_file_type', default='csv', help='File type of list of taxa. Either csv or a staight list .txt File.')
    parser.add_argument('--output_subset_name', default='selected', help='The first portion of the output file names. Best set as a cluster or a subset name. (DEFAULT: selected)')
    return parser.parse_args()

def main():
    args = parse_args()

    if args.taxa_file_type != 'csv':

        if_list(args.tree_file, args.taxa_set, args.output_subset_name, args.big_csv)

    elif args.taxa_file_type == 'csv':

        if_csv(args.tree_file, args.taxa_set, args.output_subset_name, args.big_csv)



def if_csv(tree_file, taxa_list, cluster_name, big_csv):

    names_to_keep = []
    found_names = []

    # # establish taxon namespace
    tns = dendropy.TaxonNamespace()
    #
    in_tree = dendropy.Tree.get(path=tree_file, schema='newick', taxon_namespace=tns, preserve_underscores=True)

    input_taxa_list = pd.read_csv(taxa_list)

    big_table = pd.read_csv(big_csv)

    clade_names = list(input_taxa_list['Run'])

    # print(clade_names)

    for taxon in in_tree.taxon_namespace:
        if taxon.label in clade_names:
            found_names.append(taxon.label)

    # print(found_names)
    # print(clade_names)
    print("number of taxa from the clade found in the input tree: ", len(found_names))
    print("number of taxa in the input clade: ", len(clade_names))

    found_names = set(found_names)

    mrca = in_tree.mrca(taxon_labels=found_names)
    #
    print(mrca)

    for leaf_node in mrca.leaf_iter():
        # print(dir(leaf_node))
        tip_name = str(leaf_node.taxon)

        names_to_keep.append(tip_name.strip("'"))

    # print(names_to_keep)
    print("number of names found connected to the MRCA: ", len(names_to_keep))

    updated_clade = pd.DataFrame(columns=big_table.columns)

    print(big_csv)

    for name in names_to_keep:
        print(name)
        # print(big_table['Run'])

        try:
            found_id = big_table.loc[big_table['Run'].str.contains(name, na=False)]

            updated_clade = updated_clade.append(found_id)

        except IndexError:
            print("Missing tip in chosen clade: ", name)


    updated_clade.to_csv(cluster_name + '_subset_metadata.csv')

    # print(names_to_keep)
    print("number of taxa from the clade found in the input tree: ", len(found_names))
    print("number of taxa in the input clade: ", len(clade_names))

    taxa_to_retain = set([taxon for taxon in in_tree.taxon_namespace if taxon.label in names_to_keep])

    filtered_tree = in_tree.extract_tree_with_taxa(taxa=taxa_to_retain)

    # print(filtered_tree.as_ascii_plot())

    filtered_tree.write(path=cluster_name + '_subset.tre', schema="newick")




if __name__ == '__main__':
    main()
