#!/usr/bin/env python3

import os
import argparse
import pathlib
import shutil
import subprocess
import dendropy
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(prog='Tree to clades', \
        description='Split a large tree into more manageable monophyletic clades.')
    parser.add_argument('--tree_file', default=False, help='input phylogeny option.')
    parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    return parser.parse_args()

def main():
    args = parse_args()

    # # establish taxon namespace
    tns = dendropy.TaxonNamespace()
    #
    in_tree = dendropy.Tree.get(path=args.tree_file, schema='newick', taxon_namespace=tns)

    clades = make_clades_list(in_tree)

    #Input table of metadata to pandas
    metadata_table = pd.read_csv(args.metadata_csv)
    # print(metadata_table)

    find_good_clades(clades, metadata_table)


def make_clades_list(full_tree):

    taxon_sets = []
    for node in full_tree:
        # print(node)
        if 50 < len([leaf.taxon.label for leaf in node.leaf_iter()]) < 150:
            taxon_sets.append([leaf.taxon.label for leaf in node.leaf_iter()])

    print(len(taxon_sets))
    return taxon_sets


def find_good_clades(clades_, metadata_):

    amr_gene_column_name = 'AMR genotypes core'
    location_column = 'Location'
    collection_date_column = 'Collection date'
    # run_df = metadata_['Run']

    # print(metadata_.columns)
    # print(metadata_['Run'])

    results_of_analysis = {}

    for num, clade_of_tips in enumerate(clades_):
        list_of_gene_counts = []
        num_tips_no_amr_genes = 0
        num_tips_missing_location_data = 0

        number_of_tips = len(clade_of_tips)

        for tip in clade_of_tips:
            print(tip)
            # found_tip = metadata_[(metadata_.loc[tip])]
            try:
                found_tip = metadata_[metadata_.isin([tip]).any(axis=1)]
                tips_index = found_tip.index[0]
                print(tips_index)

                tips_amr_gene_num = found_tip.at[tips_index, amr_gene_column_name].split(',')
                print(len(tips_amr_gene_num))
                list_of_gene_counts.append(len(tips_amr_gene_num))
                if len(list_of_gene_counts) == 0:
                    num_tips_no_amr_genes+=1
                    print("***")

                    date_info = found_tip.at[tips_index, location_column]
                    print(date_info)
                    if type(date_info) == float:
                        # print("####################")
                        num_tips_missing_location_data+=1

            except IndexError:
                print("error")
            # tip_amr_genes_num =

        avg_of_amr_genes = sum(list_of_gene_counts) / len(list_of_gene_counts)
        print("number of tips in clade ", number_of_tips)
        print("avg num of amr genes ", avg_of_amr_genes)
        print("num tips with no amr genes ", num_tips_no_amr_genes)
        print("num tips missing location data ", num_tips_missing_location_data)


if __name__ == '__main__':
    main()
