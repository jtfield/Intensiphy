#!/usr/bin/env python3

import os
import argparse
# import pathlib
# import shutil
# import subprocess
# import dendropy
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(prog='clade_metadata_stats.py', \
        description='calculate the average number of AMR genes for a set of pathogen detection samples metadata.')
    parser.add_argument('--metadata_csv', default=False, help='metadata csv with info on samples from the tree.')
    # parser.add_argument('--output_file', default='clade_to_keep.txt', help='output file containing taxon names separated by newlines.')
    return parser.parse_args()

def main():
    args = parse_args()

    df = pd.read_csv(args.metadata_csv)

    calculate_avg_amr(df)




def calculate_avg_amr(metadata_):

    # amr_gene_column_name = 'AMR genotypes core'
    amr_gene_column_name = 'AMR genotypes'

    total_amr_genes = 0
    total_taxa = 0

    for idx, row in metadata_.iterrows():
        # print(row[amr_gene_column_name])

        split_amr_genes = row[amr_gene_column_name].split(',')
        # print(split_amr_genes)

        gene_count = len(split_amr_genes)

        total_amr_genes+=gene_count
        total_taxa+=1

    print(total_amr_genes)
    print(total_taxa)

    average_amr_genes = total_amr_genes / total_taxa
    print(average_amr_genes)




                # tips_amr_gene_num = str(found_tip.at[tips_index, amr_gene_column_name]).split(',')



if __name__ == '__main__':
    main()
