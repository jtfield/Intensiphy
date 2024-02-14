#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
# import sys
# import json
# import dendropy

# from opentree import OT
# from opentree.annotations import *
# from annotations import write_itol_heatmap




def parse_args():
    parser = argparse.ArgumentParser(prog='amr_lister.py', \
        description='make itol annotations for ip data')
    parser.add_argument('--csv', help='csv file with AMR gene values.')
    return parser.parse_args()

def main():
    args = parse_args()

    columns = []

    csv_contents = pd.read_csv(args.csv)

    for idx, row in csv_contents.iterrows():

        amr_genes = row['AMR genotypes']

        split_row = amr_genes.split(',')
        # print(split_row)

        for gene_id in split_row:

            split_gene_info = gene_id.split('=')
            gene = split_gene_info[0]
            info = split_gene_info[1]

            if gene not in columns:
                columns.append(gene)

    print(columns)

    count_genes(columns, csv_contents)

def count_genes(genes_list, csv):
    gene_dict = {}

    for gene_name in genes_list:
        gene_dict[gene_name] = 0

    # print(gene_dict)
    for idx, row in csv.iterrows():

        amr_genes = row['AMR genotypes']

        split_row = amr_genes.split(',')
        # print(split_row)

        for gene_id in split_row:

            split_gene_info = gene_id.split('=')
            gene = split_gene_info[0]
            info = split_gene_info[1]

            gene_dict[gene]+=1
    print(gene_dict)
    print(len(gene_dict.keys()))

    make_box_plot(gene_dict)

    return gene_dict

def make_box_plot(dict):
    plt.bar(*zip(*dict.items()))
    plt.xticks(range(0, len(dict.items())), rotation='vertical')
    plt.savefig("amr_gene_prevelance.png", format="png", bbox_inches = "tight")
    plt.show()




if __name__ == '__main__':
    main()
