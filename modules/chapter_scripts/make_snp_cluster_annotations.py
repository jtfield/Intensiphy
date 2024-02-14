#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import sys
import json
import dendropy

# from opentree import OT
from opentree.annotations import *
# from annotations import write_itol_heatmap




def parse_args():
    parser = argparse.ArgumentParser(prog='make itol annotations', \
        description='make itol annotations for ip data')
    parser.add_argument('--csv', help='csv file with AMR gene values.')
    # parser.add_argument('--amr_locus', default='', help='csv file with AMR gene values.')
    # parser.add_argument('--output_file', default='converted_tree.tre', help='newick phylogeny file.')
    return parser.parse_args()

def main():
    args = parse_args()

    csv_contents = pd.read_csv(args.csv)

    cluster_info = make_cluster_dict(csv_contents)

    # print(gene_info[0])
    snp_cluster = cluster_info[1][0]
    cluster_dict = cluster_info[0]

    file_name = 'snp_cluster_' + snp_cluster + '_heatmap.txt'
    # print(file_name)
    # print(snp_cluster)
    # print(cluster_dict)

    #UNCOMMENT THESE LINES TO OUTPUT ANNOTATION FILES

    write_itol_heatmap(file_name, snp_cluster, 'presence', cluster_dict, snp_cluster)


def make_cluster_dict(csv_contents):

    columns = []
    output = {}
    wrapper_list = []
    unanimous_check_dict = {}
    genes_absent = []
    genes_ubiquitous = []

    # amr_genes = csv_contents['AMR genotypes']

    num_rows = len(csv_contents.index)
    # print(num_rows)

    # for list_of_genes in amr_genes:
    for idx, row in csv_contents.iterrows():

        cluster_id = row['SNP cluster']

        # split_row = amr_genes.split(',')
        # print(split_row)

        taxon_id = row['Run']

        internal_dict = {cluster_id : 1}

        # print(cluster_id)
        # print(taxon_id)
        # print(internal_dict)

        output[taxon_id] = internal_dict

        if cluster_id not in columns:
            columns.append(cluster_id)


    wrapper_list.append(output)
    wrapper_list.append(columns)

    # print(wrapper_list)

    return wrapper_list

        # for cluster_id in split_row:
        #
        #     split_gene_info = gene_id.split('=')
        #     gene = split_gene_info[0]
        #     info = split_gene_info[1]
        #
        #     if gene not in columns:
        #         columns.append(gene)
        #         unanimous_check_dict[gene] = 0
    # print(columns)
    # print(len(columns))

    # print(unanimous_check_dict)

    # for idx, row in csv_contents.iterrows():
    #
    #     amr_genes = row['AMR genotypes']
    #
    #     taxon_id = row['Run']
    #
    #     split_row = amr_genes.split(',')
    #     isolate_dict = {}
    #
    #     intermediate_list = []
    #     for gene_id in split_row:
    #
    #         split_gene_info = gene_id.split('=')
    #         gene = split_gene_info[0]
    #         info = split_gene_info[1]
    #
    #         isolate_dict[gene] = 1
    #         intermediate_list.append(gene)
    #         unanimous_check_dict[gene] = unanimous_check_dict[gene] + 1
    #
    #     for gene_name in columns:
    #         if gene_name not in intermediate_list:
    #             isolate_dict[gene_name] = 0
    #
    #     # print(isolate_dict)
    #     output[taxon_id] = isolate_dict
    #
    # # print(output)
    # # print(unanimous_check_dict)
    #
    # for key, value in unanimous_check_dict.items():
    #     if value == 0:
    #         genes_absent.append(key)
    #
    #     elif value == num_rows:
    #         genes_ubiquitous.append(key)
    #
    # # print(output)
    # print(genes_absent)
    # print(genes_ubiquitous)
    #
    # wrapper_list.append(output)
    # wrapper_list.append(columns)
    #
    # print(wrapper_list)
    #
    # return wrapper_list


# def check_unanimous(dict_of_isolates):






if __name__ == '__main__':
    main()
