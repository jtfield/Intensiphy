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
    parser.add_argument('--output_file', default='clade_to_keep.txt', help='output file containing taxon names separated by newlines.')
    parser.add_argument('--lower_bound', default=50, help='lowest number of taxa to consider.')
    parser.add_argument('--upper_bound', default=150, help='highest number of taxa to consider.')
    return parser.parse_args()

def main():
    args = parse_args()

    # # establish taxon namespace
    tns = dendropy.TaxonNamespace()
    #
    in_tree = dendropy.Tree.get(path=args.tree_file, schema='newick', taxon_namespace=tns, preserve_underscores=True)

    clades = make_clades_list(in_tree, args.lower_bound, args.upper_bound)
    # print(clades)

    #Input table of metadata to pandas
    metadata_table = pd.read_csv(args.metadata_csv)
    # print(metadata_table)

    processed_clades = find_good_clades(clades, metadata_table)
    # print(processed_clades)

    chosen_clade = make_clade_decision(processed_clades)
    # print(chosen_clade)

    # print(clades[chosen_clade])
    # print(chosen_clade)
    # print(clades)

    clade_metadata = pull_clade_metadata(clades[chosen_clade], metadata_table)

    clade_metadata.to_csv(args.output_file)


def make_clades_list(full_tree, lower_bound, upper_bound):

    lower_bound = int(lower_bound)
    upper_bound = int(upper_bound)
    # print(lower_bound)
    # print(upper_bound)

    taxon_sets = []
    for node in full_tree:
        # print(node)
        if lower_bound < len([leaf.taxon.label for leaf in node.leaf_iter()]) < upper_bound:
            taxon_sets.append([leaf.taxon.label for leaf in node.leaf_iter()])

    # print(taxon_sets)
    # print(len(taxon_sets))
    return taxon_sets


def find_good_clades(clades_, metadata_):

    # amr_gene_column_name = 'AMR genotypes core'
    amr_gene_column_name = 'AMR genotypes'
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
        num_missing_tips = 0

        number_of_tips = len(clade_of_tips)

        for tip in clade_of_tips:
            # print(tip)
            # found_tip = metadata_[(metadata_.loc[tip])]
            try:
                found_tip = metadata_[metadata_.isin([tip]).any(axis=1)]
                tips_index = found_tip.index[0]
                # print(tips_index)

                tips_amr_gene_num = str(found_tip.at[tips_index, amr_gene_column_name]).split(',')
                # print(len(tips_amr_gene_num))
                list_of_gene_counts.append(len(tips_amr_gene_num))
                if len(list_of_gene_counts) == 0:
                    num_tips_no_amr_genes+=1
                    # print("***")

                location_info = found_tip.at[tips_index, location_column]
                # print(location_info)
                # if type(location_info) == float:
                if str(location_info) == 'nan':
                    # print("####################")
                    num_tips_missing_location_data+=1

            except IndexError:
                print("Error: Tip missing from metadata file: ", tip)
                num_missing_tips+=1

        # print("number of tips in clade ", number_of_tips)
        # # print("avg num of amr genes ", avg_of_amr_genes)
        # print("num tips with no amr genes ", num_tips_no_amr_genes)
        # print("num tips missing location data ", num_tips_missing_location_data)
        # print("num missing tips from metadata table ", num_missing_tips)

        avg_of_amr_genes = sum(list_of_gene_counts) / len(list_of_gene_counts)

        results_of_analysis[num] = [number_of_tips, avg_of_amr_genes, num_tips_no_amr_genes, num_tips_missing_location_data, num_missing_tips]

        # print("number of tips in clade ", number_of_tips)
        # print("avg num of amr genes ", avg_of_amr_genes)
        # print("num tips with no amr genes ", num_tips_no_amr_genes)
        # print("num tips missing location data ", num_tips_missing_location_data)
        # print("num missing tips from metadata table ", num_missing_tips)

    # print(results_of_analysis)
    return results_of_analysis


def make_clade_decision(processed_clades_):
    current_best_clade = None
    best_num_tips = 0
    best_avg_amr_genes = 0
    best_num_tips_no_amr_genes = 100000000
    best_num_tips_missing_location_data = 100000000
    best_num_missing_tips = 100000000

    for key, value in processed_clades_.items():
        num_tips = value[0]
        avg_genes = value[1]
        no_amr_genes = value[2]
        missing_loc = value[3]
        missing_tips = value[4]

        if num_tips >= best_num_tips:
            # print("best number of tips")
            # print(best_num_tips, num_tips)
            if avg_genes >= best_avg_amr_genes:
                # print("best avg amr genes count")
                # print(best_avg_amr_genes, avg_genes)
                if no_amr_genes <= best_num_tips_no_amr_genes:
                    # print("best number of tips with no amr genes")
                    # print(best_num_tips_no_amr_genes, no_amr_genes)
                    if missing_loc <= best_num_tips_missing_location_data:
                        # print("best number of tips missing location")
                        # print(best_num_tips_missing_location_data, missing_loc)
                        if missing_tips <= best_num_missing_tips:
                            # print("best number of tips missing from metadata")
                            # print(best_num_missing_tips, missing_tips)

                            # print("################")
                            current_best_clade = key
                            # print("best clade ", current_best_clade)
                            best_num_tips = num_tips
                            # print("best number of tips ", best_num_tips)
                            best_avg_amr_genes = avg_genes
                            # print("best avg amr genes ", best_avg_amr_genes)
                            best_num_tips_no_amr_genes = no_amr_genes
                            # print("best num tips with no amr genes ", best_num_tips_no_amr_genes)
                            best_num_tips_missing_location_data = missing_loc
                            # print("best num tips missing loc data ", best_num_tips_missing_location_data)
                            best_num_missing_tips = missing_tips
                            # print("least number of missing tips ", best_num_missing_tips)
                            # print("!!!!!!!!!!!!!!!!!!!!!!!!!")

    print(processed_clades_[current_best_clade])
    return current_best_clade


def pull_clade_metadata(clade_, big_df_):

    #pull column names of big metadata df
    df_columns = big_df_.columns

    #initialize dataframe with columns
    clade_df = pd.DataFrame(columns=df_columns)

    big_df_drop_na_runs = big_df_.dropna(subset=['Run'])

    #Find sample ids in the big dataframe
    for id in clade_:
        try:
            # print(id)
            # print(big_df_drop_na_runs['Run'])

            found_id = big_df_drop_na_runs.loc[big_df_drop_na_runs['Run'].str.contains(id, na=False)]
            # print(found_id)

            clade_df = clade_df.append(found_id)

            # clade_df = clade_df.append(pd.Series(found_id, index=clade_df.columns[:len(found_id)]), ignore_index=True)

            # clade_df.loc[len(clade_df.index)] = found_id

        except IndexError:
            print("Missing tip in chosen clade: ", id)


    # print(clade_df)

    return clade_df


if __name__ == '__main__':
    main()
