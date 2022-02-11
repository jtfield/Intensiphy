#! /usr/bin/python3

import os
import argparse
import pathlib
import numpy as np
import pandas as pd

def check_sequence_similarities(_seq_1, _seq_2):
    """Checks similarity of two sequences"""
    # print("waffle")

    # standard_nucleotides = ['A', 'C', 'G', 'T']
    # gaps = ['-']
    # degen_nucleotides = ['W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N']

    ident_nucs = 0
    ident_gaps = 0
    ident_degens = 0
    non_ident_bases = 0

    list_seq_1 = list(_seq_1)
    list_seq_2 = list(_seq_2)

    zipped_seqs = list(zip(list_seq_1, list_seq_2))

    results = map(check_nucs, zipped_seqs)

    for i in list(results):
        if i == 1:
            ident_nucs+=1
        elif i == 2:
            ident_gaps+=1
        elif i == 3:
            ident_degens+=1
        elif i == 4:
            non_ident_bases+=1

    summed_ident = ident_nucs + ident_gaps + ident_degens
    # print(summed_ident)
    # print(summed_ident / len(list_seq_1))
    # print(summed_ident / len(list_seq_2))
    assert (summed_ident / len(list_seq_1)) == summed_ident / len(list_seq_2)
    output = summed_ident / len(list_seq_1)

    return output


def check_nucs(nuc_tuple):
    standard_nucleotides = ['A', 'C', 'G', 'T']
    gaps = ['-']
    degen_nucleotides = ['W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N']
    nuc_1 = nuc_tuple[0].upper()
    nuc_2 = nuc_tuple[1].upper()

    if nuc_1 == nuc_2:
        if nuc_1 in standard_nucleotides:
            return 1
        elif nuc_1 in gaps:
            return 2
        elif nuc_1 in degen_nucleotides:
            return 3
    else:
        return 4

def build_sim_table(taxon, taxon_list):
    """Input a list of taxon names and build an empty table matching each taxon
        pairwise"""

    taxon = [taxon]

    df = pd.DataFrame(columns=taxon, index=taxon_list)

    return df

def seq_compare(working_dir):
    """Runs comparison program to assess and record sequence similarities"""

    # Initialize seqs dir
    seqs_dir = working_dir + '/sequence_storage'

    # # Initialize taxon name list
    # taxon_names = []
    #
    # # Loop over taxon names
    # for file in os.listdir(seqs_dir):
    #
    #     # Remove the suffix frome the taxon name
    #     clean_name = file.replace(suffix,'')
    #
    #     # Append the names to the list
    #     taxon_names.append(clean_name)
    #
    # df = ''
    # if run_bool:
    #     df = build_sim_table(taxon_names)
    #
    # else:
    #     df = pd.read_csv(working_dir + '/similarity_logs/sims_database.csv')

    # print(os.listdir(seqs_dir))
    # Loop over every file
    list_of_taxa = os.listdir(seqs_dir)

    suffix = '_.fas'

    for taxon_1 in list_of_taxa:

        taxon_1_dir_path = seqs_dir + '/' + taxon_1
        table_file_name = taxon_1 + '_sim_table.csv'

        if os.path.isdir(taxon_1_dir_path):
            print(taxon_1_dir_path)

            taxon_1_seq_file = taxon_1_dir_path + '/' + taxon_1 + suffix
            df = ''
            list_taxon_1_files = os.listdir(taxon_1_dir_path)

            if table_file_name in list_taxon_1_files:

                #PROCESS TABLE
                df = pd.read_csv(taxon_1_dir_path + '/' + table_file_name, sep=',', index_col=0)

            elif table_file_name not in list_taxon_1_files:

                # BUILD TABLE
                df = build_sim_table(taxon_1, list_of_taxa)

            print(df)
            # Loop over every file that isnt the current file
            for taxon_2 in os.listdir(seqs_dir):

                taxon_2_dir_path = seqs_dir + '/' + taxon_2

                if os.path.isdir(taxon_2_dir_path):

                    taxon_2_seq_file = taxon_2_dir_path + '/' + taxon_2 + suffix

                    # print("CHECK!!!")
                    # print(file_1)
                    # print(file_2)
                    if taxon_1 != taxon_2:
                        print("CHECK!!!")
                        print(taxon_1)
                        print(taxon_2)

                        # Open each file and make the sequence comparison
                        open_file_1 = open(taxon_1_seq_file, 'r')
                        open_file_2 = open(taxon_2_seq_file, 'r')

                        # Read and split the files
                        read_file_1 = open_file_1.read().split('\n', 1)
                        read_file_2 = open_file_2.read().split('\n', 1)

                        name_1 = read_file_1[0].strip('>')
                        name_2 = read_file_2[0].strip('>')

                        seq_1 = read_file_1[1]
                        seq_2 = read_file_2[1]

                        # Use comparison function and collect the float output
                        # similarity number
                        print("Checking similarities")
                        compare_output = check_sequence_similarities(seq_1, seq_2)

                        # df.at[name_1, name_2] = compare_output
                        df.at[name_2] = compare_output

            print(df)
            # df.to_csv(taxon_1_dir_path + '/' + table_file_name, sep=',')


def build_or_update_df(working_dir, run_bool, suffix):

    # Initialize seqs dir
    seqs_dir = working_dir + '/sequence_storage'

    # Initialize taxon name list
    taxon_names = []

    # Loop over taxon names
    for file in os.listdir(seqs_dir):

        # Remove the suffix frome the taxon name
        clean_name = file.replace(suffix,'')

        # Append the names to the list
        taxon_names.append(clean_name)

    df = ''
    if run_bool == False:
        df = build_sim_table(taxon_names)

    elif run_bool:
        df = pd.read_csv(working_dir + '/similarity_logs/sims_database.csv', sep=',', index_col=0)

        df = check_duplicates(taxon_names, df)

    # Compare all files and add values to the df
    output = seq_compare(working_dir, suffix, df)

    df.to_csv(working_dir + '/similarity_logs/sims_database.csv', sep=',')

    # print(df)

def check_duplicates(taxon_list, df):
    """Check if names already exist in the df. If not, add them and return df"""
    # Get list of df columns
    columns = list(df.columns)

    # check if there are new names in the taxon list
    for name in taxon_list:

        # If name not in the dataframe already, add an empty row and column for it
        if name not in columns:

            # Add row
            df.append(pd.Series(name=name))

            # Add column
            df[name] = np.nan

    return df

def make_keep_decision(coverage_dict, taxon_list):
    """Reads the dictionary containing similarity coverage. \
    Makes decision on which sequences to keep and which to remove"""
    keep_list = []
    cov_rank_dict = {}
    coverage_list = []
    intermediate_keep_list = []

    for key, value in coverage_dict.items():
        num_covered_taxa = len(value)
        # print(num_covered_taxa)
        coverage_list.append(num_covered_taxa)
        if num_covered_taxa in cov_rank_dict.keys():
            cov_rank_dict[num_covered_taxa].append(key)
        else:

            cov_rank_dict[num_covered_taxa] = []
            cov_rank_dict[num_covered_taxa].append(key)

    set_cov = set(coverage_list)
    sorted_cov = list(set_cov)
    sorted_cov.sort()

    # print(sorted_cov)
    # print(cov_rank_dict)

    # Loop over the coverage list in reverse
    # Check if all the sequences

    # for coverage_level in sorted_cov[::-1]:
    #     current_taxa = cov_rank_dict[coverage_level]
    #     # print(current_taxa)
    #     for taxon in current_taxa:
    #         print(taxon)
    #         print(coverage_dict[taxon])
    keep_these_taxa = keep_logic_gate(sorted_cov, cov_rank_dict, coverage_dict, keep_list, taxon_list)

    return keep_these_taxa

def keep_logic_gate(sorted_coverage_list, coverage_rank_dict, cov_dict, primary_keep_list, tax_list):
    """Final logic checks on which sequences are covered by other sequences"""
    # print(sorted_coverage_list)
    # print(coverage_rank_dict)
    # print(cov_dict)
    # print(tax_list)

    added_levels = []
    checked_covering_taxa = []
    intermediate_list = []
    for coverage_level in sorted_coverage_list[::-1]:
        current_taxa = coverage_rank_dict[coverage_level]
        # print(current_taxa)
        for taxon in current_taxa:
            # print(taxon)
            # print(cov_dict[taxon])
            if coverage_level not in added_levels:
                added_levels.append(coverage_level)
                intermediate_list.append(taxon)
                checked_covering_taxa.append(taxon)
                # checked_covering_taxa.append(taxon)
                for covered_taxon in cov_dict[taxon]:
                    # if covered_taxon not in intermediate_list:
                    #     intermediate_list.append(covered_taxon)
                    if covered_taxon not in intermediate_list:
                        intermediate_list.append(covered_taxon)

    added_levels.sort()

    if added_levels == sorted_coverage_list:
        print("loop complete")
        present_taxa = 0
        for taxon in tax_list:
            if taxon in intermediate_list:
                present_taxa+=1
            else:
                print(taxon)
                checked_covering_taxa.append(taxon)
    print(checked_covering_taxa)

    return checked_covering_taxa










def check_sims_and_remove(working_dir, cutoff):
    """Assess similarities in database and decide which sequences to remove"""

    # Read similarities df
    df = pd.read_csv(working_dir + '/similarity_logs/sims_database.csv', sep=',', index_col=0)


    # Get list of taxa from the df columns
    taxa = list(df.columns)
    #
    df = df.fillna(0)
    # print(df)

    # make a copy of the df, subtracting the cutoff value from each similarity
    sub_df = df.apply(lambda x: x - cutoff)
    # print(sub_df)

    similarity_coverage_dict = check_pairing(sub_df)

    keep_seqs = make_keep_decision(similarity_coverage_dict, taxa)

    return keep_seqs


    # print(sub_df)

def check_pairing(df):
    """Function to get the taxon pairs and check if any beat the cutoff"""
    taxa = list(df.columns)

    output_dict = {}

    for taxon_1 in taxa:
        for taxon_2 in taxa:
            if taxon_1 != taxon_2:
                sim_value = df.at[taxon_1, taxon_2]
                if sim_value > 0:
                    # print(taxon_1 + " " + taxon_2)
                    if taxon_1 not in output_dict.keys():
                        output_dict[taxon_1] = []
                        output_dict[taxon_1].append(taxon_2)
                    else:
                        if taxon_2 not in output_dict[taxon_1]:
                            output_dict[taxon_1].append(taxon_2)

    return output_dict
