#! /usr/bin/python3

import os
import argparse
import pathlib
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

def build_sim_table(taxon_list):
    """Input a list of taxon names and build an empty table matching each taxon
        pairwise"""

    df = pd.DataFrame(columns=taxon_list, index=taxon_list)

    return df

def seq_compare(working_dir, suffix, df):
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

    print(os.listdir(seqs_dir))
    # Loop over every file
    for file_1 in os.listdir(seqs_dir):

        # Loop over every file that isnt the current file
        for file_2 in os.listdir(seqs_dir):
            # print("CHECK!!!")
            # print(file_1)
            # print(file_2)
            if file_1 != file_2:
                # print("CHECK!!!")
                # print(file_1)
                # print(file_2)

                # Open each file and make the sequence comparison
                open_file_1 = open(seqs_dir + '/' + file_1, 'r')
                open_file_2 = open(seqs_dir + '/' + file_2, 'r')

                # Read and split the files
                read_file_1 = open_file_1.read().split('\n', 1)
                read_file_2 = open_file_2.read().split('\n', 1)

                name_1 = read_file_1[0].strip('>')
                name_2 = read_file_2[0].strip('>')

                seq_1 = read_file_1[1]
                seq_2 = read_file_2[1]

                # Use comparison function and collect the float output
                # similarity number
                compare_output = check_sequence_similarities(seq_1, seq_2)

                df.at[name_1, name_2] = compare_output

        # return df

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
        df = pd.read_csv(working_dir + '/similarity_logs/sims_database.csv', sep=',')

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
