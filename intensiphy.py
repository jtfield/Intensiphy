#! /usr/bin/python3

import os
import argparse
import re
import csv
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(prog='Intensiphy', \
        description='Automate downloading of high-throughput sequence data and updating of alignments using Extensiphy.')
    parser.add_argument('--align_file')
    parser.add_argument('--cores')
    parser.add_argument('--accession_csv')
    parser.add_argument('--update_method')
    return parser.parse_args()

def read_csv_file(csv_file):
    """Reads the accession file provided by the users."""
    csv = pd.read_csv(csv_file)
    
    just_accessions = csv["Run"]
    # print(just_accessions)

    return just_accessions

def calulate_cores(set_cores):
    """Organizes and calulates the cores for use with Extensiphy."""
    output = []
    total_cores = int(set_cores)
    runs = 0
    cores_per_run = 0
    
    if total_cores > 2 and total_cores < 10:
        runs = 2
        cores_per_run = total_cores / 2
        
    elif total_cores >= 10:
        runs = total_cores / 2
        cores_per_run = 2
    
    assert runs * cores_per_run <= total_cores
    
    output.append(int(runs))
    output.append(int(cores_per_run))

    # test to make sure we get a list of length 2
    # that contains 2 ints
    assert type(output) == list
    assert len(output) == 2
    for num in output:
        assert type(num) == int

    return output
    
def check_duplicate_accesions(accession_db, fasta_names):
    """Checks the alignment file vs the accesion list
    and produces a list of accessions that are already in the alignment."""
    duplicate_sra_runs = []
    
    for value in accession_db:
        if value in fasta_names:
            assert type(value) == str
            assert len(value) > 1
            duplicate_sra_runs.append(value)

    return duplicate_sra_runs

def read_fasta_names(align):
    """Reads the names of the sequences from the fasta file."""
    names = []
    with open(align, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                seq_name = line.strip(">").strip('\n')
                assert len(seq_name) > 1
                assert type(seq_name) == str
                names.append(seq_name)
    
    return names



def main():
    args = parse_args()

    # calculate the core organization to pass to Extensiphy
    get_cores = calulate_cores(args.cores)
    # print(get_cores)

    # read_file = read_csv_file(args.accession_csv)
    read_fasta = read_fasta_names(args.align_file)

    #read list of accessions
    read_accessions = read_csv_file(args.accession_csv)

    #Check list of run SRA numbers vs the sequences already in the alignment to prevent duplicates.
    check_dupes = check_duplicate_accesions(read_accessions, read_fasta)
    

if __name__ == '__main__':
    main()