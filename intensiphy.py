#! /usr/bin/python3

import os
import argparse
import re
import csv
import pandas as pd
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(prog='Intensiphy', \
        description='Automate downloading of high-throughput sequence data and updating of alignments using Extensiphy.')
    parser.add_argument('--align_file')
    parser.add_argument('--cores')
    parser.add_argument('--accession_csv')
    parser.add_argument('--ep_out_dir', help='Absolute path and folder name to create for outputs')
    parser.add_argument('-b', default=False, action='store_true', help='Toggles big bactch downloading \
        of fastq files instead of continuous downloading.')
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

    print(output)
    return output
    
def check_duplicate_accesions(accession_db, fasta_names):
    """Checks the alignment file vs the accesion list
    and produces a list of accessions that are already in the alignment."""
    non_duplicate_sra_runs = []

    for num, value in enumerate(accession_db):
        if value not in fasta_names:
            assert type(value) == str
            assert len(value) > 1
            non_duplicate_sra_runs.append(value)
            
    return non_duplicate_sra_runs

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

def downloading_and_running(method, accessions, runs):
    
    # Identify if we're processing data by downloading in batches between Extesniphy runs
    # or by downloading everything all at once and running extensiphy after
    if method == False:
        # continuous gradual downloading of data has been selected/left as default.
        batches_run = batch_download(accessions, runs)

    elif method == True:
        # Bulk download all fastq files before running Extensiphy
        bulk_run = bulk_download(accessions)

        # TODO: Extensiphy run goes here
        


def batch_download(accession_list, runs_number):
    """Prepares accessions to be downloaded in batches between runs of Extensiphy."""
    batches_of_accessions = prepare_batch_accessions(accession_list, runs_number)

        for accessions in batches_of_accessions:
        
            print(accessions)
            for single_accession in accessions:
                print(single_accession)

                subprocess.run(["fasterq-dump", "--split-files", single_accession])
        
            # TODO: Extensiphy runs go here

def prepare_batch_accessions(accessions, runs):
    """Return a list of lists containing the number of files to download before each Extensiphy run"""
    chunks = [accessions[x:x+runs] for x in range(0, len(accessions), runs)]
    
    return chunks

def bulk_download(accession_list, output_folder):
    """Bulk download every sequence (unlike the batch method)"""

    subprocess.run(["mkdir", output_folder + "/bulk_reads"])

    # TODO: make sure fastq files get into the reads folder

    for single_accession in accession_list:

        subprocess.run(["fasterq-dump", "--split-files", single_accession])

def main():
    args = parse_args()

    subprocess.run(["mkdir", args.ep_out_dir])

    os.chdir(args.ep_out_dir)

    # calculate the core organization to pass to Extensiphy
    get_cores = calulate_cores(args.cores)
    # print(get_cores)

    # read_file = read_csv_file(args.accession_csv)
    read_fasta = read_fasta_names(args.align_file)

    #read list of accessions
    read_accessions = read_csv_file(args.accession_csv)

    #Check list of run SRA numbers vs the sequences already in the alignment to prevent duplicates.
    removed_dupes = check_duplicate_accesions(read_accessions, read_fasta)
    
    # Handle how we'll download SRA files: big batch or continuously while running Extensiphy
    process_data = downloading_and_running(args.b, removed_dupes, get_cores[0])

if __name__ == '__main__':
    main()