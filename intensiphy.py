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
# import sys
# sys.path.append('./modules')
# import modules.seq_similarity_assessment
# import modules.alignment_splitter
from modules.seq_similarity_assessment import *
from modules.alignment_splitter import split_alignment
from modules.fetch_and_align import *
from modules.tree_assess import *
# import tests.accessions_tests.accession_download_test
# import tests.assembly_tests.gon_phy_test

def parse_args():
    parser = argparse.ArgumentParser(prog='Intensiphy', \
        description='Automate downloading of high-throughput sequence data and updating of alignments using Extensiphy.')
    parser.add_argument('--align_file', default=False, help='input alignment option. \
        If an alignment is added using this option, samples will be added to this alignment.')
    parser.add_argument('--starting_tree', default=False, type=str, help='Path \
    to a newick tree file corresponding to an input alignment. Don\' use if not inputting an alignment.' )
    parser.add_argument('--cores')
    parser.add_argument('--accs_file', default=False, help='Accession file if accession_method is set to USER_INPUT')
    parser.add_argument('--accession_method', default="AUTO_DL", type=str, \
    help='Dictates how collecting and inputting accession numbers will be handled. \
        (OPTIONS: USER_INPUT and AUTO_DL), (DEFAULT: AUTO_DL)')
    parser.add_argument('--ep_out_dir', default='ip_output', help='Absolute path and folder name to create for outputs')
    # parser.add_argument('--path_to_ep_dir', help='Absolute path and folder name to create for outputs')
    parser.add_argument('--organism', type=str, nargs='+', help='scientific name of the organism or group of organisms you \
        would like to query SRA for and update your alignment with. Example: Neisseria gonorrhoeae[Organism] or txid482')
    parser.add_argument('--ref', default=False, type=str, help='reference sequence label (without suffix or file ending information). (Example: SRR1500345)')
    # parser.add_argument('--sim', type=float, help='sequence similarity cutoff to exclude sequences from the alignment. \
        # When two sequences surpass this cutoff, one sequence will be chosen to represent both.')
    # parser.add_argument('--ds', type=int, help='Allocated disk space Intensiphy can use.')
    # parser.add_argument('-b', default=False, action='store_true', help='Toggles big bactch downloading \
        # of fastq files instead of continuous downloading.')
    return parser.parse_args()

def main():
    args = parse_args()
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    ref_taxon = ""

    # reset_tests()

    # calculate the core organization to pass to Extensiphy
    print("Assessing allocated cores.")
    get_cores = calulate_cores(args.cores)

    # Check if output dir has been made already from a previous run
    # If output dir exists, the fundamentals of the program change to suit
    # a repeating run.
    # Returns True if output dir exists already and False if output dir doesn't exist prior to execution.
    print("Checking if this is a new run or a continuing run.")
    dir_existence = check_dir_exists(args.ep_out_dir)

    os.chdir(args.ep_out_dir)
    print(args.ep_out_dir)
    absolute_output_dir_path = os.path.abspath(os.getcwd())
    print(absolute_output_dir_path)
    #
    # # download_accessions(args.organism, args.ep_out_dir)
    # print("Working out how to handle getting accession numbers.")
    accessions = handle_accession_options(args.accession_method, args.organism, absolute_output_dir_path, args.accs_file)

    #read list of accessions
    read_accessions = ''
    if args.accession_method == 'AUTO_DL':
        read_accessions = read_csv_file(absolute_output_dir_path)

    elif args.accession_method == 'USER_INPUT':
        read_accessions = read_pathodb_csv_file(absolute_output_dir_path)


    print("Working out what kind of run we're doing.")
    if dir_existence == False:

        make_align(dir_existence, absolute_output_dir_path, args.accs_file, args.align_file)

        split_alignment(absolute_output_dir_path + '/intermediate_files/alignment.fas', absolute_output_dir_path + '/sequence_storage')

        select_ref(args.ref, absolute_output_dir_path)

        handle_starting_tree(absolute_output_dir_path, get_cores, args.starting_tree)

    else:
        print("Proceeding with new run using previous database.")

    clean_incomplete_downloads(absolute_output_dir_path)

    read_fasta = read_fasta_names(absolute_output_dir_path)

    # #read list of accessions
    # read_accessions = read_csv_file(absolute_output_dir_path)
    remove_paired_dupes = ''
    remove_single_dupes = ''
    #Check list of run SRA numbers vs the sequences already in the alignment to prevent duplicates.
    if args.accession_method == 'AUTO_DL':
        remove_paired_dupes = check_duplicate_accesions(read_accessions[0], read_fasta)
        remove_single_dupes = check_duplicate_accesions(read_accessions[1], read_fasta)
    elif args.accession_method == 'USER_INPUT':
        remove_paired_dupes = check_duplicate_accesions(read_accessions, read_fasta)

    # print(remove_paired_dupes)
    # print(remove_single_dupes)

    paired_batch_accessions = prepare_batch_accessions(remove_paired_dupes, get_cores[0])
    single_batch_accessions = prepare_batch_accessions(remove_single_dupes, get_cores[0])

    print(paired_batch_accessions)
    print(single_batch_accessions)

    write_current_run_names(absolute_output_dir_path, paired_batch_accessions)
    write_current_run_names(absolute_output_dir_path, single_batch_accessions)

    if len(paired_batch_accessions) > 0:
        print("Processing paired-end read files.")
        process_data = downloading_and_running(paired_batch_accessions, absolute_output_dir_path, get_cores, "PAIRED")

    if len(single_batch_accessions) > 0:
        print("Processing single-end read files.")
        process_data = downloading_and_running(single_batch_accessions, absolute_output_dir_path, get_cores, "SINGLE")

    # Make new alignment and perform placement into tree
    construct_align_and_place(absolute_output_dir_path)





def check_dir_exists(dir_name):
    """Check if output directory exists. If the directory exists, a previous run of Intensiphy was probably performed. \
        If the directory doesn't exist, no previous run was performed, build all necessary sub directories."""
    print("Checking if output directory already exists.")
    does_dir_exist = os.path.isdir(dir_name)
    print("Does the output specified output directory exist?: ", does_dir_exist)

    if does_dir_exist:

        assert os.path.isdir(dir_name + "/read_files")
        assert os.path.isdir(dir_name + "/run_log_files")
        assert os.path.isdir(dir_name + "/intermediate_files")
        assert os.path.isdir(dir_name + "/accession_files")
        assert os.path.isdir(dir_name + "/sequence_storage")
        # assert os.path.isdir(dir_name + "/similarity_logs")

        return does_dir_exist

    else:
        subprocess.run(["mkdir", dir_name])
        subprocess.run(["mkdir", dir_name + "/read_files"])
        subprocess.run(["mkdir", dir_name + "/run_log_files"])
        subprocess.run(["mkdir", dir_name + "/intermediate_files"])
        subprocess.run(["mkdir", dir_name + "/accession_files"])
        subprocess.run(["mkdir", dir_name + "/sequence_storage"])
        # subprocess.run(["mkdir", dir_name + "/similarity_logs"])

        assert os.path.isdir(dir_name + "/read_files")
        assert os.path.isdir(dir_name + "/run_log_files")
        assert os.path.isdir(dir_name + "/intermediate_files")
        assert os.path.isdir(dir_name + "/accession_files")
        assert os.path.isdir(dir_name + "/sequence_storage")
        # assert os.path.isdir(dir_name + "/similarity_logs")

        return does_dir_exist

if __name__ == '__main__':
    main()
