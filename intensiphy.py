#!/usr/bin/env python3

import os
import argparse
import subprocess
from typing import List
from modules.alignment_splitter import split_alignment
from modules.fetch_and_align import *
from modules.tree_assess import *

def parse_args():
    parser = argparse.ArgumentParser(prog='Intensiphy', \
        description='Automate downloading of high-throughput sequence data and updating of alignments using Extensiphy.')
    parser.add_argument('--align_file', default=False, help='input alignment option. \
        If an alignment is added using this option, samples will be added to this alignment.')
    parser.add_argument('--starting_tree', default=False, type=str, help='Path \
    to a newick tree file corresponding to an input alignment. Don\' use if not inputting an alignment.' )
    parser.add_argument('--cores', default=2, nargs='?', const=2, type=int, help='Number of cores youd like to allocate to Intensiphy')
    parser.add_argument('--accs_file', default=False, help='Accession file if accession_method is set to USER_INPUT')
    parser.add_argument('--accession_method', default="AUTO_DL", type=str, \
    help='Dictates how collecting and inputting accession numbers will be handled. \
        (OPTIONS: USER_INPUT and AUTO_DL), (DEFAULT: AUTO_DL)')
    parser.add_argument('--ip_out_dir', default='ip_output', help='Absolute path and folder name to create for outputs')
    # parser.add_argument('--path_to_ep_dir', help='Absolute path and folder name to create for outputs')
    parser.add_argument('--organism', type=str, nargs='+', default=False, help='scientific name of the organism or group of organisms you \
        would like to query SRA for and update your alignment with. Example: Neisseria gonorrhoeae[Organism] or txid482')
    parser.add_argument('--ref', default=False, type=str, help='reference sequence label (without suffix or file ending information). (Example: SRR1500345)')
    parser.add_argument('--placement', default='OFF', help='Flag controls whether to perform phylogenetic placement once sequence collection and assembly is complete. \
                        Also controls whether a starting tree is built if one isnt input by the user. \
                        (OPTIONS: ON and OFF) (DEFAULT: OFF)')
    # parser.add_argument('--sim', type=float, help='sequence similarity cutoff to exclude sequences from the alignment. \
        # When two sequences surpass this cutoff, one sequence will be chosen to represent both.')
    # parser.add_argument('--ds', type=int, help='Allocated disk space Intensiphy can use.')
    # parser.add_argument('-b', default=False, action='store_true', help='Toggles big bactch downloading \
        # of fastq files instead of continuous downloading.')
    return parser.parse_args()

def main() -> None:
    """
    The main function of the script.

    This function parses command-line arguments, calculates the core organization to pass to Extensiphy, 
    checks if the output directory has been made already from a previous run, and handles accession options.

    Returns:
        None
    """

    # Parse command-line arguments
    args: argparse.Namespace = parse_args()

    # Split the absolute path of the script into the directory and the filename
    split_path_and_name: List[str] = os.path.realpath(__file__).rsplit('/',1)

    # Initialize the reference taxon as an empty string
    ref_taxon: str = ""

    # Initialize the absolute path of the alignment file as False
    absolute_align_file_path: Optional[str] = False
    # If an alignment file has been specified
    if not args.align_file == False:
        # Get the absolute path of the alignment file
        absolute_align_file_path = os.path.abspath(args.align_file)

    # Initialize the absolute path of the accessions file as False
    absolute_accs_file_path: Optional[str] = False
    # If an accessions file has been specified
    if not args.accs_file == False:
        # Get the absolute path of the accessions file
        absolute_accs_file_path = os.path.abspath(args.accs_file)

    # Initialize the absolute path of the starting tree file as False
    absolute_tree_file_path: Optional[str] = False
    # If a starting tree file has been specified
    if not args.starting_tree == False:
        # Get the absolute path of the starting tree file
        absolute_tree_file_path = os.path.abspath(args.starting_tree)

    # Print a message indicating that the allocated cores are being assessed
    print("Assessing allocated cores.")
    # Calculate the core organization to pass to Extensiphy
    get_cores: int = calculate_cores(args.cores)

    # Check if output dir has been made already from a previous run
    # If output dir exists, the fundamentals of the program change to suit
    # a repeating run.
    # Returns True if output dir exists already and False if output dir doesn't exist prior to execution.
    print("Checking if this is a new run or a continuing run.")
    dir_existence = check_dir_exists(args.ip_out_dir)

    # Change the current working directory to the output directory
    os.chdir(args.ip_out_dir)
    # Get the absolute path of the output directory
    absolute_output_dir_path: str = os.path.abspath(os.getcwd())

# Define the path to the subprocess command log file
    sub_log_file: str = absolute_output_dir_path + '/run_log_files/subprocess_command_logs.txt'

    # Create the subprocess command log file
    subprocess.run(["touch", sub_log_file])

    # Write the names of the starting alignments to a file
    write_starting_align_names(absolute_output_dir_path, absolute_align_file_path)

    # Handle accession options and get the list of accessions
    accessions: List[str] = handle_accession_options(args.accession_method, args.organism, absolute_output_dir_path, absolute_accs_file_path, sub_log_file)

    #read list of accessions
    read_accessions = ''
    if args.accession_method == 'AUTO_DL':
        read_accessions = read_csv_file(absolute_output_dir_path)

    elif args.accession_method == 'USER_INPUT':
        read_accessions = read_pathodb_csv_file(absolute_output_dir_path)

    print("Working out what kind of run we're doing.")
    if dir_existence == False:

        make_align(dir_existence, absolute_output_dir_path, absolute_align_file_path)

        split_alignment(absolute_output_dir_path + '/intermediate_files/alignment.fas', absolute_output_dir_path + '/sequence_storage')

        # Select the reference taxon
        select_ref(args.ref, absolute_output_dir_path)
    
        # Handle the starting tree
        handle_starting_tree(absolute_output_dir_path, get_cores, absolute_tree_file_path, args.placement)

    else:
        # Print a message indicating that the program is proceeding with a new run using the previous database
        print("Proceeding with new run using previous database.")

    # Clean up any incomplete downloads
    clean_incomplete_downloads(absolute_output_dir_path)

    # Read the names of the fasta files
    read_fasta = read_fasta_names(absolute_output_dir_path)

    # Initialize the variables for removing duplicate paired and single accessions
    remove_paired_dupes = ''
    remove_single_dupes = ''

    # Check the list of run SRA numbers against the sequences already in the alignment to prevent duplicates
    if args.accession_method == 'AUTO_DL':
        remove_paired_dupes = check_duplicate_accesions(read_accessions[0], read_fasta)
        remove_single_dupes = check_duplicate_accesions(read_accessions[1], read_fasta)
    elif args.accession_method == 'USER_INPUT':
        remove_paired_dupes = check_duplicate_accesions(read_accessions, read_fasta)

    # Prepare the batch accessions for paired and single reads
    paired_batch_accessions = prepare_batch_accessions(remove_paired_dupes, get_cores[0])
    single_batch_accessions = prepare_batch_accessions(remove_single_dupes, get_cores[0])

    # Write the names of the current run accessions for paired and single reads
    write_current_run_names(absolute_output_dir_path, paired_batch_accessions)
    write_current_run_names(absolute_output_dir_path, single_batch_accessions)

    # If there are any paired batch accessions
    if len(paired_batch_accessions) > 0:
        # Print a message indicating that the program is processing paired-end read files
        print("Processing paired-end read files.")
        # Download and run the data for the paired batch accessions
        process_data = downloading_and_running(paired_batch_accessions, absolute_output_dir_path, get_cores, "PAIRED", sub_log_file)

    # If there are any single batch accessions
    if len(single_batch_accessions) > 0:
        # Print a message indicating that the program is processing single-end read files
        print("Processing single-end read files.")
        # Download and run the data for the single batch accessions
        process_data = downloading_and_running(single_batch_accessions, absolute_output_dir_path, get_cores, "SINGLE", sub_log_file)

    # Print a message indicating the end of the intensiphy process
    print('end of intensiphy')

    # If the placement option is on
    if args.placement == 'ON':
        # Print a message indicating that the program is placing new sequences into the starting tree
        print("Placing new sequences into starting tree.")
        # Make a new alignment and perform placement into the tree
        construct_align_and_place(absolute_output_dir_path)


def check_dir_exists(dir_name: str) -> bool:
    """
    Check if output directory exists. If the directory exists, a previous run of Intensiphy was probably performed. 
    If the directory doesn't exist, no previous run was performed, build all necessary sub directories.

    Args:
        dir_name (str): The name of the directory to check.

    Returns:
        bool: True if the directory exists, False otherwise.
    """

    # Print a message indicating that the program is checking if the output directory already exists
    print("Checking if output directory already exists.")

    # Check if the directory exists
    does_dir_exist: bool = os.path.isdir(dir_name)

    # Define the list of subdirectories
    ip_folders: List[str] = ['read_files', 'run_log_files', 'intermediate_files', 'accession_files', 'sequence_storage']

    # Print a message indicating whether the output directory exists
    print("Does the output specified output directory exist?: ", does_dir_exist)

    # If the directory exists
    if does_dir_exist:
        # For each subdirectory
        for sub_dir_name in ip_folders:
            # Assert that the subdirectory exists
            assert os.path.isdir(dir_name + '/' + sub_dir_name)

        # Return whether the directory exists
        return does_dir_exist

    # If the directory does not exist
    else:
        # Create the directory
        os.makedirs(dir_name)

        # For each subdirectory
        for sub_dir_name in ip_folders:
            # Create the subdirectory
            os.makedirs(dir_name + '/' + sub_dir_name)

        # Assert that the directory exists
        assert os.path.isdir(dir_name)

        # For each subdirectory
        for sub_dir_name in ip_folders:
            # Assert that the subdirectory exists
            assert os.path.isdir(dir_name + '/' + sub_dir_name)

        # Return whether the directory exists
        return does_dir_exist

if __name__ == '__main__':
    main()
