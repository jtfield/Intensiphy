#!/usr/bin/env python3

import os
# import argparse
import pathlib
import shutil
from numpy.lib.shape_base import split
import numpy as np
import pandas as pd
import subprocess
import datetime
# import dateutil
import re
from modules.alignment_splitter import split_alignment


def make_align(run_bool, output_dir_path, alignment_var):
    """
    Identify and create a symlink to the input alignment file for Intensiphy to use.

    Args:
        run_bool (bool): A boolean indicating whether this is a new run of Intensiphy.
        output_dir_path (str): The path to the output directory.
        alignment_var (str): The path to the input alignment file.

    Returns:
        None
    """

    # Check if an alignment file was provided and if this is not a new run
    if alignment_var != False and run_bool == False:
        # Print a message indicating that an input alignment was detected
        print("Input alignment detected.")

        # Get the absolute path of the input alignment file
        abs_align_path = os.path.realpath(alignment_var)
        
        # Create a Path object for the input alignment file
        symlink_file = pathlib.Path(abs_align_path)

        # Define the path for the new alignment file in the output directory
        new_align = output_dir_path + '/intermediate_files/alignment.fas'
        # Create a Path object for the new alignment file
        new_align = pathlib.Path(new_align)

        # Create a symlink to the input alignment file at the location of the new alignment file
        new_align.symlink_to(symlink_file)

    # Check if no alignment file was provided
    elif alignment_var == False:
        # Print a message indicating that an input alignment is needed
        print("You need to input a starting alignment for Intensiphy to use.")



def select_ref(ref_var, output_dir_path):
    """
    Handles the selection of a reference sequence.

    This function either uses a user-provided reference sequence or automatically selects one.
    If a reference is provided (ref_var is not False), it will be used.
    If no reference is provided (ref_var is False), the function will automatically select the first sequence it finds.

    Args:
        ref_var (str or bool): The user-provided reference sequence or False if no reference is provided.
        output_dir_path (str): The path to the output directory.

    Returns:
        None
    """

    # Define the path to the sequence_storage directory
    sep_file_path = output_dir_path + '/sequence_storage'

    # List all folders in the sequence_storage directory
    individual_seq_folders = os.listdir(sep_file_path)

    # Check if no reference sequence was provided
    if ref_var == False:
        # Select the first folder in the sequence_storage directory as the reference directory
        reference_dir = individual_seq_folders[0]
        # Define the path to the reference directory
        ref_dir_path = sep_file_path + '/' + reference_dir

        # List all files in the reference directory
        ref_dir_files = os.listdir(ref_dir_path)
        for file in ref_dir_files:
            # Check if the file is a fasta file
            if '.fas' in file:
                # Copy the fasta file to the intermediate_files directory and rename it to reference.fas
                shutil.copy(ref_dir_path + '/' + file, output_dir_path + '/intermediate_files/reference.fas')

    # Check if a reference sequence was provided
    elif ref_var != False:
        # Set the provided reference sequence as the reference name
        ref_name = ref_var
        for name in individual_seq_folders:
            # Check if the name matches the reference name
            if name == ref_name:
                ref_dir_path = sep_file_path + '/' + name

                ref_dir_files = os.listdir(ref_dir_path)
                for file in ref_dir_files:
                    # print(file)
                    if '.fas' in file:
                        shutil.copy(ref_dir_path + '/' + file, output_dir_path + '/intermediate_files/reference.fas')


def write_starting_align_names(outdir, starting_align):
    """
    Writes the list of sequence names to the file current_run_starting_taxa.txt.

    Args:
        outdir (str): Path to the output directory for this IP run.
        starting_align (str): Path to the starting alignment file.

    Returns:
        None
    """

    # Initialize an empty list to store the names
    names_list = []

    # Define the name of the file to write to
    file_name = 'current_run_starting_taxa.txt'
    # Construct the full path to the file
    path_to_file = outdir + '/' + file_name

    # Check if the file already exists
    file_exists = os.path.isfile(path_to_file)
    
    # If the file does not exist, create it
    if not file_exists:
        # Open the file in write mode
        open_record_file = open(path_to_file, 'w')

        # Open the starting alignment file
        with open(starting_align) as align_file_contents:
            # Iterate over each line in the file
            for line in align_file_contents:
                # If the line starts with '>', it is a sequence name
                if line.startswith('>'):
                    # Remove leading and trailing whitespace from the line
                    line = line.strip()
                    # Remove the '>' from the start of the line to get the name
                    name = line.strip('>')
                    # Add the name to the list of names
                    names_list.append(name)

        # Iterate over each name in the list of names
        for starting_name in names_list:
            # Write the name to the file, followed by a newline
            open_record_file.write(starting_name + '\n')

        # Close the file
        open_record_file.close()


def write_current_run_names(outdir, list_of_names):
    """
    Records the names of taxa added during the current run of IP.

    This function takes in a list of names and writes them to a file for record keeping.
    If the file already exists, the names are appended to the file.
    If the file does not exist, it is created and the names are written to it.
    The current date is also recorded in the file.

    Args:
        outdir (str): The path to the output directory.
        list_of_names (list): The list of names to be recorded.

    Returns:
        None
    """

    # Get the current date and time
    now = datetime.datetime.now()

    # Format the current date and time as a string in the format 'YYYY-MM-DD'
    current_time = now.strftime('%Y-%m-%d')

    # Define the name of the file to write to
    current_run = "current_run_taxa_added.txt"
    # Construct the full path to the file
    path_to_file = outdir + '/' + current_run

    # Check if the file already exists
    file_exists = os.path.isfile(path_to_file)

    # Initialize a variable to hold the file object
    open_record_file = ''

    # If the file exists, open it in append mode
    if file_exists:
        open_record_file = open(path_to_file, 'a+')

    # If the file does not exist, open it in write mode
    elif file_exists == False:
        open_record_file = open(path_to_file, 'w')

    # Write the current date to the file
    open_record_file.write("\ncurrent_date: " + current_time)

    # Iterate over each sublist in the list of names
    for sub_list in list_of_names:
        # Write a newline to the file
        open_record_file.write("\n")

        # Iterate over each name in the sublist
        for name in sub_list:
            # Write the name to the file
            open_record_file.write(name)

    # Close the file
    open_record_file.close()



def clean_incomplete_downloads(outdir):
    """
    Cleans up the read_files directory by removing any incomplete downloads.

    This function lists all files in the read_files directory. If the directory is empty, 
    it prints a message and returns. If there are files, it checks each one to see if it's a directory. 
    If a file is a directory, the function attempts to remove it using shutil.rmtree.

    Args:
        outdir (str): The path to the output directory.

    Returns:
        None
    """
    read_dir = outdir + '/read_files'

    # List files in read dir
    files_list = os.listdir(read_dir)

    # Check if file_dir is empty
    # if empty, no need to remove anything
    if len(files_list) == 0:
        print("Read directory is empty and ready for new reads.")

    else:

        # Loop over files and check if any are directories
        for file in files_list:
            file_path = read_dir + '/' + file
            check_dir = os.path.isdir(file_path)

            # If file is a directory, try to remove
            # throw an exception if this doesnt work
            if check_dir:
                try:
                    shutil.rmtree(file_path)
                except:
                    print("Could not delete dir :", file_path)



def handle_accession_options(accession_option, organism, folder_path, input_file, log_file):
    """
    Handles the different options for accession input.

    This function processes the accession input based on the option provided by the user.
    If the option is "AUTO_DL", it downloads accessions for the specified organism.
    If the option is "AUTO_PATHDB", it prints a message that the option is not available yet.
    If the option is "USER_INPUT", it copies the user-provided input file to the accession_files directory.

    Args:
        accession_option (str): The accession input option provided by the user.
        organism (str or bool): The organism for which to download accessions, or False if no organism is specified.
        folder_path (str): The path to the folder where the accession files should be stored.
        input_file (str): The path to the user-provided input file (only used if accession_option is "USER_INPUT").
        log_file (str): The path to the log file where output should be written.

    Returns:
        None
    """

    # Print a message indicating that the accession input option is being processed
    print("Processing accession input option.")

    # Check if the accession option is "AUTO_DL"
    if accession_option == "AUTO_DL":
        # Check if an organism was specified
        if not organism == False:
            # Download accessions for the specified organism
            download_accessions(organism, folder_path)

        # Check if no organism was specified
        elif organism == False:
            # Print a message indicating that an organism was not specified
            print("A scientific name or NCBI taxon ID was not input.")
            # Print a message indicating how to specify an organism
            print("Please use the --organism flag to specify what organism to search for.")

    # Check if the accession option is "AUTO_PATHDB"
    elif accession_option == "AUTO_PATHDB":
        # Print a message indicating that this option is not available yet
        print("option not available yet")

    # Check if the accession option is "USER_INPUT"
    elif accession_option == "USER_INPUT":
        # Open the log file in append mode
        open_log_file = open(log_file, 'a')

        # Get the current date and time
        now = datetime.datetime.now()

        # Format the current date and time as a string and use it to construct the output file name
        output_file = 'accessions_' + now.strftime('%Y-%m-%d-%H-%M-%S')

        # Construct the command to copy the input file to the accession_files directory
        cp_command = ['cp', input_file, folder_path + "/accession_files/" + output_file]

        # Run the copy command as a subprocess, redirecting stdout and stderr to the log file
        cp_run = subprocess.Popen(cp_command, stdout=open_log_file, stderr=open_log_file)

        # Wait for the subprocess to finish
        os.wait()


def download_accessions(org_name, out_dir):
    """
    Downloads accession numbers for a specified organism from the SRA database.

    This function uses the esearch and efetch utilities from the NCBI Entrez Direct package 
    to download a CSV file containing accession numbers for the specified organism. 
    The CSV file is saved in the accession_files directory within the specified output directory.

    Args:
        org_name (str): The name of the organism for which to download accession numbers.
        out_dir (str): The path to the output directory.

    Returns:
        None
    """

    # Print a message indicating that the accession number CSV file is being downloaded
    print("Downloading accession number CSV file.")

    # Get the current date and time
    now = datetime.datetime.now()

    # Change the current working directory to the accession_files directory within the output directory
    os.chdir(out_dir + "/accession_files")

    # Construct the output file name using the current date and time
    output_file = 'accessions_' + now.strftime('%Y-%m-%d-%H-%M-%S')

    # Open the output file in write mode
    accession_file = open(output_file, 'w')

    # Define the esearch command
    esearch_command = ['esearch', '-db', 'sra', '-query', org_name]
    # Define the efetch command
    efetch_command = ['efetch', '-format', 'runinfo']

    # Run the esearch command as a subprocess, capturing stdout and stderr
    esearch_accessions = subprocess.Popen(esearch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Run the efetch command as a subprocess, using the stdout of the esearch command as stdin and writing stdout to the output file
    efetch_accessions = subprocess.Popen(efetch_command, stdin=esearch_accessions.stdout, stdout=accession_file)


def read_csv_file(out_dir):
    """
    Reads a CSV file of accession numbers and filters it based on certain criteria.

    This function reads the most recent CSV file in the accession_files directory. 
    It then filters the data to include only rows where LibraryStrategy is 'WGS', 
    LibrarySource is 'GENOMIC', and Platform is 'ILLUMINA'. 
    The filtered data is further divided into two groups based on whether LibraryLayout is 'PAIRED' or 'SINGLE'. 
    The function returns a list containing the two groups of filtered data.

    Args:
        out_dir (str): The path to the output directory.

    Returns:
        list: A list containing two pandas DataFrames. The first DataFrame contains rows where LibraryLayout is 'PAIRED', 
        and the second DataFrame contains rows where LibraryLayout is 'SINGLE'.
    """

    # Print a message indicating that the CSV file is being read
    print("Reading CSV file off accession numbers.")

    # Change the current working directory to the accession_files directory within the output directory
    os.chdir(out_dir + "/accession_files")

    # Find the most recent CSV file in the accession_files directory
    current_accession = find_recent_date(out_dir + "/accession_files")

    # Read the CSV file into a pandas DataFrame
    csv = pd.read_csv(current_accession)

    # Initialize an empty list to store the filtered data
    output = []

    # Filter the DataFrame to include only rows where LibraryStrategy is 'WGS', LibrarySource is 'GENOMIC', and Platform is 'ILLUMINA'
    filtered_df = csv.query("LibraryStrategy == 'WGS' and LibrarySource == 'GENOMIC' and Platform == 'ILLUMINA'")

    # Further filter the DataFrame to include only rows where LibraryLayout is 'PAIRED'
    paired_filtered_df = filtered_df.query("LibraryLayout == 'PAIRED'")
    # Further filter the DataFrame to include only rows where LibraryLayout is 'SINGLE'
    single_filtered_df = filtered_df.query("LibraryLayout == 'SINGLE'")

    # Add the DataFrame of paired reads to the output list
    output.append(paired_filtered_df)
    # Add the DataFrame of single reads to the output list
    output.append(single_filtered_df)

    # Change the current working directory back to the output directory
    os.chdir(out_dir)


def read_pathodb_csv_file(out_dir):
    """
    Reads a CSV file of accession numbers from Pathogen DB and filters it based on certain criteria.

    This function reads the most recent CSV file in the accession_files directory. 
    It then filters the data to include only rows where the 'Run' column is not empty. 
    The function returns a DataFrame containing the filtered data.

    Args:
        out_dir (str): The path to the output directory.

    Returns:
        pandas.DataFrame: A DataFrame containing the filtered data.
    """

    # Print a message indicating that the CSV file is being read
    print("Reading CSV file off accession numbers.")

    # Change the current working directory to the accession_files directory within the output directory
    os.chdir(out_dir + "/accession_files")

    # Find the most recent CSV file in the accession_files directory
    current_accession = find_recent_date(out_dir + "/accession_files")

    # Read the CSV file into a pandas DataFrame
    csv = pd.read_csv(current_accession)

    # Initialize an empty list to store the filtered data
    output = []

    # Replace empty strings in the 'Run' column with NaN
    csv['Run'].replace('', np.nan, inplace=True)

    # Drop rows where the 'Run' column is NaN
    csv.dropna(subset=['Run'], inplace=True)

    # Change the current working directory back to the output directory
    os.chdir(out_dir)

    # Return the filtered DataFrame
    return csv


def find_recent_date(folder_of_files):
    """
    Finds the file with the most recent date in a specified folder.

    This function iterates over all files in the specified folder, extracts the date from each file name 
    using the restructure_dates function, and appends the file name and date to a list. 
    It then converts this list into a pandas DataFrame and prints the DataFrame.

    Args:
        folder_of_files (str): The path to the folder containing the files.

    Returns:
        pandas.DataFrame: A DataFrame containing the file names and their corresponding dates.
    """

    # Print a message indicating that the function is finding the file with the most recent date
    print("Finding the file with the most recent data.")

    # Get the current working directory
    cwd = os.getcwd()
    # Change the current working directory to the specified folder
    os.chdir(folder_of_files)

    # Initialize an empty list to store the file names and their corresponding dates
    list_of_files = []

    # Iterate over all files in the specified folder
    for file in os.listdir(folder_of_files):
        # Extract the date from the file name using the restructure_dates function
        file_date = restructure_dates(file)

        # Initialize an empty list to store the file name and its corresponding date
        data = []
        # Add the file name to the data list
        data.append(file)
        # Add the date to the data list
        data.append(file_date)

        # Add the data list to the list of files
        list_of_files.append(data)

    # Convert the list of files into a pandas DataFrame
    df = pd.DataFrame(list_of_files, columns=['file_name', 'creation_date'])

    # Find the most recent date in the DataFrame
    most_recent_accession_file = df['creation_date'].max()

    # Find the row in the DataFrame that corresponds to the most recent date
    most_recent_row = df.loc[df['creation_date'] == most_recent_accession_file]

    # Extract the file name from the most recent row
    current_accession_file = most_recent_row['file_name'].tolist()[0]

    # Change the current working directory back to the original directory
    os.chdir(cwd)

    # Return the path to the file with the most recent date
    return folder_of_files + '/' + current_accession_file


def restructure_dates(file_name):
    """
    Restructures the date from a file name.

    This function takes a file name, splits it to extract the date, 
    removes any hyphens from the date, and converts the date into an integer. 
    The date is assumed to be in the format "YYYY-MM-DD" and is located after the last "/" in the file name.

    Args:
        file_name (str): The name of the file, including its path.

    Returns:
        int: The date extracted from the file name, restructured as an integer in the format YYYYMMDD.
    """

    # Print a message indicating that the function is restructuring the file creation times
    print("Restructuring file creation times.")

    # Split the file name at the last "/" to separate the path and the file name
    split_path_and_file = file_name.rsplit("/", 1)

    # Extract the date from the file name
    condensed_date = split_path_and_file[0]

    # Remove any hyphens from the date
    condensed_date = split_path_and_file[0].replace("-","")

    # Split the date at "_" to separate the year, month, and day
    condensed_date = condensed_date.split("_")

    # Convert the date into an integer
    condensed_date = int(condensed_date[1])

    # Return the date
    return condensed_date


def calculate_cores(set_cores):
    """
    Calculates the number of runs and cores per run based on the total number of cores.

    Args:
        set_cores (int): The total number of cores.

    Returns:
        list: A list containing two integers. The first integer is the number of runs, and the second integer is the number of cores per run.
    """    

    # Print a message indicating that the function is calculating threads for Extensiphy
    print("Calculating threads for Extensiphy.")

    # Initialize an empty list to store the number of runs and cores per run
    output = []

    # Convert the total number of cores to an integer
    total_cores = int(set_cores)

    # Initialize variables for the number of runs and cores per run
    runs = 0
    cores_per_run = 0

    # If the total number of cores is between 2 and 10 (inclusive of 2 and exclusive of 10)
    if total_cores >= 2 and total_cores < 10:
        # Set the number of runs to 2
        runs = 2
        # Set the number of cores per run to half the total number of cores
        cores_per_run = total_cores / 2

    # If the total number of cores is 10 or more
    elif total_cores >= 10:
        # Set the number of runs to half the total number of cores
        runs = total_cores / 2
        # Set the number of cores per run to 2
        cores_per_run = 2

    # Assert that the total number of cores is not exceeded by the product of the number of runs and cores per run
    assert runs * cores_per_run <= total_cores

    # Add the number of runs to the output list
    output.append(int(runs))
    # Add the number of cores per run to the output list
    output.append(int(cores_per_run))

    # Assert that the output is a list
    assert type(output) == list
    # Assert that the output list contains exactly two items
    assert len(output) == 2
    # Assert that both items in the output list are integers
    for num in output:
        assert type(num) == int

    # Return the output list
    return output

def check_duplicate_accesions(accession_db, fasta_names):
    """
    Checks for duplicate accession numbers in the alignment.

    This function takes a DataFrame of accession numbers and a list of names from a FASTA file. 
    It converts the 'Run' column of the DataFrame to a list and checks each name in the list against the names in the FASTA file. 
    If a name from the 'Run' list is not found in the FASTA names, it is appended to a new list. 
    The function returns this new list, which contains only the names that are not already present in the FASTA file.

    Args:
        accession_db (pandas.DataFrame): A DataFrame containing accession numbers.
        fasta_names (list): A list of names from a FASTA file.

    Returns:
        list: A list of accession numbers that are not already present in the FASTA file.
    """    

    # Print a message indicating that the function is checking for accessions already present in alignment
    print("Checking for accessions already present in alignment")

    # Initialize an empty list to store the accession numbers to keep
    sras_to_keep = []

    # Print the 'Run' column of the DataFrame
    print(accession_db['Run'])

    # Convert the 'Run' column of the DataFrame to a list
    sra_ids = accession_db['Run'].tolist()

    # If the list of accession numbers is not empty
    if len(sra_ids) != 0:
        # Iterate over each accession number in the list
        for num, name in enumerate(sra_ids):
            # If the accession number is not in the list of FASTA names
            if name not in fasta_names:
                # Append the accession number to the list of accession numbers to keep
                sras_to_keep.append(name)

    # Return the list of accession numbers to keep
    return sras_to_keep


def prepare_batch_accessions(accessions, runs):
    """
    Prepares batches of accession numbers for downloading.

    This function takes a list of accession numbers and a number of runs. 
    It divides the list of accession numbers into chunks, with each chunk containing a number of accession numbers equal to the number of runs. 
    The function returns a list of these chunks.

    Args:
        accessions (list): A list of accession numbers.
        runs (int): The number of runs, which determines the size of each chunk.

    Returns:
        list: A list of lists, where each inner list is a chunk of accession numbers.
    """

    # Print a message indicating that the function is building a batch accession numbers download plan
    print("building batch accession numbers download plan.")

    # Divide the list of accession numbers into chunks, with each chunk containing a number of accession numbers equal to the number of runs
    chunks = [accessions[x:x+runs] for x in range(0, len(accessions), runs)]

    # Return the list of chunks
    return chunks


def read_fasta_names(_outdir):
    """
    Reads the names of individual sequence files from the sequence_storage directory.

    This function lists all files in the sequence_storage directory and appends each file name to a list. 
    It prints each file name as it is added to the list. The function returns the list of file names.

    Args:
        _outdir (str): The path to the output directory, which contains the sequence_storage directory.

    Returns:
        list: A list of file names from the sequence_storage directory.
    """

    print("Getting taxa labels from the current alignment file.")
    names = []

    sep_file_path = _outdir + '/sequence_storage'

    individual_seq_folders = os.listdir(sep_file_path)

    for name in individual_seq_folders:
        names.append(name)
        print(name)

    return names


def align_rename_and_move(align, outdir):
    """
    Moves and renames the alignment file to correspond to its construction time.

    This function takes the path to an alignment file and the path to an output directory. 
    It creates a symlink to the alignment file in the output directory, 
    and renames the symlink to include the current date and time. 
    The function returns the path to the new symlink.

    Args:
        align (str): The path to the alignment file.
        outdir (str): The path to the output directory.

    Returns:
        str: The path to the new symlink in the output directory.
    """

    print("Moving and renaming alignment to currespond to its construction time.")

    now = datetime.datetime.now()

    # Symlink the alignment file into the alignments folder
    abs_align_path = os.path.realpath(align)

    new_align_name = outdir + '/msa_' + now.strftime('%Y-%m-%d-%H-%M-%S')

    new_align_path = pathlib.Path(new_align_name)
    print(abs_align_path)
    print(new_align_path)

    shutil.copy(abs_align_path, new_align_name)

    return new_align_name


def fasterq_dump_reads(out_dir_, single_accession_):
    """
    Downloads reads for a single accession number using the fasterq-dump utility.

    This function changes the current working directory to the read_files directory within the specified output directory. 
    It then uses the fasterq-dump utility to download reads for the specified accession number. 
    If an error occurs during the download, the function prints an error message and sets a variable to the accession number. 
    If the download is successful, the function prints a success message and sets the variable to "DOWNLOADED". 
    The function returns the current working directory to the output directory and prints the value of the variable.

    Args:
        out_dir_ (str): The path to the output directory.
        single_accession_ (str): The accession number for which to download reads.

    Returns:
        None
    """

    os.chdir(out_dir_ + "/read_files")

    downloaded_or_not = 'UNUSED'
    print(downloaded_or_not)

    reads_dl = subprocess.Popen(["fasterq-dump", "--split-files", single_accession_], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print(reads_dl.communicate())
    errs = reads_dl.stderr.read().decode()
    out = reads_dl.stdout.read().decode()

    if "err" in errs or "err" in out:
        print("not downloaded")
        print(out)
        print(errs)
        downloaded_or_not = single_accession_

    else:
        print("downloaded")

        downloaded_or_not = "DOWNLOADED"

    print(downloaded_or_not)
    os.chdir(out_dir_)

    return downloaded_or_not


def downloading_and_running(accessions, out_dir, cores, pair_or_not_toggle, log_file):
    """
    Downloads and aligns reads for a batch of accession numbers.

    This function takes a list of batches of accession numbers, an output directory, a number of cores, 
    a toggle for whether the reads are paired or not, and a log file. 
    It sets the paths to the reference file and the alignment output file, and opens the log file in append mode. 
    For each batch of accession numbers, it prints the batch and its length, and initializes a list for accession numbers that fail to download. 
    For each accession number in the batch, it downloads the reads using the fasterq_dump_reads function. 
    If the download fails, the accession number is added to the list of failed downloads. 

    Args:
        accessions (list): A list of batches of accession numbers.
        out_dir (str): The path to the output directory.
        cores (int): The number of cores to use for alignment.
        pair_or_not_toggle (bool): A toggle for whether the reads are paired or not.
        log_file (str): The path to the log file.

    Returns:
        None
    """

    ref = out_dir +'/intermediate_files/reference.fas'
    ep_output_align = out_dir + '/intermediate_files/ep_output/RESULTS/extended.aln'
    now = datetime.datetime.now()
    
    open_log_file = open(log_file, 'a')

    for accession_batch in accessions:
        print("Batch of current accessions is ", accession_batch)
        len_of_this_batch = len(accession_batch)
        not_downloaded_this_run = []

        for num, accession in enumerate(accession_batch):
            download_status = fasterq_dump_reads(out_dir, accession)

            if download_status != "DOWNLOADED":

                not_downloaded_this_run.append(accession)
            
            # TODO: add renaming for single end read files here to accomodate extensiphy settings
            # TODO: or rework extensiphy single end read file processing

        if len(not_downloaded_this_run) <= (len_of_this_batch / 2):
            print("Batch contains at least some accessions that were downloaded.")
            print("Proceeding to EP.")
        # If the number of failed downloads is less than the number of accessions being downloaded, skip the run

            if pair_or_not_toggle == "PAIRED":

                ep_command = ["extensiphy.sh", "-a", ref, "-d", out_dir + "/read_files", "-i", "CLEAN", "-p", str(cores[0]) ,"-c", str(cores[1]), "-1", "_1.fastq", "-2", "_2.fastq", "-o", out_dir + '/intermediate_files/ep_output']

                run_ep = subprocess.Popen(ep_command, stdout=open_log_file, stderr=open_log_file)

                os.wait()

                # subprocess.run(["extensiphy.sh", "-a", ref, "-d", out_dir + "/read_files", "-i", "CLEAN", "-p", str(cores[0]) ,"-c", str(cores[1]), "-1", "_1.fastq", "-2", "_2.fastq", "-o", out_dir + '/intermediate_files/ep_output'])

            elif pair_or_not_toggle == "SINGLE":

                # TODO: understand if we're getting seg faults running extensiphy for single end reads.
                ep_command = ["extensiphy.sh", "-a", ref, "-d", out_dir + "/read_files", "-e", "SE", "-i", "CLEAN", "-p", str(cores[0]) ,"-c", str(cores[1]), "-1", "_1.fastq", "-o", out_dir + '/intermediate_files/ep_output']

                run_ep = subprocess.Popen(ep_command, stdout=open_log_file, stderr=open_log_file)

                os.wait()

                # subprocess.run(["extensiphy.sh", "-a", ref, "-d", out_dir + "/read_files", "-e", "SE", "-i", "CLEAN", "-p", str(cores[0]) ,"-c", str(cores[1]), "-1", "_1.fastq", "-o", out_dir + '/intermediate_files/ep_output'])

            print("Extensiphy runs complete.")
            print("Cleaning up files for next batch.")
            try:
                split_alignment(ep_output_align, out_dir + '/sequence_storage')

            except FileNotFoundError:
                print("EP output alignment not found. Moving dev files to logs directory.")
                ep_log_file = out_dir + '/intermediate_files/ep_output/ep_dev_log.txt'
                new_log_name = out_dir + '/run_log_files/ep_log_' + now.strftime('%Y-%m-%d-%H-%M-%S')
                shutil.copyfile(ep_log_file, new_log_name)
                shutil.rmtree(out_dir + '/intermediate_files/ep_output')

            print("Removing used read files.")
            rm_read_files(out_dir + '/read_files')

            shutil.rmtree(out_dir + '/intermediate_files/ep_output')

        elif len(not_downloaded_this_run) == len_of_this_batch:

            print("ERRORS ON DOWNLOAD: Skipping EP run due to un-downloaded accessions batch.")
            print(accession_batch)
            print("batches not downloaded this run.")
            print(not_downloaded_this_run)
            print("###")




def check_downloaded_read_files(outdir_, read_format_toggle):
    """
    Checks the downloaded read files in the read_files directory.

    This function takes an output directory and a toggle for the read format. 
    It lists all files in the read_files directory within the output directory and compiles regular expressions 
    for the accepted formats for paired SRR and ERR files. 
    The function does not return a value, but it is implied that it performs some action based on the read format toggle and the files in the directory.

    Args:
        outdir_ (str): The path to the output directory.
        read_format_toggle (bool): A toggle for the read format.

    Returns:
        None
    """

    # establish read directory path
    read_dir = outdir_ + "/read_files"

    list_of_files = os.listdir(read_dir)

    one_paired_accepted_srr_format = "SRR\w+_1.fastq"
    two_paired_accepted_srr_format = "SRR\w+_2.fastq"

    one_paired_accepted_err_format = "ERR\w+_1.fastq"
    two_paired_accepted_err_format = "ERR\w+_2.fastq"

    one_paired_compile_srr = re.compile(one_paired_accepted_srr_format)
    one_paired_compile_err = re.compile(one_paired_accepted_err_format)
    two_paired_compile_srr = re.compile(two_paired_accepted_srr_format)
    two_paired_compile_err = re.compile(two_paired_accepted_err_format)

    single_accepted_srr_format = "SRR\w+_1.fastq"
    single_accepted_err_format = "ERR\w+_1.fastq"


    single_compile_srr = re.compile(single_accepted_srr_format)
    single_compile_err = re.compile(single_accepted_err_format)

    if read_format_toggle == "PAIRED":
        for file in list_of_files:
            check_one_paired_srr = re.match(one_paired_compile_srr, file)
            check_one_paired_err = re.match(one_paired_compile_err, file)
            check_one_paired_srr = re.match(one_paired_compile_srr, file)
            check_one_paired_err = re.match(one_paired_compile_err, file)
            if check_one_paired_err:
                print("Paired ERR")
                print(check_one_paired_err[0])

                paired_read_set = check_one_paired_err[0].replace('_1','_2')
                print(paired_read_set)

def rm_read_files(read_dir):
    """
    Removes all read files in the specified directory in preparation for the next run.

    This function takes a directory and lists all files in it. 
    It then iterates over each file, constructs the path to the file, and removes the file. 
    The function does not return a value, but it prints a message indicating that it is removing the read files.

    Args:
        read_dir (str): The path to the directory containing the read files.

    Returns:
        None
    """

    print("removing read files in preparation for next run.")
    fastq_file_list = os.listdir(read_dir)
    # print(fastq_file_list)
    for fq_file in fastq_file_list:
        path_to_file = os.path.join(read_dir, fq_file)
        os.remove(path_to_file)

def bulk_download(accession_list, out_dir):
    """
    Downloads bulk data for a list of accession numbers.

    This function takes a list of accession numbers and an output directory. 
    It changes the current working directory to the read_files directory within the output directory. 
    For each accession number in the list, it prints the accession number and uses the fasterq-dump utility to download the data. 
    The function does not return a value, but it prints a message indicating that the bulk data download scheme has begun.

    Args:
        accession_list (list): A list of accession numbers for which to download data.
        out_dir (str): The path to the output directory.

    Returns:
        None
    """

    print("Bulk data download scheme begun.")
    os.chdir(out_dir + "/read_files")

    for single_accession in accession_list:

        print(single_accession)
        subprocess.run(["fasterq-dump", "--split-files", single_accession])
