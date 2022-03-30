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
import re
from modules.alignment_splitter import split_alignment


def make_align(run_bool, output_dir_path, input_accessions, alignment_var):
    """If this is the start of a run (no output folder existed prior to starting this run) \
    used gon_phyling to build an alignment based on the first 6 available samples"""
    print("Building a starting alignment.")
    if alignment_var != False:
        print("Input alignment detected.")
        print("Skipping de novo alignment production.")

        abs_align_path = os.path.realpath(alignment_var)
        #
        symlink_file = pathlib.Path(abs_align_path)

        new_align = output_dir_path + '/intermediate_files/alignment.fas'
        new_align = pathlib.Path(new_align)

        new_align.symlink_to(symlink_file)


    if alignment_var == False:
        #Test!
        # tests.assembly_tests.gon_phy_test.test_build_alignment()

        if run_bool == False:

            print("make alignment with gon_phyling")

            accessions_files = os.listdir(output_dir_path + "/accession_files")

            num_files = len(accessions_files)
            print(accessions_files)

            if num_files == 1:
                accession_file = accessions_files[0]

                accessions = pd.read_csv(output_dir_path + "/accession_files/" + accession_file)

                #Drop rows with empty values in the Run column because we only want rows with SRA numbers associated
                accessions.dropna(subset = ['Run'], inplace = True)

                indexSingleEnd = accessions[accessions['LibraryLayout'] == 'SINGLE'].index

                accessions.drop(indexSingleEnd, inplace = True)

                os.mkdir(output_dir_path + "/starting_align_files")
                os.chdir(output_dir_path + "/starting_align_files")

                #TODO: Randomize sequence selection for potentially less biased loci selection
                for num in range(0,5):
                    sra_num = accessions.loc[num,'Run']
                    subprocess.run(["fasterq-dump", "--split-files", sra_num])

                    subprocess.run(["gon_phyling.sh", "-d", output_dir_path + "/starting_align_files", "-1", "_1.fastq", "-2", "_2.fastq"])

                    os.chdir(output_dir_path + "/starting_align_files")

                    symlink_align(output_dir_path + "/starting_align_files/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas", output_dir_path)

# def get_most_recent_align(run_bool, starting_align, output_dir_path):
#     """Determines how to handle the run if the output dir already exists (indicating a continuing run) \
#         or if the output dir doesn't exist yet (indicating the starting phase of a run)."""
#     print("Processing alignment input options.")
#     if starting_align == False:
#
#         current_alignment = find_recent_date(output_dir_path + "/alignments")
#
#         return current_alignment
#
#     else:
#         if starting_align != False:
#
#             # current_alignment = symlink_align(starting_align, output_dir_path)
#             current_alignment = align_rename_and_move(starting_align, output_dir_path + "/alignments")
#
#             return current_alignment

def select_ref(ref_var, output_dir_path):
    """Small function to handle a user input reference or making the choice of reference automatically"""
    sep_file_path = output_dir_path + '/sequence_storage'

    individual_seq_folders = os.listdir(sep_file_path)

    if ref_var == False:
        reference_dir = individual_seq_folders[0]
        ref_dir_path = sep_file_path + '/' + reference_dir

        ref_dir_files = os.listdir(ref_dir_path)
        for file in ref_dir_files:
            # print(file)
            if '.fas' in file:
                shutil.copy(ref_dir_path + '/' + file, output_dir_path + '/intermediate_files/reference.fas')

    elif ref_var != False:

        ref_name = ref_var
        for name in individual_seq_folders:
            if name == ref_name:
                ref_dir_path = sep_file_path + '/' + name

                ref_dir_files = os.listdir(ref_dir_path)
                for file in ref_dir_files:
                    # print(file)
                    if '.fas' in file:
                        shutil.copy(ref_dir_path + '/' + file, output_dir_path + '/intermediate_files/reference.fas')


def write_current_run_names(outdir, list_of_names):
    """
    Takes in the list of names being added during the current run of IP.
    Writes those to a file for record keeping.
    """
    now = datetime.datetime.now()

    current_time = now.strftime('%Y-%m-%d')

    current_run = "current_run_taxa_added.txt"
    path_to_file = outdir + '/' + current_run

    file_exists = os.path.isfile(path_to_file)

    open_record_file = ''

    if file_exists:

        open_record_file = open(path_to_file, 'a+')

    elif file_exists == False:

        open_record_file = open(path_to_file, 'w')

    open_record_file.write("\ncurrent_date: " + current_time)

    for sub_list in list_of_names:

        for name in sub_list:
            open_record_file.write("\n")
            open_record_file.write(name)

    open_record_file.close()

    # elif file_exists == False:
    #
    #     open_record_file = open(path_to_file, 'w')
    #
    #     open_record_file.write("\ncurrent_date: " + current_time)
    #
    #     for name in list_of_names:
    #         open_record_file.write("\n")
    #         open_record_file.write(name)
    #
    #     open_record_file.close()






def clean_incomplete_downloads(outdir):
    """Delete folders left behind by fasterq-dump if a download is interrupted"""
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






def restructure_dates(file_name):
    """Reads a file name with a date structured as Intensiphy produces them and restructures the included date for parsing of most recent file"""
    print("Restructuring file creation times.")
    split_path_and_file = file_name.rsplit("/", 1)
    #print(split_path_and_file)
    condensed_date = split_path_and_file[0]
    #print(condensed_date)
    condensed_date = split_path_and_file[0].replace("-","")
    #print(condensed_date)
    condensed_date = condensed_date.split("_")
    #print(condensed_date)
    condensed_date = int(condensed_date[1])
    return condensed_date


def handle_accession_options(accession_option, organism, folder_path, input_file):
    """Only for use when not automatically downloading SRA numbers"""
    print("Processing accession input option.")

    if accession_option == "AUTO_DL":
        download_accessions(organism, folder_path)

    elif accession_option == "AUTO_PATHDB":
        print("option not available yet")

    elif accession_option == "USER_INPUT":
        subprocess.run(['cp', input_file, folder_path + "/accession_files/"])


def download_accessions(org_name, out_dir):
    """Download the run info file that includes run ID accession numbers and info on how the sequences were produced."""
    print("Downloading accession number CSV file.")
    # Test!
    # tests.accessions_tests.accession_download_test.test_download_accessions()

    now = datetime.datetime.now()

    os.chdir(out_dir + "/accession_files")

    # ' '.join(org_name)
    org_name = org_name[0].strip("'")

    output_file = 'accessions_' + now.strftime('%Y-%m-%d-%H-%M-%S')

    url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=' + org_name

    dl_accessions = subprocess.Popen(['wget', '-O', output_file,  url], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(dl_accessions.communicate())

    os.chdir(out_dir)


def read_csv_file(out_dir):
    """Reads the accession file provided by the users and parses compatible sequences."""
    print("Reading CSV file off accession numbers.")
    os.chdir(out_dir + "/accession_files")

    current_accession = find_recent_date(out_dir + "/accession_files")

    csv = pd.read_csv(current_accession)

    output = []

    filtered_df = csv.query("LibraryStrategy == 'WGS' and LibrarySource == 'GENOMIC' and Platform == 'ILLUMINA'")

    paired_filtered_df = filtered_df.query("LibraryLayout == 'PAIRED'")
    single_filtered_df = filtered_df.query("LibraryLayout == 'SINGLE'")

    output.append(paired_filtered_df)
    output.append(single_filtered_df)

    # print(output)

    os.chdir(out_dir)

    return output

def find_recent_date(folder_of_files):
    print("Finding the file with the most recent data.")

    cwd = os.getcwd()
    os.chdir(folder_of_files)

    list_of_files = []

    for file in os.listdir(folder_of_files):
        file_date = restructure_dates(file)
        data = []
        data.append(file)
        # split_file = file.split('_')
        # date = split_file[1]
        # parsed_date = dateutil.parser.parse(date)
        # print(parsed_date)
        # data.append(date)
        data.append(file_date)
        # created_time = os.stat(file).st_ctime
        # data.append(created_time)
        list_of_files.append(data)

    print(list_of_files)

    df = pd.DataFrame(list_of_files, columns=['file_name', 'creation_date'])
    print(df)

    # df['creation_date'] = pd.to_datetime(df['creation_date'])

    most_recent_accession_file = df['creation_date'].max()

    most_recent_row = df.loc[df['creation_date'] == most_recent_accession_file]

    current_accession_file = most_recent_row['file_name'].tolist()[0]

    os.chdir(cwd)

    return folder_of_files + '/' + current_accession_file

def calulate_cores(set_cores):
    """Organizes and calulates the cores for use with Extensiphy."""
    print("Calculating threads for Extensiphy.")
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
    print("Checking for accessions already present in alignment")
    sras_to_keep = []

    sra_ids = accession_db['Run'].tolist()
    print("Starting list of sras ", sra_ids)
    # print(fasta_names)
    if len(sra_ids) != 0:
        for num, name in enumerate(sra_ids):
            if name not in fasta_names:
                sras_to_keep.append(name)

    print("Dupes removed list of sras ", sras_to_keep)
    return sras_to_keep


    # copy_df = accession_db.copy()
    #
    # for name, value in accession_db['Run'].iteritems():
    #     if value in fasta_names:
    #         assert type(value) == str
    #         assert len(value) > 1
    #         copy_df.drop(name)
    #         print(name)
    #         print(value)
    #
    # print(copy_df)
    #
    # return copy_df


def prepare_batch_accessions(accessions, runs):
    """Return a list of lists containing the number of files to download before each Extensiphy run"""
    print("building batch accession numbers download plan.")
    # run_ids = []
    #
    # for name, value in accessions['Run'].iteritems():
    #     run_ids.append(value)

    chunks = [accessions[x:x+runs] for x in range(0, len(accessions), runs)]

    return chunks


def read_fasta_names(_outdir):
    """Reads the names of the sequences from the fasta file."""
    print("Getting taxa labels from the current alignment file.")
    names = []

    sep_file_path = _outdir + '/sequence_storage'

    individual_seq_folders = os.listdir(sep_file_path)

    for name in individual_seq_folders:
        names.append(name)

    return names


def align_rename_and_move(align, outdir):
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

# def prepare_batch_accessions(accessions, runs):
#     """Return a list of lists containing the number of files to download before each Extensiphy run"""
#     print("building batch accession numbers download plan.")
#     chunks = [accessions[x:x+runs] for x in range(0, len(accessions), runs)]
#
#     return chunks

# def download_chunk(out_dir, accessions, ds_alloc):
#     """Downloads raw read files and checks if the allocated disk space has been met \
#     if not, download files until this limit is met"""
#     total_memory_size = 0
#     accession_count = 0
#     run_ids = []
#     for name, value in accessions['Run'].iteritems():
#         run_ids.append(value)
#
#     while total_memory_size < ds_alloc and accession_count < len(accessions):
#         single_accession = run_ids[accession_count]
#         print(single_accession)




def fasterq_dump_reads(out_dir_, single_accession_):
    os.chdir(out_dir_ + "/read_files")
    # print("current working directory: ", os.getcwd())
    # print("current alignment ", current_alignment)
    reads_dl = subprocess.Popen(["fasterq-dump", "--split-files", single_accession_], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(reads_dl.communicate())
    os.chdir(out_dir_)


def downloading_and_running(accessions, out_dir, cores, pair_or_not_toggle):
    """Take batched accessions, download for each batch and run EP. Then split up the alignment and remove the original folder"""
    ref = out_dir +'/intermediate_files/reference.fas'
    ep_output_align = out_dir + '/intermediate_files/ep_output/RESULTS/extended.aln'

    for accession_batch in accessions:
        # print(accession_batch)
        print("++++")
        print("Batch of current accessions is ", accession_batch)
        for num, accession in enumerate(accession_batch):
            print("+++++++")
            # print(accession)
            fasterq_dump_reads(out_dir, accession)

        if pair_or_not_toggle == "PAIRED":

            # print("/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-a", ref, "-d", out_dir + "/read_files", "-i", "CLEAN", "-p", str(cores[0]) ,"-c", str(cores[1]), "-1", "_1.fastq", "-2", "_2.fastq", "-o", out_dir + '/intermediate_files/ep_output')
            subprocess.run(["extensiphy.sh", "-a", ref, "-d", out_dir + "/read_files", "-i", "CLEAN", "-p", str(cores[0]) ,"-c", str(cores[1]), "-1", "_1.fastq", "-2", "_2.fastq", "-o", out_dir + '/intermediate_files/ep_output'])
            # print(print("/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-a", ref, "-d", out_dir + "/read_files", "-i", "CLEAN", "-p", str(cores[0]) ,"-c", str(cores[1]), "-1", "_1.fastq", "-2", "_2.fastq", "-o", out_dir + '/intermediate_files/ep_output'))

        elif pair_or_not_toggle == "SINGLE":

            if pair_or_not_toggle == "PAIRED":

                # print("/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-a", ref, "-d", out_dir + "/read_files", "-i", "CLEAN", "-p", str(cores[0]) ,"-c", str(cores[1]), "-1", "_1.fastq", "-2", "_2.fastq", "-o", out_dir + '/intermediate_files/ep_output')
                subprocess.run(["extensiphy.sh", "-a", ref, "-d", out_dir + "/read_files", "-e", "SE", "-i", "CLEAN", "-p", str(cores[0]) ,"-c", str(cores[1]), "-1", "_1.fastq", "-o", out_dir + '/intermediate_files/ep_output'])
                # print(print("/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-a", ref, "-d", out_dir + "/read_files", "-i", "CLEAN", "-p", str(cores[0]) ,"-c", str(cores[1]), "-1", "_1.fastq", "-2", "_2.fastq", "-o", out_dir + '/intermediate_files/ep_output'))

        split_alignment(ep_output_align, out_dir + '/sequence_storage')

        rm_read_files(out_dir + '/read_files')

        shutil.rmtree(out_dir + '/intermediate_files/ep_output')

        # TEMPORARY BREAK STATEMENT
        # TODO: DONT FORGET TO REMOVE THIS!!!!
        # TESTING PURPOSES ONLY
        # break





    # # Identify if we're processing data by downloading in batches between Extesniphy runs
    # # or by downloading everything all at once and running extensiphy after
    # print("Beginning bulk downloading data and updating alignment using Extensiphy.")
    # # pull accessions from df
    # run_ids = []
    # for name, value in accessions['Run'].iteritems():
    #     run_ids.append(value)
    #
    # if method == False:
    #     # continuous gradual downloading of data has been selected/left as default.
    #     batches_run = batch_download(run_ids, run_num[0], out_dir, align)
    #
    # elif method:
    #     # Bulk download all fastq files before running Extensiphy
    #     bulk_run = bulk_download(run_ids, out_dir)
    #
    #     print("Bulk data download complete.")
    #     print("Beginning Extensiphy run to update your alignment.")
    #
    #     os.chdir(out_dir)
    #
    #     # TODO: specify cleanup of intermediate files in all runs
    #     subprocess.run(["multi_map.sh", "-a", align, "-d", out_dir + "/read_files", "-i", "CLEAN", "-p", str(run_num[0]) ,"-c", str(run_num[1]), "-1", "_1.fastq", "-2", "_2.fastq" ])

# def batch_download(accession_list, runs_number, out_dir, alignment):
#     """Prepares accessions to be downloaded in batches between runs of Extensiphy."""
#     batches_of_accessions = prepare_batch_accessions(accession_list, runs_number)
#     print("Beginning batch data download and run scheme.")
#     print("^^^^^^^^^^^^^^^^^^^")
#
#     read_dir_path = out_dir + "/read_files"
#     ep_output = out_dir + "/ep_tmp_outputs"
#     standard_align_path = ep_output + "/RESULTS/extended.aln"
#     align_outdir = out_dir + "/alignments"
#     current_alignment = alignment
#     ref_taxon = ""
#
#     os.chdir(read_dir_path)
#
# ########################################################################
# # testing EP running on a single accession instead of in a loop
#     # grouped_accession = batches_of_accessions[0]
#     # print(grouped_accession)
#     # single_accession = grouped_accession[0]
#     #
#     # fasterq_dump_reads(out_dir, single_accession)
#     #
#     # run_ep(out_dir, current_alignment, ep_output)
#     #
#     # print("DONE")
#
# ########################################################################
#
#
#
#     for num, accessions in enumerate(batches_of_accessions):
#         print("BEGINNING NEW BATCH SECTION.")
#
#         print("Currend Accession Batch ", accessions)
#
#         if num <= 2:
#
#             # os.chdir(out_dir + "/read_files")
#             print("current working directory: ", os.getcwd())
#             # print(accessions)
#             for single_accession in accessions:
#                 print("BEGINNING EP LOOP.")
#
#                 fasterq_dump_reads(out_dir, single_accession)
#
#                 run_ep(out_dir, current_alignment, ep_output)
#
#
#                 exit()
#
#                 check_ep_output = os.path.isdir(ep_output)
#                 if check_ep_output:
#
#                     print("Found EP output directory")
#                     print("current alignment ", current_alignment)
#                     print("current working directory: ", os.getcwd())
#                     # Copy and rename output alignment and put it the alignments folder for the next round
#                     align_rename_and_move(standard_align_path, align_outdir)
#
#                     print("current working directory: ", os.getcwd())
#                     current_alignment = find_recent_date(align_outdir)
#                     print("curent alignment", current_alignment)
#
#
#                     print("current working directory: ", os.getcwd())
#                     os.rename(ep_output, ep_output + "_" + str(num))
#
#                     print("current working directory: ", os.getcwd())
#                     rm_read_files(read_dir_path)
#
#                     # Remove Extensiphy output directory and contents and prepare for next phase of loop
#                     # try:
#                     #     shutil.rmtree(ep_output)
#                     # except OSError as e:
#                     #     print("Error: %s : %s" % (ep_output, e.strerror))
#
#                 elif check_ep_output == False:
#                     print("ERROR: no EP output directory was found.")
#                     exit()
#
#                 print("current working directory: ", os.getcwd())
#                 print("STARTING NEXT EP LOOP")

def fasterq_dump_reads(out_dir_, single_accession_):
    os.chdir(out_dir_ + "/read_files")
    # print("current working directory: ", os.getcwd())
    # print("current alignment ", current_alignment)
    reads_dl = subprocess.Popen(["fasterq-dump", "--split-files", single_accession_], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(reads_dl.communicate())
    os.chdir(out_dir_)

def run_ep(base_dir_, align_, outdir_):
    # print("current working directory: ", os.getcwd())
    # print("current alignment ", current_alignment)
    print("/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-a", align_, "-d", base_dir_ + "/read_files", "-1", "_1.fastq", "-2", "_2.fastq", "-o", outdir_)
    # ep_process = subprocess.Popen(["/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-a", align_, "-d", base_dir_ + "/read_files", "-1", "_1.fastq", "-2", "_2.fastq", "-o", outdir_], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print(ep_process.communicate())
    subprocess.run(["/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-a", align_, "-d", base_dir_ + "/read_files", "-1", "_1.fastq", "-2", "_2.fastq", "-o", outdir_])

    # print("current working directory: ", os.getcwd())
    # print("current alignment ", current_alignment)
    # # subprocess.run(["/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-h"])
    # print("Extensiphy run complete. BEGINNING PROCESSING")

def rm_read_files(read_dir):
    """Quick function to remove reads that have been used"""
    print("removing read files in preparation for next run.")
    fastq_file_list = os.listdir(read_dir)
    # print(fastq_file_list)
    for fq_file in fastq_file_list:
        path_to_file = os.path.join(read_dir, fq_file)
        os.remove(path_to_file)

def bulk_download(accession_list, out_dir):
    """Bulk download every sequence (unlike the batch method)"""
    print("Bulk data download scheme begun.")
    os.chdir(out_dir + "/read_files")

    for single_accession in accession_list:

        print(single_accession)
        subprocess.run(["fasterq-dump", "--split-files", single_accession])
