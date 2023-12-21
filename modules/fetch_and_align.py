#!/usr/bin/env python3

import os
import argparse
import pathlib
import shutil
from numpy.lib.shape_base import split
import numpy as np
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
    if alignment_var != False and run_bool == False:
        print("Input alignment detected.")
        print("Skipping de novo alignment production.")

        abs_align_path = os.path.realpath(alignment_var)
        #
        symlink_file = pathlib.Path(abs_align_path)

        new_align = output_dir_path + '/intermediate_files/alignment.fas'
        new_align = pathlib.Path(new_align)

        new_align.symlink_to(symlink_file)


    elif alignment_var == False:

        print("You need to input a starting alignment for Intensiphy to use.")



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


def write_starting_align_names(outdir, starting_align):
    """
    Writes the list of sequence names to the file current_run_starting_taxa.txt.

    Args:
        outdir (str): Path to the output directory for this IP run.
        starting_align (str): Path to the starting alignment file.
    """

    names_list = []

    file_name = 'current_run_starting_taxa.txt'
    path_to_file = outdir + '/' + file_name

    file_exists = os.path.isfile(path_to_file)
    
    if not file_exists:
        open_record_file = open(path_to_file, 'w')

        with open(starting_align) as align_file_contents:
            for line in align_file_contents:
                if line.startswith('>'):
                    line = line.strip()
                    name = line.strip('>')
                    names_list.append(name)

        for starting_name in names_list:
            open_record_file.write(starting_name + '\n')

        open_record_file.close()


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



def handle_accession_options(accession_option, organism, folder_path, input_file, log_file):
    """Only for use when not automatically downloading SRA numbers"""
    print("Processing accession input option.")

    if accession_option == "AUTO_DL":
        if not organism == False:
            download_accessions(organism, folder_path)

        elif organism == False:
            print("A scientific name or NCBI taxon ID was not input.")
            print("Please use the --organism flag to specify what organism to search for.")

    elif accession_option == "AUTO_PATHDB":
        print("option not available yet")

    elif accession_option == "USER_INPUT":
        open_log_file = open(log_file, 'a')

        now = datetime.datetime.now()

        output_file = 'accessions_' + now.strftime('%Y-%m-%d-%H-%M-%S')

        cp_command = ['cp', input_file, folder_path + "/accession_files/" + output_file]

        cp_run = subprocess.Popen(cp_command, stdout=open_log_file, stderr=open_log_file)

        os.wait()

        # errs = cp_run.stderr.read().decode()
        # out = cp_run.stdout.read().decode()

        # write_to_log(out, errs, log_file)

        # subprocess.run(['cp', input_file, folder_path + "/accession_files/" + output_file])


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

    accession_file = open(output_file, 'w')

    esearch_command = ['esearch', '-db', 'sra', '-query', org_name]
    efetch_command = ['efetch', '-format', 'runinfo']

    esearch_accessions = subprocess.Popen(esearch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    efetch_accessions = subprocess.Popen(efetch_command, stdin=esearch_accessions.stdout, stdout=accession_file)

    # search_errs = esearch_accessions.stderr.read().decode()
    # search_out = esearch_accessions.stdout.read().decode()

    # write_to_log(search_out, search_errs, log_file)

    # fetch_errs = efetch_accessions.stderr.read().decode()
    # fetch_out = efetch_accessions.stdout.read().decode()

    # write_to_log(fetch_out, fetch_errs, log_file)
    # print(efetch_accessions.communicate())

    # dl_accessions = subprocess.Popen(accession_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # dl_accessions = subprocess.Popen(['esearch', '-db', 'sra', '-query', org_name, '|', 'efetch', '-format', 'runinfo', '>', output_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # dl_accessions = subprocess.Popen(['wget', '-O', output_file,  url], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print(dl_accessions.communicate())

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


def read_pathodb_csv_file(out_dir):
    """Reads the accession file provided by the users. This assume the user has done some filtering to paired end only along with some other information."""
    print("Reading CSV file off accession numbers.")
    os.chdir(out_dir + "/accession_files")


    current_accession = find_recent_date(out_dir + "/accession_files")

    csv = pd.read_csv(current_accession)

    output = []

#################################################################
    # This machinary is specific to working with pathogen db csv files when no external processing has been performed
    # Leaving it here but moving forward we're assuming the user has processed any hand fed inputs to include columns/datapoints they are interested in

    # csv = csv.rename(columns={'Library layout': 'LibraryLayout'})
    #
    # filtered_df = csv.query("LibraryLayout == 'PAIRED' and Platform == 'ILLUMINA'")
    #
    # filtered_df['Location'].replace('', np.nan, inplace=True)
    #
    # filtered_df.dropna(subset=['Location'], inplace=True)
    #
    # os.chdir(out_dir)
    #
    # return filtered_df
#######################################################################

    csv['Run'].replace('', np.nan, inplace=True)

    csv.dropna(subset=['Run'], inplace=True)

    os.chdir(out_dir)

    return csv


def find_recent_date(folder_of_files):
    print("Finding the file with the most recent data.")

    cwd = os.getcwd()
    os.chdir(folder_of_files)

    list_of_files = []

    for file in os.listdir(folder_of_files):
        # print(file)
        file_date = restructure_dates(file)
        # print(file_date)
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

    print("list of files", list_of_files)

    df = pd.DataFrame(list_of_files, columns=['file_name', 'creation_date'])
    print(df)

    # df['creation_date'] = pd.to_datetime(df['creation_date'])

    most_recent_accession_file = df['creation_date'].max()

    most_recent_row = df.loc[df['creation_date'] == most_recent_accession_file]

    current_accession_file = most_recent_row['file_name'].tolist()[0]

    os.chdir(cwd)

    return folder_of_files + '/' + current_accession_file


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
    # print(condensed_date)

    return condensed_date



def calulate_cores(set_cores):
    """Organizes and calulates the cores for use with Extensiphy."""
    print("Calculating threads for Extensiphy.")
    output = []
    total_cores = int(set_cores)
    runs = 0
    cores_per_run = 0

    if total_cores >= 2 and total_cores < 10:
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

    print(accession_db['Run'])
    sra_ids = accession_db['Run'].tolist()
    # print("Starting list of sras ", sra_ids)
    # print(fasta_names)
    if len(sra_ids) != 0:
        for num, name in enumerate(sra_ids):
            if name not in fasta_names:
                sras_to_keep.append(name)

    # print("Dupes removed list of sras ", sras_to_keep)
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
        print(name)

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
    """Take batched accessions, download for each batch and run EP. Then split up the alignment and remove the original folder"""
    ref = out_dir +'/intermediate_files/reference.fas'
    ep_output_align = out_dir + '/intermediate_files/ep_output/RESULTS/extended.aln'
    now = datetime.datetime.now()
    
    open_log_file = open(log_file, 'a')

    for accession_batch in accessions:
        # print(accession_batch)
        print("++++")
        print("Batch of current accessions is ", accession_batch)
        len_of_this_batch = len(accession_batch)
        not_downloaded_this_run = []

        for num, accession in enumerate(accession_batch):
            print("+++++++")
            # print(accession)
            download_status = fasterq_dump_reads(out_dir, accession)

            if download_status != "DOWNLOADED":

                not_downloaded_this_run.append(accession)

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
    Checks if the downloaded read files look appropriate for a Extensiphy run.
    Used to ensure an EP run can be performed, otherwise it must be skipped.
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





# def run_ep(base_dir_, align_, outdir_):
#     # print("current working directory: ", os.getcwd())
#     # print("current alignment ", current_alignment)
#     print("extensiphy.sh", "-a", align_, "-d", base_dir_ + "/read_files", "-1", "_1.fastq", "-2", "_2.fastq", "-o", outdir_)
#     # ep_process = subprocess.Popen(["/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-a", align_, "-d", base_dir_ + "/read_files", "-1", "_1.fastq", "-2", "_2.fastq", "-o", outdir_], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     # print(ep_process.communicate())
#     subprocess.run(["extensiphy.sh", "-a", align_, "-d", base_dir_ + "/read_files", "-1", "_1.fastq", "-2", "_2.fastq", "-o", outdir_])

#     # print("current working directory: ", os.getcwd())
#     # print("current alignment ", current_alignment)
#     # # subprocess.run(["/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-h"])
#     # print("Extensiphy run complete. BEGINNING PROCESSING")

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
