#! /usr/bin/python3

import os
import argparse
import pathlib
import shutil
from numpy.lib.shape_base import split
import pandas as pd
import subprocess
import datetime
import tests.accessions_tests.accession_download_test
import tests.assembly_tests.gon_phy_test

def parse_args():
    parser = argparse.ArgumentParser(prog='Intensiphy', \
        description='Automate downloading of high-throughput sequence data and updating of alignments using Extensiphy.')
    parser.add_argument('--align_file', default=False, help='input alignment option. \
        If an alignment is added using this option, samples will be added to this alignment.')
    parser.add_argument('--cores')
    parser.add_argument('--accs_file', default=False, help='Accession file if accession_method is set to USER_INPUT')
    parser.add_argument('--accession_method', default="AUTO_DL", help='Dictates how collecting and inputting accession numbers will be handled. \
        (OPTIONS: USER_INPUT, AUTO_DL, AUTO_PATHDB), (DEFAULT: AUTO_DL)')
    parser.add_argument('--ep_out_dir', help='Absolute path and folder name to create for outputs')
    parser.add_argument('--organism', type=str, nargs='+', help='scientific name of the organism or group of organisms you \
         would like to query SRA for and update your alignment with. Example: Neisseria gonorrhoeae[Organism] or txid482')
    parser.add_argument('-b', default=False, action='store_true', help='Toggles big bactch downloading \
        of fastq files instead of continuous downloading.')
    return parser.parse_args()

def main():
    args = parse_args()
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    ref_taxon = ""

    reset_tests()

    # Check if output dir has been made already from a previous run
    # If output dir exists, the fundamentals of the program change to suit
    # a repeating run.
    # Returns True if output dir exists already and False if output dir doesn't exist prior to execution.
    dir_existence = check_dir_exists(args.ep_out_dir)

    os.chdir(args.ep_out_dir)
    absolute_output_dir_path = os.path.abspath(os.getcwd())

    # download_accessions(args.organism, args.ep_out_dir)
    accessions = handle_accession_options(args.accession_method, args.organism, absolute_output_dir_path, args.accs_file)

    make_align(dir_existence, absolute_output_dir_path, args.accs_file)

    current_align_file = get_most_recent_align(dir_existence, args.align_file, absolute_output_dir_path)

    # calculate the core organization to pass to Extensiphy
    get_cores = calulate_cores(args.cores)

    # download_accessions(args.organism, args.ep_out_dir)

    # read_file = read_csv_file(args.accession_csv)
    read_fasta = read_fasta_names(current_align_file, absolute_output_dir_path)

    #read list of accessions
    read_accessions = read_csv_file(args.ep_out_dir)

    #Check list of run SRA numbers vs the sequences already in the alignment to prevent duplicates.
    remove_paired_dupes = check_duplicate_accesions(read_accessions[0], read_fasta)
    remove_single_dupes = check_duplicate_accesions(read_accessions[1], read_fasta)

    # Handle how we'll download SRA files: big batch or continuously while running Extensiphy
    if not remove_paired_dupes.empty:
        process_data = downloading_and_running(args.b, remove_paired_dupes, get_cores, args.ep_out_dir, args.align_file)

def reset_tests():
    """Function to reset test outputs so this isn't done by hand each time"""
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]

    #delete files created during test gon_phyling run
    gon_phy_test_initial_dir = absolute_path + "/tests/assembly_tests/testdata/trimmed_reads"
    if os.path.exists(gon_phy_test_initial_dir):
        print("Removing previous test outputs")
        subprocess.run(["rm", "-r", absolute_path + "/tests/assembly_tests/testdata/trimmed_reads"])
        # subprocess.run(["rm", absolute_path + "/tests/assembly_tests/testdata/x"])
        subprocess.run(["rm", "-r", absolute_path + "/tests/assembly_tests/testdata/repaired_reads"])

    #delete csv downloaded during accession download test
    if os.path.exists(absolute_path + "/tests/accessions_tests/ncbi_anser_anser.csv"):
        os.remove(absolute_path + "/tests/accessions_tests/ncbi_anser_anser.csv")


def check_dir_exists(dir_name):
    """Check if output directory exists. If the directory exists, a previous run of Intensiphy was probably performed. \
        If the directory doesn't exist, no previous run was performed, build all necessary sub directories."""

    does_dir_exist = os.path.isdir(dir_name)
    print("Does the output specified output directory exist?: ", does_dir_exist)

    if does_dir_exist:

        assert os.path.isdir(dir_name + "/read_files")
        assert os.path.isdir(dir_name + "/run_log_files")
        assert os.path.isdir(dir_name + "/alignments")
        assert os.path.isdir(dir_name + "/accession_files")

        return does_dir_exist

    else:
        subprocess.run(["mkdir", dir_name])
        subprocess.run(["mkdir", dir_name + "/read_files"])
        subprocess.run(["mkdir", dir_name + "/run_log_files"])
        subprocess.run(["mkdir", dir_name + "/alignments"])
        subprocess.run(["mkdir", dir_name + "/accession_files"])

        assert os.path.isdir(dir_name + "/read_files")
        assert os.path.isdir(dir_name + "/run_log_files")
        assert os.path.isdir(dir_name + "/alignments")
        assert os.path.isdir(dir_name + "/accession_files")

        return does_dir_exist

def make_align(run_bool, output_dir_path, input_accessions):
    """If this is the start of a run (no output folder existed prior to starting this run) \
    used gon_phyling to build an alignment based on the first 6 available samples"""

    #Test!
    tests.assembly_tests.gon_phy_test.test_build_alignment()

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

            symlink_align(output_dir_path + "/starting_align_files/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/combo.fas", output_dir_path + "/alignments)


def get_most_recent_align(run_bool, starting_align, output_dir_path):
    """Determines how to handle the run if the output dir already exists (indicating a continuing run) \
        or if the output dir doesn't exist yet (indicating the starting phase of a run)."""
    if run_bool:

        current_alignment = find_recent_date(output_dir_path + "/alignments")

        return current_alignment

    else:
        if starting_align != False:

        current_alignment = symlink_align(starting_align, output_dir_path)

        return current_alignment


def handle_accession_options(accession_option, organism, folder_path, input_file):
    """Only for use when not automatically downloading SRA numbers"""
    print("tmp input")

    if accession_option == "AUTO_DL":
        download_accessions(organism, folder_path)

    elif accession_option == "AUTO_PATHDB":
        print("option not available yet")

    elif accession_option == "USER_INPUT":
        subprocess.run(['cp', input_file, folder_path + "/accession_files/"])



def download_accessions(org_name, out_dir):
    """Download the run info file that includes run ID accession numbers and info on how the sequences were produced."""

    # Test!
    tests.accessions_tests.accession_download_test.test_download_accessions()

    now = datetime.datetime.now()

    os.chdir(out_dir + "/accession_files")

    # ' '.join(org_name)
    org_name = org_name[0].strip("'")

    output_file = 'accessions_' + now.strftime('%Y-%m-%d')

    url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=' + org_name

    subprocess.run(['wget', '-O', output_file,  url])

    os.chdir(out_dir)


def read_csv_file(out_dir):
    """Reads the accession file provided by the users and parses compatible sequences."""

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

    os.chdir(folder_of_files)

    list_of_files = []

    for file in os.listdir(folder_of_files):
        data = []
        data.append(file)
        split_file = file.split('_')
        date = split_file[1]
        data.append(date)
        # created_time = os.stat(file).st_ctime
        # data.append(created_time)
        list_of_files.append(data)

    df = pd.DataFrame(list_of_files, columns=['file_name', 'creation_date'])

    df['creation_date'] = pd.to_datetime(df['creation_date'])

    most_recent_accession_file = df['creation_date'].max()

    most_recent_row = df.loc[df['creation_date'] == most_recent_accession_file]

    current_accession_file = most_recent_row['file_name'].tolist()[0]

    return current_accession_file


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

    copy_df = accession_db.copy()

    for name, value in accession_db['Run'].iteritems():
        if value in fasta_names:
            assert type(value) == str
            assert len(value) > 1
            copy_df.drop(name)
            print(name)
            print(value)

    return copy_df

def read_fasta_names(_align, _outdir):
    """Reads the names of the sequences from the fasta file."""
    names = []

    with open(_align, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                seq_name = line.strip(">").strip('\n')
                assert len(seq_name) > 1
                assert type(seq_name) == str
                names.append(seq_name)

    return names

def symlink_align(align, outdir):

    now = datetime.datetime.now()

    # Symlink the alignment file into the alignments folder
    abs_align_path = os.path.realpath(align)

    new_align_name = outdir + '/alignments/msa_' + now.strftime('%Y-%m-%d')

    symlink_file = pathlib.Path(new_align_name)

    symlink_file.symlink_to(abs_align_path)

    return new_align_name

def align_rename_and_move(align, outdir):

    now = datetime.datetime.now()

    # Symlink the alignment file into the alignments folder
    abs_align_path = os.path.realpath(align)

    new_align_name = outdir + '/alignments/msa_' + now.strftime('%Y-%m-%d-%H-%M-%S')

    new_align_path = pathlib.Path(new_align_name)

    shutil.copy(abs_align_path, new_align_name)

    return new_align_name

def prepare_batch_accessions(accessions, runs):
    """Return a list of lists containing the number of files to download before each Extensiphy run"""
    chunks = [accessions[x:x+runs] for x in range(0, len(accessions), runs)]

    return chunks


def downloading_and_running(method, accessions, run_num, out_dir, align):

    # Identify if we're processing data by downloading in batches between Extesniphy runs
    # or by downloading everything all at once and running extensiphy after

    # pull accessions from df
    run_ids = []
    for name, value in accessions['Run'].iteritems():
        run_ids.append(value)

    if method == False:
        # continuous gradual downloading of data has been selected/left as default.
        batches_run = batch_download(run_ids, run_num[0], out_dir, align)

    elif method:
        # Bulk download all fastq files before running Extensiphy
        bulk_run = bulk_download(run_ids, out_dir)

        print("Bulk data download complete.")
        print("Beginning Extensiphy run to update your alignment.")

        os.chdir(out_dir)

        # TODO: specify cleanup of intermediate files in all runs
        subprocess.run(["multi_map.sh", "-a", align, "-d", out_dir + "/read_files", "-i", "CLEAN", "-p", str(run_num[0]) ,"-c", str(run_num[1]), "-1", "_1.fastq", "-2", "_2.fastq" ])

def batch_download(accession_list, runs_number, out_dir, alignment):
    """Prepares accessions to be downloaded in batches between runs of Extensiphy."""
    batches_of_accessions = prepare_batch_accessions(accession_list, runs_number)

    read_dir_path = out_dir + "/read_files"
    ep_output = out_dir + "/ep_tmp_outputs"
    standard_align_path = ep_output + "/RESULTS/extended.aln"
    align_outdir = out_dir + "/alignments"

    os.chdir(read_dir_path)

    for accessions in batches_of_accessions:

        # print(accessions)
        for single_accession in accessions:
            # print(single_accession)
            subprocess.run(["fasterq-dump", "--split-files", single_accession])

        subprocess.run(["multi_map.sh", "-a", alignment, "-d", out_dir + "/read_files", "-i", "CLEAN", "-p", str(runs_number[0]) ,"-c", str(runs_number[1]), "-1", "_1.fastq", "-2", "_2.fastq", -o, ep_output])

        # TODO: move the output alignment and the run log file to the appropriate places
        rm_fastq_files(read_dir_path)

        # Copy and rename output alignment and put it the alignments folder for the next round
        align_rename_and_move(standard_align_path, align_outdir)

        # Remove Extensiphy output directory and contents and prepare for next phase of loop
        try:
            shutil.rmtree(ep_output)
        except OSError as e:
            print("Error: %s : %s" % (ep_output, e.strerror))


def rm_read_files(out_dir):
    """Quick function to remove reads that have been used"""
    fastq_file_list = os.listdir(read_dir_path)
    for fq_file in fastq_file_list:
        path_to_file = os.path.join(out_dir, fq_file)
        os.remove(path_to_file)

def bulk_download(accession_list, out_dir):
    """Bulk download every sequence (unlike the batch method)"""

    os.chdir(out_dir + "/read_files")

    for single_accession in accession_list:

        print(single_accession)
        subprocess.run(["fasterq-dump", "--split-files", single_accession])



if __name__ == '__main__':
    main()
