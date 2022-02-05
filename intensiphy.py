#! /usr/bin/python3

import os
import argparse
import pathlib
import shutil
from numpy.lib.shape_base import split
import pandas as pd
import subprocess
import datetime
import dateutil
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

    # reset_tests()

    # Check if output dir has been made already from a previous run
    # If output dir exists, the fundamentals of the program change to suit
    # a repeating run.
    # Returns True if output dir exists already and False if output dir doesn't exist prior to execution.
    dir_existence = check_dir_exists(args.ep_out_dir)

    os.chdir(args.ep_out_dir)
    absolute_output_dir_path = os.path.abspath(os.getcwd())

    # download_accessions(args.organism, args.ep_out_dir)
    accessions = handle_accession_options(args.accession_method, args.organism, absolute_output_dir_path, args.accs_file)

    make_align(dir_existence, absolute_output_dir_path, args.accs_file, args.align_file)

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
    print("Reseting tests!")
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
    print("Checking if output directory already exists.")
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

def make_align(run_bool, output_dir_path, input_accessions, alignment_var):
    """If this is the start of a run (no output folder existed prior to starting this run) \
    used gon_phyling to build an alignment based on the first 6 available samples"""
    print("Building a starting alignment.")
    if alignment_var != False:
        print("Input alignment detected.")
        print("Skipping de novo alignment production.")

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

def get_most_recent_align(run_bool, starting_align, output_dir_path):
    """Determines how to handle the run if the output dir already exists (indicating a continuing run) \
        or if the output dir doesn't exist yet (indicating the starting phase of a run)."""
    print("Processing alignment input options.")
    if starting_align == False:

        current_alignment = find_recent_date(output_dir_path + "/alignments")

        return current_alignment

    else:
        if starting_align != False:

            # current_alignment = symlink_align(starting_align, output_dir_path)
            current_alignment = align_rename_and_move(starting_align, output_dir_path + "/alignments")

            return current_alignment

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
    print("Getting taxa labels from the current alignment file.")
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
    print("Linking provided alignment to Intensiphy output folder.")

    now = datetime.datetime.now()

    # Symlink the alignment file into the alignments folder
    abs_align_path = os.path.realpath(align)

    new_align_name = outdir + '/alignments/msa_' + now.strftime('%Y-%m-%d')

    symlink_file = pathlib.Path(new_align_name)

    symlink_file.symlink_to(abs_align_path)

    return new_align_name

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

def prepare_batch_accessions(accessions, runs):
    """Return a list of lists containing the number of files to download before each Extensiphy run"""
    print("building batch accession numbers download plan.")
    chunks = [accessions[x:x+runs] for x in range(0, len(accessions), runs)]

    return chunks

def downloading_and_running(method, accessions, run_num, out_dir, align):

    # Identify if we're processing data by downloading in batches between Extesniphy runs
    # or by downloading everything all at once and running extensiphy after
    print("Beginning bulk downloading data and updating alignment using Extensiphy.")
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
    print("Beginning batch data download and run scheme.")
    print("^^^^^^^^^^^^^^^^^^^")

    read_dir_path = out_dir + "/read_files"
    ep_output = out_dir + "/ep_tmp_outputs"
    standard_align_path = ep_output + "/RESULTS/extended.aln"
    align_outdir = out_dir + "/alignments"
    current_alignment = alignment
    ref_taxon = ""

    os.chdir(read_dir_path)

########################################################################
# testing EP running on a single accession instead of in a loop
    # grouped_accession = batches_of_accessions[0]
    # print(grouped_accession)
    # single_accession = grouped_accession[0]
    #
    # fasterq_dump_reads(out_dir, single_accession)
    #
    # run_ep(out_dir, current_alignment, ep_output)
    #
    # print("DONE")

########################################################################



    for num, accessions in enumerate(batches_of_accessions):
        print("BEGINNING NEW BATCH SECTION.")
        print("+++++++++++++++++++++++++++++++++++++++++++++")
        print("+++++++++++++++++++++++++++++++++++++++++++++")
        print("+++++++++++++++++++++++++++++++++++++++++++++")
        print("Currend Accession Batch ", accessions)

        if num <= 2:

            # os.chdir(out_dir + "/read_files")
            print("current working directory: ", os.getcwd())
            # print(accessions)
            for single_accession in accessions:
                print("BEGINNING EP LOOP.")
                print("<>")
                print("<>")
                print("<>")
                print("<>")
                print("<>")
                print("<>")
                # print(single_accession)

                # print("Current accession", single_accession)
                # os.chdir(out_dir + "/read_files")
                # print("current working directory: ", os.getcwd())
                # print("current alignment ", current_alignment)
                # reads_dl = subprocess.Popen(["fasterq-dump", "--split-files", single_accession], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                # print(reads_dl.communicate())
                # os.chdir(out_dir)
                fasterq_dump_reads(out_dir, single_accession)


                # print("Extensiphy run options.")
                # print(current_alignment)
                # print(out_dir + "/read_files")
                # print(ep_output)
                # print("current working directory: ", os.getcwd())


                # os.chdir(out_dir + "/read_files")
                # print("current working directory: ", os.getcwd())
                # subprocess.run(["/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-a", current_alignment, "-d", out_dir + "/read_files", "-i", "CLEAN", "-1", "_1.fastq", "-2", "_2.fastq", "-o", ep_output])
                # print("current working directory: ", os.getcwd())
                # print("current alignment ", current_alignment)
                # ep_process = subprocess.Popen(["/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-a", current_alignment, "-d", out_dir + "/read_files", "-1", "_1.fastq", "-2", "_2.fastq", "-o", ep_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                # print(ep_process.communicate())
                # print("current working directory: ", os.getcwd())
                # print("current alignment ", current_alignment)
                # # subprocess.run(["/home/vortacs/tmp_git_repos/extensiphy/extensiphy.sh", "-h"])
                # print("Extensiphy run complete. BEGINNING PROCESSING")
                run_ep(out_dir, current_alignment, ep_output)

                print(">>>>")
                print(">>>>")
                print(">>>>")
                print(">>>>")
                print(">>>>")
                print(">>>>")
                exit()

                check_ep_output = os.path.isdir(ep_output)
                if check_ep_output:
                    print("000")
                    print("000")
                    print("000")
                    print("000")
                    print("000")
                    print("Found EP output directory")
                    print("current alignment ", current_alignment)
                    print("current working directory: ", os.getcwd())
                    # Copy and rename output alignment and put it the alignments folder for the next round
                    align_rename_and_move(standard_align_path, align_outdir)

                    print("current working directory: ", os.getcwd())
                    current_alignment = find_recent_date(align_outdir)
                    print("curent alignment", current_alignment)


                    print("current working directory: ", os.getcwd())
                    os.rename(ep_output, ep_output + "_" + str(num))

                    print("current working directory: ", os.getcwd())
                    rm_read_files(read_dir_path)

                    # Remove Extensiphy output directory and contents and prepare for next phase of loop
                    # try:
                    #     shutil.rmtree(ep_output)
                    # except OSError as e:
                    #     print("Error: %s : %s" % (ep_output, e.strerror))

                elif check_ep_output == False:
                    print("ERROR: no EP output directory was found.")
                    exit()

                print("current working directory: ", os.getcwd())
                print("STARTING NEXT EP LOOP")

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

def rm_read_files(out_dir):
    """Quick function to remove reads that have been used"""
    print("removing read files in preparation for next run.")
    fastq_file_list = os.listdir(out_dir)
    print(fastq_file_list)
    for fq_file in fastq_file_list:
        path_to_file = os.path.join(out_dir, fq_file)
        os.remove(path_to_file)

def bulk_download(accession_list, out_dir):
    """Bulk download every sequence (unlike the batch method)"""
    print("Bulk data download scheme begun.")
    os.chdir(out_dir + "/read_files")

    for single_accession in accession_list:

        print(single_accession)
        subprocess.run(["fasterq-dump", "--split-files", single_accession])



if __name__ == '__main__':
    main()
