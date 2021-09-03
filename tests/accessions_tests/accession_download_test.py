#! /usr/bin/python3

# run in command line to generate file to check against
# wget -O shell_run.csv "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=Anser anser[Organism]"
import subprocess
import argparse
import os
import inspect

# def parse_args():
#     parser = argparse.ArgumentParser(prog='Intensiphy', \
#         description='Automate downloading of high-throughput sequence data and updating of alignments using Extensiphy.')
#     # parser.add_argument('--accession_csv')
#     parser.add_argument('--ep_out_dir', help='Absolute path and folder name to create for outputs')
#     return parser.parse_args()

def download_accessions_test(outdir):
    # args = parse_args()

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]

    hand_checked_csv = absolute_path + '/Anser_anser_accessions.csv'
    outfile = 'ncbi_anser_anser.csv'
    url =  'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=Anser anser[Organism]'

    subprocess.run(['wget', '-O', outfile, url])

    new_file = open(outdir + '/' + outfile).readlines()

    ex_acc = set()

    for lin in new_file:
        split_line = lin.split(',')
        ex_acc.add(split_line[0])

    test_fi = open(hand_checked_csv).readlines()

    for lin in test_fi:
        split_line=lin.split(',')
        assert split_line[0] in ex_acc

# if __name__ == '__main__':
#     main()