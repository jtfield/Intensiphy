#! /usr/bin/python3

# A test to make sure we get back the expected file structure from NCBI SRA
# wget -O shell_run.csv "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=Anser anser[Organism]"
import subprocess
import os


def download_accessions_function(outfile, taxon):
    url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=' + taxon + '[Organism]' #Now generate the url from the taxon
    subprocess.run(['wget', '-O', outfile, url])


def test_download_accessions(): #Must start with test for pytest to catch it
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]
    hand_checked_csv = absolute_path + '/shell_run.csv' #can use full path here
    outfile = absolute_path + '/ncbi_anser_anser.csv' #and here
    taxon = 'Anser anser'

    download_accessions_function(outfile, taxon)

    ex_acc = set()

    examplefi = open(hand_checked_csv).readlines()

    for lin in examplefi:
        split_line = lin.split(',')
        ex_acc.add(split_line[0])

    testfi = open(outfile).readlines()

    #TODO: Seems like NCBI hase removed a few anser anser SRA entries since
    # we downloaded this hand checked CSV. Adjust later?
    for num, lin in enumerate(testfi):
        if num < 6:
            split_line=lin.split(',')
            assert split_line[0] in ex_acc

    # we don't know that the csv we include with always match the size of the accessions on NCBI
    #assert len(examplefi) == len(testfi)
