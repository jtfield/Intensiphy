#! /usr/bin/python3

# A test to make sure we get back the expected file structure from NCBI SRA
# wget -O shell_run.csv "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=Anser anser[Organism]"
import subprocess
import os


def test_build_alignment():
    print("Running Assembly Test.")
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]
    subprocess.run(["gon_phyling.sh", "-d", absolute_path + "/testdata", "-1", "_R1.fq", "-2", "_R2.fq"])

    output_dir = absolute_path + "/testdata/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing"
    output_alignment = output_dir + "/combo.fas"
    print(absolute_path)
    print(output_alignment)
    print(os.path.exists(output_alignment))
    assert os.path.exists(output_alignment) == True

    with open(output_alignment) as align_file:
        for line in align_file:
            if ">" not in line and len(line) > 1:
                assert len(line) > 100

    print("Assembly Test Complete.")
