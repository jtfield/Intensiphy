#!/usr/bin/env python3

# import dendropy
import argparse
# from dendropy.calculate import treecompare
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--orig_aln')
    parser.add_argument('--updated_aln')
    parser.add_argument('--combined_aln', default='asserted_combined_align.fas')
    return parser.parse_args()

def main():
    args = parse_args()

    updated_aln_taxa = []

    combined_output = []

    read_updated = open(args.updated_aln, 'r').read()
    split_updated = read_updated.split('>')


    read_orig = open(args.orig_aln, 'r').read()
    split_orig = read_orig.split('>')

    for up_name_and_seq in split_updated:
        if len(up_name_and_seq) > 1:
            split_up_name_and_seq = up_name_and_seq.split('\n', 1)
            upname = split_up_name_and_seq[0]
            upseq = split_up_name_and_seq[1]
            # print(upname)
            updated_aln_taxa.append(upname)

    # print(updated_aln_taxa)

    for orig_name_and_seq in split_orig:
        if len(orig_name_and_seq) > 1:
            split_orig_name_and_seq = orig_name_and_seq.split('\n', 1)
            origname = split_orig_name_and_seq[0]
            origseq = split_orig_name_and_seq[1]


            if origname not in updated_aln_taxa:
                combined_output.append(orig_name_and_seq)

    for up_name_and_seq_2 in split_updated:
        if len(up_name_and_seq_2) > 1:
            combined_output.append(up_name_and_seq_2)

    # print(combined_output[0])
    # print(len(combined_output))

    output_file = open(args.combined_aln, 'w')
    for chunk in combined_output:
        output_file.write(">")
        output_file.write(chunk)
        output_file.write("\n")

    output_file.close()
    



if __name__ == '__main__':
    main()
