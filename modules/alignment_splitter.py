#! /usr/bin/python3

import os
import argparse
import pathlib

def split_alignment(align_file, output_dir):
    """Splits an alignment file into individual sequence files."""
    print("Align splitter.")

    # Open and read alignment file
    open_align = open(align_file, 'r')
    read_align = open_align.read()

    # Split alignment file
    split_file = read_align.split('>')

    # Iterate over split file
    for chunk in split_file:

        # Split the chunk on the first newline character,
        # separating sequence name and sequence
        split_chunk = chunk.split('\n', 1)

        # If statement to make sure the chunk isnt empty
        if len(split_chunk) == 2:

            # Remove newline characters from the sequence
            name = split_chunk[0]
            seq = split_chunk[1]
            seq = seq.replace('\n','')

            # separated sequences to individual files
            output = open(output_dir + '/' + name + '_.fas', 'w')
            output.write('>' + name)
            output.write('\n')
            output.write(seq)
            output.close()
