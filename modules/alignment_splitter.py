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

            # Check if this sequence is already in the database
            dir_path = output_dir + '/' + name
            dir_exists = os.path.isdir(dir_path)

            # if no, add a folder for the associated data
            if dir_exists == False:
                print("Making new sequence repository for ", name)
                os.mkdir(dir_path)

                # separated sequences to individual files
                output = open(dir_path + '/' + name + '_.fas', 'w')
                output.write('>' + name)
                output.write('\n')
                output.write(seq)
                output.close()
