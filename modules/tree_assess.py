#!/usr/bin/env python3

import os
import pathlib
import numpy as np
import pandas as pd
import subprocess
import datetime
from typing import List, Union


def handle_starting_tree(outdir: str, threads: List[int], tree_var: Union[bool, str], placement_var: str) -> None:
    """
    Builds a starting tree if a tree isn't added as input by the user.

    Args:
        outdir (str): The path to the output directory.
        threads (List[int]): A list of two integers representing the number of threads to use.
        tree_var (Union[bool, str]): A boolean indicating whether a tree file has been input, or the path to the tree file.
        placement_var (str): A string indicating whether placement mode is on or off.

    Returns:
        None
    """

    # Define the path to the tree storage directory
    tree_storage: str = outdir + '/intermediate_files'
    # Define the path to the alignment file
    align: str = tree_storage + '/alignment.fas'
    # Define the path to the sequence storage directory
    sequence_storage: str = outdir + '/sequence_storage'
    # Convert the number of threads to a string
    threads: str = str(int(threads[0]) * int(threads[1]))

    # If placement mode is on
    if placement_var == 'ON':
        print('Placement mode has been selected.')
        print('Determining whether a tree file has been input.')
        # If no tree file has been input
        if tree_var == False:
            print('No starting tree has been provided.')
            print('Constructing a starting tree using the input starting alignment.')
            # Change the current working directory to the tree storage directory
            os.chdir(tree_storage)

            # Run RAxML to construct a starting tree
            subprocess.run(['raxmlHPC-PTHREADS', '-m', 'GTRGAMMA', '-T', threads, '-s', align, '-p', '12345', '-n', 'starting_tree.tre'])

            # Change the current working directory back to the output directory
            os.chdir(outdir)

        # If a tree file has been input
        elif tree_var != False:
            print('A starting tree has been provided.')
            print('Using the input starting tree for phylogenetic placement.')
            # Check if the tree file exists
            tree_file_exists: bool = os.path.exists(tree_var)
            if tree_file_exists:

                # Get the absolute path to the tree file
                abs_tree_path: str = os.path.realpath(tree_var)

                # Create a Path object for the tree file
                symlink_file: pathlib.Path = pathlib.Path(abs_tree_path)

                # Define the path to the new tree file
                new_tree: str = tree_storage + '/RAxML_bestTree.starting_tree.tre'
                new_tree: pathlib.Path = pathlib.Path(new_tree)

                # Create a symbolic link from the new tree file to the original tree file
                new_tree.symlink_to(symlink_file)

    # If placement mode is off
    elif placement_var == 'OFF':
        print('Placement mode has not been selected.')


def construct_align_and_place(outdir: str) -> None:
    """
    Builds an alignment using all present sequences and uses a starting tree to place new sequences in the tree (RAxML EPA).

    Args:
        outdir (str): The path to the output directory.

    Returns:
        None
    """

    # Get the current date and time
    now: datetime.datetime = datetime.datetime.now()

    # Define the path to the sequence storage directory
    seq_storage_path: str = outdir + '/sequence_storage'

    # Define the path to the starting tree file
    starting_tree_path: str = outdir + '/intermediate_files/RAxML_bestTree.starting_tree.tre'

    # Get a list of all directories in the sequence storage directory
    taxa_list: List[str] = os.listdir(seq_storage_path)

    # Define the suffix for the sequence file names
    suffix: str = '_.fas'

    # Define the name of the directory for the phylogenetic estimation
    phylo_est_dir: str = 'tree_inference_' + now.strftime('%Y-%m-%d-%H-%M-%S')

    # Define the full path to the phylogenetic estimation directory
    phylo_dir_full_path: str = outdir + '/' + phylo_est_dir

    # Define the path to the output alignment file
    output_alignment_path: str = phylo_dir_full_path + '/extended.aln'

    # Create the phylogenetic estimation directory
    os.mkdir(phylo_dir_full_path)

    # Open the output alignment file for appending
    output = open(output_alignment_path, 'a')

    # For each directory in the list of directories
    for dir in taxa_list:
        # Define the path to the directory for the taxon
        taxon_dir_path: str = seq_storage_path + '/' + dir

        # Define the path to the sequence file for the taxon
        taxon_seq_file_path: str = taxon_dir_path + '/' + dir + suffix

        # Open the sequence file for reading and read its contents
        seq_file: str = open(taxon_seq_file_path, 'r').read()

        # Write the contents of the sequence file to the output alignment file
        output.write(seq_file)

        # Write a newline character to the output alignment file
        output.write('\n')

    # Close the output alignment file
    output.close()

    # Print a message indicating that the construction of the updated alignment is complete
    print("Updated alignment construction complete.")

    # Change the current working directory to the phylogenetic estimation directory
    os.chdir(phylo_dir_full_path)

    # Define the command to run RAxML using the EPA (placement) algorithm
    place_tree_build: List[str] = ['raxmlHPC-PTHREADS', '-f', 'v', '-n', 'placement.tre', '-s', output_alignment_path, '-t', starting_tree_path, '-m', 'GTRGAMMA']

    # Print a message indicating that the placement algorithm is running to add new sequences to the tree
    print("Running placement algorithm to add new sequences to the tree.")

    # Print the command to run RAxML
    print(place_tree_build)

    # Run RAxML
    subprocess.run(place_tree_build)

    # Change the current working directory back to the output directory
    os.chdir(outdir)