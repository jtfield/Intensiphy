#!/usr/bin/env python3

import os
import argparse
import pathlib
import numpy as np
import pandas as pd
import subprocess
import datetime
# from multiprocessing import Pool, freeze_support
# import multiprocessing as mp

def handle_starting_tree(outdir, threads, tree_var):
    """Builds starting tree if a tree isnt added as input by user"""

    tree_storage = outdir + '/intermediate_files'
    align = tree_storage + '/alignment.fas'
    sequence_storage = outdir + '/sequence_storage'
    threads = str(int(threads[0]) * int(threads[1]))


    if tree_var == False:
        os.chdir(tree_storage)

        subprocess.run(['raxmlHPC-PTHREADS', '-m', 'GTRGAMMA', '-T', threads, '-s', align, '-p', '12345', '-n', 'starting_tree.tre'])

        os.chdir(outdir)

    elif tree_var != False:
        tree_file_exists = os.path.exists(tree_var)
        if tree_file_exists:

            abs_tree_path = os.path.realpath(tree_var)

            symlink_file = pathlib.Path(abs_tree_path)

            new_tree = output_dir_path + '/intermediate_files/RAxML_bestTree.starting_tree.tre'
            new_tree = pathlib.Path(new_tree)

            new_tree.symlink_to(symlink_file)



def construct_align_and_place(outdir):
    """Builds an alignment using all present sequences. \
    Uses starting tree to place new sequences in the tree (RAxML EPA)"""

    # Find current date and make folder with the date in the name
    # for record keeping purposes
    now = datetime.datetime.now()
    seq_storage_path = outdir + '/sequence_storage'
    starting_tree_path = outdir + '/intermediate_files/RAxML_bestTree.starting_tree.tre'
    taxa_list = os.listdir(seq_storage_path)
    suffix = '_.fas'

    phylo_est_dir = 'tree_inference_' + now.strftime('%Y-%m-%d-%H-%M-%S')
    phylo_dir_full_path = outdir + '/' + phylo_est_dir

    output_alignment_path = phylo_dir_full_path + '/extended.aln'

    os.mkdir(phylo_dir_full_path)

    # cat_command_start = ['cat']
    output = open(output_alignment_path, 'a')

    for dir in taxa_list:
        taxon_dir_path = seq_storage_path + '/' + dir
        taxon_seq_file_path = taxon_dir_path + '/' + dir + suffix
        # print(taxon_seq_file_path)

        # Add file path to cat command list
        # cat_command_start.append(taxon_seq_file_path)
        seq_file = open(taxon_seq_file_path, 'r').read()
        # print(type(seq_file))
        output.write(seq_file)
        output.write('\n')


    output.close()

    print("Updated alignment construction complete.")

    os.chdir(phylo_dir_full_path)
    # Run RAxML using the EPA (placement) algorithm

    # place_tree_build = ['raxmlHPC-PTHREADS', '-f', 'v', '­-n', 'placement.tre', '­-s', output_alignment_path, '-t', starting_tree_path, '-m' 'GTRGAMMA']

    place_tree_build = ['raxmlHPC-PTHREADS', '-f', 'v', '-n', 'placement.tre', '-s', output_alignment_path, '-t', starting_tree_path, '-m', 'GTRGAMMA']
    print("Running placement algorithm to add new sequences to the tree.")
    print(place_tree_build)
    # print('raxmlHPC-PTHREADS', '-s', output_alignment_path, '-t', starting_tree_path, '-m', 'GTRGAMMA', '-n', 'placement.tre', '­-f', 'v')
    # process = subprocess.Popen(['raxmlHPC-PTHREADS', '­-f', 'v', '­-s', output_alignment_path, '-t', starting_tree_path, '-m' 'GTRGAMMA', '­-n', 'placement.tre'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stdout, stderr = process.communicate()
    # print(stdout)
    # print(stderr)
    subprocess.run(place_tree_build)
    os.chdir(outdir)

    #     cat_command_start.append('\n')



    #
    # # print(cat_command_start)
    #
    # # Append the output carrot and combined alignment file name
    # cat_command_start.append('>')
    # cat_command_start.append(output_alignment_path)
    #
    # # print(cat_command_start)

    # Run cat command with subprocess
    # process = subprocess.Popen(cat_command_start, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # stdout, stderr = process.communicate()
    # # print(stdout)
    # output = open(output_alignment_path, 'w')
    # output.write(str(stdout, 'UTF-8'))
    # output.close()
    # print("Updated alignment construction complete.")
    #
    # os.chdir(phylo_dir_full_path)
    # # Run RAxML using the EPA (placement) algorithm
    # print("Running placement algorithm to add new sequences to the tree.")
    # print('raxmlHPC-PTHREADS', '­-f', 'v', '­-s', output_alignment_path, '-­t', starting_tree_path, '-­m', 'GTRCAT', '­-n', 'placement.tre')
    # # process = subprocess.Popen(['raxmlHPC-PTHREADS', '­-f', 'v', '­-s', output_alignment_path, '-­t', starting_tree_path, '-­m', 'GTRCAT', '­-n', 'placement.tre'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # # stdout, stderr = process.communicate()
    # # print(stdout)
    # # print(stderr)
    # os.chdir(outdir)
