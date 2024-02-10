#!/usr/bin/env python3

import os
import sys
import re
sys.path.insert(0, os.path.abspath(".."))
# from modules.fetch_and_align import (make_align, calculate_cores)
# from modules.fetch_and_align import *
import modules.fetch_and_align as faa
import modules.alignment_splitter as asp
import modules.tree_assess as ta

def test_calculate_cores():

    cores = 6

    core_values = faa.calculate_cores(cores)

    assert type(core_values) == list
    assert len(core_values) == 2
    for num in core_values:
        assert type(num) == int

    assert core_values[0] == 2
    assert core_values[1] == 3

    new_cores = 20

    core_values_two = faa.calculate_cores(new_cores)

    assert type(core_values_two) == list
    assert len(core_values_two) == 2
    for num in core_values_two:
        assert type(num) == int

    assert core_values_two[0] == 10
    assert core_values_two[1] == 2


def test_handle_accessions():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'
    test_accessions_file = split_path_and_name[0] + '/test_accessions'
    log_file_name = test_output_dir + '/run_log_files/subprocess_command_logs.txt'
    accessions_count = 0

    faa.handle_accession_options('USER_INPUT', False, test_output_dir, test_accessions_file, log_file_name)

    os.chdir(test_output_dir + '/accession_files')

    for ind_file in os.listdir('.'):
        if re.match('accessions_.+', ind_file):
            accessions_count+=1

    assert accessions_count >=1

    os.chdir(test_output_dir)


def test_download_accessions():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'
    organism = 'txid413501'
    accessions_count = 0

    faa.download_accessions(organism, test_output_dir)

    os.chdir(test_output_dir + '/accession_files')

    for ind_file in os.listdir('.'):
        if re.match('accessions_.+', ind_file):
            accessions_count+=1

    assert accessions_count >=1

    os.chdir(test_output_dir)


def test_make_align():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    test_align = split_path_and_name[0] + '/combo.fas'

    faa.make_align(False, test_output_dir, test_align)

    symlink_path = test_output_dir + '/intermediate_files/alignment.fas'
    
    assert os.path.islink(symlink_path)


def test_split_alignment():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'
    seq_storage = test_output_dir + '/sequence_storage'

    test_align = split_path_and_name[0] + '/combo.fas'

    asp.split_alignment(test_align, seq_storage)

    sep_seq_dirs = os.listdir(seq_storage)

    for seq_dir in sep_seq_dirs:
        # print(seq_dir)
        detected_fasta = 0
        os.chdir(seq_storage + '/' + seq_dir)
        seq_dir_files = os.listdir('.')
        for seq_file in seq_dir_files:
            if seq_file.endswith('_.fas'):
                detected_fasta+=1
        assert detected_fasta == 1


def test_unselected_ref():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    test_align = split_path_and_name[0] + '/combo.fas'

    faa.select_ref(False, test_output_dir)

    assert os.path.isfile(test_output_dir + '/intermediate_files/reference.fas')

def test_selected_ref():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    test_align = split_path_and_name[0] + '/combo.fas'

    selected_ref_align = test_output_dir + '/intermediate_files/reference.fas'

    ref_name = 'SRR7367521.ref'

    faa.select_ref(ref_name, test_output_dir)

    open_ref_file = open(selected_ref_align, 'r').readline()

    assert re.match('>' + ref_name, open_ref_file)

def test_handle_no_starting_tree():

    threads_ = [2,2]

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    ta.handle_starting_tree(test_output_dir, threads_, False, 'ON')

    assert os.path.isfile(test_output_dir + '/intermediate_files/RAxML_bestTree.starting_tree.tre')

def test_handle_present_starting_tree():

    threads_ = [2,2]

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    test_tree = split_path_and_name[0] + '/best_test_tree.tre'

    tree_symlink_path = test_output_dir + '/intermediate_files/RAxML_bestTree.starting_tree.tre'

    if os.path.isfile(tree_symlink_path):
        os.remove(tree_symlink_path)

    ta.handle_starting_tree(test_output_dir, threads_, test_tree, 'ON')

    assert os.path.islink(tree_symlink_path)

def test_clean_incomplete_downloads():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    read_dir = test_output_dir + '/read_files'

    test_dir_names = ['SRRFAKE1', 'SRRFAKE2']
    for dir_name in test_dir_names:
        os.mkdir(read_dir + '/' + dir_name)
    
    files_list = os.listdir(read_dir)

    assert len(files_list) != 0

    faa.clean_incomplete_downloads(test_output_dir)

    second_files_list = os.listdir(read_dir)

    assert len(second_files_list) == 0

def test_read_fasta_names():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    seq_dir_names = ['SRR15276367', 'SRR18324696', 'SRR26630346', 'SRR7367521.ref', 'SRR7439215']

    func_output = faa.read_fasta_names(test_output_dir)

    for seq_name in func_output:
        assert seq_name in seq_dir_names

def test_check_duplicate_accessions():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    test_accessions_file = split_path_and_name[0] + '/test_accessions'

    test_current_seqs = ['SRR15276367', 'SRR18324696', 'SRR26630346', 'SRR7367521.ref', 'SRR7439215']

    faa.handle_accession_options("USER_INPUT", False, test_output_dir, test_accessions_file, "test_log.txt")

    test_accessions_db = faa.read_pathodb_csv_file(test_output_dir)

    test_filtered_accessions = faa.check_duplicate_accesions(test_accessions_db, test_current_seqs)

    for seq_name in test_current_seqs:
        assert seq_name not in test_filtered_accessions

def test_prepare_batch_accessions():

    threads = 2

    test_list_of_ids = ['SRR01', 'SRR02', 'SRR03', 'SRR04', 'SRR05', 'SRR06', 'SRR07', 'SRR08', 'SRR09', 'SRR10']

    processed_id_list = faa.prepare_batch_accessions(test_list_of_ids, threads)

    assert len(processed_id_list) == 5
    
    for list_chunk in processed_id_list:
        assert len(list_chunk) == 2

def test_write_current_runs_names():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    test_list_of_ids = ['SRR01', 'SRR02', 'SRR03', 'SRR04', 'SRR05', 'SRR06', 'SRR07', 'SRR08', 'SRR09', 'SRR10']

    faa.write_current_run_names(test_output_dir, test_list_of_ids)

    assert os.path.isfile(test_output_dir + '/current_run_taxa_added.txt')