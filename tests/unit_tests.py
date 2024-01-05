#!/usr/bin/env python3

import os
import sys
import re
sys.path.insert(0, os.path.abspath(".."))
# from modules.fetch_and_align import (make_align, calculate_cores)
# from modules.fetch_and_align import *
import modules.fetch_and_align as faa
import modules.alignment_splitter as asp

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


def test_download_accessions():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'
    organism = 'txid413501'
    accessions_count = 0

    faa.download_accessions(organism, test_output_dir)

    for ind_file in os.listdir('.'):
        if re.match('accessions_', ind_file):
            accessions_count+=1

    assert accessions_count == 1

# TODO: make handle_accessions test

def test_make_align():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    test_align = split_path_and_name[0] + '/combo.fas'

    faa.make_align(False, test_output_dir, test_align)

    symlink_path = test_output_dir + '/intermediate_files/alignment.fas'
    
    assert os.path.islink(symlink_path)


# def alignment_splitter_test():

#     split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
#     test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

#     test_align = split_path_and_name[0] + '/combo.fas'

#     asp.split_alignment(test_align, test_output_dir)


# #TODO: Make selected_ref_test and unselected_ref_test
#     def unselected_ref_test():

#         split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
#         test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

#         test_align = split_path_and_name[0] + '/combo.fas'

#         faa.select_ref(False, test_output_dir)

#         assert os.path.isfile(test_output_dir + '/intermediate_files/reference.fas')

