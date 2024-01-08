#!/usr/bin/env python3

import os
import shutil

def main():

    print('Reseting test environment.')
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    test_output_dir = split_path_and_name[0] + '/test_ip_output_dir'

    os.chdir(test_output_dir)

    starting_dirs = ['accession_files', 'intermediate_files', 'read_files', 'run_log_files', 'sequence_storage']

    # remove files not in starting list of dirs
    for file_or_dir in os.listdir('.'):
        if file_or_dir not in starting_dirs:
            os.remove(file_or_dir)
            print('removed ', file_or_dir)

    # remove files from each directory
    for dir in os.listdir('.'):
        os.chdir(dir)
        list_of_files = os.listdir('.')
        if len(list_of_files) != 0:
            for file in list_of_files:
                try:

                    os.remove(file)
                    print('removed ', file)
                except IsADirectoryError:
                    shutil.rmtree(os.path.abspath(file))
                    print('removed dir', file)
        os.chdir('..')
    

    os.chdir(test_output_dir)
    print('Test environment reset.')






if __name__ == '__main__':
    main()