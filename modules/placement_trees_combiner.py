#!/usr/bin/env python3

import os
import argparse
import shutil
import re


def parse_args():
    parser = argparse.ArgumentParser(prog='Placement tree combiner', \
        description='Run this program to pull the all the placement trees from their respective folders into a single folder.')
    parser.add_argument('--folder_of_folders', default=False, help='Path to placement folder of folders.')
    parser.add_argument('--output_dir', default='placement_trees', help='output location to put all the placement trees.')
    return parser.parse_args()

def main():
    args = parse_args()

    absolute_dir_path = os.path.abspath(args.folder_of_folders)

    output_folder_exists = os.path.isdir(args.output_dir)

    if output_folder_exists == False:
        os.mkdir(args.output_dir)

    contents = os.listdir(absolute_dir_path)

    placement_tree_regex = 'RAxML_labelledTree\.placement_.+\.tre'

    compile_regex = re.compile(placement_tree_regex)

    for inner_folder in contents:
        print(inner_folder)
        # abs_in_folder_path = os.path.abspath(inner_folder)

        abs_in_folder_path = absolute_dir_path + '/' + inner_folder

        # print(abs_in_folder_path)
        # print(os.listdir(abs_in_folder_path))

        inner_dir_contents = os.listdir(abs_in_folder_path)

        if len(inner_dir_contents) > 0:

            for file_name in inner_dir_contents:
                # print(type(file_name))

                find_placement_tree = re.match(compile_regex, file_name)

                if find_placement_tree:
                    # print(find_placement_tree[0])
                    placement_tree = find_placement_tree[0]

                    abs_path_to_place_tree = abs_in_folder_path + '/' + find_placement_tree[0]

                    print(abs_path_to_place_tree)

                    shutil.copyfile(abs_path_to_place_tree, args.output_dir + '/' + placement_tree)









if __name__ == '__main__':
    main()
