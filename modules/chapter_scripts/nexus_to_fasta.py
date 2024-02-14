#!/usr/bin/env python3

import os
import argparse
import pathlib
import shutil
from numpy.lib.shape_base import split
import pandas as pd
import subprocess
import datetime
import dateutil
import dendropy


def parse_args():
    parser = argparse.ArgumentParser(prog='nexus to newick', \
        description='convert a tree in nexus format to a tree in newick format')
    parser.add_argument('--tree', help='Newick phylogeny file.')
    parser.add_argument('--output_file', default='converted_tree.tre', help='newick phylogeny file.')
    return parser.parse_args()

def main():
    args = parse_args()

    # tree = dendropy.Tree.get(path=args.tree_file, schema="newick", preserve_underscores=True)

    output = []

    tree = dendropy.Tree.get(path=args.tree, schema="nexus", preserve_underscores=True)

    # print(tree)

    tree.write(path=args.output_file, schema="newick")







if __name__ == '__main__':
    main()
