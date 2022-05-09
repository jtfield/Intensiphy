#!/usr/bin/env python3

import os
import argparse
import pathlib
import dendropy


def parse_args():
    parser = argparse.ArgumentParser(prog='reroot tree', \
        description='Run this program to reformat a sequence names in an alignment file and add the relevant date found in a csv.')
    parser.add_argument('--input_tree', help='Alignment file with dates added.')
    parser.add_argument('--output_tree', default='rerooted.tre', help='Alignment file with dates added.')
    return parser.parse_args()

def main():
    args = parse_args()

    tree = dendropy.Tree.get(path=args.input_tree, schema="newick")

    # mrca = tree.mrca(taxon_labels=["ERR2525602", "ERR2525603"])

    mrca = tree.mrca(taxon_labels=["SRR8171903"])

    tree.reroot_at_edge(mrca.edge, update_bipartitions=False)

    print(tree.as_ascii_plot())

    tree.write(path=args.output_tree, schema="newick")



if __name__ == '__main__':
    main()
