#!/usr/bin/env python3

import dendropy
import argparse
from dendropy.calculate import treecompare

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--true_tree')
    parser.add_argument('--new_tree')
    return parser.parse_args()

def main():
    args = parse_args()

    read_true = open(args.true_tree,'r').read()
    read_new = open(args.new_tree,'r').read()

    # establish common taxon namespace
    tns = dendropy.TaxonNamespace()

    # # ensure all trees loaded use common namespace
    # true = dendropy.Tree.get(
    #     data=read_true,
    #     schema='newick',
    #     taxon_namespace=tns, preserve_underscores=True)
    # new = dendropy.Tree.get(
    #     data=read_new,
    #     schema='newick',
    #     taxon_namespace=tns, preserve_underscores=True)
    #
    # print(treecompare.symmetric_difference(true, new))

    true = dendropy.Tree.get(
        data=read_true,
        schema='newick',
        preserve_underscores=True)
    new = dendropy.Tree.get(
        data=read_new,
        schema='newick',
        preserve_underscores=True)

    print(len(true.taxon_namespace))
    print(len(new.taxon_namespace))

    tree_1_taxa = []
    tree_2_taxa = []
    taxa_to_prune = []

    for taxon in true.leaf_node_iter():
        # print(taxon.taxon)
        tree_1_taxa.append(taxon.taxon)

    for taxon in new.leaf_node_iter():
        # print(taxon.taxon)
        tree_2_taxa.append(taxon.taxon)

    tree_1_size = len(tree_1_taxa)
    tree_2_size = len(tree_2_taxa)

    print(tree_1_size)
    print(tree_2_size)

    if tree_1_size > tree_2_size:
        print('true bigger than new')
        for tax in tree_1_taxa:
            if tax not in tree_2_taxa:
                taxa_to_prune.append(tax.label)

        print(true.as_ascii_plot())
        print(true)

        true.prune_taxa_with_labels(taxa_to_prune)
        print(true)

        print(true.as_ascii_plot())

    elif tree_2_size > tree_1_size:
        print('new bigger than true')
        for tax in tree_2_taxa:
            if tax not in tree_1_taxa:
                taxa_to_prune.append(tax.label)

        print(new.as_ascii_plot())
        print(new)

        new.prune_taxa_with_labels(taxa_to_prune)
        print(new)

        print(new.as_ascii_plot())

    # print(true)
    #     true.attach_taxon_namespace(tns)
    #
    # print(treecompare.symmetric_difference(true, new))



    ## Unweighted Robinson-Foulds distance
    # print(treecompare.symmetric_difference(true, new))

if __name__ == '__main__':
    main()
