#!/usr/bin/env python3

import dendropy
import argparse
from dendropy.calculate import treecompare

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--t1')
    parser.add_argument('--t2')
    return parser.parse_args()

def prune_trees_to_match(t1, t2):
    """requires two trees with a common taxon namespace"""
    tree_1_taxa = set()
    tree_2_taxa = set()

    for tip in t1.leaf_node_iter():
        tree_1_taxa.add(tip.taxon)

    for tip in t2.leaf_node_iter():
        tree_2_taxa.add(tip.taxon)

    shared_taxa = tree_1_taxa.intersection(tree_2_taxa)

    assert(len(shared_taxa) >= 1)
    print("These two tree have {s} shared taxa".format(s=len(shared_taxa)))

    t1.retain_taxa(shared_taxa)
    t2.retain_taxa(shared_taxa)
    return(t1, t2)

def main():
    args = parse_args()
    tns = dendropy.TaxonNamespace()
    # # ensure all trees loaded use common namespace
    t1 = dendropy.Tree.get_from_path(
         src=args.t1,
         schema='newick',
         taxon_namespace=tns, preserve_underscores=True)
    t2 = dendropy.Tree.get_from_path(
         src=args.t2,
         schema='newick',
         taxon_namespace=tns, preserve_underscores=True)

    t1, t2 = prune_trees_to_match(t1, t2)
    print("Symmetric difference is {}".format(treecompare.symmetric_difference(t1, t2)))



if __name__ == '__main__':
    main()
