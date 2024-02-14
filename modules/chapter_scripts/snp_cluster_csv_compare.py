#!/usr/bin/env python3

# import dendropy
import argparse
# from dendropy.calculate import treecompare
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv_1')
    parser.add_argument('--csv_2')
    return parser.parse_args()

def main():
    args = parse_args()

    read_csv_1 = pd.read_csv(args.csv_1)

    read_csv_2 = pd.read_csv(args.csv_2)

    csv_1_accessions = list(read_csv_1['Run'])
    len_acc_1 = len(csv_1_accessions)

    csv_2_accessions = list(read_csv_2['Run'])
    len_acc_2 = len(csv_2_accessions)

    taxa_to_prune = []

    print("length of csv_1 accesions: ", len_acc_1)
    print("length of csv_2 accesions: ", len_acc_2)

    if len_acc_1 > len_acc_2:
        print('one bigger than two')
        for tax in csv_1_accessions:
            # print(tax)
            if tax not in csv_2_accessions:
                taxa_to_prune.append(tax)

        print(taxa_to_prune)
        print(len(taxa_to_prune))

    elif len_acc_2 > len_acc_1:
        print('two bigger than one')
        for tax in csv_2_accessions:
            # print(tax)
            if tax not in csv_1_accessions:
                taxa_to_prune.append(tax)

        print(taxa_to_prune)
        print(len(taxa_to_prune))

    elif len_acc_2 == len_acc_1:
        print('datasets equal')
        for tax in csv_2_accessions:
            # print(tax)
            if tax not in csv_1_accessions:
                taxa_to_prune.append(tax)

        print(taxa_to_prune)
        print(len(taxa_to_prune))



if __name__ == '__main__':
    main()
