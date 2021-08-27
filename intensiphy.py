#! /usr/bin/python3

import os
import argparse
import re
import csv
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(prog='Intensiphy', \
        description='Automate downloading of high-throughput sequence data and updating of alignments using Extensiphy.')
    parser.add_argument('--align_file')
    parser.add_argument('--cores')
    parser.add_argument('--accession_csv')
    parser.add_argument('--update_method')
    return parser.parse_args()

def read_csv(csv_file):
    """Reads the accession file provided by the users."""
    csv = pd.read_csv(csv_file)

    return csv

def calulate_cores(set_cores):
    """Organizes and calulates the cores for use with Extensiphy."""
    output = []
    total_cores = int(set_cores)
    runs = 0
    cores_per_run = 0
    
    if total_cores > 2 and total_cores < 10:
        runs = 2
        cores_per_run = total_cores / 2
        
    elif total_cores >= 10:
        runs = total_cores / 2
        cores_per_run = 2
    
    assert runs * cores_per_run <= total_cores
    
    output.append(int(runs))
    output.append(int(cores_per_run))

    return output
        


def main():
    args = parse_args()

    # calculate the core organization to pass to Extensiphy
    get_cores = calulate_cores(args.cores)
    print(get_cores)

    

if __name__ == '__main__':
    main()