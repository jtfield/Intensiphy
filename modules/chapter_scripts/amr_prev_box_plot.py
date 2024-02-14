#!/usr/bin/env python3

import os
import argparse
import pandas as pd
# import sys
# import json
# import dendropy

# from opentree import OT
# from opentree.annotations import *
# from annotations import write_itol_heatmap




def parse_args():
    parser = argparse.ArgumentParser(prog='amr_lister.py', \
        description='make itol annotations for ip data')
    #parser.add_argument('--csv', help='csv file with AMR gene values.')
    return parser.parse_args()

def main():
    args = parse_args()

    amr_genes = ['fbar', ]




if __name__ == '__main__':
    main()
