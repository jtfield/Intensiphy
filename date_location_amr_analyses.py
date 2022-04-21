#!/usr/bin/env python3

import os
import argparse
import pathlib
import shutil
import pandas as pd
import datetime
import dateutil
import re
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(prog='Analyses on collection date, location and AMR genes', \
        description='Analyses of the tree most important surveillance date points.')
    parser.add_argument('--csv_file', help='Input CSV file containing all the relevant information')
    # parser.add_argument('--accs_file', default=False, help='Accession file if accession_method is set to USER_INPUT')
    return parser.parse_args()

def main():
    args = parse_args()

    read_csv = pd.read_csv(args.csv_file)



###################################################
# Building sample count per bin per location functions
    sort_by_year = only_years(read_csv)

    get_year_ranges = year_ranges(read_csv)

    analyze_by_region = get_regions(sort_by_year)

    make_plots(analyze_by_region)

##################################################





def only_years(df):
    year_regex = '\d\d\d\d'

    compile_year_regex = re.compile(year_regex)

    for index, row in df.iterrows():
        # print(index)
        collect_date = row['Collection date']
        find_date = re.match(compile_year_regex, collect_date)
        if find_date:
            # print(find_date)
            # print(df.at[index, 'Collection date'])
            df.at[index, 'Collection date'] = int(find_date[0])
        # print("##############")


    sorted_df = df.sort_values(by='Collection date')
    # print(sorted_df['Collection date'])
    return sorted_df


def year_ranges(df_):

    d_range = df_['Collection date'].value_counts()
    print(d_range)

    d_loc = df_['Location'].value_counts()
    print(d_loc)

def get_regions(df_):

    bins = [2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021]
    regions = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

    region_regex = 'USA.+\s(\d+)'

    compile_region_regex = re.compile(region_regex)

    output_dict = {}

    for num, date_range in enumerate(bins):
        output_dict[num] = {}
        range_dict = output_dict[num]
        for region in regions:
            output_dict[num][region] = 0

    print(output_dict)


    for num_1, date_range in enumerate(bins):
        # r_1 = 0
        # r_2 = 0
        # r_3 = 0
        # r_4 = 0
        # r_5 = 0
        # r_6 = 0
        # r_7 = 0
        # r_8 = 0
        # r_9 = 0
        # r_10 = 0

        for index, row in df_.iterrows():
            region = row['Location']
            date = int(row['Collection date'])
            # print(region)
            find_region = re.match(compile_region_regex, region)
            if find_region:
                found_region = find_region.group(1)
                # print(found_region)
                # print(type(found_region))
                if date == date_range:
                    # print(found_region)
                    # print(date)
                    output_dict[num_1][find_region.group(1)]+=1

                # for num_2, date_range in enumerate(bins):
                #     if date >= date_range[0] and date <= date_range[1]:
                #         output_dict[num][find_region.group(1)]+=1

    return output_dict





    #
    # for index, row in df_.iterrows():
    #     region = row['Location']
    #     date = int(row['Collection date'])
    #     # print(region)
    #     find_region = re.match(compile_region_regex, region)
    #     if find_region:
    #         print(find_region.group(1))
    #         found_region = find_region.group(1)
    #         for num, date_range in enumerate(bins):
    #             if date >= date_range[0] and date <= date_range[1]:
    #                 output_dict[num][find_region.group(1)]


def make_plots(dict_):

    for key, sub_dict in dict_.items():
        print(sub_dict)

        plt.bar(list(sub_dict.keys()), sub_dict.values(), color='g')
        plt.show()






if __name__ == '__main__':
    main()
