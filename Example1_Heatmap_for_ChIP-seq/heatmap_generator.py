#!/usr/bin/env python3

import os
import sys
import os.path
import subprocess
import argparse

from os import path
import pandas as pd
import numpy as np

def bin_generator(bed_file, bin_file, bin_num: int, bin_size: int):

    #check if the output file exists
    if path.exists(bin_file):
        os.remove(bin_file)

    outfh = open(bin_file, 'w')
    with open(bed_file) as bedfh:
        for line in bedfh:
            current_chr = line.strip().split('\t')[0]
            left = int(line.strip().split('\t')[1])
            right = int(line.strip().split('\t')[2])

            if bin_num % 2 == 0:
                left_point = int((left + right)/2) - (bin_num * bin_size) / 2
            else:
                left_point = left - int(bin_num / 2) * bin_size

            region_id = current_chr + '_' + str(left) + '_' + str(right)
        
        for i in range(0, bin_num):
            print(region_id + '-' + str(i), end = '\t')
            bin_left = left_point + (i*bin_size)
            bin_right = left_point + (i+1)*bin_size - 1

            bin_list = [current_chr, str(bin_left), str(bin_right)]
            
            print('\t'.join(bin_list) + '\t' + '+', file = outfh)
    outfh.close()

def mat_converter(sample_name, count_file, matrix_file):

    #get the total reads from the featureCounts output file
    summary_df = pd.read_csv(f"enhancer_bins.hg38.{sample_name}.counts.summary",\
                            skiprows=1,\
                            index_col=0,\
                            header=False)

    total_reads = np.sum(summary_df.values)

    #read the output of featureCounts and convert it to a matrix
    line1 = sys.stdin.readline().strip()
    region_id = line1.strip().split('\t')[0].split('-')[0]
    reads = list(map(int, line1.strip().split('\t')[6:]))
    cpm = np.sum(reads) / total_reads * 1e6

    print(region_id + '\t' + str(cpm), end = '')

    outfh = open(matrix_file, 'w')

    with open(count_file) as cf:
        for line in cf:
            current_region_id = line.strip().split('\t')[0].split('-')[0]
            _chr_name = line.strip().split('\t')[1]
            reads = list(map(int, line.strip().split('\t')[6]))

            cpm = np.sum(reads) / total_reads * 1e6
            
            if current_region_id == region_id:
                print ('\t', end ='', file=outfh)
                print (cpm, end = '', file=outfh)
            else:
                print ('\n', end ='', file=outfh)
                print (current_region_id + '\t' + str(cpm), end = '', file=outfh)
                region_id = current_region_id

    print('\n', end='', file=outfh)
    outfh.close()

def plot_heatmap(sample_name, bam_file, bin_num, bin_size):
    #define the output file names
    bin_file = f"Region_bins.{sample_name}.{str(bin_num)}X{str(bin_size)}.saf"
    count_file = f"Region_bins.{sample_name}.{str(bin_num)}X{str(bin_size)}.counts"
    matrix_file = f"Region_bins.{sample_name}.{str(bin_num)}X{str(bin_size)}.matrix.csv"

    # Generate the bins
    bin_generator(args.region, bin_file, args.bin_num, args.bin_size)

    # Count the reads within each bin by using featureCounts
    command = f"featureCounts -a {bin_file}"\
              f"-F SAF"\
              f"-o {count_file}"\
              f"-O"\
              f"{bam_file}"
                
    rc = subprocess.call(command, shell=True)

    # Reshape the results, calculate the CPM and output the matrix
    mat_converter(sample_name, count_file, matrix_file)

parser = argparse.ArgumentParser()
parser.add_argument('--region', help='region info in bed format.')
parser.add_argument('--sample_name', help='sample_name for the file.')
parser.add_argument('--bam', help='bam file input.')
parser.add_argument('--bin_num', help='Bin numbers')
parser.add_argument('--bin_size', help='Bin size')
args = parser.parse_args()

plot_heatmap(args.sample_name, args.bam, args.bin_num, args.bin_size)