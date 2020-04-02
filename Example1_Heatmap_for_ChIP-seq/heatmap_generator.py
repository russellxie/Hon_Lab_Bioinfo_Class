#!/usr/bin/env python3

import os
import sys
import os.path
import subprocess
import argparse

from os import path
import pandas as pd
import numpy as np

def bin_generator(bed_file, sample_name, bin_num: int, bin_size: int):
    #generate the bin file name
    bin_file = "Region_bins_%s_%d_%d.saf"%(sample_name, bin_num, bin_size)

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

def mat_generator(marker):

    #get the total reads from the featureCounts output file
    summary_df = pd.read_csv("enhancer_bins.hg38.%s.counts.summary"%(marker), skiprows=1, index_col=0, header=False)
    total_reads = np.sum(summary_df.values)

    #read the output of featureCounts and convert it to a matrix
    line1 = sys.stdin.readline().strip()
    region_id = line1.strip().split('\t')[0].split('-')[0]
    reads = list(map(int, line1.strip().split('\t')[6:]))
    cpm = np.sum(reads) / total_reads * 1e6

    print(region_id + '\t' + str(cpm), end = '')

    for line in sys.stdin:
        current_region_id = line.strip().split('\t')[0].split('-')[0]
        chr_name = line.strip().split('\t')[1]
        reads = list(map(int, line.strip().split('\t')[6:]))

        cpm = np.sum(reads) / total_reads * 1e6
        
        if current_region_id == region_id:
            print ('\t', end ='')
            print (cpm, end = '')
        else:
            print ('\n', end ='')
            print (current_region_id + '\t' + str(cpm), end = '')
            region_id = current_region_id

    print('\n', end='')

def plot_heatmap():
    ###
    # Generate the bins
    bin_generator(args.region, args.sample_name, args.bin_num, args.bin_size)

    # Count the reads within each bin by using featureCounts
    command = """featureCounts -a {0}\
                -F SAF\
                -o enhancer_bins.hg38.{1}.counts\
                -O\
                {2}""".format(saf_file, args.marker, args.bam)

    rc = subprocess.call(command, shell=True)

    # Reshape the results, calculate the CPM and output the matrix
    ###

parser = argparse.ArgumentParser()
parser.add_argument('--region', help='region info in bed format.')
parser.add_argument('--sample_name', help='sample_name for the file.')
parser.add_argument('--bam', help='bam file input.')
parser.add_argument('--bin_num', help='Bin numbers')
parser.add_argument('--bin_size', help='Bin size')
args = parser.parse_args()