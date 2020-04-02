#!/usr/bin/env python3
import sys
import subprocess
import argparse

import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--fq1', help='fastq1 - mRNA reads')
args = parser.parse_args()

def bin_generator():
    for line in sys.stdin:
        current_chr = line.strip().split('\t')[0]
        left = int(line.strip().split('\t')[1])
        right = int(line.strip().split('\t')[2])

        left_point = round((left + right)/2) - 2000

        region_id = current_chr + '_' + str(left) + '_' + str(right)
        
        for i in range(0, 80):
            print(region_id + '-' + str(i), end = '\t')
            bin_left = left_point + (i*50)
            bin_right = left_point + (i+1)*50 - 1

            bin_list = [current_chr, str(bin_left), str(bin_right)]
            
            print('\t'.join(bin_list) + '\t' + '+')

def mat_generator():
    total_reads = int(sys.argv[1])

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
    # Count the reads within each bin by using featureCounts
    command = """featureCounts -a {0}\
                -F SAF\
                -o enhancer_bins.hg38.{1}.counts\
                -O\
                {2}""".format(saf_file, marker, bam_file)

    rc = subprocess.call(command, shell=True)

    # Reshape the results, calculate the CPM and output the matrix
    ###