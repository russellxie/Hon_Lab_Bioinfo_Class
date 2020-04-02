#!/usr/bin/env python3

import sys
import numpy as np

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



