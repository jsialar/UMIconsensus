#!/usr/bin/env python
# Filter out the low occurence consensus sequences and produce 
# 1. hap.csv: Sequences of the dominant haplotype sequences
# 2. quals.csv: Highest quality score for each base in the dominant haplotype sequences

import csv
import dnaio
from glob import glob
import os
import argparse
import pandas as pd

def main(sample, noisethreshold, region, fastqpath, tablepath ):

    consensustable_detailed = pd.read_csv(tablepath, names=['hap', 'umifamilysize', 'count', 'region'])
    consensustable = consensustable_detailed.groupby('hap').agg({'count':'sum'}).reset_index()
    frac = consensustable['count'] / consensustable['count'].max()
    consensustable_filtered = consensustable[frac >= noisethreshold]

    haplist = consensustable_filtered['hap'].to_list()
    hapcount = len(haplist)
    haplength = len(haplist[0])

    score = [[0] * haplength] * hapcount

    reader = dnaio.open(fastqpath)    

    # Find max quality score for every bases
    for record in reader:
        for idx in range(hapcount):
            if record.sequence == haplist[idx]:
                for i in range(haplength):
                    newqual = ord(record.qualities[i]) - 33
                    score[idx][i] = newqual if score[idx][i] < newqual else score[idx][i]
                break

    with open('qual.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        for idx, row in enumerate(score):
            writer.writerow([f"{region}_{idx}"] + row)

    with open('hap.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        for idx, row in enumerate(haplist):
            writer.writerow([f"{region}_{idx}"] + [row])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process consensus FASTQ and table.")
    parser.add_argument('--sample', type=str, required=True, help='Sample name')
    parser.add_argument('--noisethreshold', type=float, required=True, help='Noise threshold')
    parser.add_argument('--target', type=str, required=True, help='Target name')
    parser.add_argument('--fastqpath', type=str, required=True, help='Path to consensus FASTQ file')
    parser.add_argument('--tablepath', type=str, required=True, help='Path to consensus table CSV')    
    
    args = parser.parse_args()
    main(args.sample, args.noisethreshold, args.target, args.fastqpath, args.tablepath)