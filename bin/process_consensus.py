#!/usr/bin/env python
# Filter out the low occurence consensus sequences and produce 
# 1. hap.csv: Sequences of the dominant haplotype sequences
# 2. quals.csv: Highest quality score for each base in the dominant haplotype sequences

import csv
import dnaio
from glob import glob
import os
import argparse

def main(sample, region, fastqpath, tablepath ):

    consensustable = []
    with open(tablepath, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if len(row) == 3 and row[2] == region:
                consensustable.append({'hap': row[0], 'counts': int(row[1]), 'region': row[2]})

    if not consensustable:
        raise ValueError(f"No records found for region {region}")


    max_count = max(row['counts'] for row in consensustable)
    threshold = max_count * 0.05
    filtered = [row for row in consensustable if row['counts'] > threshold]

    haplist = [row['hap'] for row in filtered]
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
    parser.add_argument('--target', type=str, required=True, help='Target name')
    parser.add_argument('--fastqpath', type=str, required=True, help='Path to consensus FASTQ file')
    parser.add_argument('--tablepath', type=str, required=True, help='Path to consensus table CSV')    
    
    args = parser.parse_args()
    main(args.sample, args.target, args.fastqpath, args.tablepath)