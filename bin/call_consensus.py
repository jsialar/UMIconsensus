#!/usr/bin/env python
import argparse
import pysam
from glob import glob
import os
from collections import Counter
import gzip
import csv

class Fastq:
    def __init__(self, name):
        self.name = '@' + name
        self.seq = []
        self.qual = []
        self.umifamilysize = None
        
    def __repr__(self):
        return self.fastqstring()
    
    def __len__(self):
        return len(self.seq)
    
    def add(self, nuc, ascii):
        self.seq.append(nuc)
        self.qual.append(ascii)

    def seqstring(self):
        seq_string = ''.join(self.seq)
        return seq_string

    def fastqstring(self):
        seq_string = self.seqstring()
        qual_string = ''.join(self.qual)
        fastq_string = f'{self.name} {self.umifamilysize}\n{seq_string}\n+\n{qual_string}'
        return fastq_string


def main(inputpath, nonpolysnp, samplename, keep):
    fastq_list = []
    bamfiles = glob(os.path.join(inputpath, '*.bam'))
    bamfiles.sort()
    region = os.path.basename(bamfiles[0]).split('.')[0].split('_')[0]

    with open(nonpolysnp, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        posdf = [row for row in reader if row['Name'] == region]

    chrno = posdf[0]['Chromosome']
    pos_list = [row['End'] for row in posdf]

    for bamfile in bamfiles:
        idx = os.path.basename(bamfile).split('.')[0].split('_')[1]                
        name = f"{region}_{idx}"
        fastq = Fastq(name)
        bam = pysam.AlignmentFile(bamfile, "rb")
        num_alignments = bam.count()

        for pos in pos_list: 
            # get consensus for 1 position each time
            called = pysam.consensus("-r", f"{chrno}:{pos}-{pos}", "--show-del", "yes",
                                     "--format", "fastq", bamfile, catch_stdout=True)
            calledsplit = called.split()
            #TODO: Minihaps with Indels are thrown out because we assume the positions have been prescreened for no indels based on HGDP_1K.
            # However, a rare indel that is not present in HGDP_1K may be found in a new sample. 
            if len(calledsplit) < 4:
                # nothing emitted at all: position not covered by this family
                nuc, qual = ('U', '!')
            else:
                _, nuc, _, qual = calledsplit
                if len(nuc) > 1:
                    # samtools anchors an insertion to the preceding base, so a
                    # 1bp window returning >1 base means an insertion after pos
                    nuc, qual = ('I', '!')
                elif nuc == '*':
                    nuc, qual = ('D', '!')

            if not keep:
                if nuc in ('N', 'D', 'I', 'U'):
                    break
                
            fastq.add(nuc, qual)
        else:
            fastq.umifamilysize = num_alignments
            fastq_list.append(fastq)

    seqcount = Counter([(i.seqstring(), i.umifamilysize) for i in fastq_list])

    # Write to CSV
    with open(f'{region}.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for (seq, umifamilysize), count in seqcount.items():
            writer.writerow([seq, umifamilysize, count, region])

    # Write to fastq
    with gzip.open(f"{samplename}.{region}.fastq.gz", "wt") as fh:
        for fastq in fastq_list:
            fh.write(fastq.fastqstring() + "\n")

if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Call consensus sequences from BAM files.")
    parser.add_argument("--inputpath", type=str, help="Path to input BAM files. All BAM files should be for the same region. ")
    parser.add_argument("--nonpolysnp", type=str, help="Path to input snp positions files.")
    parser.add_argument("--samplename", type=str, help="Sample name")
    parser.add_argument("--keep", help="Keep sequences that contain non-substitution calls: "
                                       "N (ambiguous), D (deletion), I (insertion), U (no coverage)",
                        action='store_true')
    args = parser.parse_args()

    main(args.inputpath, args.nonpolysnp, args.samplename, args.keep)