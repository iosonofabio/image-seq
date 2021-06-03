import os
import sys
import gzip
import argparse
import numpy as np
from collections import Counter
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--maxreads', type=int, default=100000)
    args = pa.parse_args()

    maxreads = args.maxreads

    muxed_fdn = '../../data/fastq/20210601_MiSeq/muxed/'
    fn_i1 = f'{muxed_fdn}Undetermined_S0_I1_001.fastq.gz'
    fn_i2 = f'{muxed_fdn}Undetermined_S0_I2_001.fastq.gz'

    counts = {
        'i1': Counter(),
        'i2': Counter(),
        'combo': Counter(),
    }

    with gzip.open(fn_i1, 'rt') as f1, gzip.open(fn_i2, 'rt') as f2:
        readiter1 = FGI(f1)
        readiter2 = FGI(f2)
        for ir, (idx1, idx2) in enumerate(zip(readiter1, readiter2)):
            if ir == maxreads:
                break

            # FIXME: add quality filter?
            name1, seq1, qual1 = idx1
            name2, seq2, qual2 = idx2

            counts['i1'][seq1] += 1
            counts['i2'][seq2] += 1
            counts['combo'][seq1+'+'+seq2] += 1


