import os
import sys
import gzip
import argparse
import numpy as np
import pandas as pd
from collections import Counter
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI
from Bio.Seq import reverse_complement


def read_index_map():
    idx_fn = '../../data/nextera_indices/miseq_indices_map.tsv'
    i7_namedict = {
        'Shruti i7': 'Shruthi i7_2.2',
        'FZ i7': 'i7_test',
    }
    combos = {}
    with open(idx_fn, 'rt') as f:
        i7_idx = f.readline().rstrip('\n').split('\t')[1:]
        col_idx = f.readline().rstrip('\n').split('\t')[1:]
        for i in range(11):
            fields = f.readline().rstrip('\n').split('\t')
            row = fields[0]
            for ic, i5name in enumerate(fields[1:]):
                if i5name == '':
                    continue

                i5name = 'FZi5_'+str(int(i5name[-2:]))
                i7name = i7_namedict[i7_idx[ic]]
                col = col_idx[ic]
                combos[(i7name, i5name)] = row+col

    return combos


def read_designed_indices():
    '''Read index sequences as designed

    NOTE: i7 ends up in the I1 fastq, i5 in the I2 fastq.
    '''
    idx_fn = '../../data/nextera_indices/miseq_indices.csv'
    idx_dict = {
        'i7': {},
        'i5': {},
    }
    with open(idx_fn, 'rt') as f:
        header = f.readline().rstrip('\n').split(',')
        for line in f:
            fields = line.rstrip('\n').split(',')
            name = fields[0]
            idx_seq = fields[3]

            if 'i5' in name:
                idx_dict['i5'][name] = idx_seq
            elif 'i7_' in name:
                if 'Shruthi' in line:
                    # This is an 8-base barcode, so fill the rest
                    seqlong = reverse_complement(fields[4])
                    tmp = seqlong.find(idx_seq)
                    idx_seq = seqlong[tmp: tmp+12]
                idx_dict['i7'][name] = idx_seq


    idx_map = read_index_map()
    idx_dict['combo'] = {}
    for (i7name, i5name), well in idx_map.items():
        # NOTE: i5 needs to be RC
        idx_dict['combo'][(i7name, i5name)] = \
                idx_dict['i7'][i7name]+'+'+reverse_complement(idx_dict['i5'][i5name])

    return idx_dict


def get_read_lengths():
    muxed_fdn = '../../data/fastq/20210601_MiSeq/muxed/'
    fnd = {
        'i1': f'{muxed_fdn}Undetermined_S0_I1_001.fastq.gz',
        'i2': f'{muxed_fdn}Undetermined_S0_I2_001.fastq.gz',
        'r1': f'{muxed_fdn}Undetermined_S0_R1_001.fastq.gz',
        'r2': f'{muxed_fdn}Undetermined_S0_R2_001.fastq.gz',
        }
    read_lengths = {}
    for key, fn in fnd.items():
        with gzip.open(fn, 'rt') as f:
            readiter = FGI(f)
            for name, seq, qual in readiter:
                read_lengths[key] = len(seq)
                break
    return read_lengths


def fill_samplesheet():
    idx_map = read_index_map()
    idx_dict = read_designed_indices()['combo']

    fn_empty = '../../data/raw/20210601_MiSeq/SampleSheet.csv'
    fn_filled = '../../data/raw/20210601_MiSeq/SampleSheet_filled.csv'
    with open(fn_empty, 'rt') as fe, open(fn_filled, 'wt') as ff:
        for line in fe:
            ff.write(line)

        for (i7name, i5name), seq_combo in idx_dict.items():
            seq_i7, seq_i5 = seq_combo.split('+')
            sampleID = idx_map[(i7name, i5name)]
            samplename = sampleID
            samplePlate = ''
            sampleWell = sampleID
            fields = [
                sampleID, samplename, samplePlate, sampleWell,
                i7name, seq_i7, i5name, seq_i5,
                # Project, description
                'image-seq_pilot1', '',
            ]
            line = ','.join(fields)+'\n'
            ff.write(line)


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--maxreads', type=int, default=100000)
    args = pa.parse_args()

    read_length_dict = get_read_lengths()

    idx_dict = read_designed_indices()

    fill_samplesheet()
    sys.exit()

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

    for key, val in counts.items():
        counts[key] = pd.Series(val).sort_values(ascending=False).to_frame(name='n_reads')
        counts[key]['designed'] = ''

    idx_corr = {
        #'i1': 'i7',
        #'i2': 'i5',
        'combo': 'combo',
    }
    for key, key2 in idx_corr.items():
        count = counts[key]
        idx_des = idx_dict[key2]
        for name, seq in idx_des.items():
            if seq in count.index:
                count.at[seq, 'designed'] = name
                continue
            seqrc = reverse_complement(seq)
            if seqrc in count.index:
                count.at[seqrc, 'designed'] = name+' (RC)'
                continue

            

    # Index 2: all good, mostly our i7 but ~25% Shruti's, extended by 4 bases
    # Index 1: all good except for 1 barcode (TCTTTCCCTACA), let's figure out if it's
    # at close Hamming distance from another barcode
    if False:
        i1mat = np.array(counts['i2'].index[:33].str.split('').str.slice(1, -1).tolist())
        dmat = np.zeros((33, 33), np.int32)
        for i, seq in enumerate(i1mat):
            for j, seq2 in enumerate(i1mat[:i]):
                dmat[i, j] = dmat[j, i] = (seq != seq2).sum()
        dmat = pd.DataFrame(dmat, index=counts['i2'].index[:33], columns=counts['i2'].index[:33])
        # The intruder index is not close to anything... ??? not a biggie, but what the heck?

    
