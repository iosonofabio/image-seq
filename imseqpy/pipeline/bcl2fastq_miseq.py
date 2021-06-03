import os
import sys
import argparse
import subprocess as sp


if __name__ == '__main__':

    pa = argparse.ArgumentParser()
    pa.add_argument('--demux', action='store_true')
    args = pa.parse_args()
    
    call = '''
    bcl2fastq
    --runfolder-dir ../../data/raw/20210601_MiSeq/210531_M02082_0067_000000000-JPL2D/
    --loading-threads 16
    --processing-threads 16
    --writing-threads 0
    --no-lane-splitting
    '''
    call = call.replace('\n', ' ').replace('  ', ' ')

    if args.demux:
        call += ' --sample-sheet ../../data/raw/20210601_MiSeq/SampleSheet_filled.csv'
        call += ' --output-dir ../../data/fastq/20210601_MiSeq/demuxed/'
    else:
        call += ' --create-fastq-for-index-reads'
        call += ' --output-dir ../../data/fastq/20210601_MiSeq/muxed/'


    sp.run(call, shell=True, check=True)
