import os
import sys
import subprocess as sp


if __name__ == '__main__':

    call = '''
    bcl2fastq
    --runfolder-dir ../../data/raw/20210601_MiSeq/210531_M02082_0067_000000000-JPL2D/
    --output-dir ../../data/fastq/20210601_MiSeq/muxed/
    --loading-threads 16
    --processing-threads 16
    --writing-threads 0
    --create-fastq-for-index-reads
    --no-lane-splitting
    '''
    call = call.replace('\n', ' ').replace('  ', ' ')

    sp.run(call, shell=True, check=True)
