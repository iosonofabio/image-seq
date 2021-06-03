import os
import sys
import numpy as np
import pandas as pd
import subprocess as sp


def get_samplefiles():
    from collections import defaultdict

    fn_samplesheet = '../../data/raw/20210601_MiSeq/SampleSheet_filled.csv'
    df = pd.read_csv(fn_samplesheet, header=None).iloc[-32:].set_index(0)
    samplenames = list(df.index)
    
    fdn_fastq = '../../data/fastq/20210601_MiSeq/demuxed/image-seq_pilot1'
    # All fastqs in one folder
    sample_dict = defaultdict(dict)
    fns = os.listdir(fdn_fastq)
    for fn in fns:
        sn = fn.split('_')[0]
        if '_R1_001' in fn:
            key = 'R1'
        else:
            key = 'R2'
        sample_dict[sn][key] = fdn_fastq+'/'+fn

    return sample_dict


if __name__ == '__main__':

    fdn_ref = '../../data/genome_reference/STAR_DIR/'

    sample_dict = get_samplefiles()
    
    for sn, dic in sample_dict.items():
        fn_r1 = dic['R1']
        fn_r2 = dic['R2']
        fdn_out = f'../../data/bam/20210601_MiSeq/{sn}/'

        print(f'Sample {sn}')
        os.makedirs(fdn_out, exist_ok=True)

        call = '''
        STAR
        --runMode alignReads
        --runThreadN 16
        --genomeDir '''+fdn_ref+'''
        --readFilesCommand zcat
        --outFileNamePrefix '''+fdn_out+'''
        --outSAMtype BAM Unsorted
        --outSAMunmapped None
        --outFilterMismatchNmax 10
        --quantMode GeneCounts
        --readFilesIn '''+fn_r1+' '+fn_r2
        # this was 100 when generated and it's not worth making a
        # new STAR_DIR just for this
        #--sjdbOverhang 73

        call = call.replace('\n', ' ').replace('  ', ' ')
        print(call)
        if False:
            sp.run(call, shell=True, check=True)

    counts = {}
    for sn in sample_dict:
        fn_counts = f'../../data/bam/20210601_MiSeq/{sn}/ReadsPerGene.out.tab'
        # The second column is the unstranded one
        count = pd.read_csv(
            fn_counts, sep='\t', index_col=0, usecols=[0, 1],
            squeeze=True, header=None)
        counts[sn] = count
    counts = pd.DataFrame(counts)
    counts.index.name = 'GeneName'

    fdn_counts = '../../data/counts/20210601_MiSeq/'
    os.makedirs(fdn_counts, exist_ok=True)
    counts.to_csv(fdn_counts+'counts.tsv', sep='\t', index=True)



    
