# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/06/21
content:    Make dataset file from counts and using gene names.
'''
import os
import sys
import numpy as np
import pandas as pd

sys.path.append('/home/fabio/university/postdoc/singlet')
import singlet


if __name__ == '__main__':

    print('Read counts and metadata from file')
    fdn_data = '../../data/counts/20210601_MiSeq/'
    fn_dataset = f'{fdn_data}/raw.h5ad'
    ds = singlet.Dataset(
        dataset={'path': fn_dataset},
    )

    dsm = ds.query_features_by_name(ds.featurenames[:-4])
    dsm.obs['coverage'] = dsm.counts.sum(axis=0)
    dsm.counts.normalize(inplace=True)

    avg_exp = dsm.counts.mean(axis=1)
