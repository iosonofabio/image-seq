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

import matplotlib.pyplot as plt
import seaborn as sns

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

    cellnames_imaged = [
        'L5', 'E5', 'F5', 'C5', 'D5', 'G5', 'K5', 'I5', 'H5', 'J5',
    ]
    cellnames_not = [x for x in dsm.samplenames if x not in cellnames_imaged]
    dsmi = dsm.query_samples_by_name(cellnames_imaged)
    dsmn = dsm.query_samples_by_name(cellnames_not)

    from scipy.stats import ks_2samp
    xi = dsmi.obs['coverage'].values
    xn = dsmn.obs['coverage'].values
    res = ks_2samp(xi, xn)
    pval = res[1]

    print('Plot coverage')
    fig, ax = plt.subplots(figsize=(3, 3))
    # all cells
    x = np.sort(dsm.obs['coverage'])
    y = 1.0 - np.linspace(0, 1, len(x))
    ax.plot(x, y, 'o-', color='k', lw=2, label='all')

    # imaged
    x = np.sort(dsmi.obs['coverage'])
    y = 1.0 - np.linspace(0, 1, len(x))
    ax.plot(x, y, 'o-', color='tomato', lw=2, label='imaged')

    # imaged
    x = np.sort(dsmn.obs['coverage'])
    y = 1.0 - np.linspace(0, 1, len(x))
    ax.plot(x, y, 'o-', color='grey', lw=2, label='not imaged')

    ax.set_xlabel('N. uniquely mapped reads')
    ax.set_ylabel('Fraction of cells > x')
    ax.set_xscale('log')
    ax.grid(True)
    ax.legend(loc='lower left')

    ax.text(0.88, 0.93,
            '$P_{KS} = '+'{:.2f}'.format(pval)+'$',
            ha='right', va='top',
            transform=ax.transAxes,
            )

    fig.tight_layout()


    print('Check MT genes')
    mt_genes = dsm.featurenames[dsm.featurenames.str.startswith('MT-')]

    fig, ax = plt.subplots(figsize=(3.3, 3))
    # all cells
    percent_mito = dsm.counts.loc[mt_genes].sum(axis=0) * 1e-4
    x = np.sort(percent_mito)
    y = 1.0 - np.linspace(0, 1, len(x))
    ax.plot(x, y, 'o-', color='k', lw=2, label='all')

    # imaged
    percent_mito = dsmi.counts.loc[mt_genes].sum(axis=0) * 1e-4
    x = np.sort(percent_mito)
    y = 1.0 - np.linspace(0, 1, len(x))
    ax.plot(x, y, 'o-', color='tomato', lw=2, label='imaged')

    # imaged
    percent_mito = dsmn.counts.loc[mt_genes].sum(axis=0) * 1e-4
    x = np.sort(percent_mito)
    y = 1.0 - np.linspace(0, 1, len(x))
    ax.plot(x, y, 'o-', color='grey', lw=2, label='not imaged')

    ax.set_xlabel('Mitochonrdial reads [%]')
    ax.set_ylabel('Fraction of cells > x')
    ax.set_xscale('log')
    ax.set_xlim(0.1, 100)
    ax.grid(True)
    ax.legend(loc='lower left')

    fig.tight_layout()
