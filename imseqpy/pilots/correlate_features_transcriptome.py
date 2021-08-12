# vim: fdm=indent
'''
author:     Fabio Zanini
date:       11/06/21
content:    Correlate image features and transcriptome.
'''
import os
import sys
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
import singlet



def load_features():
    import pickle
    fdn_images = '../../data/hyperspectral/20210601_MiSeq/'
    fn_feas = fdn_images+'features.pkl'
    with open(fn_feas, 'rb') as f:
        features = pickle.load(f)

    return features


if __name__ == '__main__':

    features = load_features()
    # FIXME: bad definition ;-)
    features['eccentricity'] -= 1

    wls = features['wavelengths'].iloc[0]
    wlx, wly = list(zip(*wls))
    wlx = list(np.unique(wlx))
    wly = list(np.unique(wly))

    fig, axs = plt.subplots(2, 5, figsize=(15, 3.6), sharex=True, sharey=True)
    axs = axs.ravel()
    vmax = np.vstack(features['spectrum']).max()
    for i, ax in enumerate(axs):
        feas = features.iloc[i]
        spectrum = feas['spectrum']
        ax.set_title(feas['well'])
        for j, ((wlxj, wlyj), val) in enumerate(zip(wls, spectrum)):
            x = wlx.index(wlxj)
            y = wly.index(wlyj)
            val_norm = (val + 1) / (vmax + 1)
            s = 0.05 + 0.45 * val_norm
            alpha = 0.3 + 0.7 * val_norm
            h = plt.Circle((x, y), s, color='red', alpha=alpha)
            ax.add_artist(h)

        ax.set_yticks(np.arange(len(wly)))
        ax.set_xticks(np.arange(len(wlx)))
        ax.set_xlim(-0.6, len(wlx) - 0.4)
        ax.set_ylim(-0.6, len(wly) - 0.4)
        if (i % 5) == 0:
            ax.set_yticklabels(wly)
        if i >= 5:
            ax.set_xticklabels(wlx, rotation=90)
    fig.text(0.52, 0.02, 'Excitation [nm]', ha='center')
    fig.text(0.02, 0.52, 'Emission [nm]', rotation=90, va='center')
    fig.suptitle('Spectra')
    fig.tight_layout(rect=(0.04, 0.05, 1, 0.995))
    fig.savefig('../../figures/MiSeq_pilot/spectra.png')

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

    print('Restrict to imaged cells')
    dsim = dsm.query_samples_by_name(features['well'])
    features.set_index('well', inplace=True, drop=False)
    for col in features.columns:
        dsim.obs[col] = features.loc[dsim.samplenames, col]

    print('Correlate with some simple features')
    feas = ['area', 'eccentricity']
    labels = ['Area [$px^2$]', 'Eccentricity (y/x - 1)']
    corr = dsim.correlation.correlate_features_phenotypes(feas, fillna=0)
    n_genes = 3
    colors = sns.color_palette('husl', n_colors=n_genes * 2)
    colors = [colors[:n_genes], colors[n_genes:]]
    fig, axs = plt.subplots(2, 2, figsize=(6, 6))
    for i, (fea, ax_row) in enumerate(zip(feas, axs)):
        genes_both = [
            corr.nlargest(n_genes, fea).index.tolist(),
            corr.nsmallest(n_genes, fea).index.tolist(),
            ]
        for j, (ax, genes) in enumerate(zip(ax_row, genes_both)):
            for ig, gene in enumerate(genes):
                x = dsim.obs[fea]
                y = dsim.counts.loc[gene] + 0.1
                idx = np.argsort(x)
                ax.plot(
                    x[idx], y[idx], 'o-', alpha=0.8, label=gene, lw=2,
                    color=colors[j][ig],
                    )
            ax.legend(fontsize=8)
            ax.set_xlabel(labels[i])
            if j == 0:
                ax.set_ylabel('Gene exp [cpm]')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.grid(True)
    axs[0, 0].set_title('Positive correlation')
    axs[0, 1].set_title('Negative correlation')
    fig.tight_layout()
    fig.savefig('../../figures/MiSeq_pilot/simple_corrs.png')
