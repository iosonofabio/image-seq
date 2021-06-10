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
import anndata


if __name__ == '__main__':


    print('Load raw count data')
    fn_counts = '../../data/counts/20210601_MiSeq/counts.tsv'
    counts_raw = pd.read_csv(fn_counts, index_col=0, sep='\t').T.astype(np.float32)
    n_samples = counts_raw.shape[0]

    print('Load biomart conversion table')
    fn_biomart = '../../data/ensembl/mart_export.txt'
    gene_dict = pd.read_csv(fn_biomart, sep='\t', index_col=0)
    gene_dict['Gene name'] = gene_dict['Gene name'].fillna('')
    gene_dict = gene_dict.loc[gene_dict['Gene name'] != '']

    print('Compute intersection')
    idx = np.intersect1d(gene_dict.index, counts_raw.columns)
    gene_dict = gene_dict.loc[idx]
    counts_raw = counts_raw[list(counts_raw.columns[:4]) + list(gene_dict.index)]

    print('Extract gene names')
    genes = gene_dict['Gene name'].unique()
    genes = list(genes) + list(counts_raw.columns[:4])
    n_genes = len(genes)

    print('Prepare empty cout table')
    counts = pd.DataFrame(
        np.zeros((n_samples, n_genes), np.float32),
        index=counts_raw.index,
        columns=genes,
        )
    print('Add other features first')
    for fea in genes[-4:]:
        counts[fea] = counts_raw[fea]

    print('Add human genes')
    for eid, row in counts_raw.T.iloc[4:].iterrows():
        gene = gene_dict.at[eid, 'Gene name']
        counts.loc[:, gene] += row


    print('Complement feature sheet')
    for gene in genes[-4:]:
        gene_dict.loc[gene] = ''
        gene_dict.at[gene, 'Gene name'] = gene
    gene_dict.rename(columns={'Chromosome/scaffold name': 'Chromosome'}, inplace=True)
    var = gene_dict.drop_duplicates('Gene name').set_index('Gene name').loc[counts.columns]

    print('Store to file')
    fn_counts = '../../data/counts/20210601_MiSeq/raw_with_gene_names.tsv'
    counts.to_csv(fn_counts, sep='\t', index=True)
    fn_varmeta = '../../data/counts/20210601_MiSeq/raw_with_gene_names_var.tsv'
    var.to_csv(fn_varmeta, sep='\t', index=True)

    if True:
        X = counts.values

        print('Make anndata')
        adata = anndata.AnnData(
            X=X,
            var=var,
            obs=pd.DataFrame([], index=counts.index),
        )

        print('Write to file')
        fn_counts = '../../data/counts/20210601_MiSeq/raw.h5ad'
        adata.write(fn_counts)
