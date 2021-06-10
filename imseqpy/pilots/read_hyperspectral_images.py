# vim: fdm=indent
'''
author:     Fabio Zanini
date:       10/06/21
content:    Parse hyperspectral H5 from Abbas.
'''
import os
import sys
import numpy as np
import scipy as sp
import pandas as pd
import h5py

import matplotlib.pyplot as plt


if __name__ == '__main__':


    fdn_images = '../../data/hyperspectral/20210601_MiSeq/'
    fns = [
        'MCF7_C3.mat',
        'MCF7_C4.mat',
        'MCF7_N4.mat',
    ]

    fn = fdn_images+fns[0]
    with h5py.File(fn) as f:
        print('Read metadata about excitation/emission lambdas')
        keys = ['excitationWavelength', 'emission']
        wls = {key: [] for key in keys}
        for key in keys:
            tmp = f['HAC_Image']['imageStruct']['protocol']['channel'][key]
            n_colors = tmp.size
            for i in range(n_colors):
                tmpi = f[f[tmp[i, 0]][0, 0]][:, 0]
                wl = int(tmpi.astype(dtype=np.uint8).tobytes().decode())
                wls[key].append(wl)
        wls['combo'] = []
        for i in range(n_colors):
            wls['combo'].append(
                (wls['excitationWavelength'][i], wls['emission'][i]),
                )

        print('Read image data')
        img = f['HAC_Image']['imageStruct']['data'][:, :, :]

    cell_bboxes = {
        # xlim, ylim, angle to rotate (counterclockwise, in 360 degrees), path (after bbox and rotation)
        fn+' #1': [
            (1664.9294140231955, 2024.5584277188732), (3436.625717828795, 3229.996766977331), 0,
            [],
            ],
        fn+' #2': [
            (968.9291281436347, 1332.594225971346), (2812.239789502408, 2566.1778327632596), 30,
            [],
            ],
        fn+' #3': [
            (1140.8564210130933, 1245.1407344471465), (1830.815551197974, 1448.929332988765), 92,
            [],
            ],
    }

    def plot_image(img, bbox=None, rotate=True):
        from scipy.ndimage import rotate as rotfun
        if bbox is not None:
            (x0, x1), (y0, y1), theta, path = bbox
            img = img[:, int(y1): int(y0) + 1, int(x0): int(x1) + 1]
            if rotate:
                img = rotfun(img, theta, axes=(1, 2), reshape=True)

        fig, axs = plt.subplots(3, 5, figsize=(16, 10), sharex=True, sharey=True)
        axs = axs.ravel()
        for i in range(n_colors):
            ax = axs[i]
            imgi = img[i]
            ax.imshow(imgi, interpolation='nearest')
            ax.set_title(str(wls['combo'][i]))
            ax.set_xticks([])
            ax.set_yticks([])
        fig.tight_layout()
        return fig

    for cname, bbox in cell_bboxes.items():
        fig = plot_image(img, bbox=bbox)
        fig.suptitle(cname)
        input()
        plt.close(fig)

