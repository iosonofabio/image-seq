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


# FIXME: not great as a global, but it's a pilot script so oh, well
wls = {'combo': [
    (325, 414),
    (343, 414),
    (366, 414),
    (343, 451),
    (366, 451),
    (373, 451),
    (343, 575),
    (393, 575),
    (406, 575),
    (441, 575),
    (400, 594),
    (406, 594),
    (431, 594),
    (480, 594),
    (339, 575),
]}


def read_image(fn):
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
    return {
        'data': img,
        'wavelengths': wls['combo'],
        'image': os.path.basename(fn).split('.')[0],
    }

if __name__ == '__main__':


    fdn_images = '../../data/hyperspectral/20210601_MiSeq/'
    fns = [
        'MCF7_C3.mat',
        'MCF7_C4.mat',
        'MCF7_N4.mat',
    ]
    fn = fdn_images+fns[2]

    cell_bboxes = {
        # xlim, ylim, angle to rotate (counterclockwise, in 360 degrees), path (after bbox and rotation),
        # sequencing well
        'MCF7_C3 #1': [
            (1664.9294140231955, 2024.5584277188732), (3436.625717828795, 3229.996766977331), 0,
            [(3, 90), (59, 157), (120, 196), (216, 201), (296, 174), (333, 151),
             (356, 114), (279, 31), (141, 14), (130, 15)],
            'L5',
            ],
        'MCF7_C3 #2': [
            (968.9291281436347, 1332.594225971346), (2812.239789502408, 2566.1778327632596), 30,
            [(56, 220), (149, 291), (186, 293), (238, 288), (309, 265), (379, 228), (413, 192),
             (392, 154), (317, 112), (257, 106), (174, 121), (103, 170), (75, 193)],
            None,
            ],
        'MCF7_C3 #3': [
            (1140.8564210130933, 1245.1407344471465), (1830.815551197974, 1448.929332988765), 92,
            [(36, 53), (114, 76), (204, 85), (267, 78), (306, 71), (373, 47),
             (326, 25), (219, 14), (132, 25), (67, 37)],
            None,
            ],
        'MCF7_C4 #1': [
            (344.893470825001, 752.7833823903588), (3190.5434957298794, 2795.8824461612357), 0,
            [(35, 195), (75, 279), (151, 333), (193, 353), (275, 323),
             (349, 247), (375, 167), (367, 89), (241, 43), (149, 81),
             (107, 109), (49, 159)],
            'E5',
            ],
        'MCF7_C4 #2': [
            (2836.2139524080576, 3173.7336675186716), (526.5080939704727, 186.96730272147158), 0,
            [(50, 173), (65, 246), (155, 307), (174, 305), (246, 276), (291, 231),
             (310, 128), (196, 38), (141, 42), (81, 77), (54, 147)],
            'F5',
            ],
        'MCF7_C4 #3': [
            (1714.0685468503084, 2115.344027417307), (621.0561497506573, 281.5153585016579), -10,
            [(54, 229), (101, 327), (211, 311), (322, 275), (406, 236), (391, 206),
             (345, 152), (250, 72), (196, 51), (116, 105), (80, 164), (54, 203)],
            'C5',
            ],
        # The following two are actually sister cells during mitosis
        'MCF7_C4 #4': [
            (3126.160615231529, 3450.8189367342256), (1773.2836872857968, 1406.733969460171), 30,
            [(141, 265), (259, 336), (330, 352), (385, 346), (364, 150), (236, 143),
             (149, 184), (134, 226)],
            'D5',
            ],
        'MCF7_C4 #5': [
            (3329.654019240843, 3621.7867126825813), (1843.283230297852, 1546.537915380718), 30,
            [(75, 209), (107, 314), (251, 297), (307, 287), (337, 221), (317, 177),
             (193, 99), (114, 104)],
            'K5',
            ],
        'MCF7_N4 #1': [
            (124.10896933879481, 476.5395876553565), (2317.670230100273, 1973.4836028554435), -30,
            [(45, 239), (84, 274), (158, 350), (230, 396), (254, 401), (333, 366),
             (410, 279), (412, 252), (347, 149), (267, 104), (196, 125), (145, 162),
             (92, 210)],
            'I5',
            ],
        'MCF7_N4 #2': [
            (2672.422431695017, 3192.8243680981273), (2927.1475175609808, 2507.21922234461), 30,
            [(93, 332), (235, 431), (313, 479), (381 ,470), (468, 423), (545, 354),
             (623, 298), (524, 255), (403, 185), (300, 134), (231, 190), (166, 276),
             (149, 298)],
            'G5',
            ],
        'MCF7_N4 #3': [
            (3399.087017715191, 3798.701979528641), (1366.6388899318329, 1019.9849471539005), -5,
            [(66, 187), (105, 259), (179, 314), (288, 300), (339, 261), (380, 222),
             (178, 164), (343, 113), (255, 76), (123, 87), (70, 136), (68, 145)],
            'H5',
            ],
        'MCF7_N4 #4': [
            (3592.8415025306913, 4008.446854928823), (1615.8381342167213, 1278.1562882601022), 40,
            [(108, 282), (162, 345), (293, 382), (356, 359), (444, 268), (387, 188),
             (270, 129), (207, 140), (159, 191), (131, 142)],
            'J5',
            ],
    }

    def segment_image(img, bbox):
        from scipy.ndimage import rotate as rotfun
        from matplotlib.path import Path
        (x0, x1), (y0, y1), theta, poly_verts, _ = bbox
        img = img[:, int(y1): int(y0) + 1, int(x0): int(x1) + 1]
        img = rotfun(img, theta, axes=(1, 2), reshape=True)

        if poly_verts:
            _, ny, nx = img.shape
            x, y = np.meshgrid(np.arange(nx), np.arange(ny))
            x, y = x.flatten(), y.flatten()
            points = np.vstack((x, y)).T
            path = Path(poly_verts)
            grid = path.contains_points(points)
            grid = grid.reshape((ny, nx))
            img[:, ~grid] = 0
        return img

    def extract_features(img):
        img_max = img.max(axis=0)
        img_bin = img_max > 0
        area = (img_max > 0).sum()
        # Horizontal length
        length = img_bin.any(axis=0).nonzero()[0]
        length = length[-1] - length[0]
        # Vertical width (waist line)
        width = img_bin.any(axis=1).nonzero()[0]
        width = width[-1] - width[0]

        # Eccentricity
        ecc = length / width

        # Spectrum
        spectrum = img.sum(axis=1).sum(axis=1)

        feas = {
            'area': area,
            'length': length,
            'width': width,
            'eccentricity': ecc,
            'spectrum': spectrum,
        }

        return feas

    def plot_image(img, bbox=None, rotate=True):
        from scipy.ndimage import rotate as rotfun
        from matplotlib.path import Path
        if bbox is not None:
            (x0, x1), (y0, y1), theta, poly_verts, _ = bbox
            img = img[:, int(y1): int(y0) + 1, int(x0): int(x1) + 1]
            if rotate:
                img = rotfun(img, theta, axes=(1, 2), reshape=True)

            if poly_verts:
                _, ny, nx = img.shape
                x, y = np.meshgrid(np.arange(nx), np.arange(ny))
                x, y = x.flatten(), y.flatten()
                points = np.vstack((x, y)).T
                path = Path(poly_verts)
                grid = path.contains_points(points)
                grid = grid.reshape((ny, nx))
                img[:, ~grid] = 0

        fig, axs = plt.subplots(3, 5, figsize=(7.87, 5.69), sharex=True, sharey=True)
        axs = axs.ravel()
        n_colors = img.shape[0]
        for i in range(n_colors):
            ax = axs[i]
            imgi = img[i]
            ax.imshow(imgi, interpolation='nearest')
            ax.set_title(str(wls['combo'][i]))
            ax.set_xticks([])
            ax.set_yticks([])
        fig.tight_layout()
        return fig

    plt.ion()

    if False:
        print('Plot whole image')
        plot_image(img)
        sys.exit()

        print('Plot single cell ROIs')
        for cname, bbox in cell_bboxes.items():
            pattern = cname.split()[0]
            if pattern not in fn:
                continue
            fig = plot_image(img, bbox=bbox)
            fig.suptitle(cname)
            input()

    if False:
        print('Extract image features')
        features_all = []
        for fn in fns:
            print(fn)
            fn = fdn_images+fn
            dic = read_image(fn)
            img = dic['data']
            wavelengths = dic['wavelengths']
            image = dic['image']
            for cname, bbox in cell_bboxes.items():
                imagei = cname.split()[0]
                if image != imagei:
                    continue
                if bbox[-1] is None:
                    continue
                print(cname)
                img_seg = segment_image(img, bbox)
                features = extract_features(img_seg)
                features['image'] = image
                features['well'] = bbox[-1]
                features['wavelengths'] = wavelengths
                features_all.append(features)
        features_all = pd.DataFrame(features_all)

        import pickle
        fn_out = fdn_images+'features.pkl'
        with open(fn_out, 'wb') as f:
            pickle.dump(features_all, f)

    if True:
        print('Load feature cache')
        import pickle
        fn_out = fdn_images+'features.pkl'
        with open(fn_out, 'rb') as f:
            features_all = pickle.load(f)

    def plot_bubble_spectrum(spectrum, wls):
        wls_ex = list(np.unique([x[0] for x in wls]))
        wls_em = list(np.unique([x[1] for x in wls]))
        nex = len(wls_ex)
        nem = len(wls_em)
        smax = spectrum.max()
        snorm = spectrum / smax

        fig, ax = plt.subplots(figsize=(0.5 + 0.4 * nem, 0.1 + 0.4 * nex))
        for (wex, wem), val in zip(wls, snorm):
            x = wls_em.index(wem)
            y = wls_ex.index(wex)
            r = (3 + 97 * val) / 200.
            alpha = 0.2 + 0.8 * val
            h = plt.Circle(
                (x, y), r, ec='none', fc='black', alpha=alpha,
                )
            ax.add_artist(h)
        ax.set_xlim(-0.6, nem - 0.4)
        ax.set_ylim(-0.6, nex - 0.4)
        ax.set_xticks(np.arange(nem))
        ax.set_yticks(np.arange(nex))
        ax.set_xticklabels([str(x) for x in wls_em])
        ax.set_yticklabels([str(x) for x in wls_ex])
        ax.set_xlabel('$\lambda_{emiss}$ [nm]')
        ax.set_ylabel('$\lambda_{exit}$ [nm]')

        fig.tight_layout()
        return {
            'fig': fig,
            'ax': ax,
        }



    if True:
        print('Plot images for a single cell as an example: G5')
        cname = 'MCF7_N4 #2'
        bbox = cell_bboxes[cname]
        fn = fdn_images+cname.split()[0]+'.mat'
        dic = read_image(fn)
        img = dic['data']

        fig = plot_image(img, bbox=bbox)
        fig.suptitle(cname)

        features = features_all.iloc[-3]

        pdic = plot_bubble_spectrum(features['spectrum'], features['wavelengths'])
        fig = pdic['fig']
        fig.suptitle('G5')
        fig.tight_layout()
