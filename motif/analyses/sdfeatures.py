#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : sdmatrix.py
#   Author       : Yufeng Liu
#   Date         : 2023-01-30
#   Description  : 
#
#================================================================
import os
import sys
import numpy as np
import pandas as pd
import pickle
from sklearn.decomposition import PCA

sys.path.append('../../common_lib')
from common_utils import load_type_from_excel, stype2struct, plot_sd_matrix, struct_dict, CorticalLayers, PstypesToShow

FEAT_NAMES = ['PathLength', 'EucLength', 'OrientX', 'OrientY', 'OrientZ', 'Volume']

class MainTractFeatures(object):
    def __init__(self, path_coords):
        self.coords = path_coords

    def plength(self):
        pts1 = self.coords[:-1]
        pts2 = self.coords[1:]
        length = np.linalg.norm(pts2 - pts1, axis=1).sum()
        return length

    def elength(self):
        pt1 = self.coords[0]
        pt2 = self.coords[-1]
        length = np.linalg.norm(pt2 - pt1)
        return length

    def orientation(self):
        pt1 = self.coords[0] # termini
        pt2 = self.coords[-1]
        v = pt1 - pt2
        v = v / np.linalg.norm(v)
        return v

    def volume(self):
        pca = PCA(3)
        new = pca.fit_transform(self.coords)
        diff = new.max(axis=0) - new.min(axis=0)
        return diff.prod()

    def run(self, include_coord=False):
        pl = self.plength()
        el = self.elength()
        orient = self.orientation()
        vo = self.volume()
        if include_coord:
            return pl, el, *orient, vo, *self.coords[-1]
        else:
            return pl, el, *orient, vo


def calc_feature(main_tract_file):
    with open(main_tract_file, 'rb') as fp:
        cdict = pickle.load(fp)

    features = []
    for region, clist in cdict.items():
        print(region)
        for neuron in clist:
            prefix = neuron[0]
            coords = neuron[1]
            mtf = MainTractFeatures(coords)
            fs = mtf.run()
            features.append([prefix, *fs])

    df = pd.DataFrame(features, columns=['Cell name', *FEAT_NAMES])
    tmp = df[FEAT_NAMES]
    df.loc[:, FEAT_NAMES] = (tmp - tmp.mean()) / (tmp.std() + 1e-10)
    return df

def reassign_classes(main_tract_file, celltype_file):
    df_ct = pd.read_csv(celltype_file, index_col=0)
    df_f = calc_feature(main_tract_file)
    df = df_f.merge(df_ct, how='inner', on='Cell name')

    # assign cortical_layer and ptype
    cstypes = []
    ptypes = []
    for stype, cl, pt in zip(df.Manually_corrected_soma_region, df.Cortical_layer, df.Subclass_or_type):
        if cl is np.NaN:
            cl = ''
        else:
            cl = f'-{cl}'
        cstypes.append(f'{stype}{cl}')

        if pt is np.NaN:
            pt = ''
        else:
            pt = pt.split('_')[-1]
            pt = f'-{pt}'
        ptypes.append(f'{stype}{pt}')
    df['cstype'] = cstypes
    df['ptype'] = ptypes

    return df

def calc_regional_features(df, out_dir):
    # s-type regional mean features
    stypes = struct_dict['CTX'] + struct_dict['TH'] + struct_dict['STR']
    df_st = df[df.Manually_corrected_soma_region.isin(stypes)][['Manually_corrected_soma_region', *FEAT_NAMES]]
    mean_st = df_st.groupby('Manually_corrected_soma_region').mean().reindex(stypes)
    std_st = df_st.groupby('Manually_corrected_soma_region').std().reindex(stypes)
    mean_st.to_csv(os.path.join(out_dir, 'stype_mean_features_motif.csv'))
    std_st.to_csv(os.path.join(out_dir, 'stype_std_features_motif.csv'))

    # p-type regional features
    ptypes = PstypesToShow['CTX'] + PstypesToShow['TH'] + PstypesToShow['STR']
    df_pt = df[df.ptype.isin(ptypes)][['ptype', *FEAT_NAMES]]
    mean_pt = df_pt.groupby('ptype').mean().reindex(ptypes)
    std_pt = df_pt.groupby('ptype').std().reindex(ptypes)
    mean_pt.to_csv(os.path.join(out_dir, 'ptype_mean_features_motif.csv'))
    std_pt.to_csv(os.path.join(out_dir, 'ptype_std_features_motif.csv'))

    # cortical layer for CTX
    cstypes = CorticalLayers
    df_cs = df[df.cstype.isin(cstypes)][['cstype', *FEAT_NAMES]]
    mean_cs = df_cs.groupby('cstype').mean().reindex(cstypes)
    std_cs = df_cs.groupby('cstype').std().reindex(cstypes)
    mean_cs.to_csv(os.path.join(out_dir, 'cstype_mean_features_motif.csv'))
    std_cs.to_csv(os.path.join(out_dir, 'cstype_std_features_motif.csv'))


if __name__ == '__main__':
    main_tract_file = 'main_tract.pkl'
    celltype_file = '../../common_lib/41586_2021_3941_MOESM4_ESM.csv'
    out_dir = '../../sd_features/data'

    df = reassign_classes(main_tract_file, celltype_file)
    calc_regional_features(df, out_dir)
    

