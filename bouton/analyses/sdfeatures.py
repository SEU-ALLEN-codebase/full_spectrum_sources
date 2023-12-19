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

sys.path.append('../common_lib')
from common_utils import load_type_from_excel, stype2struct, plot_sd_matrix, struct_dict, CorticalLayers, PstypesToShow

FEAT_NAMES = ["Bouton Number", "TEB Ratio",
              "Bouton Density", "Geodesic Distance",
              "Bouton Interval", "Project Regions"]

def load_feature(feat_file):
    df = pd.read_csv(feat_file)
    tmp = df[FEAT_NAMES]
    df.loc[:, FEAT_NAMES] = (tmp - tmp.mean()) / (tmp.std() + 1e-10)
    return df

def reassign_classes(feat_file, celltype_file):
    df_ct = pd.read_csv(celltype_file, index_col=0)
    df_f = load_feature(feat_file)
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
    mean_st.to_csv(os.path.join(out_dir, 'stype_mean_features_bouton.csv'))
    std_st.to_csv(os.path.join(out_dir, 'stype_std_features_bouton.csv'))

    # p-type regional features
    ptypes = PstypesToShow['CTX'] + PstypesToShow['TH'] + PstypesToShow['STR']
    df_pt = df[df.ptype.isin(ptypes)][['ptype', *FEAT_NAMES]]
    mean_pt = df_pt.groupby('ptype').mean().reindex(ptypes)
    std_pt = df_pt.groupby('ptype').std().reindex(ptypes)
    mean_pt.to_csv(os.path.join(out_dir, 'ptype_mean_features_bouton.csv'))
    std_pt.to_csv(os.path.join(out_dir, 'ptype_std_features_bouton.csv'))

    # cortical layer for CTX
    cstypes = CorticalLayers
    df_cs = df[df.cstype.isin(cstypes)][['cstype', *FEAT_NAMES]]
    mean_cs = df_cs.groupby('cstype').mean().reindex(cstypes)
    std_cs = df_cs.groupby('cstype').std().reindex(cstypes)
    mean_cs.to_csv(os.path.join(out_dir, 'cstype_mean_features_bouton.csv'))
    std_cs.to_csv(os.path.join(out_dir, 'cstype_std_features_bouton.csv'))


if __name__ == '__main__':
    feat_file = 'bouton_features/bouton_features.csv'
    celltype_file = '../common_lib/41586_2021_3941_MOESM4_ESM.csv'
    out_dir = '../sd_features/data'

    df = reassign_classes(feat_file, celltype_file)
    calc_regional_features(df, out_dir)
    

