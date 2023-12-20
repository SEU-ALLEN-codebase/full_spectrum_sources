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

sys.path.append('../../common_lib')
from common_utils import load_type_from_excel, stype2struct, plot_sd_matrix, struct_dict, CorticalLayers, PstypesToShow

FEAT_NAMES_DENDRITE = ['max_density', 'num_nodes', 'total_path_length', 'volume',
                        'num_branches', 'dist_to_soma', 'dist_to_soma2', 'num_hubs',
                        'variance_ratio']
FEAT_NAMES_AXON = [ 'max_density', 'num_nodes', 'total_path_length', 'volume',
                    'num_branches', 'dist_to_soma', 'dist_to_soma2', 
                    'num_hubs', 'variance_ratio']

def load_feature(feat_files):
    df_ax = pd.read_csv(feat_files[0])
    df_ba = pd.read_csv(feat_files[1])
    df = df_ax.merge(df_ba, how='inner', on='prefix')

    include_apical = len(feat_files) == 3
    if include_apical:
        print('Include apical...')
        df_ap = pd.read_csv(feat_files[2])
        cidx = df.shape[1]
        apical_fnames = [f'{fn}_apical' for fn in FEAT_NAMES_DENDRITE]
        nf = len(apical_fnames)
        df.loc[:, apical_fnames] = np.zeros((df.shape[0], nf))
        df_ap_prefixs = set(df_ap['prefix'].to_numpy().tolist())
        for prefix in df['prefix']:
            if prefix in df_ap_prefixs:
                idx = np.nonzero((df['prefix'] == prefix).to_numpy())[0][0]
                idx2 = np.nonzero((df_ap['prefix'] == prefix).to_numpy())[0][0]
                df.iloc[idx, cidx:] = df_ap.iloc[idx2, -nf:].to_numpy()

    df.set_index('region_x', inplace=True)
    df.drop(['Unnamed: 0_x', "Unnamed: 0_y", 'region_y'], axis=1, inplace=True)
    FEAT_ALL = [cn for cn in df.columns if cn != 'prefix']

    tmp = df[FEAT_ALL]
    df.loc[:, FEAT_ALL] = (tmp - tmp.mean()) / (tmp.std() + 1e-10)
    return df, FEAT_ALL

def reassign_classes(feat_files, celltype_file):
    df_ct = pd.read_csv(celltype_file, index_col=0)
    df_f, FEAT_ALL = load_feature(feat_files)
    df = df_f.merge(df_ct, how='inner', left_on='prefix', right_on='Cell name')

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

    return df, FEAT_ALL

def calc_regional_features(df, out_dir, FEAT_ALL):
    # s-type regional mean features
    stypes = struct_dict['CTX'] + struct_dict['TH'] + struct_dict['STR']
    df_st = df[df.Manually_corrected_soma_region.isin(stypes)][['Manually_corrected_soma_region', *FEAT_ALL]]
    mean_st = df_st.groupby('Manually_corrected_soma_region').mean().reindex(stypes)
    std_st = df_st.groupby('Manually_corrected_soma_region').std().reindex(stypes)
    mean_st.to_csv(os.path.join(out_dir, 'stype_mean_features_arbor.csv'))
    std_st.to_csv(os.path.join(out_dir, 'stype_std_features_arbor.csv'))

    # p-type regional features
    ptypes = PstypesToShow['CTX'] + PstypesToShow['TH'] + PstypesToShow['STR']
    df_pt = df[df.ptype.isin(ptypes)][['ptype', *FEAT_ALL]]
    mean_pt = df_pt.groupby('ptype').mean().reindex(ptypes)
    std_pt = df_pt.groupby('ptype').std().reindex(ptypes)
    mean_pt.to_csv(os.path.join(out_dir, 'ptype_mean_features_arbor.csv'))
    std_pt.to_csv(os.path.join(out_dir, 'ptype_std_features_arbor.csv'))

    # cortical layer for CTX
    cstypes = CorticalLayers
    df_cs = df[df.cstype.isin(cstypes)][['cstype', *FEAT_ALL]]
    mean_cs = df_cs.groupby('cstype').mean().reindex(cstypes)
    std_cs = df_cs.groupby('cstype').std().reindex(cstypes)
    mean_cs.to_csv(os.path.join(out_dir, 'cstype_mean_features_arbor.csv'))
    std_cs.to_csv(os.path.join(out_dir, 'cstype_std_features_arbor.csv'))


if __name__ == '__main__':
    feat_files = ['min_num_neurons9_l2/features_r2_somaTypes_axonal.csv',
                  'min_num_neurons9_l2/features_r2_somaTypes_dendrite.csv',
                  ]
    celltype_file = '../../common_lib/41586_2021_3941_MOESM4_ESM.csv'
    out_dir = '../../sd_features/data'

    df, FEAT_ALL = reassign_classes(feat_files, celltype_file)
    calc_regional_features(df, out_dir, FEAT_ALL)
    

