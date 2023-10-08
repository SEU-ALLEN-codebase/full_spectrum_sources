#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : predict_types.py
#   Author       : Yufeng Liu
#   Date         : 2023-02-05
#   Description  : 
#
#================================================================
import sys
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

sys.path.append('../../common_lib')
from common_utils import struct_dict, normalize_features, assign_subtypes

__CMP_FEAT_NAMES__ = ['Stems', 'Bifurcations', 'Branches', 'Tips', 'OverallWidth',
                      'OverallHeight', 'OverallDepth', 'Length', 'Volume', 'MaxEuclideanDistance',
                      'MaxPathDistance', 'MaxBranchOrder', 'AverageContraction',
                      'AverageFragmentation', 'AverageParent-daughterRatio', 
                      'AverageBifurcationAngleLocal', 'AverageBifurcationAngleRemote', 
                      'HausdorffDimension']

def load_features(gs_file, me_file, celltype_file, nodes_range=(500,1500)):
    df_gs = pd.read_csv(gs_file, index_col=0)
    df_me = pd.read_csv(me_file, index_col=0)
    df_ct = pd.read_csv(celltype_file, index_col=0)
    
    # filter by node number and region name
    #nmin, nmax = nodes_range
    #df_me = df_me[(df_me['Nodes'] >= nmin) & (df_me['Nodes'] <= nmax)]
    # remove reconstructions with nan
    #df_me = df_me[~(df_me.isna().sum(axis=1).astype(np.bool))]
    # by regions
    regions = []
    for rs in struct_dict.values():
        regions.extend(rs)
    df_me = df_me[df_me.region_name_r316.isin(regions)]
    df_gs = df_gs[df_gs.region_name_r316.isin(regions)]

    # merge celltype and gs
    df_gs = df_gs.merge(df_ct, how='inner', left_on=df_gs.index, right_on='Cell name')
    print(df_gs.shape, df_me.shape)

    # normalize
    normalize_features(df_me, __CMP_FEAT_NAMES__, inplace=True)
    normalize_features(df_gs, __CMP_FEAT_NAMES__, inplace=True)

    # get the subtypes
    assign_subtypes(df_gs, inplace=True)

    return df_gs, df_me

def load_features_microenviron(gs_file, me_file, cell_type, nodes_range=(500,1500)):
    df_gs = pd.read_csv(gs_file, index_col=0)
    df_me = pd.read_csv(me_file, index_col=0)
    df_ct = pd.read_csv(celltype_file, index_col=0)
    
    # filter by node number and region name
    nmin, nmax = nodes_range
    df_me = df_me[(df_me['Nodes'] >= nmin) & (df_me['Nodes'] <= nmax)]
    # remove reconstructions with nan
    df_me = df_me[~(df_me.isna().sum(axis=1).astype(np.bool))]
    # get the features of me
    me_keys = [f'{fn}_mean' for fn in __CMP_FEAT_NAMES__]
    df_me = df_me[['region_name_r316', *me_keys]]
    mapper = dict(zip(me_keys, __CMP_FEAT_NAMES__))
    df_me.rename(columns=mapper, inplace=True)

    # merge celltype and gs
    df_gs = df_gs.merge(df_ct, how='inner', left_on=df_gs.index, right_on='Cell name')

    # normalize
    normalize_features(df_me, __CMP_FEAT_NAMES__, inplace=True)
    normalize_features(df_gs, __CMP_FEAT_NAMES__, inplace=True)

    return df_gs, df_me

def predict_ptype(df_gs, df_me):
    print('==> Predict p-type')
    gs_ptype = df_gs.loc[:, ['ptype', *__CMP_FEAT_NAMES__]].groupby('ptype').mean()
    regions = ['-'.join(x.split('-')[:max(1, len(x.split('-'))-1)]) for x in gs_ptype.index]
    gs_ptype['region'] = regions
    gs_ptype.set_index(['region', gs_ptype.index], inplace=True)
    # intialize
    df_me['ptype'] = ['' for i in range(df_me.shape[0])]
    # prediction by matching the centroid of group
    for region in np.unique(regions):
        ridx = np.nonzero((df_me.region_name_r316 == region).to_numpy())[0]
        if gs_ptype.loc[region].shape[0] > 1:
            cur_feat = df_me.iloc[ridx].loc[:,__CMP_FEAT_NAMES__]
            pdist = distance_matrix(cur_feat, gs_ptype.loc[region])
            cur_ptypes = gs_ptype.loc[region].index[pdist.argmin(axis=1)]
            df_me.iloc[ridx, df_me.shape[-1]-1] = cur_ptypes
        else:
            df_me.iloc[ridx, df_me.shape[-1]-1] = gs_ptype.loc[region].index.values[0]

def predict_cstype(df_gs, df_me):
    print('==> Predict cs-type')
    gs_cstype = df_gs.loc[:, ['cstype', *__CMP_FEAT_NAMES__]].groupby('cstype').mean()
    regions = ['-'.join(x.split('-')[:max(1, len(x.split('-'))-1)]) for x in gs_cstype.index]
    gs_cstype['region'] = regions
    gs_cstype.set_index(['region', gs_cstype.index], inplace=True)
    # intialize
    df_me['cstype'] = ['' for i in range(df_me.shape[0])]
    # prediction by matching the centroid of group
    for region in np.unique(regions):
        ridx = np.nonzero((df_me.region_name_r316 == region).to_numpy())[0]
        cur_feat = df_me.iloc[ridx].loc[:,__CMP_FEAT_NAMES__]
        pdist = distance_matrix(cur_feat, gs_cstype.loc[region])
        cur_cstypes = gs_cstype.loc[region].index[pdist.argmin(axis=1)]
        df_me.iloc[ridx, df_me.shape[-1]-1] = cur_cstypes

def predict_microenviron_gs(df_gs, df_me):
    print(f'--> Predict ME features for GS')
    gs_features = df_gs[__CMP_FEAT_NAMES__]
    me_features = df_me[__CMP_FEAT_NAMES__]
    pdist = distance_matrix(gs_features, me_features)
    min_idx = np.argmin(pdist, axis=1)
    df_gs.loc[:, __CMP_FEAT_NAMES__] = me_features.iloc[min_idx].to_numpy()

    df_gs.to_csv('../me_map_new20230510/data/gold_standard_me.csv')
    
    

if __name__ == '__main__':
    # Surface-related features are inconsistent in Gold standards
    gs_file = '../gs_local/src/lm_gs_dendrite.csv'
    me_file = '../me_map_new20230510/data/lm_features_d22_15441.csv'
    celltype_file = '../../common_lib/41586_2021_3941_MOESM4_ESM.csv'
    nodes_range = (500, 1500)
    
    if 0:
        df_gs, df_me = load_features(gs_file, me_file, celltype_file, nodes_range)
        predict_ptype(df_gs, df_me)
        predict_cstype(df_gs, df_me)
        df_me.to_csv('../me_map_new20230510/data/lm_features_d22_15441_with_ptype_cstype.csv')
    

    if 1:
        me_file = '../me_map_new20230510/data/micro_env_features_nodes300-1500_statis.csv'
        df_gs, df_me = load_features_microenviron(gs_file, me_file, nodes_range)
        predict_microenviron_gs(df_gs, df_me)

