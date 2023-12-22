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
from common_utils import load_type_from_excel, stype2struct, plot_sd_matrix, struct_dict, CorticalLayers, PstypesToShow, normalize_features

sys.path.append('../generation')
from config import __FEAT_NAMES__, __FEAT_ALL__

__FEAT_MEAN__ = [f'{fn}_mean' for fn in __FEAT_NAMES__]
__FEAT_ALL__ = [f'{fn}_mean' for fn in __FEAT_NAMES__] + \
               [f'{fn}_median' for fn in __FEAT_NAMES__] + \
               [f'{fn}_std' for fn in __FEAT_NAMES__]

def load_feature(feat_file):
    df = pd.read_csv(feat_file, index_col=1)
    # make sure they are normalized
    #tmp = df[__FEAT_ALL__]
    #df.loc[:, __FEAT_ALL__] = (tmp - tmp.mean()) / (tmp.std() + 1e-10)
    return df

def calc_regional_features(df, out_dir):
    # s-type regional mean features
    rn = 'region_name_r316'
    stypes = struct_dict['CTX'] + struct_dict['TH'] + struct_dict['STR']
    df_st = df[df[rn].isin(stypes)][[rn, *__FEAT_MEAN__]]
    mean_st = df_st.groupby(rn).mean().reindex(stypes)
    std_st = df_st.groupby(rn).std().reindex(stypes)
    mean_st.to_csv(os.path.join(out_dir, 'stype_mean_features_microenviron.csv'))
    std_st.to_csv(os.path.join(out_dir, 'stype_std_features_microenviron.csv'))

    # p-type features
    ptypes = PstypesToShow['CTX'] + PstypesToShow['TH'] + PstypesToShow['STR']
    df_pt = df[df.ptype.isin(ptypes)][['ptype', *__FEAT_MEAN__]]
    mean_pt = df_pt.groupby('ptype').mean().reindex(ptypes)
    std_pt = df_pt.groupby('ptype').std().reindex(ptypes)
    mean_pt.to_csv(os.path.join(out_dir, 'ptype_mean_features_microenviron.csv'))
    std_pt.to_csv(os.path.join(out_dir, 'ptype_std_features_microenviron.csv'))

    # cortical layer for CTX
    cstypes = CorticalLayers
    df_cs = df[df.cstype.isin(cstypes)][['cstype', *__FEAT_MEAN__]]
    mean_cs = df_cs.groupby('cstype').mean().reindex(cstypes)
    std_cs = df_cs.groupby('cstype').std().reindex(cstypes)
    mean_cs.to_csv(os.path.join(out_dir, 'cstype_mean_features_microenviron.csv'))
    std_cs.to_csv(os.path.join(out_dir, 'cstype_std_features_microenviron.csv'))


if __name__ == '__main__':
    feat_file = '../me_map_new20230510/data/micro_env_features_nodes300-1500_statis_with_ptype_cstype.csv'
    out_dir = '../../sd_features/data'

    df = load_feature(feat_file)
    calc_regional_features(df, out_dir)
    

