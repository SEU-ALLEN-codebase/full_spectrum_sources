#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : sdmatrix.py
#   Author       : Yufeng Liu
#   Date         : 2023-02-02
#   Description  : 
#
#================================================================
import os
import string
import itertools
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

import sys
sys.path.append('../../common_lib')
from common_utils import plot_sd_matrix, struct_dict, PstypesToShow, CorticalLayers
sys.path.append('../generation')
from config import __FEAT_ALL__


# Say, "the default sans-serif font is COMIC SANS"
matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"
#matplotlib.rcParams['font.weight'] = 'bold'


if 0:
    # DS matrix for stypes (separate and all)
    mefeature_file = '../me_map_new20230510/data/micro_env_features_nodes300-1500_statis_with_ptype_cstype.csv'
    outdir = '../../sd_matrix/levels'
    normalize = True
    
    structures = [key for key in struct_dict.keys()] + ['all']
    regions_list = [value for value in struct_dict.values()]
    regions_list = regions_list + [list(itertools.chain(*regions_list))]
    df_raw = pd.read_csv(mefeature_file, index_col=0)
    for structure, regions in zip(structures, regions_list):
        print(f'==> Processing for structure: {structure}, including {len(regions)} regions')
        df = df_raw[df_raw['region_name_r316'].isin(regions)]
        outname = os.path.join(outdir, f'sdmatrix_microenviron_stype_{structure.lower()}')

        if normalize:
            tmp = df.loc[:, __FEAT_ALL__]
            df_feat = (tmp - tmp.mean()) / (tmp.std() + 1e-10)
        else:
            df_feat = df[__FEAT_ALL__]
        print(df_feat.mean().mean())
        corr = df_feat.transpose().corr()
        structs = df['region_name_r316']
        plot_sd_matrix(structs, regions, corr, outname, '', vmin=-0.4, vmax=0.8, annot=False)


if 0:
    # DS matrix for ptypes
    mefeature_file = f'../me_map_new20230510/data/micro_env_features_nodes300-1500_statis_with_ptype_cstype.csv'
    outdir = '../../sd_matrix/levels'
    normalize = True
    
    structures = [key for key in PstypesToShow.keys()] + ['all']
    regions_list = [value for value in PstypesToShow.values()]
    regions_list = regions_list + [list(itertools.chain(*regions_list))]
    df_raw = pd.read_csv(mefeature_file, index_col=0)
    for structure, regions in zip(structures, regions_list):
        print(f'==> Processing for structure: {structure}, including {len(regions)} regions')
        df = df_raw[df_raw['ptype'].isin(regions)]
        outname = os.path.join(outdir, f'sdmatrix_microenviron_ptype_{structure.lower()}')

        if normalize:
            tmp = df.loc[:, __FEAT_ALL__]
            df_feat = (tmp - tmp.mean()) / (tmp.std() + 1e-10)
        else:
            df_feat = df[__FEAT_ALL__]
        print(df_feat.mean().mean())
        corr = df_feat.transpose().corr()
        structs = df['ptype']
        plot_sd_matrix(structs, regions, corr, outname, '', vmin=-0.4, vmax=0.8, annot=False)

if 1:
    # DS matrix for cstypes
    mefeature_file = f'../me_map_new20230510/data/micro_env_features_nodes300-1500_statis_with_ptype_cstype.csv'
    outdir = '../../sd_matrix/levels'
    normalize = True
    
    df_raw = pd.read_csv(mefeature_file, index_col=0)
    structure = 'clayer'
    regions = CorticalLayers
    print(f'==> Processing for structure: {structure}, including {len(regions)} regions')
    df = df_raw[df_raw['cstype'].isin(regions)]
    outname = os.path.join(outdir, f'sdmatrix_microenviron_cstype_all')

    if normalize:
        tmp = df.loc[:, __FEAT_ALL__]
        df_feat = (tmp - tmp.mean()) / (tmp.std() + 1e-10)
    else:
        df_feat = df[__FEAT_ALL__]
    print(df_feat.mean().mean())
    corr = df_feat.transpose().corr()
    structs = df['cstype']
    plot_sd_matrix(structs, regions, corr, outname, '', vmin=-0.4, vmax=0.8, annot=False)

