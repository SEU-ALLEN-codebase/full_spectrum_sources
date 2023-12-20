#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : ds_matrix_condense.py
#   Author       : Yufeng Liu
#   Date         : 2023-02-08
#   Description  : 
#
#================================================================
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

import sys
sys.path.append('../../common_lib')
from common_utils import struct_dict, CorticalLayers, PstypesToShow, plot_sd_matrix

sns.set_theme(style="whitegrid", rc={'legend.labelspacing': 1.0, 'font.weight': 'bold'})

def summary_class_level_ds(df, cols=(14,9,3)):
    if type(df) is str:
        df = pd.read_csv(df, index_col=0)

    icols = [0]
    icol = 0
    for ncol in cols:
        icol += ncol
        icols.append(icol)

    intras = []
    inters = []
    diffs = []
    for i in range(len(icols)-1):
        ix = icols[i]
        iy = icols[i+1]
        #print(len(df.columns[ix:iy]), df.columns[ix:iy])
        intra = df.iloc[ix:iy, ix:iy].sum().sum()
        inter = df.iloc[ix:iy].sum().sum() - intra
        intra /= (iy - ix)**2
        inter /= ((df.shape[0] - iy + ix) * (iy - ix))
        intras.append(intra)
        inters.append(inter)
        diffs.append(intra - inter)
    intras, inters, diffs = np.round(np.array(intras), 4), np.round(np.array(inters), 4), np.round(np.array(diffs), 4)
    return intras, inters, diffs

def plot_dsmatrix_illustration():
    import seaborn as sns
    import matplotlib.pyplot as plt

    data = pd.DataFrame([[0.6, -0.3],[-0.3, 0.8]], 
                        columns=['Cell type1', 'Cell type2'], 
                        index=['Cell type1', 'Cell type2'])
    labels = [['intra-type', 'inter-type'],
            ['inter-type', 'intra-type']]
    fig, axes = plt.subplots(figsize=(8,6))
    g = sns.heatmap(data=data, cmap='vlag', annot=labels, fmt='', cbar=True, ax=axes, 
                    annot_kws={"size": 25, 'color':'black'}, 
                    vmin=-.5, vmax=1)
    plt.setp(axes.xaxis.get_majorticklabels(), fontsize=25)
    plt.setp(axes.yaxis.get_majorticklabels(), fontsize=25)
    g.figure.axes[-1].set_yticks([-0.5,0,0.5,1.0], [-0.5,0,0.5,1.0], fontsize=20)

    plt.savefig('example_dsmatrix.png', dpi=300)
    plt.close('all')

def load_cross_scale_features(feat_files, celltype_file):
    df_ct = pd.read_csv(celltype_file, index_col=0)
    for lname in ['microenviron', 'fullMorpho', 'arbor', 'motif', 'bouton']:
        level = feat_files[lname]
        print(f'--> Loading data: {lname}')
        path = level['path']
        neuron = level['neuron']
        fn = level['feat_names']
        if type(path) == str:
            if lname == 'fullMorpho' or lname == 'bouton':
                df_tmp = pd.read_csv(path)
            else:
                df_tmp = pd.read_csv(path, index_col=0)
            
            if neuron == 'index':
                df_tmp.reset_index(inplace=True)
            if fn is not None:
                df_tmp = df_tmp[[neuron, *fn]]
            else:
                if level['drop_key'] is not None:
                    df_tmp = df_tmp.drop(level['drop_key'], axis=1, inplace=True)
        else:
            df_ax = pd.read_csv(path[0])
            df_ba = pd.read_csv(path[1])
            df_tmp = df_ax.merge(df_ba, how='inner', on='prefix')

            include_apical = len(feat_files) == 3
            if include_apical:
                print('Include apical...')
                df_ap = pd.read_csv(path[2])
                fnames = ['max_density', 'num_nodes', 'total_path_length', 'volume',
                          'num_branches', 'dist_to_soma', 'dist_to_soma2', 'num_hubs',
                          'variance_ratio']
                cidx = df_tmp.shape[1]
                df_tmp.loc[:, fnames] = np.zeros((df_tmp.shape[0], len(fnames)))
                df_ap_prefixs = set(df_ap['prefix'].to_numpy().tolist())
                for prefix in df['prefix']:
                    if prefix in df_ap_prefixs:
                        idx = np.nonzero((df_tmp['prefix'] == prefix).to_numpy())[0][0]
                        idx2 = np.nonzero((df_ap['prefix'] == prefix).to_numpy())[0][0]
                        df_tmp.iloc[idx, cidx:] = df_ap.iloc[idx2, -len(fnames):]

            df_tmp.set_index('region_x', inplace=True)
            df_tmp.drop(['Unnamed: 0_x', 'region_y', "Unnamed: 0_y"], axis=1, inplace=True)
         
        if lname == 'microenviron':
            df = df_tmp
        else:
            df = df.merge(df_tmp, how='inner', left_on='Cell name', right_on=neuron)
    
    df = df.merge(df_ct, how='inner', on='Cell name')
    # generate stype, sptype and sctype
    sctypes, sptypes = [], []
    for stype, cl, pt in zip(df.Manually_corrected_soma_region, df.Cortical_layer, df.Subclass_or_type):
        if cl is np.NaN:
            cl = ''
        else:
            cl = f'-{cl}'
        sctypes.append(f'{stype}{cl}')

        if pt is np.NaN:
            pt = ''
        else:
            pt = pt.split('_')[-1]
            pt = f'-{pt}'
        sptypes.append(f'{stype}{pt}')
    df['sctype'] = sctypes
    df['sptype'] = sptypes
    df.rename(columns={'Manually_corrected_soma_region': 'stype'}, inplace=True)

    # remove redundent columns
    df.drop(['formatted name', 'Soma_x', 'Soma_y', 'Soma_z', 'Registered_soma_region', 
            'Transgenic_line', 'Brain_id', 'prefix', 'index', 
            'Cortical_layer', 'Subclass_or_type'], axis=1, inplace=True)
    df.set_index('Cell name', inplace=True)

    return df

def calc_dsmatrix_stype(df, figname):
    regions = struct_dict['CTX'] + struct_dict['TH'] + struct_dict['STR']
    df = df[df.stype.isin(regions)]
    stypes = df.stype
    df.drop(['stype', 'sctype', 'sptype'], axis=1, inplace=True)
    # normalize
    df = (df - df.mean()) / (df.std() + 1e-10)
    corr = df.transpose().corr()
    print(df.shape, corr.shape, stypes, regions)
    plot_sd_matrix(stypes, regions, corr, figname, '', annot=False, vmin=-0.5, vmax=0.5)

def calc_dsmatrix_sptype(df, figname):
    regions = PstypesToShow['CTX'] + PstypesToShow['TH'] + PstypesToShow['STR']
    df = df[df.sptype.isin(regions)]
    sptypes = df.sptype
    df.drop(['stype', 'sctype', 'sptype'], axis=1, inplace=True)
    # normalize
    df = (df - df.mean()) / (df.std() + 1e-10)
    corr = df.transpose().corr()
    print(df.shape, corr.shape, sptypes, regions)
    plot_sd_matrix(sptypes, regions, corr, figname, '', annot=False, vmin=-0.5, vmax=0.5)


def calc_dsmatrix_sctype(df, figname):
    regions = CorticalLayers
    df = df[df.sctype.isin(regions)]
    sctypes = df.sctype
    df.drop(['stype', 'sctype', 'sptype'], axis=1, inplace=True)
    # normalize
    df = (df - df.mean()) / (df.std() + 1e-10)
    corr = df.transpose().corr()
    print(df.shape, corr.shape, sctypes, regions)
    plot_sd_matrix(sctypes, regions, corr, figname, '', annot=False, vmin=-0.5, vmax=0.5)

def plot_dsmatrix_relplot(matfile, type_str='stype'):
    df = pd.read_csv(matfile, index_col=0)
    
    height = 12
    ds_value = 'DS value'
    if type_str == 'sptype':
        keys = PstypesToShow['CTX']
        df = df.loc[keys, keys]
    elif type_str == 'stype':
        keys = struct_dict['STR'] + struct_dict['TH'] + struct_dict['CTX']
        df = df.loc[keys, keys]

    df_stack = df.stack().reset_index(name=ds_value)
    norm = (-0.2, 0.6)
    g = sns.relplot(
        data = df_stack,
        x="level_0", y="level_1", hue=ds_value, size=ds_value,
        palette="icefire", edgecolor="1.", #hue_norm=norm,    # icefire
        height=height, sizes=(30, 500), #size_norm=norm,
    )
    g.set(xlabel="", ylabel="", aspect="equal")
    g.despine(left=True, bottom=True)
    g.ax.margins(x=0.02, y=0.02)
    for label in g.ax.get_xticklabels():
        label.set_rotation(90)
    
    plt.xticks(fontsize=21)
    plt.yticks(fontsize=21)

    for artist in g.legend.legendHandles:
        artist.set_color("1.")
    plt.setp(g._legend.get_texts(), fontsize=16)
    plt.setp(g._legend.get_title(), fontsize=16)
    g._legend.set_bbox_to_anchor([1.01,0.3])

    plt.subplots_adjust(left=0.14, bottom=0.15)

    #plt.tight_layout()
    plt.savefig(f'dsmatrix_multi-scale_{type_str}.png', dpi=300)
    plt.close('all')
    

if __name__ == '__main__':
    import sys; sys.path.append('../../microenviron/generation')
    from config import __FEAT_NAMES__

    celltype_file = '../../common_lib/41586_2021_3941_MOESM4_ESM.csv'
    feat_files = {
        'microenviron': {
            'path': '../../microenviron/me_map_new20230510/data/gold_standard_me.csv',
            'feat_names': __FEAT_NAMES__,
            'neuron': 'Cell name'
        },
        'fullMorpho': {
            'path': '../../full_morpho/data/feature.txt',
            'feat_names': __FEAT_NAMES__,
            'neuron': 'Cell name'
        },
        'arbor': {
            'path': ('../../arbors/src/min_num_neurons10_l2/features_r2_somaTypes_axonal.csv',
                  '../../arbors/src/min_num_neurons10_l2/features_r2_somaTypes_dendrite.csv'),
            'feat_names': None,
            'drop_key': ['region'],
            'neuron': 'prefix',
        },
        'motif': {
            'path': '../../sd_matrix/motif_features.csv',
            'feat_names': None,
            'neuron': 'index',
            'drop_key': None,
        },
        'bouton': {
            'path': '../../bouton/bouton_features.csv',
            'feat_names': ("Bouton Number", "TEB Ratio",
                          "Bouton Density", "Geodesic Distance",
                          "Bouton Interval", "Project Regions"),
            'neuron': 'Cell name'
        }
    }

 
    if 0:   
        figname = 'sdmatrix_heatmap_stype_all'
        df = load_cross_scale_features(feat_files, celltype_file)
        calc_dsmatrix_stype(df, figname)
        figname = 'sdmatrix_heatmap_sctype_all'
        calc_dsmatrix_sctype(df, figname)
        figname = 'sdmatrix_heatmap_sptype_all'
        calc_dsmatrix_sptype(df, figname)
    
    
    if 1:
        for type_str in ['stype', 'sctype', 'sptype']:
            if type_str != 'stype':
                continue
            matfile = f'../multi-scale/corr_regionLevel_sdmatrix_heatmap_{type_str}_all.csv'
            plot_dsmatrix_relplot(matfile, type_str=type_str)
    

