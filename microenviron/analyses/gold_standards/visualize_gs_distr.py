#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : visualize_gs_distr.py
#   Author       : Yufeng Liu
#   Date         : 2023-01-24
#   Description  : 
#
#================================================================
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from anatomy.anatomy_core import parse_ana_tree

sys.path.append('../../src')
from config import __FEAT_NAMES__, plot_sd_matrix

def calc_regional_mean_features(lmfile, outfile, min_counts=10):
    rnkey = 'region_name_r316'
    print(rnkey)
    df = pd.read_csv(lmfile, index_col=0)
    print(df.shape)
    df = df[~df[rnkey].isna()]
    print(df.shape)

    # annotation dict
    ana_dict = parse_ana_tree(keyname='name')
    struct_dict = {
        688: 'CTX',
        623: 'CNU',
        512: 'CB',
        343: 'BS'
    }
    
    regions, counts = np.unique(df[rnkey], return_counts=True)
    output = []
    for region, count in zip(regions, counts):
        print(region, count)
        if count < min_counts:
            continue
        region_index = df.index[df[rnkey] == region]
        # find the brain structure
        for pid in ana_dict[region]['structure_id_path']:
            if pid in struct_dict:
                struct_name = struct_dict[pid]
                break
        else:
            print(f'!! {new_rname}')
            struct_name = np.NaN

        feat = df.loc[region_index, __FEAT_NAMES__]
        fmean = feat.mean().to_numpy().tolist()
        output.append([region, struct_name, len(region_index), *fmean])
    columns = [rnkey, 'brain_structures', 'NumRecons', *__FEAT_NAMES__]
    rmf = pd.DataFrame(output, columns=columns)
    # normalization required
    temp = rmf.loc[:, __FEAT_NAMES__]
    rmf.loc[:, __FEAT_NAMES__] = (temp - temp.mean()) / (temp.std() + 1e-10)

    rmf.to_csv(outfile, float_format='%g', index=False)

def visualize_distr(feat_file, region='CP'):
    #fnames = __FEAT_NAMES__
    fnames = ['Stems', 'Bifurcations', 'Branches', 'Tips', 'Length', 'Volume', 
              'MaxEuclideanDistance', 'MaxPathDistance', 'MaxBranchOrder', 
              'AverageContraction', 'AverageFragmentation', 'AverageBifurcationAngleLocal',
              'AverageBifurcationAngleRemote']
    df = pd.read_csv(feat_file, index_col=0)
    rcoords = pd.DataFrame(list(zip(range(len(fnames)), df.loc[region,fnames])), columns=['Region', 'feature'])
    
    fs = list(zip(np.repeat(df.index.to_list(), df.shape[0]),
                fnames * df.shape[0],
                df.loc[:,fnames].to_numpy().reshape(-1)))
    df_f = pd.DataFrame(data=fs, columns=['Region', 'feature_name', 'feature'])
    sns.boxplot(data=df_f, x='feature_name', y='feature', showfliers=False)

    sns.scatterplot(data=rcoords, x='Region', y='feature', marker='o', color='red', label=region)
    plt.xticks(ticks=range(len(fnames)), labels=fnames, rotation=90, fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('')
    plt.ylabel('mean-features', fontsize=18)
    #plt.ylim(-3, 2.5)
    plt.legend(loc='upper left')
    plt.tight_layout()

    plt.savefig(f'{region}_feature_distr_gs_local.png', dpi=300)
    plt.close('all')

def visualize_sd_matrix(region_feature_file):
    df = pd.read_csv(region_feature_file, index_col=0)
    brain_structures = df.pop('brain_structures')
    nr = df.pop('NumRecons')
    corr = df.transpose().corr()
    plot_sd_matrix(brain_structures, corr, 'gold_standard_sd.png', 'gold_standard-SD-matrix')
    
if __name__ == '__main__':
    lmfile = './lm_gs_local.csv'
    regional_feature_file = 'regional_features_gs.csv'

    #calc_regional_mean_features(lmfile, regional_feature_file, min_counts=10)
    #visualize_distr(regional_feature_file)
    visualize_sd_matrix(regional_feature_file)


