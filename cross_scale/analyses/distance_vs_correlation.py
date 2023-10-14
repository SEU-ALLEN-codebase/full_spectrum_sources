#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : distance_vs_correlation.py
#   Author       : Yufeng Liu
#   Date         : 2023-02-15
#   Description  : 
#
#================================================================
import sys
sys.path.append(r'D:\code\support_lib\pylib')
import os
import sys
import random
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

from file_io import load_image
from anatomy.anatomy_config import MASK_CCF25_FILE
from anatomy.anatomy_core import parse_ana_tree

sys.path.append('../../common_lib')
from common_utils import struct_dict

__NUM_BINS__ = 25
__LABEL_FONTS__ = 20

def load_spos_corr(spos_file, df_corr, is_spos=True, zc=228):
    df_spos = pd.read_csv(spos_file, index_col=0)
    if is_spos:
        # all neurons mapped to left-hemisphere
        right_h = df_spos.z_pos >= zc
        nzi = np.nonzero(right_h.to_numpy())[0]
        z_loc = df_spos.columns.get_loc('z_pos')
        df_spos.iloc[nzi, z_loc] = zc * 2 - df_spos.iloc[nzi, z_loc]

        # to 1um resolution
        df_spos.loc[:, ['x_pos', 'y_pos', 'z_pos']] *= 0.04
    else:
        df_spos.loc[:] *= 0.001

    # keep only the target neurons
    df_spos = df_spos.loc[df_corr.index]
    
    pdists = distance_matrix(df_spos, df_spos)
    indices = np.triu_indices_from(pdists, k=1)
    dists = pdists[indices]
    corrs = df_corr.to_numpy()[indices]

    return dists, corrs

def plot_func(dists, corrs, xname, yname, nsample, figname, xmax, bs=5):
    sample_indices = random.sample(range(len(dists)), nsample)
    dists_sample = dists[sample_indices]
    corrs_sample = corrs[sample_indices]

    data = np.array([dists_sample, corrs_sample]).transpose().astype(float)
    df = pd.DataFrame(data, columns=[xname, yname])
    glm = sns.lmplot(data=df, x=xname, y=yname,
               scatter_kws={'s':.5, 'color':'green'}, fit_reg=False
    )

    # show mean average
    bin_means, bin_edges, binnumber = stats.binned_statistic(dists, corrs,
        statistic='mean', bins=__NUM_BINS__)
    bin_stds, bin_edges, binnumber = stats.binned_statistic(dists, corrs,
        statistic='std', bins=__NUM_BINS__)
    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width/2
    fit_color = 'orangered'
    plt.plot(bin_centers, bin_means, color=fit_color)
    plt.fill_between(bin_centers, bin_means-bin_stds, bin_means+bin_stds, color=fit_color, alpha=0.2)
    
    # cc = 0 anchor line
    plt.axhline(y=0, linewidth=2, linestyle='-', color='gray', clip_on=False, alpha=1.)

    plt.xlim(0, xmax)
    plt.xticks(np.arange(0,xmax+1,bs), np.arange(0,xmax+1,bs), fontsize=__LABEL_FONTS__-5)
    yticks = np.round(np.arange(-0.8,1,0.4),1)
    plt.yticks(yticks, fontsize=__LABEL_FONTS__-5)
    plt.xlabel(f'{xname} (mm)', fontsize=__LABEL_FONTS__)
    plt.ylabel(f'{yname}', fontsize=__LABEL_FONTS__)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.gca().xaxis.set_tick_params(width=2, direction='in')
    plt.gca().yaxis.set_tick_params(width=2, direction='in')
    plt.gca().yaxis.set_label_coords(-0.13, 0.42)
    plt.tight_layout()

    plt.savefig(figname, dpi=300)
    plt.close('all')
    

def soma_dist_vs_corr(spos_file, corr_file):
    df_corr = pd.read_csv(corr_file, index_col=0)
    df_corr.drop(['type'], axis=1, inplace=True)
    dists, corrs = load_spos_corr(spos_file, df_corr)

    nsample = 10000
    xname = 'Soma-soma distance'
    yname = 'Correlation'
    figname = 'soma_distance_vs_correlation.png'
    plot_func(dists, corrs, xname, yname, nsample, figname, xmax=10, bs=5)
    
 
def proj_dist_vs_corr(proj_file, corr_file):
    df_corr = pd.read_csv(corr_file, index_col=0)
    df_corr.drop(['type'], axis=1, inplace=True)
    dists, corrs = load_spos_corr(proj_file, df_corr, is_spos=False)

    nsample = 10000
    xname = 'Axon-axon "distance"'
    yname = 'Correlation'
    figname = 'axon_distance_vs_correlation.png'
    plot_func(dists, corrs, xname, yname, nsample, figname, xmax=60, bs=6)

def soma_proj_dist_vs_corr(spos_file, proj_file, corr_file):
    df_corr = pd.read_csv(corr_file, index_col=0)
    df_corr.drop(['type'], axis=1, inplace=True)
    spos_dists, _ = load_spos_corr(spos_file, df_corr, is_spos=True)
    proj_dists, corrs = load_spos_corr(proj_file, df_corr, is_spos=False)

    dists = np.sqrt(spos_dists * proj_dists)    # reweight for balance

    nsample = 10000
    xname = 'Neuron-neuron distance'
    yname = 'Correlation'
    figname = 'neuron_distance_vs_correlation.png'
    plot_func(dists, corrs, xname, yname, nsample, figname, xmax=25)

def calc_region_size_from_pdists(regions, df):
    means = []
    stds = []
    for region in regions:
        print(region)
        cur_proj = df.loc[region]
        pdists = distance_matrix(cur_proj, cur_proj)
        indices = np.triu_indices_from(pdists, k=1)
        dists = pdists[indices]
        means.append(dists.mean())
        stds.append(dists.std())
    data = np.array([regions, means, stds]).transpose()
    df = pd.DataFrame(data, columns=['region', 'mean', 'std'])
    df = df.astype({'mean': float, 'std': float})
    return df

def estimate_region_size(regions=None, scale=25., region_file='region_size.csv'):
    if os.path.exists(region_file):
        df = pd.read_csv(region_file, index_col=0)
    else:
        mask = load_image(MASK_CCF25_FILE)
        # we should zeroing the right-hemisphere
        sz, sy, sx = mask.shape
        mask[sz//2:] = 0

        ana_dict = parse_ana_tree(keyname='name')
        #get the indices of ana_dict
        if regions is None:
            regions = struct_dict['CTX'] + struct_dict['TH'] + struct_dict['STR']

        majors = []
        minors = []
        for region in regions:
            print(region)
            rids = ana_dict[region]['orig_ids']
            for i, rid in enumerate(rids):
                if i == 0:
                    cmask = mask == rid
                else:
                    cmask = cmask | (mask == rid)
            cs = np.vstack(np.nonzero(cmask)).transpose()
            pca = PCA()
            cs_t = pca.fit_transform(cs)
            cs_max = np.max(cs_t, axis=0)
            cs_min = np.min(cs_t, axis=0)
            dim = (cs_max - cs_min) * scale
            minors.append(np.min(dim))
            majors.append(np.max(dim))
            
        data = np.array([regions, majors, minors]).transpose()
        df = pd.DataFrame(data, columns=['region', 'major', 'minor'])
        df = df.astype({'major': float, 'minor': float})
        df.to_csv('region_size.csv')
        
    # ploting
    # scale to millimeter
    fsize_tick = 9
    fsize_label = 17
    pal = {
        'major': 'blue',
        'minor': 'magenta'
    }
    df.loc[:, ['major', 'minor']] = df[['major', 'minor']] / 1000.
    
    plt.rcParams['figure.figsize'] = [4,4]
    dfs = df.set_index('region').stack().reset_index()
    dfs.rename(columns={0:'value'}, inplace=True)
    dfs['INDEX'] = np.repeat(np.arange(df.shape[0]), 2)
    g = sns.scatterplot(dfs, x='INDEX', y='value', hue='level_1', s=45, palette=pal)
    mean_major = df.major.mean()
    mean_minor = df.minor.mean()
    plt.axhline(y=mean_major, linewidth=2, linestyle='--', color=pal['major'], alpha=0.8)
    plt.axhline(y=mean_minor, linewidth=2, linestyle='--', color=pal['minor'], alpha=0.8)
    plt.xticks(np.arange(df.shape[0]), df.region, fontsize=fsize_tick, rotation=90, ha='center')
    plt.yticks(fontsize=fsize_tick)
    plt.xlabel("")
    plt.ylabel('Length (mm)', fontsize=fsize_label)
    plt.legend(title='Axes', loc='upper right', bbox_to_anchor=(0.9,0.97))

    ax = plt.gca()
    ax.text(1.04, 0.5, f"avg_major={mean_major:.2f}; avg_minor={mean_minor:.2f}",
        horizontalalignment='center',
        verticalalignment='center',
        rotation=-90,
        transform=ax.transAxes,
        fontsize=12)
    ax.xaxis.set_tick_params(width=2, direction='in')
    ax.yaxis.set_tick_params(width=2, direction='in')
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    plt.tight_layout()

    plt.savefig('region_size.png', dpi=300)
    plt.close('all')

def estimate_region_projection_size(proj_file, celltype_file, regions=None, region_file='region_proj_size.csv'):
    if os.path.exists(region_file):
        df = pd.read_csv(region_file, index_col=0)
    else:
        df_proj = pd.read_csv(proj_file, index_col=0)
        df_ct = pd.read_csv(celltype_file, index_col=0).set_index('Cell name')
        df_proj['region'] = df_ct['Manually_corrected_soma_region']
        df_proj.set_index('region', inplace=True)

        if regions is None:
            regions = struct_dict['CTX'] + struct_dict['TH'] + struct_dict['STR']

        df = calc_region_size_from_pdists(regions, df_proj)
        df.to_csv(region_file)
        
    # ploting
    # scale to millimeter
    fsize_tick = 9
    fsize_label = 17
    df.loc[df.index, ['mean', 'std']] = df[['mean', 'std']] / 1000.
    
    plt.rcParams['figure.figsize'] = [4,4]
    g = sns.scatterplot(df.reset_index(), x='index', y='mean', color='darkorange', s=45)
    mean_radius = df['mean'].mean()
    plt.axhline(y=mean_radius, linewidth=2, linestyle='--', color='darkorange', alpha=0.8)
    plt.xticks(np.arange(df.shape[0]), df.region, fontsize=fsize_tick, rotation=90, ha='center')
    plt.yticks(fontsize=fsize_tick)
    plt.xlabel("")
    plt.ylabel('Mean intra-dist (mm)', fontsize=fsize_label)

    ax = plt.gca()
    ax.text(1.04, 0.5, f"average={mean_radius:.2f}",
        horizontalalignment='center',
        verticalalignment='center',
        rotation=-90,
        transform=ax.transAxes,
        fontsize=12)
    ax.xaxis.set_tick_params(width=2, direction='in')
    ax.yaxis.set_tick_params(width=2, direction='in')
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

    plt.subplots_adjust(bottom=0.15, left=0.15)

    plt.savefig('region_proj_size.png', dpi=300)
    plt.close('all')

if __name__ == '__main__':
    spos_file = '../../common_lib/misc/soma_pos_1891_v20230110.csv'
    proj_file = '../../common_lib/41586_2021_3941_MOESM4_ESM_proj.csv'
    corr_file = '../multi-scale/corr_neuronLevel_sdmatrix_heatmap_stype_all.csv'
    celltype_file = '../../common_lib/41586_2021_3941_MOESM4_ESM.csv'

    soma_dist_vs_corr(spos_file, corr_file)
    proj_dist_vs_corr(proj_file, corr_file)
    soma_proj_dist_vs_corr(spos_file, proj_file, corr_file)
    
    #estimate_region_size()
    #estimate_region_projection_size(proj_file, celltype_file)
    


