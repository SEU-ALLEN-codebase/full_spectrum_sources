#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : multi-levels_analyzer.py
#   Author       : Yufeng Liu
#   Date         : 2023-02-07
#   Description  : 
#
#================================================================
import os
import glob
import numpy as np
from scipy.spatial import distance_matrix
from scipy.stats import skew, kurtosis
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage

from ds_matrix_condense import summary_class_level_ds

import sys
sys.path.append('../../common_lib')
from common_utils import struct_dict, CorticalLayers, PstypesToShow

__LEVELS__ = ('microenviron', 'fullMorpho', 'arbor', 'motif', 'bouton')
sns.set_theme(style="whitegrid", rc={'legend.labelspacing': 1.0, 'font.weight': 'normal'})

def calc_ds_similarity_among_levels(data_dir='./levels', type_str='stype', struct='all'):
    vecs = []
    for level in __LEVELS__:
        print(level)
        csvfile = os.path.join(data_dir, f'corr_regionLevel_sdmatrix_{level}_{type_str}_{struct}.csv')
        matrix = pd.read_csv(csvfile, index_col=0)
        vec = matrix.to_numpy().reshape(-1)
        #vec = vec / np.linalg.norm(vec)
        vecs.append(vec)
    vecs = np.array(vecs)
    #pdists = np.matmul(vecs, vecs.transpose())
    #pdists = pd.DataFrame(pdists, index=__LEVELS__, columns=__LEVELS__)
    pdists = pd.DataFrame(vecs.transpose(), columns=__LEVELS__)
    pdists.rename(columns={'bouton': 'varicosity'}, inplace=True)
    pdists = pdists.corr()

    rpdists = 1 - pdists

    ds_value = 'distance'
    pdists = pdists.stack().reset_index(name=ds_value)
    tmp = pdists[ds_value]
    tmp.clip(0, 1)
    pdists.loc[:, ds_value] = 1 - tmp
    
    g = sns.relplot(
        data = pdists, 
        x="level_0", y="level_1", hue=ds_value, size=ds_value,
        palette="icefire", edgecolor="0.7", hue_norm=(0, 0.8),
        height=4, sizes=(10, 500), size_norm=(0, 0.8)
    )
    g.set(xlabel="", ylabel="", aspect="equal")
    g.despine(left=True, bottom=True)
    g.ax.margins(0.1)
    for label in g.ax.get_xticklabels():
        label.set_rotation(90)
        label.set_fontsize(15)
    for label in g.ax.get_yticklabels():
        label.set_fontsize(15)

    for artist in g.legend.legendHandles:
        artist.set_color("1.")

    plt.subplots_adjust(left=0.28, right=0.8, bottom=0.35)

    plt.savefig(f'ds_similarity_{type_str}_{struct}.png', dpi=300)
    plt.close('all')

def load_multiLevel_features(cache_file, data_dir, type_str, cnames, vname):
    if os.path.exists(cache_file):
        print(f'--> Loading cache file {cache_file}, make sure it is what you want!')
        df = pd.read_csv(cache_file, index_col=0)
        return df

    data = []
    for level in __LEVELS__:
        print(f'--> Processing for level: {level}')
        feat_file = os.path.join(data_dir, f'corr_neuronLevel_sdmatrix_{level}_{type_str}_all.csv')
        df = pd.read_csv(feat_file, index_col=0)

        c_cnt = 0
        for cname in cnames:
            print(f'    --> class: {cname}')
            indices = np.nonzero((df.type == cname).to_numpy())[0]
            dfs = df.iloc[indices, indices].to_numpy()
            triu_indices = np.triu_indices_from(dfs, k=1)
            vds = dfs[triu_indices[0], triu_indices[1]]
            levels = np.tile(level, vds.shape[0])
            names = np.tile(cname, vds.shape[0])
            data.append(np.vstack((vds, levels, names)).transpose())
            
            c_cnt += 1
            #if c_cnt == 5:
            #    break
        
    data = np.vstack(data)
    df = pd.DataFrame(data, columns=[vname, 'level', 'type'])
    df = df.astype({vname: float})
    print(f'Total shape of data: {df.shape}')
    # save file
    df.to_csv(cache_file)

    return df

def calc_intraregional_stereotypy(data_dir='./levels', type_str='stype'):
    if type_str == 'stype':
        cols = (14,9,3)
        cnames = ('CTX', 'TH', 'STR')
    elif type_str == 'ptype':
        cols = (10, 9, 5, 3, 1, 1, 1, 1, 1)
        cnames = ('CTX-ET', 'CTX-IT', 'TH-core', 'TH-matrix', 'TH-RT', 'CP-GPe', 'CP-SNr', 'STR-ACB', 'STR-OT')
    elif type_str == 'cstype':
        cols = (6, 5, 11, 2)
        cnames = ('L2/3', 'L4', 'L5', 'L6')


    data = []
    data_reg = []
    for level in __LEVELS__:
        feat_file = os.path.join(data_dir, f'corr_regionLevel_sdmatrix_{level}_{type_str}_all.csv')
        df = pd.read_csv(feat_file, index_col=0)
        names = df.columns

        intras, inters, diffs = summary_class_level_ds(df, cols)
        print(f'{level:<15} {intras} {inters} {diffs}')
        data.append(intras)

        intras_reg, inters_reg, diffs_reg = summary_class_level_ds(df, np.ones(df.shape[1], dtype=np.int32))
        #print(f'{level:<15} {intras_reg} {inters_reg} {diffs_reg}')
        data_reg.append(intras_reg)

    

    # region-level
    if type_str == 'ptype':
        height = 15
    else:
        height = 12

    df_reg = pd.DataFrame(data_reg, columns=names, index=__LEVELS__)
    ds_value = 'intra-type\nDS value'
    df_reg_stack = df_reg.transpose().stack().reset_index(name=ds_value)
    orders = [0, 0.2, 0.4, 0.6, 0.8]
    g = sns.relplot(
        data = df_reg_stack, 
        x="level_0", y="level_1", hue=ds_value, size=ds_value,
        palette="icefire", edgecolor="1.", hue_norm=(0.2, 1.0),    # icefire
        height=height, sizes=(30, 500), size_norm=(0.2, 1.0),
    )
    g.set(xlabel="", ylabel="", aspect="equal")
    g.despine(left=True, bottom=True)
    g.ax.margins(x=0.02, y=0.15)
    for label in g.ax.get_xticklabels():
        label.set_rotation(90)
    for artist in g.legend.legendHandles:
        artist.set_edgecolor("1.")
    
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    # adjust the legend
    plt.setp(g._legend.get_texts(), fontsize=16)
    plt.setp(g._legend.get_title(), fontsize=16)

    if type_str == 'ptype':
        print(f'Adjust figure layout for {type_str}')
        plt.subplots_adjust(left=0.12, right=0.92)
    else:
        plt.subplots_adjust(left=0.15, right=0.9)

    plt.savefig(f'intra-region_stereotypy_{type_str}.png', dpi=300)
    plt.close('all')

    
def calc_intraregional_stereotypy_ridge_plot(data_dir='./levels', type_str='stype'):
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    vname = 'DS value'
    fontsize = 30
    pal = {
        'microenviron': 'orangered',
        'fullMorpho': 'darkviolet',
        'arbor': 'blue',
        'motif': 'peru',
        'bouton': 'forestgreen'
    }
    type_names = {
        'stype': 's-type',
        'ptype': 'sp-type',
        'cstype': 'sl-type'
    }

    if type_str == 'stype':
        cnames = struct_dict['CTX'] + struct_dict['TH'] + struct_dict['STR']
        hspace = -0.9
        height = 1
    elif type_str == 'ptype':
        cnames = PstypesToShow['CTX']# + PstypesToShow['TH'] + PstypesToShow['STR']
        hspace = -0.8
        height = 1
    elif type_str == 'cstype':
        cnames = CorticalLayers
        hspace = -0.8
        height = 1
    
    cache_file = f'data/intraregional_dsmatrix_{type_str}.csv'
    df = load_multiLevel_features(cache_file, data_dir, type_str, cnames, vname)
    import ipdb; ipdb.set_trace()
    
    # Initialize the facegrid
    g = sns.FacetGrid(df, row='type', col='level', hue='level', aspect=2.5, 
                      height=height, palette=pal)
    
    # Draw the densities in a few steps
    print('Draw kdeplot...')
    bw_adjust = 1.
    g.map(sns.kdeplot, vname,
          bw_adjust=bw_adjust, clip=(-1,1),
          fill=True, alpha=0.3, linewidth=1.5)
    g.map(sns.kdeplot, vname, clip=(-1,1), lw=2, bw_adjust=bw_adjust, alpha=0.7)

    # passing color=None to refline() uses the hue mapping
    #g.refline(y=0, linewidth=2, linestyle="-", color='k', clip_on=False)
    #g.map(plt.axhline, y=0, linewidth=2, linestyle="-", color='k', clip_on=False)

    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        if label == 'microenviron':
            ax = plt.gca()
            ax.text(-.5, (1+hspace)*0.4, x.iloc[0],# fontweight="bold", 
                    color='k',
                    ha="left", va="center", transform=ax.transAxes,
                    fontsize=fontsize, fontweight='normal')

    g.map(label, 'type')

    # Set the subplots to overlap
    g.figure.subplots_adjust(hspace=hspace, wspace=-0.15)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[], ylabel='')
    #g.set(xticks=np.arange(-.5,1+0.5,0.5))
    #g.set_xticklabels(np.arange(-.5,1.5,0.5), fontsize=fontsize-5, fontweight='normal')
    g.set(xticks=[], xlabel='')

    #g.add_legend()
    #sns.move_legend(g, 'lower center', ncols=5)
    #plt.setp(g._legend.get_title(), fontsize=0, fontweight='bold')
    #plt.setp(g._legend.get_texts(), fontsize=fontsize+5, fontweight='bold')
    for i, level in enumerate(__LEVELS__):
        g.axes[-1,i].set_xlabel(level, fontsize=fontsize, fontweight='normal', rotation=30)
    
    g.despine(bottom=True, left=True)
    plt.subplots_adjust(left=0.15, bottom=0.1)
    plt.ylabel('Pairwise correlation', fontsize=fontsize, fontweight='normal')

    plt.savefig(f'intra-regional_ridge_plot_{type_str}.png', dpi=300)
    plt.close('all')


def calc_distribution_statistics(data_dir='./levels', type_str='stype'):
    sns.set_theme(style='white', rc={'legend.labelspacing': 1.0, 'font.weight': 'normal'})


    vname = 'DS value'
    figname = f'statistics_{type_str}.png'
    pal = {
        'microenviron': 'orangered',
        'fullMorpho': 'darkviolet',
        'arbor': 'blue',
        'motif': 'peru',
        'bouton': 'forestgreen'
    }
    type_names = {
        'stype': 's-type',
        'ptype': 'sp-type',
        'cstype': 'sl-type'
    }

    if type_str == 'stype':
        cnames = struct_dict['CTX'] + struct_dict['TH'] + struct_dict['STR']
        hspace = -0.7
        height = 1.0
    elif type_str == 'ptype':
        cnames = PstypesToShow['CTX']# + PstypesToShow['TH'] + PstypesToShow['STR']
        hspace = -0.6
        height = 1.0
    elif type_str == 'cstype':
        cnames = CorticalLayers
        hspace = -0.5
        height = 1.0

    cache_file = f'data/intraregional_dsmatrix_{type_str}.csv'
    df = load_multiLevel_features(cache_file, data_dir, type_str, cnames, vname)

    jname = 'level-type'
    df[jname] = df['level'] + '-' + df['type']

    loaded_jnames = np.unique(df[jname])

    # resort the data
    skeys = []
    for level in __LEVELS__:
        for cname in cnames:
            skey = f'{level}-{cname}'
            if skey not in loaded_jnames:
                continue
            skeys.append(skey)
    levels = [level.split('-')[0] for level in skeys]
    types = ['-'.join(level.split('-')[1:]) for level in skeys]

    def plot_relplot(dfi, levels, types, title):
        name = dfi.columns[0]
        dfi['level'] = levels
        dfi['type'] = types
        nl = len(__LEVELS__)
        dfi['x'] = np.tile(np.arange(dfi.shape[0]//nl), nl)
        #sns.relplot(dfi, x='level', y='x', size=name, hue='level', 
        #            aspect=1, kind='scatter', palette=pal); 
        plt.savefig(f'{title}.png'); 
        plt.close('all')

    def plot_heatmap(dfi, types):
        name = dfi.columns[0]
        keys = [k if k != 'bouton' else 'varicosity' for k in __LEVELS__]
        dfi = pd.DataFrame(dfi[name].to_numpy().reshape(5,-1),
                           index=keys, 
                           columns=types[:len(types)//len(__LEVELS__)])
        g = sns.heatmap(dfi, cmap='coolwarm', xticklabels=1, square=True,
                        cbar_kws={'aspect':5, 'shrink':0.15, 'pad':0.03,
                                  })
        g.tick_params(left=True, bottom=False, labelsize=12, direction='out')
        g.tick_params(axis='x', pad=-2)
        plt.ylabel(name, fontsize=14)
        plt.subplots_adjust(left=0.22)
        plt.savefig(f'{name}_intravalues_{type_str}.png', dpi=300)
        plt.close()
        


    df_d = df.drop(['level', 'type'], axis=1)
    means = df_d.groupby(jname).mean().loc[skeys].rename(columns={'DS value': 'mean'})
    stds = df_d.groupby(jname).std().loc[skeys].rename(columns={'DS value': 'std'})
    skews = df_d.groupby(jname).skew().loc[skeys].rename(columns={'DS value': 'skew'})
    kurts = df_d.groupby(jname).apply(pd.DataFrame.kurt).loc[skeys].rename(columns={'DS value': 'kurtosis'})
    
    nc = len(cnames)
    for i in range(5):
        i1, i2 = i*nc, (i+1)*nc
        print(f'MEAN: {means[i1:i2].mean()}, {means[i1:i2].std()}')
        print(f'STD: {stds[i1:i2].mean()}, {stds[i1:i2].std()}')
        print(f'SKEW: {skews[i1:i2].mean()}, {skews[i1:i2].std()}')
        print(f'KURTOSIS: {kurts[i1:i2].mean()}, {kurts[i1:i2].std()}')

    plot_heatmap(means, types)
    plot_heatmap(stds, types)
    plot_heatmap(skews, types)
    plot_heatmap(kurts, types)

    
    
   
if __name__ == '__main__':

    for type_str in ['stype', 'ptype', 'cstype']:
        calc_ds_similarity_among_levels(type_str=type_str)
        #calc_intraregional_stereotypy(type_str=type_str)
        calc_intraregional_stereotypy_ridge_plot(type_str=type_str)
        calc_distribution_statistics(type_str=type_str)
    
