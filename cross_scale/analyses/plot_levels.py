#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : plot.py
#   Author       : Yufeng Liu
#   Date         : 2023-01-30
#   Description  : 
#
#================================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def calc_regional_ds(df):
    if type(df) is str:
        df = pd.read_csv(df, index_col=0)
    r_intra = np.diag(df).mean()
    r_inter = df.to_numpy()[np.triu_indices_from(df, k=1)].mean()
    r_diff = r_intra - r_inter
    return r_intra, r_inter, r_diff

def calc_ctx_diff_ptype(df, icols=(0, 3, 13, 22)):
    if type(df) is str:
        df = pd.read_csv(df, index_col=0)

    c_intra = 0
    n_intra = 0
    for i in range(len(icols)-1):
        for j in range(len(icols)-1):
            if i == j:
                sub = df.iloc[icols[i]:icols[i+1], icols[j]:icols[j+1]]
                n_intra += sub.size
                c_intra += sub.sum().sum()
    c_inter = (df.sum().sum() - c_intra) / df.size
    c_intra /= n_intra
    return c_intra, c_inter, c_intra - c_inter


if 1:
    # regional level
    # brain-structral level
    fdir = './levels'
    levels = ['MicroEnv', 'FullMorpho', 'Arbor', 'ProjectMotif', 'Bouton']
    #Metrics = ['Intra-region', 'Inter-region', 'Diff = Intra - Inter']
    Metrics = ['Intra-region', 'Inter-region']
    ylabels = {
        'ctx': 'Avg DS score of \nCortical types',
        'str': 'Avg DS score of \nStriatum types',
        'th': 'Avg DS score of \nStriatum types'
    }
    
    for struct in ['ctx', 'str', 'th']:
        print(f'==> Processing for brain structure: {struct}')
        ds_values = []
        for lstr in ['microenviron', 'fullMorpho', 'arbor', 'motif', 'bouton']:
            fn = os.path.join(fdir, f'corr_regionLevel_sdmatrix_{lstr}_{struct}.csv')
            ds_intra, ds_inter, ds_diff = calc_regional_ds(fn)
            
            ds_values.append([ds_intra, ds_inter])
        ds_values = np.array(ds_values)
        ds_values = ds_values.transpose().reshape(-1)
        df = pd.DataFrame(ds_values, columns=['DS_score'])
            
        df['Level'] = levels * len(Metrics)
        df['Metric'] = np.repeat(Metrics, len(levels))

        sns.set_style("darkgrid")
        sns.lineplot(df, x='Level', y='DS_score', hue='Metric', markers=True, style='Metric', markersize=10)
        plt.xticks(levels, fontsize=17)
        #plt.yticks([-0.1, 0, 0.1, 0.2, 0.3], fontsize=12)
        plt.yticks(fontsize=15)
        if struct == 'ctx':
            plt.ylim(-0.1, 0.4)
        elif struct == 'th':
            plt.ylim(-0.1, 0.5)
        elif struct == 'str':
            plt.ylim(-0.1, 0.4)
        else:
            raise ValueError

        plt.legend(fontsize=14)#, loc='upper left')
        plt.xlabel('', fontsize=0)
        plt.ylabel(ylabels[struct], fontsize=20)
        plt.tight_layout()
        plt.savefig(f'sdmatrix_comparison_{struct}.png', dpi=300)
        plt.close('all')

if 1:
    
    def plot_evolve(ctype, brain_structure='ctx'):
        if ctype == 'CorticalLayer':
            keyword = 'Layer'
            label_key = 'Layer'
            cols = (0,6,11,22,24)
        elif ctype == 'Ptype':
            keyword = 'Ptype'
            label_key = 'Projection'
            cols = (0,3,13,22)
        else:
            raise ValueError
        
        levels = ['fullMorpho', 'arbor', 'motif', 'bouton']
        Metrics = ['Intra-region', 'Inter-region', 'Intra-class', 'Inter-class']
        values = []
        for level in levels:
            csvfile = f'with{ctype}/corr_regionLevel_sdmatrix_{level}_{brain_structure}_with{keyword}.csv'
            data = pd.read_csv(csvfile, index_col=0)
            r_intra, r_inter = calc_regional_ds(data)[:2]
            c_intra, c_inter = calc_ctx_diff_ptype(data, icols=cols)[:2]
            values.append([r_intra, r_inter, c_intra, c_inter])
        values = np.array(values)
        values = values.transpose().reshape(-1)
        
        df = pd.DataFrame(values, columns=['DS_score'])
        df['Level'] = levels * len(Metrics)
        df['Metric'] = np.repeat(Metrics, len(levels))

        sns.set_style("darkgrid")
        sns.lineplot(df, x='Level', y='DS_score', hue='Metric', markers=True, style='Metric', markersize=10)
        plt.xticks(levels, fontsize=20)
        plt.yticks(fontsize=15)
        plt.ylim(-0.1, 0.4)
        plt.legend(fontsize=14)#, loc='upper left')
        plt.xlabel('', fontsize=0)
        plt.ylabel(f'DS score of Cortical\n{label_key} types', fontsize=20)
        plt.tight_layout()
        plt.savefig(f'comparison_{brain_structure}_{ctype}.png', dpi=300)
        plt.close('all')

    plot_evolve('Ptype', 'ctx')
    plot_evolve('CorticalLayer', 'ctx')


