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
import itertools
import numpy as np
import pandas as pd

sys.path.append('../../common_lib')
sys.path.append('../../microenviron/generation')
from common_utils import load_type_from_excel, stype2struct, plot_sd_matrix, struct_dict, CorticalLayers, PstypesToShow
from config import __FEAT_NAMES__

def load_data(feat_file, celltype_file, normalize=True, min_num_neurons=10):
    df_ct = pd.read_csv(celltype_file)
    df_f = pd.read_csv(feat_file)
    df = df_f.merge(df_ct, how='inner', on='Cell name')

    # remove number of neurons less than `min_num_neurons`
    regions, counts = np.unique(df.Manually_corrected_soma_region, return_counts=True)
    regions = regions[counts >= min_num_neurons]

    df_sel = df[df.Manually_corrected_soma_region.isin(regions)]
    return df_sel

def load_data_with_cortical_layer(feat_file, celltype_file, normalize=True, min_num_neurons=10):
    df_ct = pd.read_csv(celltype_file)
    df_f = pd.read_csv(feat_file)
    df = df_f.merge(df_ct, how='inner', on='Cell name')

    # remove number of neurons less than `min_num_neurons`
    regions, counts = np.unique(df.Manually_corrected_soma_region, return_counts=True)
    regions = regions[counts >= min_num_neurons]
    df_sel = df[df.Manually_corrected_soma_region.isin(regions)]

    cstype = []
    for stype, cl in zip(df_sel.Manually_corrected_soma_region, df_sel.Cortical_layer):
        if cl is np.NaN:
            cl = ''
        else:
            cl = f'-{cl}'
        cstype.append(f'{stype}{cl}')
    df_sel['cstype'] = cstype

    regions = CorticalLayers
    df_sel = df_sel[df_sel.cstype.isin(regions)]

    return df_sel

def load_data_with_ptype(feat_file, celltype_file, target_ptypes, normalize=True, min_num_neurons=10):
    df_ct = pd.read_csv(celltype_file)
    df_f = pd.read_csv(feat_file)
    df = df_f.merge(df_ct, how='inner', on='Cell name')

    # remove number of neurons less than `min_num_neurons`
    regions, counts = np.unique(df.Manually_corrected_soma_region, return_counts=True)
    regions = regions[counts >= min_num_neurons]
    df_sel = df[df.Manually_corrected_soma_region.isin(regions)]

    ptypes = []
    for stype, pt in zip(df_sel.Manually_corrected_soma_region, df_sel.Subclass_or_type):
        if pt is np.NaN:
            pt = ''
        else:
            pt = pt.split('_')[-1]
            pt = f'-{pt}'
        ptypes.append(f'{stype}{pt}')
    df_sel['ptype'] = ptypes

    df_sel = df_sel[df_sel.ptype.isin(target_ptypes)]

    return df_sel

def plot_all(feat_file, celltype_file, normalize=True, min_num_neurons=10):
    df_sel = load_data(feat_file, celltype_file, normalize, min_num_neurons)
    structs = [stype2struct[stype] for stype in df_sel.Manually_corrected_soma_region]
    df_sel = df_sel[__FEAT_NAMES__]
    if normalize:
        df_sel = (df_sel - df_sel.mean()) / (df_sel.std() + 1e-10)
    corr = df_sel.transpose().corr()
    plot_sd_matrix(structs, corr, 'sdmatrix_fullMorph', '')

def plot_struct(feat_file, celltype_file, figname, normalize=True, min_num_neurons=10, regions=None, vmin=-0.4, vmax=0.8, annot=False):
    df_sel = load_data(feat_file, celltype_file, normalize, min_num_neurons)
    if regions is not None:
        df_sel = df_sel[df_sel.Manually_corrected_soma_region.isin(regions)]
    structs = [stype for stype in df_sel.Manually_corrected_soma_region]
    df_sel = df_sel[__FEAT_NAMES__]
    if normalize:
        df_sel = (df_sel - df_sel.mean()) / (df_sel.std() + 1e-10)
    corr = df_sel.transpose().corr()
    plot_sd_matrix(structs, regions, corr, figname, '', vmin=vmin, vmax=vmax, annot=annot)

def plot_struct_with_cortical_layer(feat_file, celltype_file, figname, normalize=True, min_num_neurons=10, regions=None, vmin=-0.4, vmax=0.8, annot=False):
    df_sel = load_data_with_cortical_layer(feat_file, celltype_file, normalize, min_num_neurons)
    structs = [stype for stype in df_sel.cstype]
    df_sel = df_sel[__FEAT_NAMES__]
    if normalize:
        df_sel = (df_sel - df_sel.mean()) / (df_sel.std() + 1e-10)
    corr = df_sel.transpose().corr()
    regions = CorticalLayers
    plot_sd_matrix(structs, regions, corr, figname, '', vmin=vmin, vmax=vmax, annot=annot)

def plot_struct_with_ptype(feat_file, celltype_file, figname, normalize=True, min_num_neurons=10, regions=None, vmin=-0.4, vmax=0.8, annot=False):
    df_sel = load_data_with_ptype(feat_file, celltype_file, regions, normalize, min_num_neurons)
    structs = [ptype for ptype in df_sel.ptype]
    df_sel = df_sel[__FEAT_NAMES__]
    if normalize:
        df_sel = (df_sel - df_sel.mean()) / (df_sel.std() + 1e-10)
    corr = df_sel.transpose().corr()
    plot_sd_matrix(structs, regions, corr, figname, '', vmin=vmin, vmax=vmax, annot=annot)

if __name__ == '__main__':
    feat_file = '../data/feature.txt'
    celltype_file = '../../common_lib/41586_2021_3941_MOESM4_ESM.csv'
    outdir = '../data/sd_matrix/levels'

    if 1:
        # stypes
        structures = [key for key in struct_dict.keys()] + ['all']
        regions_list = [value for value in struct_dict.values()]
        regions_list = regions_list + [list(itertools.chain(*regions_list))]
        for structure, regions in zip(structures, regions_list):
            figname = os.path.join(outdir, f'sdmatrix_fullMorpho_stype_{structure.lower()}')
            plot_struct(feat_file, celltype_file, figname, regions=regions, vmin=-0.4, vmax=0.8, annot=False)

    if 1:
        # ptypes
        structures = [key for key in PstypesToShow.keys()] + ['all']
        regions_list = [value for value in PstypesToShow.values()]
        regions_list = regions_list + [list(itertools.chain(*regions_list))]
        for structure, regions in zip(structures, regions_list):
            figname = os.path.join(outdir, f'sdmatrix_fullMorpho_ptype_{structure.lower()}')
            plot_struct_with_ptype(feat_file, celltype_file, figname, regions=regions, vmin=-0.4, vmax=0.8, annot=False)

    if 1:
        figname = os.path.join(outdir, f'sdmatrix_fullMorpho_cstype_all')
        plot_struct_with_cortical_layer(feat_file, celltype_file, figname, regions=CorticalLayers, vmin=-0.4, vmax=0.8, annot=False)


