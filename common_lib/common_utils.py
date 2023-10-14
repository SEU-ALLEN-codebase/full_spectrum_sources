#!/usr/bin/env python

#================================================================
#   Copyright (C) 2022 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : common_utils.py
#   Author       : Yufeng Liu
#   Date         : 2022-09-12
#   Description  : 
#
#================================================================
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

stype2struct = {
    'AId': 'CTX',
    'MOs': 'CTX',
    'MOp': 'CTX',
    'CLA': 'CTX',
    'SSs': 'CTX',
    'VISp': 'CTX',
    'RSPv': 'CTX',
    'CP': 'STR',
    'VPL': 'TH',
    'MG': 'TH',
    'VPM': 'TH',
    'LGd': 'TH',
    'LD': 'TH',
    'LP': 'TH',
    'SSp-m': 'CTX',
    'SSp-ll': 'CTX',
    'SSp-bfd': 'CTX',
    'VPLpc': 'TH',  # remove it
    'SMT': 'TH',
    'SSp-n': 'CTX',
    'SSp-un': 'CTX',
    'SSp-ul': 'CTX',
    'VM': 'TH',
    'OT': 'STR',
    'ACB': 'STR',
    'RT': 'TH',
    'VISrl': 'CTX'
}

struct_dict = {
    'CTX': ['AId', 'CLA', 'MOp', 'MOs', 'RSPv', 'SSp-bfd', 'SSp-ll', 'SSp-m', 
            'SSp-n', 'SSp-ul', 'SSp-un', 'SSs', 'VISp', 'VISrl'],
    'TH': ['LGd', 'MG', 'SMT', 'VPL', 'VPM', 'LD', 'LP', 'VM', 'RT'],
    'STR': ['ACB', 'CP', 'OT']
}

CorticalLayers = [
    'MOp-2/3', 'MOs-2/3', 'SSp-bfd-2/3', 'SSs-2/3', 'VISp-2/3', 'VISrl-2/3',
    'SSp-bfd-4', 'SSp-m-4', 'SSp-n-4', 'SSs-4', 'VISp-4', 'MOp-5', 'MOs-5',
    'RSPv-5', 'SSp-bfd-5', 'SSp-ll-5', 'SSp-m-5', 'SSp-n-5', 'SSp-ul-5',
    'SSp-un-5', 'SSs-5', 'VISp-5', 'AId-6', 'SSs-6'
]

PstypesToShow = {
    #'CTX': ['AId-Car3', 'CLA-Car3', 'SSs-Car3', 'MOp-ET', 'MOs-ET', 'RSPv-ET', 
    'CTX': ['MOp-ET', 'MOs-ET', 'RSPv-ET', 
            'SSp-bfd-ET', 'SSp-ll-ET', 'SSp-m-ET', 'SSp-n-ET', 'SSp-ul-ET', 
            'SSp-un-ET', 'SSs-ET', 'MOp-IT', 'MOs-IT', 'SSp-bfd-IT', 'SSp-m-IT', 
            'SSp-n-IT', 'SSp-ul-IT', 'SSs-IT', 'VISp-IT', 'VISrl-IT'],
    'TH': ['LGd-core', 'MG-core', 'SMT-core', 'VPL-core', 'VPM-core', 'LD-matrix',
           'LP-matrix', 'VM-matrix', 'RT'],
    'STR': ['CP-GPe', 'CP-SNr', 'ACB', 'OT']
}

def load_celltypes(celltype_file, column_name='Subclass_or_type', soma_type_merge=True):
    data = pd.read_csv(celltype_file)
    ptypes = data['Subclass_or_type']
    stypes = data['Manually_corrected_soma_region']

    is_soma_type = column_name == 'Manually_corrected_soma_region'
    if is_soma_type:
        ctypes = stypes
    else:
        ctypes = ptypes

    prefixs = data['Cell name']
    ctype_dict = {}
    if is_soma_type and (not soma_type_merge):
        for name, ptype, ctype in zip(prefixs, ptypes, ctypes):
            if ctype is np.nan or ptype is np.nan:
                continue
            key = f'{ptype}-{ctype}'
            if key not in ctype_dict:
                ctype_dict[key] = [name]
            else:
                ctype_dict[key].append(name)
    else:
        for name, ctype in zip(prefixs, ctypes):
            if ctype is np.nan:
                continue
            if ctype not in ctype_dict:
                ctype_dict[ctype] = [name]
            else:
                ctype_dict[ctype].append(name)

    rev_dict = {}
    for key, value in ctype_dict.items():
        for v in value:
            rev_dict[v] = key

    # load correspondence between different cell type level
    p2stypes = {}
    for ptype, stype in zip(ptypes, stypes):
        if ptype is np.nan or stype is np.nan:
            continue
        if ptype not in p2stypes:
            p2stypes[ptype] = set([stype])
        else:
            p2stypes[ptype].add(stype)
        
    return ctype_dict, rev_dict, p2stypes

def load_type_from_excel(celltype_file, column_name='Manually_corrected_soma_region', use_abstract_ptype=False, keep_Car3=False):
    if celltype_file.endswith('xslx'):
        data = pd.read_excel(celltype_file, skiprows=1)
    elif celltype_file.endswith('csv'):
        data = pd.read_csv(celltype_file, index_col=0)
    ctypes = data[column_name]

    prefixs = data['Cell name']
    ctype_dict = {}
    for name, ctype in zip(prefixs, ctypes):
        if ctype is np.nan:
            continue
        if ctype not in ctype_dict:
            ctype_dict[ctype] = [name]
        else:
            ctype_dict[ctype].append(name)

    rev_dict = {}
    for key, value in ctype_dict.items():
        for v in value:
            rev_dict[v] = key

    # load correspondence between different cell type level
    stypes = data['Manually_corrected_soma_region']
    ptypes = data['Subclass_or_type']
    if use_abstract_ptype:
        nptypes = []
        for i in range(len(ptypes)):
            ptname = ptypes[i]
            if ptname is np.nan:
                nptypes.append(np.nan)
                continue
            if not keep_Car3 and ptname == 'Car3':
                nptypes.append(np.nan)
                continue

            nptypes.append(ptname.split('_')[0])
        ptypes = nptypes

    p2stypes = {}
    for ptype, stype in zip(ptypes, stypes):
        if ptype is np.nan or stype is np.nan:
            continue
        if ptype not in p2stypes:
            p2stypes[ptype] = set([stype])
        else:
            p2stypes[ptype].add(stype)
       
    return ctype_dict, rev_dict, p2stypes


def load_pstype_from_excel(celltype_file, sname='Manually_corrected_soma_region', pname='Subclass_or_type', use_abstract_ptype=False, keep_Car3=False):
    if celltype_file.endswith('xslx'):
        data = pd.read_excel(celltype_file, skiprows=1)
    elif celltype_file.endswith('csv'):
        data = pd.read_csv(celltype_file, index_col=0)
    stypes = data[sname]
    ptypes = data[pname]
    if use_abstract_ptype:
        nptypes = []
        for i in range(len(ptypes)):
            ptname = ptypes[i]
            if ptname is np.nan:
                nptypes.append(np.nan)
                continue
            if not keep_Car3 and ptname == 'Car3':
                nptypes.append(np.nan)
                continue
            nptypes.append(ptname.split('_')[0])
        ptypes = nptypes

    prefixs = data['Cell name']
    ctype_dict = {}
    for name, ptype, stype in zip(prefixs, ptypes, stypes):
        if stype is np.nan or ptype is np.nan:
            continue
        key = f'{ptype}-{stype}'
        if key not in ctype_dict:
            ctype_dict[key] = [name]
        else:
            ctype_dict[key].append(name)

    rev_dict = {}
    for key, value in ctype_dict.items():
        for v in value:
            rev_dict[v] = key

    # load correspondence between different cell type level
    p2stypes = {}
    for ptype, stype in zip(ptypes, stypes):
        if ptype is np.nan or stype is np.nan:
            continue
        if ptype not in p2stypes:
            p2stypes[ptype] = set([stype])
        else:
            p2stypes[ptype].add(stype)
       
    return ctype_dict, rev_dict, p2stypes

def get_structures_from_regions(region_ids, ana_dict, struct_dict=None, return_name=True):
    if struct_dict is None:
        struct_dict = {
            688: 'CTX',
            623: 'CNU',
            512: 'CB',
            343: 'BS'
        }
    structures = []
    for idx in region_ids:
        id_path = ana_dict[idx]['structure_id_path']
        for pid in id_path:
            if pid in struct_dict:
                if return_name:
                    structures.append(struct_dict[pid])
                else:
                    structures.append(pid)
                break
        else:
            structures.append(np.NaN)
    return np.array(structures)

def normalize_features(df, feat_names=None, inplace=False):
    if feat_names is None:
        feat_names = df.columns
    tmp = df.loc[:, feat_names]
    if inplace:
        df.loc[:, feat_names] = (tmp - tmp.mean()) / (tmp.std() + 1e-10)
        return df
    else:
        df_new = df.copy()
        df_new.loc[:, feat_names] = (tmp - tmp.mean()) / (tmp.std() + 1e-10)
        return df_new

def assign_subtypes(df, reg_key='Manually_corrected_soma_region', 
                    cortical_layer='Cortical_layer', subclass='Subclass_or_type', 
                    inplace=True):
    cstypes, ptypes = [], []
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
    
    if inplace:
        df['cstype'] = cstypes
        df['ptype'] = ptypes
        return df
    else:
        df_new = df.copy()
        df_new['cstype'] = cstypes
        df_new['ptype'] = ptypes
        return df_new
    

def plot_sd_matrix(brain_structures, structs, corr_raw, figname, title, annot=True, vmin=-0.5, vmax=0.9):
    
    def plot_single_matrix(matrix, structs, figname, vmin, vmax, annot):
        nstructs = len(structs)

        matrix = np.round(matrix, 2)
        df_sd = pd.DataFrame(matrix, columns=structs, index=structs)
        outdir, outname = os.path.split(figname)
        df_sd.to_csv(os.path.join(outdir, f'corr_regionLevel_{outname}.csv'), float_format='%.4f')

        sd_self = np.diag(df_sd).mean()
        sd_inter = df_sd.to_numpy()[np.triu_indices_from(df_sd, k=1)].mean()
        sd_diff = sd_self - sd_inter
        print(f'Mean SD for self, inter and diff are: {sd_self:.3f}, {sd_inter:.3f} and {sd_diff:.3f}')

        if nstructs <=5:
            fs0 = 25
            fs1 = 18
            fs2 = 25
        elif nstructs > 10:
            fs0 = 8
            fs1 = 12
            fs2 = 14
        else:
            fs0 = 10
            fs1 = 16
            fs2 = 17

        fig, ax_sd = plt.subplots(figsize=(6,6))
        sns.heatmap(data=df_sd, ax=ax_sd, cmap='coolwarm', annot=annot,
                    annot_kws={"size": fs0}, cbar=False, vmin=vmin, vmax=vmax)
        #ax_sd.set_title(title, fontsize=30)
        ax_sd.set_xlabel('', fontdict={'fontsize': fs1})
        ax_sd.set_ylabel('', fontdict={'fontsize': fs1})
        plt.setp(ax_sd.xaxis.get_majorticklabels(), fontsize=fs2)
        plt.setp(ax_sd.yaxis.get_majorticklabels(), fontsize=fs2)
        plt.tight_layout()
        plt.savefig(f'{figname}.png', dpi=300)
        plt.close('all')


    brain_structures = np.array(brain_structures)
    #structs = np.unique(brain_structures)
    nstructs = len(structs)
    sd_matrix = np.zeros((nstructs, nstructs))
    std_sd_matrix = sd_matrix.copy()
    print(sd_matrix.shape, corr_raw.shape)
    for i in range(nstructs):
        struct1 = corr_raw.index[brain_structures == structs[i]]
        for j in range(i, nstructs):
            struct2 = corr_raw.index[brain_structures == structs[j]]
            cc = corr_raw.loc[struct1, struct2].values.mean()
            cc_std = corr_raw.loc[struct1, struct2].values.std()
            sd_matrix[i][j] = cc
            sd_matrix[j][i] = cc
            std_sd_matrix[i][j] = cc_std
            std_sd_matrix[j][i] = cc_std

    corr = corr_raw.copy()
    corr['type'] = brain_structures
    outdir, outname = os.path.split(figname)
    corr.to_csv(os.path.join(outdir, f'corr_neuronLevel_{outname}.csv'), float_format='%.4f')
    
    plot_single_matrix(sd_matrix, structs, figname, vmin, vmax, annot)

    


