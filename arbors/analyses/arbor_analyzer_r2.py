#!/usr/bin/env python

#================================================================
#   Copyright (C) 2022 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : arbor_analyzer.py
#   Author       : Yufeng Liu
#   Date         : 2022-10-13
#   Description  : 
#
#================================================================

import os
import glob
import sys
import pandas as pd
from itertools import permutations
import pickle

import numpy as np
from scipy.spatial import distance_matrix
from sklearn import decomposition
from skimage import morphology
import skimage.io, skimage.transform

import cc3d

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('../../common_lib')
from common_utils import load_pstype_from_excel, load_type_from_excel


def load_stat_file(stat_file):
    stats = []
    with open(stat_file) as fp:
        line = fp.readline()
        for line in fp.readlines():
            aid, ncount, cx, cy, cz = list(map(float, line.split()))
            aid = int(aid)
            ncount = int(ncount)
            stats.append([aid, ncount, cx, cy, cz])
    return stats

def load_swc_as_array(swcfile):
    data = np.loadtxt(swcfile, usecols=(0,1,2,3,4,5,6), skiprows=1)
    return data

def load_soma_pos(soma_file):
    spos_dict = {}
    with open(soma_file) as fp:
        for line in fp.readlines():
            line = line.strip()
            if not line: continue
            if line[0] == '#': continue
            prefix, cx, cy, cz = line.split()
            spos_dict[prefix] = np.array(list(map(float, [cx, cy, cz])))
    return spos_dict

def sort_arbors(features_of_arbors):
    min_score = 1000000
    #return sorted(features_of_arbors, key=lambda x: -x['max_density'] + x['dist_to_soma'])
    return sorted(features_of_arbors, key=lambda x: x['dist_to_soma2'])

def get_arbor_correspondence(subject_features, reference_features):
    # we use max_density and dist_to_soma to determinate the correspondence jointly
    subj_f = np.array([[f['max_density'], f['dist_to_soma2']] for f in subject_features])
    ref_f = np.array([[f['max_density'], f['dist_to_soma2']] for f in reference_features])

    min_score = 1000000
    n = len(subj_f)
    for indices in permutations(range(n), n):
        curr_f = [subj_f[i] for i in indices]
        score = np.fabs(curr_f - ref_f).sum()
        if score < min_score:
            min_score = score
            min_indices = indices
    return list(min_indices)

class ArborFeatures(object):
    def __init__(self, swcfile, arbor_identifier=None, density_dist_thresh=20., ndigits=4):
        if type(swcfile) is str:
            self.swc_array = load_swc_as_array(swcfile)
        else:
            self.swc_array = np.array(swcfile)
        self.density_dist_thresh = density_dist_thresh
        self.ndigits = ndigits
        self.arbor_identifier = arbor_identifier # debug only

    def calc_center(self):
        coords = self.swc_array[:, 2:5]
        center = coords.mean(axis=0)
        return center

    def calc_density(self, return_max_info=True, density_dist_th=None):
        if density_dist_th is None:
            density_dist_th = self.density_dist_thresh

        coords = self.swc_array[:, 2:5]
        pdists = distance_matrix(coords, coords)
        density = (pdists < density_dist_th).sum(axis=0)
        if return_max_info:
            self.max_density = round(density.max(), self.ndigits)
            max_id = np.argmax(density)
            self.max_density_coord = coords[max_id]
        return density

    def max_density(self):
        if not hasattr(self, 'max_density_coord'):
            self.calc_density()
        return self.max_density

    def dist_to_soma2(self):
        if not hasattr(self, 'max_density_coord'):
            self.calc_density()
        dist = np.linalg.norm(self.max_density_coord - self.soma_pos)
        return dist

    def num_nodes(self):
        return self.swc_array.shape[0]

    def set_center(self, center):
        self.center = np.array(center)
    
    def set_soma_pos(self, soma_pos):
        self.soma_pos = np.array(soma_pos)

    def get_center(self):
        return self.center

    def get_soma_pos(self):
        return self.soma_pos

    def dist_to_soma(self):
        dist = np.linalg.norm(self.soma_pos - self.center)
        return round(dist, self.ndigits)

    def detect_hubs(self):
        """
        A hub is defined as:
            - connective regions high density nodes
        Algorithms:
            1. convert the swc into 3D image space
            2. get all the high density nodes
            3. Expand each node
        """
        density_dist_th = 20
        density_th = 50
        density = self.calc_density(return_max_info=False, density_dist_th=density_dist_th)
        high_density = density > density_th
        if high_density.sum() == 0:
            return 0
        print(f'High density ratio: {1.0*high_density.sum()/len(density)}')
        
        coords = np.round(self.swc_array[:, 2:5]).astype(np.int)[:, ::-1]
        cmin = coords.min(axis=0)
        cmax = coords.max(axis=0)
        coords = coords - cmin
        bbox = cmax - cmin + 1
        # initialize an image3D with bbox
        img = np.zeros(bbox, dtype=np.uint8)
        coords_high_density = coords[high_density]
        
        img[tuple(zip(*coords_high_density))] = 1    # masking the points
        # dilation or expansion
        k = 15
        img_dil = morphology.dilation(img, np.ones((k,k,k)))
        _, N = cc3d.connected_components(img_dil, return_N=True)

        # for visualization
        visualize = False
        if visualize and self.arbor_identifier is not None:
            vis = np.zeros(bbox, dtype=np.uint8)
            vis[img_dil == 1] = 45
            vis[tuple(zip(*coords))] = 128
            vis[tuple(zip(*coords_high_density))] = 255
            vis = vis.max(axis=0)
            vis = vis[::-1] #vertical flip
            skimage.io.imsave(f"{self.arbor_identifier}.png", skimage.transform.rescale(vis, 4, order=0))

        return N

    def num_hubs(self):
        return self.detect_hubs()
        

    def total_path_length(self):
        pos_dict = {}
        for node in self.swc_array:
            pos_dict[int(node[0])] = node

        coords1, coords2 = [], []
        for node in self.swc_array:
            idx, pid = int(node[0]), int(node[6])
            if pid in pos_dict:
                coords1.append(node[2:5])
                coords2.append(pos_dict[pid][2:5])
        coords1 = np.array(coords1)
        coords2 = np.array(coords2)

        diff = coords2 - coords1
        tpl = np.linalg.norm(diff, axis=1).sum()
        return round(tpl, self.ndigits)

    def calc_bbox3d_r(self):
        coords = self.swc_array[:, 2:5]
        pca = decomposition.PCA()
        coords_transformed = pca.fit_transform(coords)
        bbox3d_r = coords_transformed.max(axis=0) - coords_transformed.min(axis=0)
        self.bbox3d_r = bbox3d_r
        return bbox3d_r

    def variance_ratio(self):
        if self.swc_array.shape[0] < 10:
            return 0

        if not hasattr(self, 'bbox3d_r'):
            self.calc_bbox3d_r()
        pc1, pc2, pc3 = np.sort(self.bbox3d_r)
        return pc3 / (pc1 + pc2 + pc3)

    def volume(self):
        if not hasattr(self, 'bbox3d_r'):
            self.calc_bbox3d_r()
        volume = np.prod(self.bbox3d_r)
        return round(volume, self.ndigits)

    def num_branches(self):
        cn_dict = {}
        for node in self.swc_array:
            idx, pid = int(node[0]), int(node[6])
            try:
                cn_dict[pid] += 1
            except KeyError:
                cn_dict[pid] = 1

        num_branches = 0
        for key, nc in cn_dict.items():
            if nc > 1:
                num_branches += 1

        return num_branches

    def run(self):
        features = {}
        feature_keys = ['max_density', 'num_nodes', 'total_path_length', 'volume', 'num_branches', 'dist_to_soma', 'dist_to_soma2', 'num_hubs', 'variance_ratio']
        #feature_keys = ['max_density', 'volume', 'num_branches', 'dist_to_soma2', 'num_hubs']
        for key in feature_keys:
            func = getattr(self, key)
            feature = func()
            #if (type(feature) is list) or (type(feature) is tuple):
            features[key] = feature   
                
        return features

class NeuriteArbors(object):
    def __init__(self, stat_file, arbor_lt_file, soma_pos):
        self.stats = load_stat_file(stat_file)
        self.load_arbors(arbor_lt_file, soma_pos)

    def num_arbors(self):
        return len(self.arbors)

    def load_arbors(self, arbor_lt_file, soma_pos):
        prefix = os.path.split(arbor_lt_file)[-1].replace('_axon.swc._m3_lt.eswc', '')

        swc_array = load_swc_as_array(arbor_lt_file)
        arbors = []
        for i in range(len(self.stats)):
            flag = swc_array[:,1] == i
            af = ArborFeatures(swc_array[flag], arbor_identifier=f'{prefix}_{i}')
            af.set_center(self.stats[i][2:5])
            af.set_soma_pos(soma_pos)
            arbors.append(af)
        
        self.arbors = arbors

class AxonalArbors(NeuriteArbors):
    def __init__(self, stat_file, arbor_file, soma_pos):
        super(AxonalArbors, self).__init__(stat_file, arbor_file, soma_pos)

def calc_basal_features(soma_types, arbor_dir, spos_dict, min_num_neurons=5, median=True, neurite_type='basal'):
    fdict = {}   
    for stype, prefixs in soma_types.items():
        print(f'\n<---------------- Soma type: {stype} ---------------->')

        all_features = []
        for prefix in prefixs:
            print(f'====> neuron: {prefix}')
            swcfile = os.path.join(arbor_dir, f'{prefix}_{neurite_type}den.swc')
            if not os.path.exists(swcfile):
                continue
         
            features_list = []
            af = ArborFeatures(swcfile, arbor_identifier=f'{prefix}_{neurite_type}', density_dist_thresh=20., ndigits=4)
            af.set_center(af.calc_center())
            af.set_soma_pos(spos_dict[prefix])
            features = [af.run()]
            all_features.append(features)
            print(f'--> features: {features}, center: {af.center}')

        if len(all_features) < min_num_neurons:
            continue

        # get the median
        mfeatures = []
        keys = list(all_features[0][0].keys())
        num_neurons = len(all_features)
        num_arbors = len(all_features[0])
        for i in range(num_arbors):
            f_i = {}
            for key in keys:
                values = []
                for j in range(num_neurons):
                    values.append(all_features[j][i][key])
                if median:
                    vm = np.median(values)
                else:
                    vm = np.mean(values)
                std = np.std(values)

                f_i[key] = (vm, std)
            mfeatures.append(f_i)
        # check and re-order according to distance to soma
        mfeatures.sort(key=lambda x: x['dist_to_soma2'][0])
        print(mfeatures)

        fdict[stype] = mfeatures

    return fdict

def calc_axon_features(soma_types, arbor_dir, spos_dict, min_num_neurons=5, median=True, neurite_type='axonal'):
    fdict = {}
    if neurite_type == 'axonal':
        keyword = 'axon'
    elif neurite_type == 'apical':
        keyword = 'apicalden'

    for stype, prefixs in soma_types.items():
        print(f'\n<---------------- Soma type: {stype} ---------------->')

        all_features = []
        for prefix in prefixs:
            print(f'====> neuron: {prefix}')
            stat_file = os.path.join(arbor_dir, f'{prefix}_{keyword}.swc.autoarbor_m3.arborstat.txt')
            #print(stat_file)
            if not os.path.exists(stat_file):
                continue
         
            arbor_file = os.path.join(arbor_dir, f'{prefix}_{keyword}.swc._m3_lt.eswc')
            aa = AxonalArbors(stat_file, arbor_file, spos_dict[prefix])
            print(f'    arbor num: {aa.num_arbors()}')
            features_list = []
            for af in aa.arbors:
                features = af.run()
                features_list.append(features)
                print(f'--> features: {features}, center: {af.center}')

            # find correspondence and estimate the median feature for each class
            if len(all_features) == 0:
                # sort by distance to soma
                features_list = sort_arbors(features_list)
                all_features.append(features_list)
                reference_features = features_list
            else:
                min_indices = get_arbor_correspondence(features_list, reference_features)
                #print(min_indices)
                sorted_features = [features_list[i] for i in min_indices]
                all_features.append(sorted_features)
            print('\n')

        if len(all_features) < min_num_neurons:
            continue

        # get the median
        mfeatures = []
        keys = list(all_features[0][0].keys())
        num_neurons = len(all_features)
        num_arbors = len(all_features[0])
        for i in range(num_arbors):
            f_i = {}
            for key in keys:
                values = []
                for j in range(num_neurons):
                    values.append(all_features[j][i][key])
                if median:
                    vm = np.median(values)
                else:
                    vm = np.mean(values)
                std = np.std(values)

                f_i[key] = (vm, std)
            mfeatures.append(f_i)
        # check and re-order according to distance to soma
        mfeatures.sort(key=lambda x: x['dist_to_soma2'][0])
        print(mfeatures)

        fdict[stype] = mfeatures

    return fdict

def plot_stypes_features(fdict, p2stypes, neurite_type='axonal', plot_value='mean'):
    """
    fdict structure:
        - key: cell_type
        - value: list of dict
            - item: dict
                - key: feat_name
                - value: list[float], [feat_mean, feat_std]
    """
    if plot_value == 'mean':
        fidx = 0
    elif plot_value == 'std':
        fidx = 1
    else:
        raise ValueError

    # each feature a separate subfig
    cnames = []
    for ptype, stypes in p2stypes.items():
        for stype in sorted(stypes):
            if stype in fdict:
                cnames.append(stype)
    
    if neurite_type == 'axonal':
        keys = ['max_density', 'volume', 'num_branches', 'dist_to_soma2']
    elif neurite_type == 'basal':
        keys = ['volume', 'num_branches']
    elif neurite_type == 'apical':
        keys = ['max_density', 'volume', 'num_branches', 'dist_to_soma2']
    m = len(cnames)
    n = 0
    for ctype in cnames:
        flist = fdict[ctype]
        n = max(len(flist), n)

    tmat = np.zeros((n, m), dtype=np.int)
    mi = 0
    dist_th = 30.0
    for ctype in cnames:
        flist = fdict[ctype]
        for fi, fs in enumerate(flist):
            d2s = fs['dist_to_soma2'][0]
            if d2s < dist_th:
                tmat[fi, mi] = 1
            else:
                tmat[fi, mi] = 2 
        mi += 1
    arbor_type = []
    tdict = {0:' ', 1:'Local', 2:'Distal'}
    for tmi in tmat.reshape(-1):
        arbor_type.append(tdict[tmi])

    # concat all data into a unit dataframe
    features = []
    for ikey, key in enumerate(keys):
        fmat = np.zeros((n, m), dtype=np.float)
        mi = 0
        for ctype in cnames:
            flist = fdict[ctype]
            for ifeat, feat in enumerate(flist):
                #fmat[ifeat, mi] = feat[key][fidx]
                if fidx == 0:
                    fmat[ifeat, mi] = feat[key][fidx]
                else:
                    fmat[ifeat, mi] = feat[key][1] / (feat[key][0]+1e-10)
            mi += 1
        if key != 'variance_ratio' and fidx == 0:
            fmat = (fmat - fmat.min()) / (fmat.max() - fmat.min() + 1e-10)
        if (fidx == 1) and (neurite_type in ['apical', 'basal']):
            fmat = np.clip(fmat, 0, 2.5)  # the dist_to_soma for apical & basal is sensitive
        
        feature_value_name = 'scale'
        df = pd.DataFrame(data=fmat, columns=cnames)
        df_mat0 = df.stack().reset_index(name=feature_value_name)
        df_mat0['arbor_type'] = arbor_type
        df_mat0['feature'] = [key for i in range(n*m)]
        

        if ikey == 0:
            df_mat = df_mat0
        else:
            df_mat = pd.concat([df_mat, df_mat0], ignore_index=True)

    # Draw each cell as a scatter point with varying size and color
    if n == 3:
        #palette = ['w', 'b', 'r']
        palette = {
                ' ': 'w',
                'Local': 'b',
                'Distal': 'r'
            }
    elif n == 2:
        palette = ['w', 'b']
    elif n == 1:
        if neurite_type == 'apical':
            palette = ['m']
        elif neurite_type == 'basal':
            palette = ['g']
    
    
    base_size = 15
    plt.rcParams['axes.labelsize'] = base_size * 1.5
    g = sns.relplot(
        data=df_mat,
        x="level_0", y="level_1", hue="arbor_type", size=feature_value_name, col="feature",
        palette=palette, edgecolor="1.",
        height=10,
        aspect=0.2,
        sizes=(30, 300),
    )
    g.set_titles(col_template="{col_name}")
    g.fig.subplots_adjust(wspace=0.1)
    g.sharex = True
    plt.setp(g._legend.get_texts(), fontsize=base_size)

    # Tweak the figure to finalize
    #g.set(xlabel="Arbors", ylabel="Cell types")
    g.set(xlabel="", ylabel='Cell types')
    if fidx == 0:
        fname = 'normalized mean'
    else:
        fname = 'mean-normalized std'
    g.fig.suptitle(f'Features of {neurite_type} arbors', x=0.5, y=1.00, va='bottom', fontsize=int(base_size*1.8))
    for i, ax in enumerate(g.axes[0]):
        ax.set_xlim(-1,n)
        ax.set_xticks(list(range(n)))
        ax.tick_params(axis='x', colors=(0,0,0,0), labelsize=0)
        ax.tick_params(axis='y', labelsize=base_size)
        ax.grid(alpha=0.5)
        ax.margins(0.02)
        kname = keys[i]
        if kname == 'num_branches':
            kname = '#branches'
        elif kname == 'dist_to_soma2':
            kname = 'dist2soma'
        ax.set_title(kname, fontsize=18)
        
    plt.savefig(f'{neurite_type}_arbors_f{plot_value}.png', bbox_inches='tight', pad_inches=0.1, dpi=300)


def plot_ptypes_features(fdict, p2stypes, neurite_type='axonal', plot_value='mean'):
    """
    fdict structure:
        - key: cell_type
        - value: list of dict
            - item: dict
                - key: feat_name
                - value: list[float], [feat_mean, feat_std]
    """
    if plot_value == 'mean':
        fidx = 0
    elif plot_value == 'std':
        fidx = 1
    else:
        raise ValueError

    # each feature a separate subfig
    cnames = sorted(fdict.keys())

    # get only some types
    keep_stypes = set(['MOp', 'MOs', 'SSp-bfd', 'SSp-m'])
    keep_ptypes = []
    for ctype in cnames:
        ctype_name = '-'.join(ctype.split('-')[1:])
        if ctype_name in keep_stypes:
            keep_ptypes.append(ctype)

    
    if neurite_type == 'axonal':
        keys = ['max_density', 'volume', 'num_branches', 'dist_to_soma2']
    elif neurite_type == 'basal':
        keys = ['volume', 'num_branches']
    elif neurite_type == 'apical':
        keys = ['volume', 'num_branches']
    m = len(keep_ptypes)
    n = 0
    for ctype in keep_ptypes:
        flist = fdict[ctype]
        n = max(len(flist), n)

    tmat = np.zeros((n, m), dtype=np.int)
    mi = 0
    dist_th = 30.0
    for ctype in keep_ptypes:
        flist = fdict[ctype]
        for fi, fs in enumerate(flist):
            d2s = fs['dist_to_soma2'][0]
            if d2s < dist_th:
                tmat[fi, mi] = 1
            else:
                tmat[fi, mi] = 2 
        mi += 1

    arbor_type = []
    tdict = {0:' ', 1:'Local', 2:'Distal'}
    for tmi in tmat.reshape(-1):
        arbor_type.append(tdict[tmi])

    # concat all data into a unit dataframe
    features = []
    for ikey, key in enumerate(keys):
        fmat = np.zeros((n, m), dtype=np.float)
        mi = 0
        for ctype in keep_ptypes:
            flist = fdict[ctype]
            for ifeat, feat in enumerate(flist):
                #fmat[ifeat, mi] = feat[key][fidx]
                if fidx == 0:
                    fmat[ifeat, mi] = feat[key][fidx]
                else:
                    fmat[ifeat, mi] = feat[key][1] / (feat[key][0]+1e-10)
            mi += 1
        if key != 'variance_ratio' and fidx == 0:
            fmat = (fmat - fmat.min()) / (fmat.max() - fmat.min() + 1e-10)
        if (fidx == 1) and (neurite_type in ['apical', 'basal']):
            fmat = np.clip(fmat, 0, 2.5)  # the dist_to_soma for apical & basal is sensitive
        
        feature_value_name = 'scale'
        print(fmat.shape, len(keep_ptypes))
        df = pd.DataFrame(data=fmat, columns=keep_ptypes)
        df_mat0 = df.stack().reset_index(name=feature_value_name)
        df_mat0['arbor_type'] = arbor_type
        df_mat0['feature'] = [key for i in range(m*n)]
        

        if ikey == 0:
            df_mat = df_mat0
        else:
            df_mat = pd.concat([df_mat, df_mat0], ignore_index=True)

    # Draw each cell as a scatter point with varying size and color
    if n == 3:
        #palette = ['w', 'b', 'r']
        palette = {
                ' ': 'w',
                'Local': 'b',
                'Distal': 'r'
            }
    elif n == 2:
        palette = {
                'Local': 'b',
                'Distal': 'r'
        }
    elif n == 1:
        if neurite_type == 'apical':
            palette = ['m']
        elif neurite_type == 'basal':
            palette = ['g']
    
    if neurite_type == 'basal' or neurite_type == 'apical':
        aspect = 0.5
    else:
        aspect = 0.3    

    base_size = 8
    plt.rcParams['axes.labelsize'] = base_size * 1.5
    g = sns.relplot(
        data=df_mat,
        x="level_0", y="level_1", hue="arbor_type", size=feature_value_name, col="feature",
        palette=palette, edgecolor="1.",
        height=4,
        aspect=aspect,
        sizes=(30, 300),
    )
    g.set_titles(col_template="{col_name}")
    g.fig.subplots_adjust(wspace=0.1)
    g.sharex = True
    plt.setp(g._legend.get_texts(), fontsize=base_size)

    # Tweak the figure to finalize
    #g.set(xlabel="Arbors", ylabel="Cell types")
    g.set(xlabel="", ylabel='Cell types')
    if fidx == 0:
        fname = 'normalized mean'
    else:
        fname = 'mean-normalized std'
    g.fig.suptitle(f'Features of {neurite_type} arbors', x=0.5, y=1.00, va='bottom', fontsize=int(base_size*1.8))
    for i, ax in enumerate(g.axes[0]):
        ax.set_xlim(-1,n)
        ax.set_xticks(list(range(n)))
        ax.tick_params(axis='x', colors=(0,0,0,0), labelsize=0)
        ax.tick_params(axis='y', labelsize=base_size)
        ax.grid(alpha=0.5)
        ax.margins(0.1)
        kname = keys[i]
        if kname == 'num_branches':
            kname = '#branches'
        elif kname == 'dist_to_soma2':
            kname = 'dist2soma'
        ax.set_title(kname, fontsize=base_size*1.2)
        
    plt.savefig(f'{neurite_type}_arbors_f{plot_value}.png', bbox_inches='tight', pad_inches=0.1, dpi=300)


def plot_ptypes_features_dendrite(fdicts, p2stypes, neurite_type='dendrite', plot_value='mean'):
    if plot_value == 'mean':
        fidx = 0
    elif plot_value == 'std':
        fidx = 1
    else:
        raise ValueError

    # each feature a separate subfig
    cnames = sorted(fdicts[0].keys())

    # get only some types
    keep_stypes = set(['MOp', 'MOs', 'SSp-bfd', 'SSp-m'])
    keep_ptypes = []
    for ctype in cnames:
        ctype_name = '-'.join(ctype.split('-')[1:])
        if ctype_name in keep_stypes:
            keep_ptypes.append(ctype)

    
    keys = ['volume', 'num_branches']
    m = len(keep_ptypes)
    n = 2

    # do not consider local and nonlocal
    tmat = np.ones((n, m), dtype=np.int)
    tmat[0] = 1
    tmat[1] = 2

    arbor_type = []
    tdict = {0:' ', 1:'basal', 2:'apical'}
    for tmi in tmat.reshape(-1):
        arbor_type.append(tdict[tmi])

    # concat all data into a unit dataframe
    features = []
    for ikey, key in enumerate(keys):
        fmat = np.zeros((n, m), dtype=np.float)
        mi = 0
        for ctype in keep_ptypes:
            flist1 = fdicts[0][ctype]
            flist2 = fdicts[1][ctype]
            for ifeat, feat in enumerate(flist1):
                if fidx == 0:
                    fmat[ifeat*2, mi] = feat[key][fidx]
                    fmat[ifeat*2+1, mi] = flist2[ifeat][key][fidx]
                else:
                    fmat[ifeat*2, mi] = feat[key][1] / (feat[key][0]+1e-10)
                    fmat[ifeat*2+1, mi] = flist2[ifeat][key][1] / (flist2[ifeat][key][0]+1e-10)
            mi += 1
        if key != 'variance_ratio' and fidx == 0:
            fmat = (fmat - fmat.min()) / (fmat.max() - fmat.min() + 1e-10)
        if (fidx == 1) and (neurite_type in ['apical', 'basal']):
            fmat = np.clip(fmat, 0, 2.5)  # the dist_to_soma for apical & basal is sensitive
        
        feature_value_name = 'scale'
        print(fmat.shape, len(keep_ptypes))
        df = pd.DataFrame(data=fmat, columns=keep_ptypes)
        df_mat0 = df.stack().reset_index(name=feature_value_name)
        df_mat0['arbor_type'] = arbor_type
        df_mat0['feature'] = [key for i in range(m*n)]
        

        if ikey == 0:
            df_mat = df_mat0
        else:
            df_mat = pd.concat([df_mat, df_mat0], ignore_index=True)

    # Draw each cell as a scatter point with varying size and color
    palette = {
        'basal': 'g',
        'apical': 'm'
    }
    
    base_size = 8
    plt.rcParams['axes.labelsize'] = base_size * 1.5
    g = sns.relplot(
        data=df_mat,
        x="level_0", y="level_1", hue="arbor_type", size=feature_value_name, col="feature",
        palette=palette, edgecolor="1.",
        height=4,
        aspect=0.5,
        sizes=(30, 300),
    )
    g.set_titles(col_template="{col_name}")
    g.fig.subplots_adjust(wspace=0.1)
    g.sharex = True
    plt.setp(g._legend.get_texts(), fontsize=base_size)

    # Tweak the figure to finalize
    #g.set(xlabel="Arbors", ylabel="Cell types")
    g.set(xlabel="", ylabel='Cell types')
    if fidx == 0:
        fname = 'normalized mean'
    else:
        fname = 'mean-normalized std'
    g.fig.suptitle(f'Features of {neurite_type} arbors', x=0.5, y=1.00, va='bottom', fontsize=int(base_size*1.8))
    for i, ax in enumerate(g.axes[0]):
        ax.set_xlim(-1,n)
        ax.set_xticks(list(range(n)))
        ax.tick_params(axis='x', colors=(0,0,0,0), labelsize=0)
        ax.tick_params(axis='y', labelsize=base_size)
        ax.grid(alpha=0.5)
        ax.margins(0.1)
        kname = keys[i]
        if kname == 'num_branches':
            kname = '#branches'
        elif kname == 'dist_to_soma2':
            kname = 'dist2soma'
        ax.set_title(kname, fontsize=base_size*1.2)
        
    plt.savefig(f'{neurite_type}_arbors_f{plot_value}.png', bbox_inches='tight', pad_inches=0.1, dpi=300)




if __name__ == '__main__':
    #celltype_file = '/PBshare/SEU-ALLEN/Users/yfliu/transtation/seu_mouse/swc/1741_Celltype.csv'
    celltype_file = '../../common_lib/41586_2021_3941_MOESM4_ESM.csv'
    neurite_type = 'axonal'

    arbor_dir_dict = {
        'axonal': '../data/axon_arbors_round2_ln',
        'basal': '../data/basal80_sort',
        'apical': '../data/apical80_sort',
        'dendrite': ['../data/basal80_sort', '../data/apical80_sort']
    }

    soma_file = '../data/1891_somalist.txt'
    soma_type_merge = True
    use_abstract_ptype = True
    min_num_neurons = 10
    out_dir = f'min_num_neurons{min_num_neurons}_l2'
    
    arbor_dir = arbor_dir_dict[neurite_type]
    
    if soma_type_merge:
        feat_file = f'{out_dir}/features_r2_somaTypes_{neurite_type}.pkl'
        soma_types, _, p2stypes = load_type_from_excel(celltype_file, use_abstract_ptype=use_abstract_ptype)
    else:
        feat_file = f'{out_dir}/features_r2_projAndSomaTypes_{neurite_type}.pkl'
        soma_types, _, p2stypes = load_pstype_from_excel(celltype_file, use_abstract_ptype=use_abstract_ptype)

    spos_dict = load_soma_pos(soma_file)
    print(p2stypes)
 
    if 1:
        if neurite_type == 'axonal':
            fdict = calc_axon_features(soma_types, arbor_dir, spos_dict, min_num_neurons=min_num_neurons, median=False, neurite_type=neurite_type)
        elif neurite_type == 'basal' or neurite_type == 'apical':
            fdict = calc_basal_features(soma_types, arbor_dir, spos_dict, min_num_neurons=min_num_neurons, median=False, neurite_type=neurite_type)
        
        # save to file
        with open(feat_file, 'wb') as fp:
            pickle.dump(fdict, fp)
    
 
    if 0:   
        with open(feat_file, 'rb') as fp:
            fdict = pickle.load(fp)
        for plot_value in ['mean', 'std']:
            if soma_type_merge:
                plot_stypes_features(fdict, p2stypes, neurite_type, plot_value=plot_value)
            else:
                plot_ptypes_features(fdict, p2stypes, neurite_type, plot_value=plot_value)
    if 0:
        # plot ptype for both basal and apical
        assert(not soma_type_merge)
        fdicts = []
        with open(f'{out_dir}/features_r2_projAndSomaTypes_basal.pkl', 'rb') as fp:
            fdicts.append(pickle.load(fp))
        with open(f'{out_dir}/features_r2_projAndSomaTypes_apical.pkl', 'rb') as fp:
            fdicts.append(pickle.load(fp))
        
        for plot_value in ['mean', 'std']:
            plot_ptypes_features_dendrite(fdicts, p2stypes, neurite_type, plot_value=plot_value)

