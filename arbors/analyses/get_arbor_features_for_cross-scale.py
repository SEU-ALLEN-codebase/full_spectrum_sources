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
    def __init__(self, swcfile, arbor_identifier=None, density_dist_thresh=20., ndigits=4, scale=25.):
        if type(swcfile) is str:
            self.swc_array = load_swc_as_array(swcfile) / scale
        else:
            self.swc_array = np.array(swcfile) / scale

        self.scale = scale
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
        self.soma_pos = np.array(soma_pos) / self.scale

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
        
        coords = np.round(self.swc_array[:, 2:5]).astype(np.int32)[:, ::-1]
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

def calc_basal_features(soma_types, arbor_dir, spos_dict, min_num_neurons=5, median=True, neurite_type='basal', outfile=None, scale=25.):
    regions = []
    feat_names = []
    features_all = []
    for stype, prefixs in soma_types.items():
        print(f'\n<---------------- Soma type: {stype} ---------------->')

        fnames = []
        all_features = []
        for prefix in prefixs:
            print(f'====> neuron: {prefix}')
            #swcfile = os.path.join(arbor_dir, f'{prefix}_{neurite_type}den.swc')
            swcfile = os.path.join(arbor_dir, f'{prefix}.swc')
            if neurite_type == 'basal':
                swcfile = os.path.join(arbor_dir, f'{prefix}_basal.swc')
            if not os.path.exists(swcfile):
                continue
         
            features_list = []
            #af = arborfeatures(swcfile, arbor_identifier=f'{prefix}_{neurite_type}', density_dist_thresh=20., ndigits=4)
            af = ArborFeatures(swcfile, arbor_identifier=f'{prefix}_{neurite_type}', density_dist_thresh=20., ndigits=4)
            
            af.set_center(af.calc_center())
            af.set_soma_pos(spos_dict[prefix])
            features = af.run()
            all_features.append(features)
            print(f'--> features: {features}, center: {af.center}')
            fnames.append(prefix)

        if len(all_features) < min_num_neurons:
            continue

        for curr_features, fname in zip(all_features, fnames):
            feat_names = [f'{fn}_{neurite_type}' for fn in curr_features.keys()]
            feat_values = list(curr_features.values())
            features_all.append([fname, stype, *feat_values])
    
    df = pd.DataFrame(features_all, columns=['prefix', 'region', *feat_names])
    df.to_csv(outfile)

def calc_axon_features(soma_types, arbor_dir, spos_dict, min_num_neurons=5, median=True, neurite_type='axonal', outfile=None):
    #keyword = '_axon'
    keyword = ''
    features_all = []
    for stype, prefixs in soma_types.items():
        print(f'\n<---------------- Soma type: {stype} ---------------->')

        all_features = []
        fnames = []
        for prefix in prefixs:
            print(f'====> neuron: {prefix}')
            stat_file = os.path.join(arbor_dir, f'{prefix}{keyword}.swc.autoarbor_m3.arborstat.txt')
            #print(stat_file)
            if not os.path.exists(stat_file):
                continue
         
            arbor_file = os.path.join(arbor_dir, f'{prefix}{keyword}.swc._m3_lt.eswc')
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
            fnames.append(prefix)
            print('\n')

        if len(all_features) < min_num_neurons:
            continue

        for curr_features, fname in zip(all_features, fnames):
            feat_names = [f'{fn}_axonal1' for fn in curr_features[0].keys()] + \
                         [f'{fn}_axonal2' for fn in curr_features[1].keys()]
            feat_values = list(curr_features[0].values()) + list(curr_features[1].values())
            features_all.append([fname, stype, *feat_values])
    
    df = pd.DataFrame(features_all, columns=['prefix', 'region', *feat_names])
    df.to_csv(outfile)




if __name__ == '__main__':
    celltype_file = '../../common_lib/41586_2021_3941_MOESM4_ESM.csv'
    neurite_type = 'dendrite'

    arbor_dir_dict = {
        'axonal': '/PBshare/SEU-ALLEN/Users/Sujun/2023_12_11_1886data/axon_arbors_l2',
        'basal': '/PBshare/SEU-ALLEN/Users/Sujun/2023_12_11_1886data/basal80_sort',
        'apical': '/PBshare/SEU-ALLEN/Users/Sujun/2023_12_11_1886data/apical80_sort',
        'dendrite': '/PBshare/SEU-ALLEN/Users/Sujun/2023_12_11_1886data/den_80_recolor'
    }

    soma_file = '../data/1886_somalist.txt'
    soma_type_merge = True
    use_abstract_ptype = True
    min_num_neurons = 9
    out_dir = f'min_num_neurons{min_num_neurons}_l2'
    
    arbor_dir = arbor_dir_dict[neurite_type]
    
    feat_file = f'{out_dir}/features_r2_somaTypes_{neurite_type}.csv'
    soma_types, _, p2stypes = load_type_from_excel(celltype_file, use_abstract_ptype=use_abstract_ptype)

    spos_dict = load_soma_pos(soma_file)
 
    if neurite_type == 'axonal':
        calc_axon_features(soma_types, arbor_dir, spos_dict, min_num_neurons=min_num_neurons, median=False, neurite_type=neurite_type, outfile=feat_file)
    elif neurite_type == 'basal' or neurite_type == 'apical' or neurite_type == 'dendrite':
        calc_basal_features(soma_types, arbor_dir, spos_dict, min_num_neurons=min_num_neurons, median=False, neurite_type=neurite_type, outfile=feat_file)
        
    
 

