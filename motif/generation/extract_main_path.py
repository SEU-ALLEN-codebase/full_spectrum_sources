#!/usr/bin/env python

#================================================================
#   Copyright (C) 2022 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : main_projection_analyzer.py
#   Author       : Yufeng Liu
#   Date         : 2022-08-03
#   Description  : 
#
#================================================================
import os
import sys
sys.path.append(r'..\..\pylib')
import glob
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from swc_handler import parse_swc, write_swc, get_specific_neurite, NEURITE_TYPES, scale_swc, flip_swc, find_soma_index
from morph_topo.morphology import Morphology

sys.path.append('../../common_lib')
from common_utils import load_pstype_from_excel, load_type_from_excel


class MainPathValidator(object):
    def __init__(self):
        pass

    def analyze_path_lengths(self, plen_file):
        with open(plen_file) as fp:
            for line in fp.readlines():
                line = line.strip()
                if not line: continue
                ctxt = line.split()
                class_name = ctxt[0]
                lens = np.array(list(map(float, ctxt[1:]))).reshape(-1,2)
                # plot the length distribution
                orig_lens = lens[:,0]
                pruned_lens = lens[:,1]
                
                fig, axes = plt.subplots(2, 2)
                plt.suptitle(class_name)
                
                xmin, xmax = lens.min(), lens.max()
                axes[0,0].set_title('Max path length')
                axes[0,0].hist(orig_lens, bins=20, range=(xmin, xmax), label=f'std:{orig_lens.std():.2f}')

                axes[0,1].set_title('Geodesic length of main path')
                axes[0,1].hist(pruned_lens, bins=20, range=(xmin, xmax), label=f'std:{pruned_lens.std():.2f}')

                ratios = pruned_lens / (orig_lens + 0.1)
                axes[1,0].set_title('Ratio main path length to max')
                axes[1,0].hist(ratios, bins=20, range=(ratios.min(), ratios.max()))
                # legend
                for i in range(2):
                    if i == 1: continue
                    for j in range(2):
                        axes[i,j].legend()
                
                figname = f'{class_name}_path_lengths.png'
                print(f'Saving figure for: {class_name}...')
                fig.tight_layout()
                plt.savefig(figname, dpi=150)
                plt.close()
                


class MainPathProjection(object):
    def __init__(self, swcfile, scale_factor=1.0, flip_to_left=True):
        tree = parse_swc(swcfile)
        if flip_to_left:
            soma_index = find_soma_index(tree)
            zs = tree[soma_index][4]
            if zs > 228 * 25.:
                tree = flip_swc(tree, axis='z', dim=456*25.)

        if scale_factor != 1:
            tree = scale_swc(tree, scale_factor)
        self.morph = Morphology(tree)
        # calculate some feature in advance
        self.tip_paths_dict = self.morph.get_all_paths()
        
    def calc_seg_lengths(self):
        """
        Before pruning short segments while extracting main path, we should know the 
        segment length distribution so that the pruning step is guarenteed. 
        """
        lengths, frag_lengths_dict = self.morph.calc_frag_lengths()
        topo_tree, seg_dict = self.morph.convert_to_topology_tree()
        seg_lengths_dict = self.morph.calc_seg_path_lengths(seg_dict, frag_lengths_dict)
        tip_path_length_dict = self.morph.get_path_len_dict(self.tip_paths_dict, lengths)
        return seg_lengths_dict, tip_path_length_dict


    def extract_main_tract(self, outswc):
        seg_lengths_dict, tip_path_length_dict = self.calc_seg_lengths()
        # select the longest axonal path
        max_length = 0
        max_idx = -1
        for tip, plen in tip_path_length_dict.items():
            if self.morph.pos_dict[tip][1] in NEURITE_TYPES['axon']:
                if plen > max_length:
                    max_length = plen
                    max_idx = tip
        longest_path = self.tip_paths_dict[max_idx]

        # get the all axonal seg_lengths, so as to adaptively prune termini
        axon_seg_lengths = []
        for node_id in longest_path:
            if (node_id in seg_lengths_dict) and (self.morph.pos_dict[node_id][1] in NEURITE_TYPES['axon']):
                axon_seg_lengths.append(seg_lengths_dict[node_id])
        axon_seg_lengths.sort()
        #print(f'Stat of segment length: ')
        #print(f'   max: {max(axon_seg_lengths)}, min: {min(axon_seg_lengths)}')

        # we can remove only the short segments, and do not care whether some remaining
        k = max(len(axon_seg_lengths) - 2, 0)
        length_thr1 = axon_seg_lengths[k]  # the threshold could be refined
        #print(f'length_thr: {length_thr1}')
        # iterative pruning 
        pruned_path = []
        pruned_len = max_length
        for i, idx in enumerate(longest_path):
            #print(i, idx)
            if idx in seg_lengths_dict:
                if seg_lengths_dict[idx] >= length_thr1:
                    #print(seg_lengths_dict[idx])
                    pruned_path = longest_path[i:]
                    break
                else:
                    pruned_len -= seg_lengths_dict[idx]
                    #print(seg_lengths_dict[idx], pruned_len)
        #print(f'Max axonal path length and pruned path length are: {max_length} / {pruned_len}')
        
        pruned_tree = []
        for idx in pruned_path:
            _, type_, x, y, z, r, p = self.morph.pos_dict[idx][:7]
            if (idx in self.morph.bifurcation) or (idx in self.morph.multifurcation):
                type_ = 3
                r = 5
            elif idx == self.morph.idx_soma:
                r = 10
            else:
                type_ = 2
            pruned_tree.append((idx,type_,x,y,z,r,p))

        write_swc(pruned_tree, outswc)

        # save the original longest path
        outswc_l = f'{outswc}_longest.swc'
        tree_l = []
        for idx in longest_path:
            tree_l.append(self.morph.pos_dict[idx][:7])
        write_swc(tree_l, outswc_l)
        
        return max_length, pruned_len

def wrapper(swcfile, scale_factor, outswc):
    print(outswc)
    mpp = MainPathProjection(swcfile, scale_factor=scale_factor)
    # check
    mlen, plen = mpp.extract_main_tract(outswc)

if __name__ == '__main__':
    swc_dir = '../data/registration/S3_registered_ccf'
    out_dir = '../main_tracts_types'
    ctype_file = '../../common_lib/41586_2021_3941_MOESM4_ESM.csv'
    min_files = 0
    scale_factor = 1.
  
    ptypes, rev_ptypes, p2stypes = load_pstype_from_excel(ctype_file)
 
    args_list = []   
    for swcfile in glob.glob(os.path.join(swc_dir, '*')):
        prefix = os.path.split(swcfile)[-1].split('.swc')[0]
        if prefix not in rev_ptypes:
            ctype = 'unk'
        else:
            ctype = rev_ptypes[prefix]

        if ctype is np.nan:
            ctype = 'unk'

        print(prefix)
        brain_id = prefix.split('_')[0]
        outswc = os.path.join(out_dir, f'{ctype}_{prefix}_axonal_tract.swc')
        args_list.append((swcfile, scale_factor, outswc))

    # multiprocessing
    from multiprocessing import Pool
    pool = Pool(processes=24)
    pool.starmap(wrapper, args_list)
    pool.close()
    pool.join()
        

        
    

