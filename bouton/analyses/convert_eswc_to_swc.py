#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : convert_eswc_to_swc.py
#   Author       : Yufeng Liu
#   Date         : 2023-01-03
#   Description  : 
#
#================================================================

import os
import glob
import numpy as np
import pandas as pd
import pickle

import sys
sys.path.append('../../pylib')
from swc_handler import write_swc, parse_swc

def to_swc(eswc, swc):
    #COL_NAMES = ['type', 'x', 'y', 'z', 'radius', 'parent', 'seg_id',
    #             'level', 'mode', 'timestamp', 'teraflyindex', 'x_ccf',
    #             'y_ccf', 'z_ccf', 'flag', 'bou_radius', 'bra_radius_mean',
    #             'bra_radius_std', 'bou_intensity', 'bra_intensity_mean',
    #             'bra_intensity_std', 'bou_density', 'pdist_to_soma',
    #             'unknown', 'edist_to_soma']
    COL_NAMES = ['id', 'type', 'x', 'y', 'z', 'radius', 'parent', 'seg_id',
                 'level', 'mode', 'timestamp', 'teraflyindex', 'x_ccf',
                 'y_ccf', 'z_ccf', 'radius2', 'flag', 'u2', 'u3', 'u4', 'u5', 'u6',
                 'u7']
    data = pd.read_csv(eswc, sep=' ', header=2, names=COL_NAMES, index_col=0)
    tree = []
    for irow, row in data.iterrows():
        node = (irow, int(row['type']), row['x_ccf']/25., row['y_ccf']/25., row['z_ccf']/25., 
                row['flag'], int(row['parent']))
        tree.append(node)
    #print(len(tree), tree[0])
    write_swc(tree, swc)

def swc_to_pickle(swc_file, pickle_file):
    print(swc_file)
    tree = parse_swc(swc_file)
    with open(pickle_file, 'wb') as fp:
        pickle.dump(tree, fp)

    
if __name__ == '__main__':
    eswc_dir = '/data/lyf/data/fullNeurons/V2023_12_11/beswc_featured'
    swc_dir = 'bouton_v20231211_swc'
    pickle_dir = 'bouton_v20231211_swc_pickle'

    
    args_list = []
    for eswc_file in glob.glob(os.path.join(eswc_dir, '*eswc')):
        fn = os.path.splitext(os.path.split(eswc_file)[-1])[0]
        print(fn)
        swc_file = os.path.join(swc_dir, f'{fn}.swc')
        to_swc(eswc_file, swc_file)
        if os.path.exists(swc_file):
            continue
        args_list.append((eswc_file, swc_file))

    from multiprocessing import Pool
    pool = Pool(processes=24)
    pool.starmap(to_swc, args_list)
    pool.close()
    pool.join()
    '''

    
    args_list = []
    for swc_file in glob.glob(os.path.join(swc_dir, '*swc')):
        fn = os.path.split(swc_file)[-1]
        pkl_file = os.path.join(pickle_dir, f'{fn}.pkl')
        args_list.append((swc_file, pkl_file))
        
    from multiprocessing import Pool
    pool = Pool(processes=24)
    pool.starmap(swc_to_pickle, args_list)
    pool.close()
    pool.join()
    '''
