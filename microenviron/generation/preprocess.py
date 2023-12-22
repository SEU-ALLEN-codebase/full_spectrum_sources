#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : preprocess.py
#   Author       : Yufeng Liu
#   Date         : 2023-01-17
#   Description  : 
#
#================================================================
import sys
sys.path.append('../../pylib')
import os
import glob
import re
import numpy as np
import pandas as pd
from sklearn import decomposition

from file_io import load_image
from anatomy.anatomy_config import MASK_CCF25_FILE
from anatomy.anatomy_core import parse_ana_tree

import sys
sys.path.append('../../common_lib')
from common_func import load_regions, get_region_mapper

def add_r316_regions(df, ana_dict, mask_regions):
    rids_set = load_regions(level=0)
    orig_ids = mask_regions

    rc_dict = get_region_mapper(rids_set, orig_ids, ana_dict, reverse_mapper=True)[1]
    
    new_rids = []
    new_rnames = []
    brain_structures = []
    struct_dict = {
        688: 'CTX',
        623: 'CNU',
        512: 'CB',
        343: 'BS'
    }
    for idx in df['region_id_r671']:
        if idx in rc_dict:
            new_id = rc_dict[idx]
            new_id = new_id
            new_rids.append(new_id)
            new_rname = ana_dict[new_id]['acronym']
            new_rnames.append(new_rname)
            
            # find the brain structure
            for pid in ana_dict[new_id]['structure_id_path']:
                if pid in struct_dict:
                    struct_name = struct_dict[pid]
                    brain_structures.append(struct_name)
                    break
            else:
                print(f'!! {new_rname}')
                brain_structures.append(np.NaN)
                
        else:
            new_rids.append(0)
            new_rnames.append(np.NaN)
            brain_structures.append(np.NaN)

    df['region_id_r316'] = new_rids
    df['region_name_r316'] = new_rnames
    df['brain_structure'] = brain_structures
    
    return df
        

def aggregate_information(feature_dir, swc_dir):
    """
       Aggregate more information besides the 22d features, and combine all features into an
     unified file. Additional information including: 
       - soma_pos[x,y,z] in CCFv3 space
       - registered regions in the finest granuity, that is 671
       - name of regions
       - brain id
       - dataset name: 150k or 50k
       - and original name and features
    """
    DEBUG = False

    # soma positions
    swcnames = []
    file_infos = []
    for ds in glob.glob(os.path.join(swc_dir, '*')):
        dname = os.path.split(ds)[-1]
        for brain_dir in glob.glob(os.path.join(ds, '[1-9]*[0-9]')):
            brain_id = os.path.split(brain_dir)[-1]
            if DEBUG:
                if brain_id != '17051': continue
            print(f'--> Loading soma position for {dname}/{brain_id}')
            for swcfile in glob.glob(os.path.join(brain_dir, '*.swc')):
                # find the soma postion from file with decoding
                with open(swcfile) as fp:
                    soma_str = re.search('.* -1\n', fp.read()).group()
                # get the orientation of morphology
                coords = np.genfromtxt(swcfile, usecols=(2,3,4))
                pca = decomposition.PCA()
                pca.fit(coords)

                spos = soma_str.split()[2:5]
                swcnames.append(os.path.split(swcfile)[-1])
                file_infos.append((*spos, dname, brain_id, *pca.components_[0], *pca.explained_variance_ratio_))
    
    df = pd.DataFrame(file_infos, 
        columns=['soma_x', 'soma_y', 'soma_z', 'dataset_name', 'brain_id', 
                 'pc11', 'pc12', 'pc13', 'pca_vr1', 'pca_vr2', 'pca_vr3'],
        index=swcnames
    ).astype({'soma_x': np.float64, 'soma_y': np.float64, 'soma_z': np.float64, 'pc11': np.float64,
            'pc12': np.float64, 'pc13': np.float64, 'pca_vr1': np.float64, 'pca_vr2': np.float64,
            'pca_vr3': np.float64})
    print(df.dtypes)
    
    # get the region
    print('>> Extract the region id and name')
    ana_dict = parse_ana_tree()
    mask = load_image(MASK_CCF25_FILE)  # z,y,x order!
    mask_regions = [idx for idx in np.unique(mask) if idx != 0]

    sposs = np.round(df[['soma_x', 'soma_y', 'soma_z']].to_numpy()).astype(int)
    # clip
    np.clip(sposs, 0, np.array(mask.shape[::-1])-1, out=sposs)
    # region accessing in vectorized style
    region_ids = mask[sposs[:,2], sposs[:,1], sposs[:,0]]
    region_names = [ana_dict[idx]['acronym'] if idx != 0 else np.NaN for idx in region_ids]

    df['region_id_r671'] = region_ids
    df['region_name_r671'] = region_names
    # also get the r316 
    add_r316_regions(df, ana_dict, mask_regions)

    # also merge the L-measure features
    print('===> Merge with the features')
    fdim = 22
    features = []
    #for ds in glob.glob(os.path.join(feature_dir, '*')):
    #    dname = os.path.split(ds)[-1]
    #    for brain_file in glob.glob(os.path.join(ds, '*.txt')):
    for brain_file in glob.glob(os.path.join(feature_dir, '*csv')):
        prefix = os.path.split(brain_file)[-1][:-4]
        brain_id = prefix.split('_')[-1]
        dname = '_'.join(prefix.split('_')[-3:-1])

        if DEBUG:
            if brain_id != '17051': continue
        print(f'======> parsing feature for {dname}/{brain_id}')
        brain_features = pd.read_csv(brain_file, index_col=0)
        features.append(brain_features)
    features = pd.concat(features, join='inner').astype(np.float64)
    # concat with the overall information
    df = df.join(features, how='inner')

    df.to_csv('lm_features_d22_all.csv', float_format='%g')
    
    
    


if __name__ == '__main__':
    # The directory `swc_dir` contains all the reconstructed morphologies in the format of SWC. The morphologies are registered to the CCFv3 space.
    swc_dir = '/PBshare/SEU-ALLEN/Users/yfliu/transtation/Research/platform/micro_environ/data/improved_reg/42k_local_morphology_new20230510_gcoord_final'
    # Driectory containing the global features for each neuron
    feature_dir = './lm_features'
    aggregate_information(feature_dir, swc_dir)

    
