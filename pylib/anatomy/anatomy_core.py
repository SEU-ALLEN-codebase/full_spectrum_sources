#!/usr/bin/env python

#================================================================
#   Copyright (C) 2022 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : anatomy_core.py
#   Author       : Yufeng Liu
#   Date         : 2022-11-10
#   Description  : 
#
#================================================================

import json
import numpy as np
import pandas as pd

from anatomy.anatomy_config import *
from file_io import load_image

def parse_ana_tree(tree_file=None, map_file=None, keyname='id'):
    if tree_file is None:
        tree_file = ANATOMY_TREE_FILE

    with open(tree_file) as f1:
        tree = json.load(f1)

    mapper, rev_mapper = parse_id_map(map_file)   

    ids_dict = {}
    for reg in tree:
        name = reg['acronym']
        idx = reg['id']
        id_path = reg['structure_id_path']
        for _idx in id_path:
            if _idx in ids_dict:
                ids_dict[_idx].append(idx)
            else:
                ids_dict[_idx] = [idx]
    # to set
    for key, value in ids_dict.items():
        ids_dict[key] = list(set(value))
 
    ana_dict = {}
    for reg in tree:
        name = reg['acronym']
        idx = reg['id']
        reg['mapped_id'] = rev_mapper[idx]
        reg['orig_ids'] = ids_dict[idx]
        mapped_ids = []
        for _idx in ids_dict[idx]:
            mapped_ids.extend(rev_mapper[_idx])
        reg['mapped_ids'] = mapped_ids

        if keyname == 'name':
            ana_dict[name] = reg
        elif keyname == 'id':
            ana_dict[idx] = reg
        else:
            raise ValueError

    return ana_dict

def parse_regions316(regions_file=None):
    if regions_file is None:
        regions_file = REGIONS316_FILE

    regions = pd.read_excel(regions_file)
    return regions
        
        
def parse_id_map(map_file=None):
    if map_file is None:
        map_file = ID_MAP_RES25_FILE

    with open(map_file) as f2:
        mapper = json.load(f2)
        # to string
        new_mapper = {}
        keys = list(mapper.keys())
        for key in keys:
            new_mapper[int(key)] = mapper[key]

        rev_mapper = {}
        for id1, id2 in new_mapper.items():
            try:
                rev_mapper[id2].append(id1)
            except KeyError:
                rev_mapper[id2] = [id1]
    
    return new_mapper, rev_mapper

def get_regional_neighbors(mask_file=None, num_voxels=5):
    """
    Find neighbor region ids for every region in given mask brain
    """
    import skimage.morphology

    if mask_file is None:
        mask_file = MASK_CCF25_FILE
    mask = load_image(mask_file)
    if mask.ndim == 4:
        mask = mask[0]
    values = np.unique(mask)

    rn_dict = {}
    processed = 0
    for v in values:
        processed += 1
        if v == 0:  # out-of-brain voxel
            continue

        mask_i = mask == v
        mask_i_dil = skimage.morphology.binary_dilation(mask_i, skimage.morphology.ball(radius=num_voxels))
        neighbors = np.unique(mask[mask_i_dil])
        # remove self
        neighbors = [ni for ni in neighbors if ni != v]
        rn_dict[v] = neighbors

        print(f'--> Processed: {processed} / {len(values)}')

    return rn_dict

def get_regional_neighbors_cuda(mask_file=None, radius=5):
    """
       Pytorch-CUDA version of the function `get_regional_neighbors`. This version if orders of magnitude
    faster than the previous skimage implementation.
    """
    import skimage.morphology
    import torch
    import torch.nn.functional as F

    if mask_file is None:
        mask_file = MASK_CCF25_FILE
    mask = load_image(mask_file)
    if mask.ndim == 4:
        mask = mask[0]
    
    # use cuda
    values = np.unique(mask)
    values = values[values > 0]
    mask = torch.from_numpy(mask.astype(np.int64)).cuda()
    print(f'Size of mask: {len(values)}')

    weights = torch.from_numpy(skimage.morphology.ball(radius=radius).astype(np.float32)).cuda().unsqueeze(0).unsqueeze(0)

    rn_dict = {}
    processed = 0
    cnt = 0
    for v in values:
        processed += 1
        mask_i = (mask == v).float().unsqueeze(0).unsqueeze(0)
        mask_i_dil = (F.conv3d(mask_i, weights, stride=1, padding=radius) > 0)[0,0]
        neighbors = torch.unique(mask[mask_i_dil]).cpu()
        # remove self
        neighbors = [ni.item() for ni in neighbors if (ni != v) and (ni != 0)]
        rn_dict[v] = neighbors
        cnt += 1
        if cnt % 5 == 0:
            print(f'[{processed}/{len(values)}]: {mask_i.sum().item()}->{mask_i_dil.sum().item()}-->{len(neighbors)}')

    return rn_dict

if __name__ == '__main__':
    import pickle
    
    radius = 5
    rn_dict = get_regional_neighbors_cuda(radius=radius)
    with open('./resources/regional_neighbors_res25_radius5.pkl', 'wb') as fp:
        pickle.dump(rn_dict, fp)

