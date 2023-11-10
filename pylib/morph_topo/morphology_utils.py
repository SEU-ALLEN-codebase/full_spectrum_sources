#!/usr/bin/env python

#================================================================
#   Copyright (C) 2022 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : morphology_utils.py
#   Author       : Yufeng Liu
#   Date         : 2022-04-18
#   Description  : 
#
#================================================================

import numpy as np

def get_outside_soma_mask(morph, dist_thresh):
    sc = np.array(morph.pos_dict[morph.idx_soma][2:5])
    coords = np.array([morph.pos_dict[node[0]][2:5] for node in morph.tree])
    indices = [node[0] for node in morph.tree]
        
    
    vec = coords - sc #self.spacing
    dists = np.linalg.norm(vec, axis=1)

    out_mask = dists > dist_thresh
    out_dict = {}
    for idx, m in zip(indices, out_mask):
        out_dict[idx] = m

    return out_dict

