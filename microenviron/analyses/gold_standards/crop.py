#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : crop.py
#   Author       : Yufeng Liu
#   Date         : 2023-01-24
#   Description  : 
#
#================================================================
import sys
sys.path.append(r'D:\code\support_lib\pylib')
import os
import glob
import sys
import numpy as np
from swc_handler import parse_swc, write_swc, shift_swc, find_soma_index, trim_swc, get_specific_neurite, NEURITE_TYPES

def crop_swc(swcfile, outfile, crop_size=(308,308,512), only_dendrite=False):
    tree = parse_swc(swcfile)
    soma_index = find_soma_index(tree, p_soma=-1)
    sx, sy, sz = tree[soma_index][2:5]
    fx = sx - crop_size[0]/2
    fy = sy - crop_size[1]/2
    fz = sz - crop_size[2]/2

    # move the tree to center_crop
    tree = shift_swc(tree, fx, fy, fz)
    if only_dendrite:
        tree = get_specific_neurite(NEURITE_TYPES['dendrite'])
    tree = trim_swc(tree, crop_size[::-1], bfs=True)
    tree = shift_swc(tree, -fx, -fy, -fz)

    write_swc(tree, outfile)

def crop_all(swc_dir, out_dir, crop_size=(308,308,512)):
    i = 0
    for swcfile in glob.glob(os.path.join(swc_dir, '*.swc')):
        swcname = os.path.split(swcfile)[-1]
        print(f'[{i}]: {swcname}')
        outfile = os.path.join(out_dir, swcname)
        crop_swc(swcfile, outfile, crop_size)

        i += 1
    

if __name__ == '__main__':
    swc_dir = '../data/registration/S3_registered_ccf'
    out_dir = '../crop_dendrite'
    crop_all(swc_dir, out_dir)

