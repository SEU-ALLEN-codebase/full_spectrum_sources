#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : compare_CTX_SSp-m.py
#   Author       : Yufeng Liu
#   Date         : 2023-03-25
#   Description  : 
#
#================================================================

import os
import glob
import numpy as np

def get_the_subtypes(subtype_dir='../CTX_ET-SSp-m-subclasses'):
    ctxfiles = sorted(glob.glob(os.path.join(subtype_dir, '*swc')))
    ctx_classes = {}
    for cf in ctxfiles:
        fn = os.path.split(cf)[-1]
        cls = int(fn[0])
        prefix = fn[15:-17]
        if cls in ctx_classes:
            ctx_classes[cls].append(prefix)
        else:
            ctx_classes[cls] = [prefix]
    return ctx_classes

def montage_image(mip_dir, files, sw, sh, prefix):
    swh = sw * sh
    for i in range(0, len(files), swh):
        subset = [os.path.join(mip_dir, f'{prefix}.png') for prefix in files[i:i+swh]]
        args_str = f'montage -mode concatenate {" ".join(subset)} -tile {sw}x{sh} montage_{prefix}_{i:04d}.png'
        os.system(args_str)

if __name__ == '__main__':
    mip_dir = '../../neurite_arbors/axon_on_ccf'
    
    ctx_classes = get_the_subtypes()
    for cls, files in ctx_classes.items():
        print(cls, len(files))
        montage_image(mip_dir, files, 1, 4, cls)
    

