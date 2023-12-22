#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : match_gs_recon.py
#   Author       : Yufeng Liu
#   Date         : 2023-04-04
#   Description  : 
#
#================================================================

import os
import glob
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

import sys
sys.path.append('../../../../pylib')
from swc_handler import get_soma_from_swc

gs_dir = '/data/lyf/data/fullNeurons/V2023_12_11/profiled_final'
crop_dir = '/home/lyf/Research/publication/full-spectrum/micro_environ/benchmark/gs_crop_dendrite_tmd'
recon_dir = '/home/lyf/Research/publication/full-spectrum/micro_environ/benchmark_old/recon1891_weak1854_tmd'
outfile = 'file_mapping1854_tmd.csv'

def _make_gen(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024*1024)

def count_lines(filename):
    with open(filename, 'rb') as fp:
        f_gen = _make_gen(fp.raw.read)
        return sum(buf.count(b'\n') for buf in f_gen)

# get all gold standard soma position
gs = []
cnt = 0
for swcfile in glob.glob(os.path.join(gs_dir, '*swc')):
    fn = os.path.split(swcfile)[-1]
    prefix = os.path.splitext(fn)[0]
    crop_file = os.path.join(crop_dir, f'{prefix}.swc')
    if fn.startswith('pre'):
        brain_id = fn.split('_')[1]
    else:
        brain_id = fn.split('_')[0]
    soma = get_soma_from_swc(swcfile, match_str='.* -1 -1 .*\n')
    spos = list(map(float, soma[2:5]))
    nodes = count_lines(swcfile) - 1
    gs.append([prefix, brain_id, *spos, crop_file, nodes])
    
    cnt += 1
    if cnt % 100 == 0:
        print(cnt)
gs = pd.DataFrame(gs, columns=['filename', 'brain_id', 'xpos', 'ypos', 'zpos', 'path', 'nodes'])

# load all recon
res = []
scale = 2.
'''
for bdir in glob.glob(os.path.join(recon_dir, '[1-9]*')):
    brain_id = os.path.split(bdir)[-1]
    for swcfile in glob.glob(os.path.join(bdir, '*swc')):
        fn = os.path.split(swcfile)[-1]
        xyz = np.array(list(map(float, fn[:-4].split('_')))) * scale
        nodes = count_lines(swcfile) - 1
        res.append([fn, brain_id, *xyz, swcfile, nodes])
'''
for swcfile in glob.glob(os.path.join(recon_dir, '*swc')):
    fn = os.path.split(swcfile)[-1]
    brain_id, x, y, z = list(map(float, fn[:-4].split('_')))
    nodes = count_lines(swcfile) - 1
    brain_id = str(int(brain_id))
    x *= 2
    y *= 2
    z *= 2
    res.append([fn[:-4], brain_id, x, y, z, swcfile, nodes])

res = pd.DataFrame(res, columns=gs.columns)
 

dthr = 2*1.732
brains = np.unique(res.brain_id)
mindices = []
for brain in brains:
    print(brain)
    res_mask = res.brain_id == brain
    ridx = np.nonzero(res_mask.to_numpy())[0]
    res_i = res[res_mask]
    if brain in ['182712', '18467', '18469', '18866']:
        continue
    gs_mask = gs.brain_id == brain
    gidx = np.nonzero(gs_mask.to_numpy())[0]
    gs_i = gs[gs_mask]
    dm_i = distance_matrix(res_i[['xpos', 'ypos', 'zpos']], gs_i[['xpos', 'ypos', 'zpos']])
    print(dm_i.shape)
    
    imin = dm_i.argmin(axis=1)
    vmin = dm_i.min(axis=1)
    # double check the size
    dflag = vmin < dthr
    fidx = np.nonzero(dflag)[0]
    
    pairs = [[ri, gi] for (ri, gi) in zip(ridx, gidx[imin])]
    mindices.extend(pairs)

mindices = np.array(mindices)
mres = res.iloc[mindices[:,0]].reset_index(drop=True)
mgs = gs.iloc[mindices[:,1]].reset_index(drop=True)
merged = pd.merge(mres, mgs, left_index=True, right_index=True)
merged.to_csv(outfile, index=False)


