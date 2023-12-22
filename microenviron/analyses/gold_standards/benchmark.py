#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : benchmark.py
#   Author       : Yufeng Liu
#   Date         : 2023-04-04
#   Description  : 
#
#================================================================
import os, glob
import subprocess
import math
import random
import numpy as np
import pandas as pd
import tmd
import cv2
import matplotlib.pyplot as plt

import sys
sys.path.append('../../../pylib')
from swc_handler import parse_swc, get_specific_neurite, NEURITE_TYPES, flip_swc
from neuron_quality.metrics import DistanceEvaluation
from morph_topo.morphology import Morphology

np.set_printoptions(precision=4)

def evaluate(match_file, dsa_thr=2., outfile='temp.csv', only_dendrite=False, flip_y=True):
    de = DistanceEvaluation(dsa_thr=dsa_thr, resample2=False)
    df = pd.read_csv(match_file)

    dmatrix = []
    for idx, row in df.iterrows():
        swc1 = row['path_x']
        swc2 = row['path_y']
        tree1 = parse_swc(swc1)
        if flip_y:
            tree1 = flip_swc(tree1, axis='y', dim=512)
        tree2 = parse_swc(swc2)

        nodes1 = len(tree1)
        if only_dendrite:
            dendrites = get_specific_neurite(tree2, NEURITE_TYPES['dendrite'])
            nodes2 = len(dendrites)
            ds = de.run(tree1, dendrites)
        else:
            nodes2 = len(tree2)
            ds = de.run(tree1, tree2)
        # path length
        morph1 = Morphology(tree1)
        morph2 = Morphology(tree2)
        pl1 = morph1.calc_total_length()
        if only_dendrite:
            seg_lengths, lengths_dict = morph2.calc_frag_lengths()
            is_dendrite = np.array([node[1] in [3,4] for node in tree2])
            pl2 = seg_lengths[is_dendrite].sum()
        else:
            pl2 = morph2.calc_total_length()
        
        if idx % 10 == 0:
            print(f'[{idx}] #nodes1={nodes1}, nodes2={nodes2}, metrics=\n{ds}\n')

        dmatrix.append([*ds[2], nodes1, nodes2, pl1, pl2])
    dmatrix = pd.DataFrame(dmatrix, columns=['pds12', 'pds21', 'pds', 'nodes1', 'nodes2', 'path_length1', 'path_length2'])
    dmatrix.to_csv(outfile, index=False)

    return dmatrix

def sort_swc(swc_in, swc_out=None, vaa3d='/opt/Vaa3D_x.1.1.4_ubuntu/Vaa3D-x'):
    cmd_str = f'xvfb-run -a -s "-screen 0 640x480x16" {vaa3d} -x sort_neuron_swc -f sort_swc -i {swc_in} -o {swc_out}'
    #p = subprocess.check_output(cmd_str, shell=True)

    # retype
    df = pd.read_csv(swc_out, sep=' ', names=('#id', 'type', 'x', 'y', 'z', 'r', 'p'), comment='#', index_col=False)
    df['type'] = 3
    df.loc[0, 'type'] = 1
    df.to_csv(swc_out, sep=' ', index=False)

    return True

def sort_swcs(in_dir, out_dir):
    args_list = []
    for in_swc in glob.glob(os.path.join(in_dir, '*swc')):
        in_name = os.path.split(in_swc)[-1]
        print(in_name)
        out_swc = os.path.join(out_dir, in_name)

        args_list.append((in_swc, out_swc))

    # multiprocessing
    from multiprocessing import Pool
    pool = Pool(processes=24)
    pool.starmap(sort_swc, args_list)
    pool.close()
    pool.join()


def benchmark_with_tmd(match_file, ctype_file, ref_region='CP'):

    def get_ph(swc1, swc2):
        if swc1 == swc2:
            return None, None
        
        neuron1 = tmd.io.load_neuron(swc1)
        neuron2 = tmd.io.load_neuron(swc2)
        try:
            ph1 = tmd.methods.get_ph_neuron(neuron1)
            ph2 = tmd.methods.get_ph_neuron(neuron2, neurite_type="basal_dendrite")
            return ph1, ph2
        except:
            print(idx)
            return None, None


    df = pd.read_csv(match_file)
    mnames, mcnts = np.unique(df.filename_y, return_counts=True)
    mnames = mnames[mcnts == 1]
    df = df[df.filename_y.isin(mnames)]
    df_m = df.copy()

    # find all neurons from the reference region
    df_ct = pd.read_csv(ctype_file, index_col=0)
    if ref_region == 'CP':
        neurons = df_ct[df_ct.Manually_corrected_soma_region == ref_region]['Cell name']
    elif ref_region == 'CP_GPe':
        neurons = df_ct[df_ct.Subclass_or_type == ref_region]['Cell name']
    else:
        raise ValueError
    

    phs1, phs2 = [], []
    names = []
    for idx, row in df.iterrows():
        swc1 = row['path_x']
        swc2 = row['path_y']
        print(idx, swc1, swc2)
        ph1, ph2 = get_ph(swc1, swc2)
        if ph1 is None:
            continue

        phs1.append(ph1)
        phs2.append(ph2)
        names.append(row.filename_x)

        #if idx >= 20:   # for debug
        #    break

    # we should estimate TMD distance for neurons from the same region as a reference
    # reconstructed neurons in reference region
    ref_neurons = df_m[df_m.filename_y.isin(neurons)].path_y.to_numpy().tolist()
    #print(df_m[df_m.filename_y.isin(neurons)].filename_y); sys.exit()
    
    min_pairs = 1000
    nrep = int(math.ceil(min_pairs * 1.0 / len(ref_neurons)))
    
    random.seed(2023)   # for reproducibility
    ns1 = random.sample(ref_neurons*nrep, min_pairs)
    ns2 = random.sample(ref_neurons*nrep, min_pairs)
    ref_phs1, ref_phs2 = [], []
    icnt = 0
    for nname1, nname2 in zip(ns1, ns2):
        ph1, ph2 = get_ph(nname1, nname2)
        if ph1 is None:
            continue
        ref_phs1.append(ph1)
        ref_phs2.append(ph2)

        icnt += 1
        if icnt % 20 == 0:
            print(icnt)


    # Normalize the limits
    xlim, ylim = tmd.analysis.get_limits(phs1 + phs2 + ref_phs1 + ref_phs2)
    print(f'xlim and ylim are {xlim}, {ylim}')

    # estimate the TMD distance for the reference pairs
    ref_dists = []
    for rph1, rph2 in zip(ref_phs1, ref_phs2):
        try:
            img1 = tmd.analysis.get_persistence_image_data(rph1, xlim=xlim, ylim=ylim)
            img2 = tmd.analysis.get_persistence_image_data(rph2, xlim=xlim, ylim=ylim)
        except:
            continue

        dimg = tmd.analysis.get_image_diff_data(img1, img2)
        dist = np.sum(np.abs(dimg))
        ref_dists.append(dist)
    ref_dists = np.array(ref_dists)
    print(f'Number of reference pairs: {ref_dists.shape}')
    print(f'Statis for reference TMD pairs: max={ref_dists.max():.2f}, min={ref_dists.min():.2f}, mean={ref_dists.mean():.2f}, std={ref_dists.std():.2f}')


    # Get the color map by name:
    cm = plt.get_cmap('gist_rainbow')
    df_tmd = []
    for i, name, ph1, ph2 in zip(range(len(names)), names, phs1, phs2):
        try:
            img1 = tmd.analysis.get_persistence_image_data(ph1, xlim=xlim, ylim=ylim)
            img2 = tmd.analysis.get_persistence_image_data(ph2, xlim=xlim, ylim=ylim)
        except:
            print(i, name)
            continue

        #img = np.hstack((img1, img2))
        # Apply the colormap like a function to any array:
        #cimg = (cm(img)[:,:,:3] * 255).astype(np.uint8)
        #cv2.imwrite(f'TMD_{name}.png', cimg)
        
        dimg = tmd.analysis.get_image_diff_data(img1, img2)
        dist = np.sum(np.abs(dimg))
        df_tmd.append([name, dist])
    df_tmd = pd.DataFrame(df_tmd, columns=['filename', 'TMD dist'])
    df_tmd.to_csv('tmd.csv')


if __name__ == '__main__':
    # the `match_file` contains the correspondence between the gold standards and reconstructed morphologies
    match_file = './utils/file_mapping1854.csv'
    only_dendrite = True    # Using only the dendrites of the gold standards


    if 0:
        # sort all swc to meet the strict requirement for TMD estimation
        in_dir = '../../benchmark_old/recon1891_weak1854'
        out_dir = '../../benchmark_old/recon1891_weak1854_tmd'
        #sort_swc_for_tmd(in_dir, out_dir)
        sort_swcs(in_dir, out_dir)

    if 1:
        # Part 3
        match_file = './utils/file_mapping1854_tmd.csv'
        ctype_file = '../../../common_lib/41586_2021_3941_MOESM4_ESM.csv'
        ref_region = 'CP_GPe'
        benchmark_with_tmd(match_file, ctype_file, ref_region=ref_region)
    

