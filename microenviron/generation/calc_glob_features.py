#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : calc_glob_features.py
#   Author       : Yufeng Liu
#   Date         : 2023-01-24
#   Description  : 
#
#================================================================

import os, glob
import sys
import time
import subprocess
import pandas as pd
import numpy as np

sys.path.append('../src')
from config import __FEAT_NAMES22__

def calc_global_features(swc_file, vaa3d='D:/Vaa3D/Vaa3D-x.1.1.2_Windows_64bit/Vaa3D-x.exe'):
    cmd_str = f'{vaa3d} /x global_neuron_feature /f compute_feature /i {swc_file}'
    p = subprocess.check_output(cmd_str, shell=True)
    output = p.decode().splitlines()[15:-1]
    info_dict = {}
    for s in output:
        it1, it2 = s.split(':')
        it1 = it1.strip()
        it2 = it2.strip()
        info_dict[it1] = float(it2)

    # extract the target result
    #print(info_dict)
    features = []
    features.append(int(info_dict['N_node']))
    features.append(info_dict['Soma_surface'])
    features.append(int(info_dict['N_stem']))
    features.append(int(info_dict['Number of Bifurcatons']))
    features.append(int(info_dict['Number of Branches']))
    features.append(int(info_dict['Number of Tips']))
    features.append(info_dict['Overall Width'])
    features.append(info_dict['Overall Height'])
    features.append(info_dict['Overall Depth'])
    features.append(info_dict['Average Diameter'])
    features.append(info_dict['Total Length'])
    features.append(info_dict['Total Surface'])
    features.append(info_dict['Total Volume'])
    features.append(info_dict['Max Euclidean Distance'])
    features.append(info_dict['Max Path Distance'])
    features.append(info_dict['Max Branch Order'])
    features.append(info_dict['Average Contraction'])
    features.append(info_dict['Average Fragmentation'])
    features.append(info_dict['Average Parent-daughter Ratio'])
    features.append(info_dict['Average Bifurcation Angle Local'])
    features.append(info_dict['Average Bifurcation Angle Remote'])
    features.append(info_dict['Hausdorff Dimension'])

    return features

def calc_global_features_all(swc_dir, outfile):

    features_all = []
    iswc = 0
    t0 = time.time()
    for swcfile in glob.glob(os.path.join(swc_dir, '*swc')):
        print(swcfile)
        #prefix = os.path.splitext(os.path.split(swcfile)[-1])[0]
        prefix = os.path.split(swcfile)[-1]
        features = calc_global_features(swcfile)
        features_all.append([prefix, *features])
        
        iswc += 1
        if iswc % 10 ==  0:
            print(f'--> {iswc} in {time.time() - t0:.2f} s')

    df = pd.DataFrame(features_all, columns=['',  *__FEAT_NAMES22__])
    df.to_csv(outfile, float_format='%g', index=False)

if __name__ == '__main__':
    #swc_dir = 'D:\\mygitcode\\full_spectrum_sources\\microenviron\\generation\\data\\improved_reg\\42k_local_morphology_new20230510_gcoord_final'
    swc_dir = '.\\data\\improved_reg\\42k_local_morphology_new20230510_gcoord_final'
    out_dir  = './lm_features'

    
    args_list = []
    for brain_dir in glob.glob(os.path.join(swc_dir, '*/*')):
        pdset, brain_id = os.path.split(brain_dir)
        dset = os.path.split(pdset)[-1]
        outfile = f'lm_d22_{dset}_{brain_id}.csv'
        outfile = os.path.join(out_dir, outfile)
        print(dset, outfile)
        if os.path.exists(outfile):
            continue
        args_list.append((brain_dir, outfile))
        #calc_global_features_all(brain_dir, outfile)


    from multiprocessing import Pool
    pool = Pool(processes=24)
    pool.starmap(calc_global_features_all, args_list)
    pool.close()
    pool.join()
    '''
    for args in args_list:
        calc_global_features_all(*args)
        break
    '''
