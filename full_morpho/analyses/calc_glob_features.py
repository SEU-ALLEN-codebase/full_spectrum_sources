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

sys.path.append('../../microenviron/generation')
from config import __FEAT_NAMES22__

def calc_global_features(swc_file, vaa3d='/opt/Vaa3D_x.1.1.4_ubuntu/Vaa3D-x'):
    cmd_str = f'xvfb-run -a -s "-screen 0 640x480x16" {vaa3d} -x global_neuron_feature -f compute_feature -i {swc_file}'
    p = subprocess.check_output(cmd_str, shell=True)
    output = p.decode().splitlines()[37:-2] # Vaa3D-x [37:-2], Vaa3D: [32:-2]
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
    t0 = time.time()
    iswc = 0
    for swcfile in glob.glob(os.path.join(swc_dir, '*swc')):
        prefix = os.path.split(swcfile)[-1][:-4]

        features = calc_global_features(swcfile)
        features_all.append([prefix, *features])
        
        iswc += 1
        if iswc % 5 ==  0:
            print(f'--> {iswc} in {time.time() - t0:.2f} s')

    df = pd.DataFrame(features_all, columns=['Cell name',  *__FEAT_NAMES22__])
    df.to_csv(outfile, float_format='%g', index=False)

if __name__ == '__main__':
    swc_dir = '/data/lyf/data/fullNeurons/V2023_12_11/S3_1um_final'
    outfile = '../data/feature.txt'
    
    calc_global_features_all(swc_dir, outfile)

    
    
    
