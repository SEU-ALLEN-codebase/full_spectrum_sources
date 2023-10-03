#!/usr/bin/env python

#================================================================
#   Copyright (C) 2022 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : run_arbor_r2.py
#   Author       : Yufeng Liu
#   Date         : 2022-10-14
#   Description  : 
#
#================================================================

import os
import glob
import pandas as pd
from multiprocessing.pool import Pool

import sys
sys.path.append('../../common_lib')
from common_utils import load_celltypes


def run_autoarbor(*args):
    swc_file, low, high = args
    os.system(f'python autoarbor_v1_yf.py --filename {swc_file} --L {low} --H {high}')


parser = argparse.ArgumentParser()
parser.add_argument('--r', help='round of autoarborization', type=str)
parser.add_argument('--L', help='Lower Bound', type=int)
parser.add_argument('--H', help='Higher Bound', type=int)
parser.add_argument('--swc_dir', help='the folder containing SWCs',
                    type=str, default='../data/axon80_sort')
args = parser.parse_args()


celltype_file = '../../common_lib/41586_2021_3941_MOESM4_ESM.csv'
soma_types, soma_types_r, p2stypes = load_celltypes(celltype_file, column_name='Soma_region')
low, high = args.L, args.H
r = args.r  # `r=1` is the first round of autoarbor, using L=${low} and H=${high}
            # `r=2` is the second round of autoarbor, using the consensus arbor number for each class


if r == 2:
    # `params_file` is the file containing the median arbor for each class
    params_file = 'log/axonal_arbor_params_round2.csv'
    pdict = {}
    params = pd.read_csv(params_file)
    for param in params.iterrows():
        print(param[1])
        ptype = param[1]['proj_type']
        stype = param[1]['soma_type']
        narbor = param[1]['median_arbor_num']
        pdict[(ptype, stype)] = narbor
    

args_list = []
for stype in soma_types:
    print(f'<------------------- soma type: {stype} --------------------->')
    prefixs = soma_types[stype]
    if r == 2:
        na = pdict[(ptype, stype)]
        low, high = na, na

    for prefix in prefixs:
        print(f'===> {prefix}')
        swc_file = os.path.join(swc_dir, f'{prefix}_axon.swc')
        if os.path.exists(f'{swc_file}._m3_l.eswc'):
            continue
        args_list.append((swc_file, na, na))
    
# do multiprocessing processing
nprocessors = 10
pt = Pool(nprocessors)
pt.starmap(run_autoarbor, args_list)
pt.close()
pt.join()
