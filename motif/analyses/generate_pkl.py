#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : generate_pkl.py
#   Author       : Yufeng Liu
#   Date         : 2023-05-11
#   Description  : 
#
#================================================================
import os, glob
import numpy as np
import pickle

import sys
sys.path.append('../../pylib')
from swc_handler import parse_swc

def load_and_dump(tract_dir='', key='stype'):
    tract_dict = {}
    n = 0
    for tract_file in glob.glob(os.path.join(tract_dir, '*tract.swc')):
        filename = os.path.split(tract_file)[-1]
        if filename.startswith('unk'):
            ptype = 'unk'
            stype = filename[4:-17].split('_')[0]
            prefix = filename[4+len(stype)+1:-17]
        else:
            ptype = filename.split('-')[0]
            stype = '-'.join(filename.split('-')[1:]).split('_')[0]
            prefix = '_'.join('-'.join(filename.split('-')[1:]).split('_')[1:])[:-17]

        tree = parse_swc(tract_file)
        coords = np.array([node[2:5] for node in tree])

        if key == 'stype':
            kname = stype
            fout = 'main_tract.pkl'
        elif key == 'ptype':
            kname = f'{ptype}-{stype}'
            fout = 'main_tract_ptype.pkl'
    
        if kname in tract_dict:
            tract_dict[kname].append((prefix, coords))
        else:
            tract_dict[kname] = [(prefix, coords)]

        n+= 1
        if n % 10 == 0:
            print(n, ptype, stype, prefix)
            #break

    with open(fout, 'wb') as fp:
        pickle.dump(tract_dict, fp)


if __name__ == '__main__':
    tract_dir = '../main_tracts_types'
    key = 'stype'
    load_and_dump(tract_dir, key=key)

