#!/usr/bin/env python

#================================================================
#   Copyright (C) 2021 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : image_io.py
#   Author       : Yufeng Liu
#   Date         : 2021-05-17
#   Description  : 
#
#================================================================
import os
import glob
import SimpleITK as sitk
import pickle
from v3d.io import *
from pathlib import Path


def load_image(img_file, flip_tif=True):
    img_file = Path(img_file)
    if img_file.suffix in ['.v3draw', '.V3DRAW', '.raw', '.RAW']:
        return load_v3draw(img_file)
    if img_file.suffix in ['.v3dpbd', '.V3DPBD']:
        return PBD().load_image(img_file)
    img = sitk.GetArrayFromImage(sitk.ReadImage(str(img_file)))
    if flip_tif and img_file.suffix in ['.TIF', '.TIFF', '.tif', '.tiff']:
        img = np.flip(img, axis=-2)
    return img


def save_image(outfile, img: np.ndarray, flip_tif=True):
    outfile = Path(outfile)
    if outfile.suffix in ['.v3draw', '.V3DRAW']:
        save_v3draw(img, outfile)
    elif outfile.suffix in ['.TIF', '.TIFF', '.tif', '.tiff']:
        if flip_tif:
            img = np.flip(img, axis=-2)
        sitk.WriteImage(sitk.GetImageFromArray(img), str(outfile))
    else:
        sitk.WriteImage(sitk.GetImageFromArray(img), str(outfile))
    return True


def load_pickle(pkl_file):
    with open(pkl_file, 'rb') as fp:
        data = pickle.load(fp)
    return data


def save_pickle(obj, outfile):
    with open(outfile, 'wb') as fp:
        pickle.dump(obj, outfile)


def save_markers(outfile, markers, radius=0, shape=0, name='', comment='', c=(0,0,255)):
    with open(outfile, 'w') as fp:
        fp.write('##x,y,z,radius,shape,name,comment, color_r,color_g,color_b\n')
        for marker in markers:
            x, y, z = marker
            fp.write(f'{x:3f}, {y:.3f}, {z:.3f}, {radius},{shape}, {name}, {comment},0,0,255\n')

def get_tera_res_path(tera_dir, res_ids=None, bracket_escape=True):
    '''
    res_ids: int or tuple
    - if int: it represents the resolution level, starting from the lowest resolution. The value -1
        means the highest resolution
    - if tuple: multiple levels 
    '''
    
    resfiles = list(glob.glob(os.path.join(tera_dir, 'RES*')))
    if len(resfiles) == 0:
        print(f'Error: the brain {os.path.split(tera_dir)[-1]} is not found!')
        return 
    
    # sort by resolutions
    ress = sorted(resfiles, key=lambda x: int(x.split('/')[-1][4:-1].split('x')[0]))
    if type(res_ids) is int:
        res_path = ress[res_ids]
        if bracket_escape:
            res_path = res_path.replace('(', '\(')
            res_path = res_path.replace(')', '\)')
        return res_path
    elif type(res_ids) is tuple:
        res_pathes = []
        for idx in res_ids:
            res_path = ress[idx]
            if bracket_escape:
                res_path = res_path.replace('(', '\(')
                res_path = res_path.replace(')', '\)')
            res_pathes.append(res_path)
        return res_pathes

        

