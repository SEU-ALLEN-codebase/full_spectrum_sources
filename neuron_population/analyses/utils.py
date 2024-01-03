#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : utils.py
#   Author       : Yufeng Liu
#   Date         : 2023-03-02
#   Description  : 
#
#================================================================
import os,glob
import json
import pickle
import copy
import math
import numpy as np
import pandas as pd
import seaborn as sns
from swc_handler import parse_swc

from file_io import load_image

def write_json(write_file,des_file):
    dict_object = json.dumps(write_file)
    with open(des_file, "w") as f:
        f.write(dict_object)
        
def load_json(des_file):
    with open(des_file, "r") as f:
        load_file = json.loads(f.read())
    return load_file

def load_json_key2int(dict_file):
    with open(dict_file, "r") as f:
        load_dict = json.loads(f.read())
    new_dict = {}
    for k,v in load_dict.items():
        new_dict[int(k)] = v
    return new_dict

def load_json_keys2int(dict_file):
    with open(dict_file, "r") as f:
        load_dict = json.loads(f.read())
    new_dict = {}
    for k,v in load_dict.items():
        new_dict[k] = {}
        for kk,vv in v.items():
            new_dict[k][int(kk)] = vv
    return new_dict

def get_tree():
    with open('../data/Fig2/assets/tree_yzx.json','r') as f:
        tree = json.loads(f.read())
    return tree

def get_tree_from(info,hold=[]):
    info_tree = {}
    tree = get_tree()
    for t in tree:
        if (not hold) or (t[info] in hold):
            info_tree[t[info]] = t  
    return info_tree

def load_level_u32csv():
    level_file = '../data/Fig2/assets/level.csv'
    level_df = pd.read_csv(level_file)
    return level_df

def check_region_u32list(u32list,n):
    if n not in [671,316]: raise
    u32_list = np.array(u32list)    
    u32_list[u32_list==997]=0
    u32_list[u32_list==1024]=0
    u32_list[u32_list==1009]=0
    u32_list[u32_list==567]=0
    u32_list[u32_list==73]=0
    u32_list[u32_list==304325711]=0
    level_df = load_level_u32csv()
    new_list = []
    for i,u32 in enumerate(u32_list):
        if u32==0: new_list.append(0); continue
        n4 = level_df[level_df['n1327_u32_id']==u32]['n4_u32_id'].values[0]
        n13 = level_df[level_df['n1327_u32_id']==u32]['n13_u32_id'].values[0]
        if 0 in [n4,n13]: new_list.append(0); continue
        new_list.append(int(u32))
    return new_list

def get_n1327_n4_u32dict():
    level_df = load_level_u32csv()
    n1327_n5_u32dict = dict(zip(level_df['n1327_u32_id'].tolist(),level_df['n4_u32_id'].tolist()))
    return n1327_n5_u32dict                          

def get_u32_u16_id_dict():
    u32_u16_id_dict_file = f'../data/Fig2/assets/u32_u16_dict.json'
    u32_u16_id_dict = load_json_key2int(u32_u16_id_dict_file)
    return u32_u16_id_dict

def get_u16_u32_id_dict():
    u16_u32_id_dict_file = '../data/Fig2/assets/u16_u32_dict.json'
    u16_u32_id_dict = load_json_key2int(u16_u32_id_dict_file)
    return u16_u32_id_dict

def get_u32_from_u16_list(u16id_list):
    u32id_list = []
    u16_u32_id_dict = get_u16_u32_id_dict()
    for u16id in u16id_list:
        u32id = u16_u32_id_dict[u16id]
        u32id_list.append(u32id)
    return u32id_list

def get_u16_from_u32_list(u32id_list,left_list):
    u16id_list = []
    u32_u16_id_dict = get_u32_u16_id_dict()
    for u32id,left in zip(u32id_list,left_list):
        if int(u32id)==0: u16id = 0
        elif left: u16id = int(u32_u16_id_dict[u32id]['l'])
        else: u16id = int(u32_u16_id_dict[u32id]['r'])
        u16id_list.append(u16id)
    return u16id_list

def get_level_id_list(id_list,level):
    level_id_list = []
    level_df = load_level_u32csv()
    for n1327id in id_list:
        if n1327id==0:
            level_id_list.append(0)
            continue
        nid = level_df[level_df['n1327_u32_id']==n1327id][f'n{level}_u32_id'].values[0]
        level_id_list.append(nid)
    return level_id_list   
    
def get_voxel_center_neighbor_dict(level,uint):
    dict_file = f'../data/Fig2/assets/n{level}_u{uint}_voxel_center_neighbor_dict.json'
    load_dict = load_json_keys2int(dict_file)
    return load_dict
#nrrd u25 voxel, nrrd u25 center, nrrd neighbor 
def get_dict(level,uint,key=['voxel','center','neighbor'][0]):
    voxel_center_neighbor_dict = get_voxel_center_neighbor_dict(level,uint)
    get_dict = voxel_center_neighbor_dict[f'n{level}_u{uint}_{key}_dict']
    return get_dict
  
def get_u25voxel_dict_from_nrrd(level,uint,imgfile='',allow0=False):    
    if imgfile=='': 
        imgfile = f'../data/Fig2/assets/n{level}_u{uint}.nrrd'
    dtype = np.uint32 if uint==32 else np.uint16
    image = load_image(imgfile).astype(dtype)
    vs,cs = np.unique(image,return_counts=True)
    voxel_dict = {int(v):int(c) for v,c in zip(vs,cs) if (v!=0 and ~allow0)}
    return voxel_dict
  
def get_region_name_list(region_unit_id_list,uint):
    u32_tree = get_tree_from('id')
    u16_u32_id_dict = get_u16_u32_id_dict()
    region_name_list = []
    for region_unit_id in region_unit_id_list:
        if uint==32:
            try:
                region_name = u32_tree[region_unit_id]['acronym']
            except:
                region_name = ''
        else:
            try:
                ends = ['_l','_r'][region_unit_id%2]
                region_u32_id = u16_u32_id_dict[region_unit_id]
                region_name = u32_tree[region_u32_id]['acronym'] + ends 
            except:
                region_name = ''
        region_name_list.append(region_name)
    return region_name_list

def get_region_u32id_list(region_name_list):
    acronym_tree = get_tree_from('acronym')
    region_u32id_list = []
    for region_name in region_name_list:
        try:
            region_u32id = int(acronym_tree[region_name]['id'])
        except:
            region_u32id = 0
        region_u32id_list.append(region_u32id)
    return region_u32id_list

def get_region_rgb_list(region_uint_list,uint=32):
    region_rgb_list = []
    u32_id_tree = get_tree_from('id')
    u16_u32_id_dict = get_u16_u32_id_dict()
    for region_uint in region_uint_list:
        if uint==16:
            region_uint_id = u16_u32_id_dict[region_unit_id]
        try:
            region_rgb_list.append(u32_id_tree[region_uint]['rgb_triplet'])
        except:
            region_rgb_list.append([255,255,255])
    return region_rgb_list

def color_map(map_dict,colors):
    values = np.unique(list(map_dict.values())).tolist()
    colors = (colors*10)[:len(values)]
    value_color_dict = dict(zip(values,colors))
    new_dict = dict([(k,value_color_dict[v]) for k,v in map_dict.items()])
    return new_dict,value_color_dict

def get_colorcode16_list(value_list,nc=10,cmap='spring',reverse=True):        
    rgblist = (np.array(sns.color_palette(cmap,nc+1))*255).astype(int)
    if reverse: 
        rgblist = rgblist[nc::-1,:]
    code16_list = ["".join([hex(c)[2:].rjust(2,'0').upper() for c in rgb]) for rgb in rgblist]
    value_code16_list = [code16_list[v] for v in value_list]
    return value_code16_list
