#!/usr/bin/env python

#================================================================
#   Copyright (C) 2021 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : morphology_angles.py
#   Author       : Yufeng Liu
#   Date         : 2021-07-19
#   Description  : 
#
#================================================================

import numpy as np

from math_utils import calc_included_angles_from_coords
from swc_handler import NEURITE_TYPES

class MorphAngles(object):
    def __init__(self):
        pass

    def calc_outgrowth_angles(self, morph, spacing=None, indices_set=None, epsilon=1e-7):
        """
        calculate the outgrowth angles for given morphology object.
        Arguments are:
        - morph: Morphology object
        - spacing: resolution in micrometer, in (x,y,z) order
        - indices_set: set of indices that want to calculate
        """
        if indices_set is None:
            indices_set = set(morph.pos_dict.keys())

        angles = []
        cur_coords = []
        par_coords = []
        chi_coords = []
        for idx in indices_set:
            p_idx = morph.pos_dict[idx][6]
            if p_idx not in morph.pos_dict:
                continue
            if idx not in morph.child_dict:
                continue

            cur_coord = morph.pos_dict[idx][2:5]
            par_coord = morph.pos_dict[p_idx][2:5]
            for c_idx in morph.child_dict[idx]:
                chi_coord = morph.pos_dict[c_idx][2:5]
                cur_coords.append(cur_coord) 
                par_coords.append(par_coord)
                chi_coords.append(chi_coord)
        
        angs = calc_included_angles_from_coords(cur_coords, par_coords, chi_coords, spacing=spacing)
        return angs

class MorphCurvature(object):
    def __init__(self, morph, neurite_type='all', spacing=None):
        self.morph = morph
        self.neurite_type = neurite_type
        if spacing is not None:
            self.spacing = np.array(spacing)
        else:
            self.spacing = None
        # get all paths
        self.paths = self.morph.get_all_paths()

    def estimate_coplanarity(self, discard_multifurcate=True, ignore_thresh=0.05):
        npt = 4 # number of points to test copalanarity
        
        cop_dict = {}
        multifurcation = self.morph.multifurcation | self.morph.bifurcation
        for tip, path in self.paths.items():
            coords = np.array([self.morph.pos_dict[idx][2:5] for idx in path])
            mflags = [idx in multifurcation for idx in path]
            if self.neurite_type != 'all':
                nflags = [self.morph.pos_dict[idx][1] in NEURITE_TYPES[self.neurite_type] for idx in path]
            # path with insufficient nodes
            if len(path) < npt: continue

            # primary direction
            pd1 = coords[1:] - coords[:-1]
            if self.spacing is not None:
                pd1 *= self.spacing.reshape(1, -1)
            
            # normalize
            pd1 = pd1 / (np.linalg.norm(pd1, axis=1).reshape(-1,1) + 1e-7)

            for ii in range(len(path)-npt-1):
                # check multifurcate
                if discard_multifurcate:
                    mflag = sum(mflags[ii:ii+npt])
                    if mflag:
                        continue
                else:
                    mflag = sum(mflags[ii:ii+npt])
                    if not mflag:
                        continue

                # check neurite types
                if self.neurite_type != 'all':
                    nflag = sum(nflags[ii:ii+npt])
                    if nflag:
                        continue
                        
                cur_pd = pd1[ii:ii+npt-1]
                '''# get the normal vector
                nv = np.cross(cur_pd[0], cur_pd[1])
                v = cur_pd[2]
                nvl = np.linalg.norm(nv)
                if nvl < ignore_thresh:
                    # try another pair
                    nv = np.cross(cur_pd[1], cur_pd[2])
                    nvl = np.linalg.norm(nv)
                    if nvl < ignore_thresh:
                        continue
                    else:
                        v = cur_pd[0]

                # check angle
                cos_ang = v.dot(nv) / nvl
                if cos_ang < 0:
                    cos_ang = -cos_ang
                ang = np.pi/2. - np.arccos(cos_ang)
                '''
                nv1 = np.cross(cur_pd[0], cur_pd[1])
                nv2 = np.cross(cur_pd[1], cur_pd[2])
                nvl1, nvl2 = np.linalg.norm(nv1), np.linalg.norm(nv2)
                if nvl1 < ignore_thresh or nvl2 < ignore_thresh:
                    continue
                cos_ang = nv1.dot(nv2) / nvl1 / nvl2
                ang = np.arccos(cos_ang)
                ang = min(ang, np.pi - ang)

                cop_dict[path[ii+1]] = ang
        return cop_dict

    def estimate_angular_dependence(self, discard_multifurcate=True):
        npt = 4 # at least four points to estimate the dependence
        angd_dict = {}
        multifurcation = self.morph.multifurcation | self.morph.bifurcation
        for tip, path in self.paths.items():
            coords = np.array([self.morph.pos_dict[idx][2:5] for idx in path])
            if discard_multifurcate:
                mflags = [idx in multifurcation for idx in path]
            if self.neurite_type != 'all':
                nflags = [self.morph.pos_dict[idx][1] in NEURITE_TYPES[self.neurite_type] for idx in path]
            # path with insufficient nodes
            if len(path) < npt: continue

            # primary direction
            pd1 = coords[1:] - coords[:-1]
            # normalize
            pd1 = pd1 / (np.linalg.norm(pd1, axis=1).reshape(-1,1) + 1e-7)

            for ii in range(len(path)-npt-1):
                # check multifurcate
                if discard_multifurcate:
                    mflag = sum(mflags[ii:ii+npt])
                    if mflag:
                        continue
                # check neurite types
                if self.neurite_type != 'all':
                    nflag = sum(nflags[ii:ii+npt])
                    if nflag:
                        continue

                cur_pd = pd1[ii:ii+npt-1]
                cos_ang1 = (-cur_pd[0]).dot(cur_pd[1])
                cos_ang2 = (-cur_pd[1]).dot(cur_pd[2])
                ang1 = np.arccos(cos_ang1)
                ang2 = np.arccos(cos_ang2)
                angd_dict[path[ii+1]] = (ang1, ang2)
        return angd_dict
                
        

if __name__ == '__main__':
    import time
    import os, glob
    import pickle
    import matplotlib.pyplot as plt

    from morphology import Morphology
    from swc_handler import parse_swc, load_spacings

    swc_dir = '/media/lyf/storage/seu_mouse/swc/xy1z1'
    spacing_file = '/media/lyf/storage/seu_mouse/swc/AllbrainResolutionInfo.csv'
    
    # load the spacings
    spacing_dict = load_spacings(spacing_file)

    
    '''
    ########## angular distribution for swc file ##########
    neurite_type = 'dendrite'
    npd = 0
    all_angles = []
    t0 = time.time()
    np.random.seed(2048)
    for swcfile in glob.glob(os.path.join(swc_dir, '*.swc')):
        #if np.random.random() > 0.8: continue   # subset test

        swcname = os.path.split(swcfile)[-1]
        brain_id = int(swcname.split('_')[0])
        try:
            spacing = spacing_dict[brain_id]
        except KeyError:
            #spacing = (0.23,0.23,1.0)
            continue
        
        tree = parse_swc(swcfile)
        morph = Morphology(tree)
        ma = MorphAngles()
        if neurite_type == 'all':
            filter_set = None
        else:
            filter_set = morph.get_nodes_by_types(neurite_type)
        angles = ma.calc_outgrowth_angles(morph, spacing=spacing, indices_set=filter_set)
        all_angles.extend(angles.tolist())
        #print(f'Number of angles: {len(angles)}')

        npd += 1
        if npd % 10 == 0:
            print(f'--> {len(all_angles)} angles found in {npd} files. Calculation takes {time.time() - t0:.2f} seconds')
    
    all_angles = np.array(all_angles)
    with open(f'{neurite_type}_angles.pkl', 'wb') as fp:
        pickle.dump(all_angles, fp)
    
    plt.hist(all_angles, bins=100)
    plt.savefig(f'{neurite_type}.png', dpi=300)
    #########################################################
    '''
    


    
    ############### coplanarity estimation ##################
    discard_multifurcate = True
    ignore_normal_thresh = 0.1
    neurite_type = 'dendrite'
    total_process_file = 50
    isotropic = False
    
    t0 = time.time()
    all_angles = []
    npd = 0
    for swcfile in glob.glob(os.path.join(swc_dir, '*.swc')):
        swcname = os.path.split(swcfile)[-1]
        brain_id = int(swcname.split('_')[0])
        if isotropic:
            try:
                spacing = spacing_dict[brain_id]
            except KeyError:
                #spacing = (0.23,0.23,1.0)
                continue
        else:
            spacing = None

        tree = parse_swc(swcfile)
        morph = Morphology(tree)
        mc = MorphCurvature(morph, neurite_type=neurite_type, spacing=spacing)
        cop_dict = mc.estimate_coplanarity(
                        discard_multifurcate=discard_multifurcate, 
                        ignore_thresh=ignore_normal_thresh)
        all_angles.extend(list(cop_dict.values()))
        
        npd += 1
        if npd % 2 == 0:
            print(f'--> {len(all_angles)} angles found in {npd} files. Calculation takes {time.time() - t0:.2f} seconds')
        if npd >= total_process_file:   # use subset for statistics
            break

    all_angles = np.array(all_angles)
    with open(f'coplanarity_{neurite_type}_swc{total_process_file}.pkl', 'wb') as fp:
        pickle.dump(all_angles, fp)
    print(f'#angles: {len(all_angles)}')
    plt.hist(all_angles, bins=100, range=(0, 1.57), density=True)
    plt.savefig(f'coplanarity_{neurite_type}_swc{total_process_file}.png', dpi=300)
    


    '''    
    ################# angular dependence #####################
    discard_multifurcate = True
    neurite_type = 'dendrite'
    total_process_file = 50
    
    
    t0 = time.time()
    angle_pairs = []
    npd = 0
    for swcfile in glob.glob(os.path.join(swc_dir, '*.swc')):
        tree = parse_swc(swcfile)
        morph = Morphology(tree)
        mc = MorphCurvature(morph, neurite_type=neurite_type)
        angd_dict = mc.estimate_angular_dependence(
                        discard_multifurcate=discard_multifurcate)
        angle_pairs.extend(list(angd_dict.values()))
        
        npd += 1
        if npd % 2 == 0:
            print(f'--> {len(angle_pairs)} angles found in {npd} files. Calculation takes {time.time() - t0:.2f} seconds')
        if npd >= total_process_file:   # use subset for statistics
            break

    angle_pairs = np.array(angle_pairs)
    with open(f'angd_{neurite_type}_swc{total_process_file}.pkl', 'wb') as fp:
        pickle.dump(angle_pairs, fp)
    
    #with open(f'../results/curvature/angd_{neurite_type}_swc{total_process_file}.pkl', 'rb') as fp:
    #    angle_pairs = pickle.load(fp)

    print(f'#angles: {len(angle_pairs)}')
    plt.figure(figsize=(8,6))
    #plt.scatter(angle_pairs[:,0], angle_pairs[:,1], s=1.0, c='b', alpha=0.8)
    plt.hist2d(angle_pairs[:,0], angle_pairs[:,1], bins=(100,100), cmap=plt.cm.nipy_spectral, density=True)
    plt.colorbar()
    plt.xlabel('angle1 (rad)')
    plt.ylabel('angle2 (rad)')
    #plt.legend()
    plt.xlim(2.09, 3.15)
    plt.ylim(2.09, 3.15)
    #plt.show()
    plt.savefig(f'angd_{neurite_type}_swc{total_process_file}.png', dpi=300)
    '''

