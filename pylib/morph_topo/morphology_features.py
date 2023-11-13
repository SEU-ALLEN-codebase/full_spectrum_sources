#!/usr/bin/env python

#================================================================
#   Copyright (C) 2022 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : morphology_features.py
#   Author       : Yufeng Liu
#   Date         : 2022-04-15
#   Description  : 
#
#================================================================

import os, sys, glob
import subprocess
import numpy as np

from swc_handler import parse_swc, write_swc, scale_swc
from file_io import load_image
from math_utils import calc_included_angles_from_coords, calc_included_angles_from_vectors
from morph_topo.morphology import Morphology, Topology
from neuron_quality.find_break_crossing import find_point_by_distance, BreakFinder, CrossingFinder


class TopoFeatures(object):
    """
    Calculate the morphological features according from swc
    """

    def __init__(self, swcfile, line_length=8.0, z_factor=1/0.23):
        if isinstance(swcfile, str):
            tree = parse_swc(swcfile)
        elif isinstance(swcfile, list):
            tree = swcfile
        else:
            raise NotImplementedError(f"Parameter swcfile must be either str or list, but not {type(swcfile)}")

        tree_iso = scale_swc(tree, scale=(1.0,1.0,z_factor))

        self.morph = Morphology(tree_iso)
        topo_tree, seg_dict = self.morph.convert_to_topology_tree()
        self.topo = Topology(topo_tree)
        self.seg_dict = seg_dict
        self.line_length = line_length

    def dists_to_soma(self, morph_lengths_dict):
        def dist_dfs(idx, idx_soma, pd, child_dict, pos_dict, morph_lengths_dict):
            if idx == idx_soma:
                pd[idx] = 0
            else:
                par_id = pos_dict[idx][6]
                pd[idx] = pd[par_id] + morph_lengths_dict[par_id]

            if idx in child_dict:
                for child_id in child_dict[idx]:
                    dist_dfs(child_id, idx_soma, pd, child_dict, pos_dict, morph_lengths_dict)
            
        
        # path distance        
        path_dists = {}
        idx = self.morph.idx_soma
        dist_dfs(idx, self.morph.idx_soma, path_dists, self.morph.child_dict, self.morph.pos_dict, morph_lengths_dict)
        # spatial distance
        spatial_dists = {}
        for node, dist in zip(self.morph.tree, self.morph.get_distances_to_soma()):
            spatial_dists[node[0]] = dist
        return path_dists, spatial_dists
        
    def dists_to_parent_seg(self, morph_lengths_dict, topo_lengths_dict):
        path_dists = self.morph.calc_seg_path_lengths(self.seg_dict, morph_lengths_dict)

        return path_dists, topo_lengths_dict
  
    def get_angles(self):
        # local angle
        local_angs = {}
        ang_default = np.pi
        local_angs[self.morph.idx_soma] = ang_default
        global_angs = {}
        global_angs[self.morph.idx_soma] = ang_default
        soma_pos = np.array(self.morph.pos_dict[self.morph.idx_soma][2:5])
        for seg_id, seg_nodes in self.seg_dict.items():
            if seg_id == self.topo.idx_soma: continue
            
            par_topo = self.topo.pos_dict[seg_id][6]
            if par_topo == self.morph.idx_soma:
                local_angs[seg_id] = ang_default
                global_angs[seg_id] = ang_default
            else:
                c = np.array(self.topo.pos_dict[par_topo][2:5])   # xyz coordinates
                if len(seg_nodes) == 0:
                    start_node = seg_id
                else:
                    start_node = seg_nodes[-1]
                c1 = find_point_by_distance(c, start_node, False, self.morph, self.line_length, False)
                seg_nodes2 = self.seg_dict[par_topo]
                if len(seg_nodes2) == 0:
                    start_node2 = self.topo.pos_dict[par_topo][6]
                else:
                    start_node2 = seg_nodes2[0]
                c2 = find_point_by_distance(c, start_node2, True, self.morph, self.line_length, False)
                local_ang = calc_included_angles_from_coords(c, c1, c2, return_rad=True)[0]
                local_angs[seg_id] = local_ang

                global_ang = calc_included_angles_from_coords(c, c1, soma_pos, return_rad=True)[0]
                global_angs[seg_id] = global_ang
        
        return local_angs, global_angs

    def get_num_childs(self):
        nchilds_dict = {}
        for idx in self.topo.pos_dict:
            par_id = self.topo.pos_dict[idx][6]
            if idx == self.topo.idx_soma:
                nchilds_dict[idx] = (1, 0)    # nchilds of parent, current
            else:
                if idx not in self.topo.child_dict:
                    nchilds_dict[idx] = (len(self.topo.child_dict[par_id]), 0)
                else:
                    nchilds_dict[idx] = (len(self.topo.child_dict[par_id]), len(self.topo.child_dict[idx]))
        return nchilds_dict

    def get_relative_coords(self):
        soma_pos = np.array(self.morph.pos_dict[self.morph.idx_soma][2:5])
        rcoords = {}
        for seg_id, seg_nodes in self.seg_dict.items():
            coord = np.array(self.morph.pos_dict[seg_id][2:5])
            rcoord = coord - soma_pos
            rcoords[seg_id] = rcoord
        return rcoords

    def calc_all_features(self):
        if len(self.morph.tree) == 0:
            return False

        _, morph_lengths_dict = self.morph.calc_frag_lengths()
        _, topo_lengths_dict = self.topo.calc_frag_lengths()
        
        pdists_soma, sdists_soma = self.dists_to_soma(morph_lengths_dict)
        pdists_seg, sdists_seg = self.dists_to_parent_seg(morph_lengths_dict, topo_lengths_dict)
        local_angs, global_angs = self.get_angles()
        nchilds_dict = self.get_num_childs()
        
        topo_features = {}
        topo_features['pdists_soma'] = pdists_soma  # idx: dist
        topo_features['sdists_soma'] = sdists_soma  # idx: dist
        topo_features['pdists_seg'] = pdists_seg    # idx: dist
        topo_features['sdists_seg'] = sdists_seg    # idx: dist
        topo_features['local_angs'] = local_angs    # idx: ang
        topo_features['global_angs'] = global_angs  # idx: ang
        topo_features['nchilds_dict'] = nchilds_dict    # idx: (#par, #cur)
        topo_features['order_dict'] = self.topo.order_dict  # idx: order(int)

        return topo_features
     

class TopoImFeatures(object):
    def __init__(self, swcfile, imgfile, radii_cache_dir='./radii_cache', vaa3d_path='/home/lyf/Softwares/installation/Vaa3D/v3d_external/bin/vaa3d'):
        self.radii_cache_dir = radii_cache_dir
        self.vaa3d_path = vaa3d_path
        self.imgfile = imgfile
        
        if isinstance(swcfile, str):
            tree = parse_swc(swcfile)
        elif isinstance(swcfile, list):
            tree = swcfile
        else:
            raise NotImplementedError(f"Parameter swcfile must be either str or list, but not {type(swcfile)}")

        # get anisotropic object
        self.morph = Morphology(tree)
        topo_tree, seg_dict = self.morph.convert_to_topology_tree()
        self.topo = Topology(topo_tree)
        self.seg_dict = seg_dict
        self.swcfile = swcfile
        
        if isinstance(imgfile, str):
            self.img = load_image(imgfile)
        elif isinstance(imgfile, np.ndarray):
            self.img = imgfile
        else:
            raise NotImplementedError(f"Parameter imgfile must be either str or np.ndarray, not {type(imgfile)}")

    def get_node_intensities(self):
        int_dict = {}
        yshape = self.img.shape[1]
        for node in self.morph.tree:
            idx, type_, x, y, z = node[:5]
            x, y, z = map(lambda t: int(round(t)), [x, y, z])
            int_dict[idx] = self.img[z,yshape-y-1,x]  # the y-axis is reverted
        return int_dict

    def seg_intensities(self):
        int_dict = self.get_node_intensities()
        seg_int_dict = {}
        # the seg_id, that is tip or branch point is not calculated, as it may be error!
        for seg_id, seg_nodes in self.seg_dict.items():
            ints = [int_dict[seg_node] for seg_node in seg_nodes]
            if len(ints) == 0:
                seg_int_dict[seg_id] = (-1., -1., -1., -1.)
            else:
                int_median = np.median(ints)
                int_max = max(ints)
                int_min = min(ints)
                int_mean = np.mean(ints)
                seg_int_dict[seg_id] = (int_max, int_min, int_mean, int_median)
        
        return seg_int_dict

    def seg_radii(self):
        cache_dir = self.radii_cache_dir
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)

        swcname = os.path.split(self.swcfile)[-1]
        outswcfile = os.path.join(cache_dir, f'radius_{swcname}')
        if not os.path.exists(outswcfile):
            subprocess.check_output(f'{self.vaa3d_path} -x neuron_radius -f neuron_radius -i {self.imgfile} {self.swcfile} -o {outswcfile} -p AUTO 1', shell=True)
        # parse the new generated swc with estimated radius
        rad_tree = parse_swc(outswcfile)
        # the idx are modified while saving with vaa3d, correct it!
        for i, onode, cnode in zip(range(len(rad_tree)), self.morph.tree, rad_tree):
            rad_tree[i] = (onode[0], *cnode[1:6], *onode[6:])
        radius_dict = {}
        for node in rad_tree:
            idx, r = node[0], node[5]
            radius_dict[idx] = r

        # seg-level estimation
        rad_dict = {}
        for seg_id, seg_nodes in self.seg_dict.items():
            rads = [radius_dict[seg_node] for seg_node in seg_nodes]
            if len(rads) == 0:
                rad_dict[seg_id] = (-1., -1., -1., -1.)
            else:
                rad_median = np.median(rads)
                rad_max = max(rads)
                rad_min = min(rads)
                rad_mean = np.mean(rads)
                rad_dict[seg_id] = (rad_max, rad_min, rad_mean, rad_median)

        return rad_dict
        
    def calc_all_features(self):
        if len(self.morph.tree) == 0:
            return False
        
        # image features
        tif = {}
        tif['intensity'] = self.seg_intensities()   # idx: (max, min, mean, median)
        tif['radii'] = self.seg_radii()             # idx: (max, min, mean, median)

        return tif

if __name__ == '__main__':
    swcfile = '/media/lyf/storage/seu_mouse/crop_data/processed/dendriteImageMaxR/app2/18869/10048_8350_4580.swc'
    imgfile = '/media/lyf/storage/seu_mouse/crop_data/processed/dendriteImageMaxR/tiff/18869/10048_8350_4580.tiff'
    radius_dir = '/media/lyf/storage/seu_mouse/crop_data/processed/dendriteImageMaxR/app2_radii'
    line_length = 8
    z_factor = 1/0.23

    tf = TopoFeatures(swcfile, line_length=line_length, z_factor=z_factor)
    tf_feat = tf.calc_all_features()
    tif = TopoImFeatures(swcfile, imgfile, radii_cache_dir=radius_dir)
    tif_feat = tif.calc_all_features()
    feats = tf_feat
    for k, v in tif_feat.items():
        feats[k] = v

    for key, feat in feats.items():
        print(key, len(feat))
    print(feats['intensity'], feats['radii'])
    
