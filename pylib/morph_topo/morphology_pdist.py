#!/usr/bin/env python

#================================================================
#   Copyright (C) 2021 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : morphology_pdist.py
#   Author       : Yufeng Liu
#   Date         : 2021-07-22
#   Description  : 
#
#================================================================

import numpy as np
import time
from sklearn.metrics import pairwise_distances
from swc_handler import parse_swc
from morphology import Morphology


class PDist(object):
    """
    This is an optimized version of previous pairwise_distance.py, mainly by 
    full leveraging of parallelization to speed up calculation
    """

    def __init__(self, ignore_radius_from_soma=50, offspring_thresh=10):
        self.ignore_radius_from_soma = ignore_radius_from_soma
        self.offspring_thresh = offspring_thresh

    def set_morph(self, morph):
        self.morph = morph
        tmp = [node[2:5] for node in self.morph.tree]
        #import ipdb; ipdb.set_trace()
        self.coords = np.array([node[2:5] for node in self.morph.tree])
        self.idxs = np.array([node[0] for node in self.morph.tree])
        self.indexs = np.arange(len(self.morph.tree))

    def get_soma_nearby_nodes(self):
        soma_coord = self.coords[self.morph.index_soma]
        dists = np.linalg.norm(self.coords - soma_coord, axis=1)

        mask = dists < self.ignore_radius_from_soma
        nodes_near_soma = self.idxs[mask]
        nodes_away_soma = self.idxs[~mask]
        
        return nodes_near_soma, nodes_away_soma

    def get_linkages_with_thresh(self):
        # get parents within thresh
        parent_dict = {}
        for idx in self.idxs:
            leaf = self.morph.pos_dict[idx]
            # exclude parent
            os_id = 0 
            cur_set = []
            while os_id < self.offspring_thresh:
                try:
                    p_leaf = self.morph.pos_dict[leaf[-1]]
                    cur_set.append(p_leaf[0])

                    leaf = p_leaf # update leaf
                    os_id += 1
                except KeyError:
                    break
            parent_dict[idx] = set(cur_set)

        offspring_dict = {}
        for ofs, parents in parent_dict.items():
            for p_idx in parents:
                try:
                    offspring_dict[p_idx].append(ofs)
                except KeyError:
                    offspring_dict[p_idx] = [ofs]
        # convert to set
        for key, value in offspring_dict.items():
            offspring_dict[key] = set(value)

        return parent_dict, offspring_dict

    def find_crossing_pairs(self, crossing_thresh=3.0):
        # auxilary functions
        def filter_by_common_parent(parent_dict, cur_idx, crossing_idxs, crossing_coords):
            ncrossing = len(crossing_idxs)
            for i in range(ncrossing-1, -1, -1):
                idx = crossing_idxs[i]
                has_common_parent = False
                for pr0 in parent_dict[cur_idx]:
                    for pr1 in parent_dict[idx]:
                        if pr0 == pr1:
                            has_common_parent = True
                            break
                    if has_common_parent: break
                if has_common_parent:
                    # pop those points
                    crossing_idxs.pop(i)
                    crossing_coords.pop(i)


        nodes_away_soma = self.get_soma_nearby_nodes()[-1]
        print(f'Total number of nodes: {len(nodes_away_soma)}')
        tt = crossing_thresh * 1.5  # threshold for duplicated pair detection
        
        # get the linkages with thresh
        parent_dict, offspring_dict = self.get_linkages_with_thresh()

        t0 = time.time()
        nc = 0
        crossing_dict = {}
        for idx in nodes_away_soma[:-1]:
            if nc % 1000 == 0:
                print(f'--> {nc / len(nodes_away_soma):.2%} finished in {time.time() - t0}s')

            # get the nodes
            leaf = self.morph.pos_dict[idx]
            pts = set(nodes_away_soma[idx:]) - parent_dict[idx] - set([idx])
            try:
                pts = pts - offspring_dict[idx]
            except KeyError:
                pass

            # all curr_distances
            cur_pos = self.coords[self.morph.index_dict[idx]]
            node_idxs = list(pts)
            node_indexs = [self.morph.index_dict[ii] for ii in node_idxs]
            node_coords = self.coords[node_indexs]
            
            offset = node_coords - cur_pos
            dists = np.linalg.norm(offset, axis=1)
            # crossing points
            crossing_mask = dists < crossing_thresh
            crossing_is = np.nonzero(crossing_mask)[0]
            crossing_idxs = [node_idxs[ii] for ii in crossing_is]
            crossing_indexs = [node_indexs[ii] for ii in crossing_is]
            crossing_coords = list(node_coords[crossing_is])
           
            # filter by common parent in place
            filter_by_common_parent(parent_dict, idx, crossing_idxs, crossing_coords)
            for icoord, iidx in zip(crossing_coords, crossing_idxs):
                # check if two close to existing pairs
                has_near_pair = False
                for pidxs, pcs in crossing_dict.items():
                    if (np.linalg.norm(cur_pos - pcs[0]) < tt and \
                        np.linalg.norm(icoord - pcs[1]) < tt) or \
                       (np.linalg.norm(cur_pos - pcs[1]) < tt and \
                        np.linalg.norm(icoord - pcs[0]) < tt):
                        has_near_pair = True
                        break
                if not has_near_pair:
                    crossing_dict[(idx, iidx)] = (cur_pos, icoord)
                
            nc += 1

        print(f"{len(crossing_dict)} crossing points within {crossing_thresh} voxels")
        #print(f"{crossing_dict}")
            
        return crossing_dict

    @staticmethod
    def get_crossing_point(pts0, pts1, sampling_step=0.02):
        """
        Find out the crossing point from given anchor points. 
        Arguments are:
        - pts0: anchor points set1, with first (pts0[0]) is the nodes nearest to 
                pts1[0]. The other points in pts0 are nodes directly connected to 
                pts0[0].
        - pts1: as pts0
        - sampling_step: sampling step size for point detection
        """
        def get_finer_anchors(pts, grids):
            pt0 = np.array(pts[0])
            anchors_all = []
            for i in range(1, len(pts)):
                v = np.array(pts[i]) - pt0
                anchors = pt0 + grids.reshape(-1,1) * v
                anchors_all.append(anchors)
            anchors_all = np.vstack(anchors_all)
            return anchors_all

        # grids for finer anchors interpolation
        grids = np.arange(0, 1+1e-7, sampling_step)
        anchors0 = get_finer_anchors(pts0, grids)
        anchors1 = get_finer_anchors(pts1, grids)
        pdists = pairwise_distances(anchors0, anchors1)
        dmin = pdists.min()
        px, py = np.unravel_index(pdists.argmin(), pdists.shape)
        # find out the position
        pi = px // len(grids) + 1
        pj = py // len(grids) + 1
        return pi, pj, dmin, anchors0[px], anchors1[py]


    def get_crossing_points(self, cur_crossing_dict):
        def get_points(idx, pos_dict, child_dict):
            # find out the direct connecting points
            idxs = [idx, pos_dict[idx][6]]
            if idx in child_dict:
                for ii in child_dict[idx]:
                    idxs.append(ii)
            coords = np.array([pos_dict[idx][2:5] for idx in idxs])
            return coords, idxs
        
    
        crossing_points = []
        for (idx0, idx1), (c0, c1) in cur_crossing_dict.items():
            try:
                pts0, idxs0 = get_points(idx0, self.morph.pos_dict, self.morph.child_dict)
                pts1, idxs1 = get_points(idx1, self.morph.pos_dict, self.morph.child_dict)
            except KeyError:
                print(f"Error node found, mostly probably the node with parent id -1")
                continue
            # 
            pi, pj, dmin, c_near0, c_near1 = self.get_crossing_point(pts0, pts1)
            
            # The meaning of pi (pj as well) is the index on pts0, which is in order: 
            # [anchor point, parent point, child points], thus 1 means parent, and >=2 
            # denotes a child node
            # change to parent to child order for each fiber
            #print(pi, pj, pts0, pts1)
            if pi == 1: # parent node
                i00, i01 = idxs0[pi], idx0
                c00, c01 = pts0[pi], c0
            elif pi > 1:
                i00, i01 = idx0, idxs0[pi]
                c00, c01 = c0, pts0[pi]
            else:
                raise ValueError(f"Incorrect index {pi} for pi")
            if pj == 1:
                i10, i11 = idxs1[pj], idx1
                c10, c11 = pts1[pj], c1
            elif pj > 1:
                i10, i11 = idx1, idxs1[pj]
                c10, c11 = c1, pts1[pj]
            else:
                raise ValueError(f"Incorrect index {pj} for pj")

            parent_to_cur0 = pi == 1
            parent_to_cur1 = pj == 1
            crossing_points.append(((i00, i01, i10, i11), (c00, c01, c10, c11), (parent_to_cur0, parent_to_cur1), (c_near0, c_near1), dmin))

        return crossing_points

    def get_isometric_anchors(self, crossing_points, anchor_dist, return_center_point=True):
        """
        find out the anchors on both crossing fibers, which are isometric to anchor points
        """
        def find_point_with_distance(pt, anchor_idx, is_parent, morph, dist, return_center_point=True, epsilon=1e-7):
            """
            Find the point of exact `dist` to the start pt on tree structure. args are:
            - pt: the start point
            - anchor_idx: the first node on swc tree to trace. 
            - is_parent: whether the anchor_idx is the parent of `pt`, otherwise child. 
                         if the node has several child, a random one is selected
            - morph: Morphology object for current tree
            - dist: distance threshold
            - return_center_point: whether to return the point with exact distance or
                         geometric point of all traced nodes
            - epsilon: small float to avoid zero-division error 
            """

            d = 0
            ci = pt
            pts = [pt]
            while d < dist:
                try:
                    cc = np.array(morph.pos_dict[anchor_idx][2:5])
                except KeyError:
                    print(f"Parent/Child node not found within distance: {dist}")
                    break
                d0 = np.linalg.norm(ci - cc)
                d += d0
                if d < dist:
                    ci = cc  # update coordinates
                    pts.append(cc)

                    if is_parent:
                        anchor_idx = morph.pos_dict[anchor_idx][6]
                    else:
                        if anchor_idx not in morph.child_dict:
                            break
                        else:
                            anchor_idxs = morph.child_dict[anchor_idx]
                            anchor_idx = anchor_idxs[np.random.randint(0, len(anchor_idxs))]

            # interpolate to find the exact point
            dd = d - dist
            if dd < 0:
                pt_a = cc
            else:
                dcur = np.linalg.norm(cc - ci)
                assert(dcur - dd >= 0)
                pt_a = ci + (cc - ci) * (dcur - dd) / (dcur + epsilon)
                pts.append(pt_a)

            if return_center_point:
                pt_a = np.mean(pts, axis=0)

            return pt_a


        all_pts = []
        for cpt in crossing_points:
            idxs, coords, pcs, c_nears, dmin = cpt
            anchor_pts = []
            #print(len(idxs), idxs)
            #print(len(coords), coords)
            for i, pc, c_near in zip(range(2), pcs, c_nears):
                idxs0 = idxs[i*2: 2*i+2]
                coords0 = coords[i*2: 2*i+2]
                #print(i, idxs0, coords0)
                #import ipdb; ipdb.set_trace()
                pt0 = find_point_with_distance(c_near, idxs0[0], True, self.morph, anchor_dist)
                pt1 = find_point_with_distance(c_near, idxs0[1], False, self.morph, anchor_dist)
                # test only:
                #print(f'Direct distance for {i}: {np.linalg.norm(pt0 - c_near)}, {np.linalg.norm(pt1 - c_near)}, {np.linalg.norm(pt0 - pt1)}')
                anchor_pts.extend([c_near, pt0, pt1])
                
            all_pts.append(anchor_pts)

        return all_pts

    def get_crossing_points_for_file(self, swcfile, cur_dict, filter_thresh=None):
        tree = parse_swc(swcfile)
        morph = Morphology(tree)
        self.set_morph(morph)
        cpts = self.get_crossing_points(cur_dict)
        if filter_thresh is not None:
            # do filtering with pairwise point
            for i in range(len(cpts)-1, -1, -1):
                cpt = cpts[i]
                if cpt[-1] > filter_thresh:
                    cpts.pop(i)
        
        return cpts


if __name__ == '__main__':
    import os, glob, sys
    import pickle
    import cProfile, pstats
    import io
    from morphology import Morphology
    from swc_handler import parse_swc

    swc_dir = '/media/lyf/storage/seu_mouse/swc/xy1z1'    
    ignore_radius_from_soma = 50.
    offspring_thresh = 10
    crossing_thresh = 5.0

    pd = PDist(ignore_radius_from_soma, offspring_thresh)
    """
    #Generatation of crossing_dict
    
    crossing_dict = {}
    n_p = 0
    for swcfile in glob.glob(os.path.join(swc_dir, '*swc')):
        print(swcfile)
        swc_name = os.path.split(swcfile)[-1]
        tree = parse_swc(swcfile)
        morph = Morphology(tree)
        pd.set_morph(morph)
        cur_dict = pd.find_crossing_pairs(crossing_thresh)
        crossing_dict[swc_name] = cur_dict
    
        n_p += 1
        print(f'===> Finished {n_p} files...')

    with open('crossing_5voxels_maxRes.pkl', 'wb') as fp:
        pickle.dump(crossing_dict, fp)
    """

    anchor_dist = 5.0
    with open('crossing_5voxels_maxRes.pkl', 'rb') as fp:
        crossing_dict = pickle.load(fp)

    n_bi = 0
    n_uni = 0
    n_bi_tot = 0
    n_uni_tot = 0
    nprocessed = 0
    for swc_name, cur_dict in crossing_dict.items():
        if swc_name != '18868_4944_x7089_y8035.swc': continue
        print(f'==> Processing for swc: {swc_name}')

        swcfile = os.path.join(swc_dir, swc_name)
        cpts = pd.get_crossing_points_for_file(swcfile, cur_dict, filter_thresh=3)
        pd.morph.get_critical_points()
        print(f'Nuber of unifurcation and bifurcation nodes: {len(pd.morph.unifurcation)} and {len(pd.morph.bifurcation)}')
        n_bi_tot += len(pd.morph.bifurcation)
        n_uni_tot += len(pd.morph.unifurcation)
        for cpt in cpts:
            idx0 = cpt[0][1] if cpt[2][0] else cpt[0][0]
            idx1 = cpt[0][3] if cpt[2][1] else cpt[0][2]
            if idx0 in pd.morph.bifurcation:
                n_bi += 1
            else:
                n_uni += 1
        print(f'Unifurcation and bifurcation number: {n_uni} and {n_bi}')

        # get the crossing structure
        crossing_structures = pd.get_isometric_anchors(cpts, anchor_dist)
        #print(crossing_structures)

        nprocessed += 1
        if nprocessed % 20 == 0:
            print(f'----> Processed {nprocessed} files')
    print(f'Total uni- and bi-furcation are: {n_uni_tot}, {n_bi_tot}')

    
