#!/usr/bin/env python

#================================================================
#   Copyright (C) 2021 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : neurite_shape.py
#   Author       : Yufeng Liu
#   Date         : 2021-07-30
#   Description  : size, intensity, and direction
#
#================================================================

import numpy as np
import copy
import time
import subprocess
from morphology import Morphology, Topology
from image_io import load_image, save_image
from swc_handler import parse_swc, write_swc, trim_swc, NEURITE_TYPES

class AbstractNeuriteShape(object):
    def __init__(self, morph):
        self.morph = morph
        if not hasattr(self.morph, 'tips'):
            self.morph.get_critical_points()

    def get_branch_dict(self):
        # firstly, find out the number of unifurcate nodes for each branch
        branch_dict = {}
        nodes_with_parents = self.morph.tips | self.morph.multifurcation
        for midx in nodes_with_parents:  # does not contain soma
            up_nodes = []
            idx = midx
            while idx in self.morph.pos_dict:
                pidx = self.morph.pos_dict[idx][-1]
                up_nodes.append(pidx)
                if pidx not in self.morph.unifurcation:
                    break

                idx = pidx  # update
            branch_dict[midx] = up_nodes
        return branch_dict

    def resample_tree(self):
        """
        Resample branches with nodes less than $min_nodes
        """
        def resample_branch(idx, pidx, new_pos_dict, nidx):
            coord = np.array(new_pos_dict[idx][2:5])
            pcoord =  np.array(new_pos_dict[pidx][2:5])
            c = (coord + pcoord) / 2

            node = new_pos_dict[idx]
            new_pos_dict[idx] = (*node[:-1], nidx)
            new_pos_dict[nidx] = (nidx, node[1], *c, node[5], pidx)

    
        branch_dict = self.get_branch_dict()
        new_pos_dict = copy.deepcopy(self.morph.pos_dict)
        # get the maximal idx for current tree. Note we could not assume
        # the idx are continuous as it may editted from standard swc.
        idx_max = max(self.morph.pos_dict.keys())

        for idx, ups in branch_dict.items():
            if len(ups) == 1:
                # resampling required
                resample_branch(idx, ups[0], new_pos_dict, idx_max)
                idx_max += 1
        # convert the new tree
        resampled_tree = [it[-1] for it in sorted(new_pos_dict.items(), key=lambda x:x[0])]
        print(f'#{len(resampled_tree)} nodes from #{len(self.morph.tree)}')
        return resampled_tree

class NeuriteShapeSingle(object):
    def __init__(self, swcfile, imgfile, use_local_maximal=True, lsize=2, normalize_image=True):
        self.use_local_maximal = use_local_maximal
        self.lsize = lsize
        self.imgfile = imgfile

        img = load_image(imgfile)
        if normalize_image:
            dtype = img.dtype
            img = img.astype(np.float32)
            img = (img - img.min()) / (img.max() - img.min() + 1e-7)
            img = (img * 255.).astype(dtype)
        self.img = img[:,::-1,:]    # flip the image to be consistent with swc
        self.imgshape = img.shape
        self.resampled_tree = self.get_resampled_tree(swcfile)

    def get_resampled_tree(self, tree):
        tree = parse_swc(swcfile)
        trimed_tree = trim_swc(tree, self.imgshape, keep_candidate_points=False)
        morph = Morphology(trimed_tree)
        ns = AbstractNeuriteShape(morph)
        resampled_tree = ns.resample_tree()
        if self.use_local_maximal:
            s = self.lsize
            # re-locate the nodes according to intensity
            for i, node in enumerate(resampled_tree):
                xi,yi,zi = map(int, map(round, node[2:5]))
                block = self.img[max(zi-s,0):zi+s+1, 
                                 max(yi-s,0):yi+s+1, 
                                 max(xi-s,0):xi+s+1]
                mindices = np.unravel_index(block.argmax(), block.shape)
                zn = mindices[0] + max(zi-s, 0)
                yn = mindices[1] + max(yi-s, 0)
                xn = mindices[2] + max(xi-s, 0)
                resampled_tree[i] = (*node[:2],xn,yn,zn,*node[5:])

        return resampled_tree

    def get_branch_intensity_dict(self, debug_vis=False):
        morph = Morphology(self.resampled_tree)
        morph.get_critical_points()
        topo_tree = morph.convert_to_topology_tree()
        topo = Topology(topo_tree)
        ns = AbstractNeuriteShape(morph)
        branch_dict = ns.get_branch_dict()

        ins_dict = {}
        ins_std_dict = {}
        if debug_vis:
            imgout = np.zeros(self.img.shape, dtype=img.dtype)
        
        s = self.lsize
        for idx, pidxs in branch_dict.items():
            if morph.idx_soma in pidxs:
                continue    # discard soma connecting branch
            # get the positions
            ins = []
            for pi in pidxs[:-1]:
                xi,yi,zi = np.round(morph.pos_dict[pi][2:5]).astype(int)
                ic = self.img[zi,yi,xi]
                if debug_vis:
                    imgout[max(zi-s,0):zi+s+1, 
                           max(yi-s,0):yi+s+1, 
                           max(xi-s,0):xi+s+1] = \
                          self.img[max(zi-s,0):zi+s+1, 
                             max(yi-s,0):yi+s+1, 
                             max(xi-s,0):xi+s+1]
                
                ins.append(ic)
            # use the median intensity so as to get stable result
            ins.sort()
            ic = ins[len(ins)//2]
            ins_dict[idx] = ic
            
            if len(ins) > 2:
                ins = np.array(ins)
                ins_std_dict[idx] = ins.std()# / ins.mean()
        
        if debug_vis:
            save_image('debug.tiff', imgout)

        # push the morph and topo object to self
        self.morph = morph
        self.topo = topo
        self.ins_dict = ins_dict
        self.ins_std_dict = ins_std_dict

        return ins_dict, ins_std_dict

    def get_branch_radius_dict(self, vaa3d_path):
        def estimate_radius(vaa3d_path, tree, imgfile):
            current_time = time.localtime()
            time_str = time.strftime('%H-%M-%S', current_time)
            inswcfile = f'temp_in_{time_str}.swc'
            outswcfile = f'{inswcfile}.out.swc'
            write_swc(tree, inswcfile)
            # call vaa3d to estimate
            subprocess.check_output(f'{vaa3d_path} -x neuron_radius -f neuron_radius -i {imgfile} {inswcfile} -o {outswcfile} -p AUTO 1', shell=True)
            # parse the new generate swc with estimated radius
            tree = parse_swc(inswcfile)
            new_tree = parse_swc(outswcfile)
            # the idx are modified while saving with vaa3d, correct it
            for i, onode, cnode in zip(range(len(tree)), tree, new_tree):
                new_tree[i] = (onode[0], *cnode[1:-1], onode[-1])

            radius_dict = {}
            for node in new_tree:
                idx, r = node[0], node[-2]
                radius_dict[idx] = r

            #write image for debuging
            #os.system(f'cp {imgfile} temp_img_{time_str}.tiff')
            # remove the temporary files
            os.system(f'rm -f {inswcfile} {outswcfile}')
            return radius_dict
            

        if hasattr(self, 'morph'):
            morph = self.morph
            topo = self.topo
        else:
            morph = Morphology(self.resampled_tree)
            morph.get_critical_points()
            topo_tree = morph.convert_to_topology_tree()
            topo = Topology(topo_tree)
        ns = AbstractNeuriteShape(morph)
        branch_dict = ns.get_branch_dict()
        radius_dict = estimate_radius(vaa3d_path, self.resampled_tree, self.imgfile)

        rad_dict = {}
        rad_std_dict = {}
        for idx, pidxs in branch_dict.items():
            if morph.idx_soma in pidxs:
                continue    # discard soma connecting branch
            # get the positions
            rads = []
            for pi in pidxs[:-1]:
                rad = radius_dict[pi]
                rads.append(rad)
            # use the median intensity so as to get stable result
            rads.sort()
            rad = rads[len(rads)//2]
            rad_dict[idx] = rad
            
            if len(rads) > 2:
                rads = np.array(rads)
                rad_std_dict[idx] = rads.std()# / ins.mean()
        
        # push the morph and topo object to self
        if not hasattr(self, 'morph'):
            self.morph = morph
            self.topo = topo
        self.rad_dict = rad_dict
        self.rad_std_dict = rad_std_dict

        return rad_dict, rad_std_dict

    def calc_global_shape(self, shape_dict, neurite_type):
        # pre-calculate the path length to soma
        seg_lengths = self.morph.calc_seg_lengths()
        path_dict = self.morph.get_path_idx_dict()
        plen_dict = self.morph.get_path_len_dict(path_dict, seg_lengths)

        order_shape = []
        plen_shape = []
        for idx, ss in shape_dict.items():
            order = self.topo.order_dict[idx]
            plen = plen_dict[idx]
            type_ = self.morph.pos_dict[idx][1]
            if (neurite_type == 'axon') and (type_ not in NEURITE_TYPES['axon']):
                continue
            elif (neurite_type == 'dendrite') and (type_ not in NEURITE_TYPES['dendrite']):
                continue
            order_shape.append([order, ss])
            plen_shape.append([plen, ss])
        return order_shape, plen_shape


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    import os, glob
    import pickle
    import time

    data_dir = '/media/lyf/storage/seu_mouse/crop_data/processed/dendriteImageMaxR'
    brain_id = '18455'
    use_local_maximal = True
    lsize = 2 # local path size
    shape_type = 'radius' # radius or intensity
    neurite_type = 'dendrite'   # axon, dendrite, overall
    vaa3d_path = '/home/lyf/Softwares/installation/Vaa3D/v3d_external/bin/vaa3d'

    
    ################### intensity/size statistics ##################
    order_ins_all = []
    plen_ins_all = []
    num = 0
    t0 = time.time()
    for swcfile in glob.glob(f'{data_dir}/swc_fitted_zhongye/{brain_id}/*swc'):
        prefix = os.path.split(swcfile)[-1][:-4]

        print(f'---> {prefix}')
        imgfile = f'{data_dir}/tiff/{brain_id}/{prefix}.tiff'
        nss = NeuriteShapeSingle(swcfile, imgfile, use_local_maximal, lsize, normalize_image=True)
        if shape_type == 'intensity':
            nss.get_branch_intensity_dict()
            order_ins, plen_ins = nss.calc_global_shape(nss.ins_dict, neurite_type)
        else:
            nss.get_branch_radius_dict(vaa3d_path)
            order_ins, plen_ins = nss.calc_global_shape(nss.rad_dict, neurite_type)
        order_ins_all.extend(order_ins)
        plen_ins_all.extend(plen_ins)

        num += 1
        if num % 50 == 0:
            print(f'#{num} files finished in {time.time() - t0:.2f}')
            break

    order_ins_all = np.array(order_ins_all)
    plen_ins_all = np.array(plen_ins_all)

    with open(f'{shape_type}_brain{brain_id}_{neurite_type}.pkl', 'wb') as fp:
        pickle.dump([order_ins_all, plen_ins_all], fp)
    
    #with open(f'intensity_brain{brain_id}.pkl', 'rb') as fp:
    #    order_ins_all, plen_ins_all = pickle.load(fp)
    
    plt.scatter(plen_ins_all[:,0], plen_ins_all[:,1], s=1.0, c='r', alpha=0.5)
    plt.xlabel('path length')
    plt.ylabel(f'{shape_type}')
    plt.legend()
    plt.savefig(f'path_vs_{shape_type}_{brain_id}_50images_{neurite_type}.png', dpi=150)
    plt.close()
    plt.scatter(order_ins_all[:,0], order_ins_all[:,1], s=1.0, c='b', alpha=0.5)
    plt.xlabel('topology order')
    plt.ylabel(f'{shape_type}')
    plt.legend()
    plt.savefig(f'order_vs_{shape_type}_{brain_id}_50images_{neurite_type}.png', dpi=150)
    
    
    """
    ################## intensity/radius decreasing assumption ##############
    num_true, num_false = 0, 0
    nprocessed = 0
    stds = []
    for swcfile in glob.glob(f'{data_dir}/swc_fitted_zhongye/{brain_id}/*swc'):
        prefix = os.path.split(swcfile)[-1][:-4]

        print(f'---> {prefix}')
        imgfile = f'{data_dir}/tiff/{brain_id}/{prefix}.tiff'
        nss = NeuriteShapeSingle(swcfile, imgfile, use_local_maximal, \
                                 lsize, normalize_image=True)
        if shape_type == 'intensity':
            nss.get_branch_intensity_dict()
        elif shape_type == 'radius':
            nss.get_branch_radius_dict(vaa3d_path)

        for node in nss.topo.tree:
            idx, pidx = node[0], node[-1]
            if idx == 1: 
                continue   # skip soma
            if shape_type == 'intensity':
                if pidx not in nss.ins_dict:
                    continue
                i0 = nss.ins_dict[idx]
                pi0 = nss.ins_dict[pidx]
            elif shape_type == 'radius':
                if pidx not in nss.rad_dict:
                    continue
                i0 = nss.rad_dict[idx]
                pi0 = nss.rad_dict[pidx]
            
            # skip points do not match neurite type
            type_ = nss.morph.pos_dict[idx][1]
            if (neurite_type == 'axon') and (type_ not in NEURITE_TYPES['axon']):
                continue
            elif (neurite_type == 'dendrite') and (type_ not in NEURITE_TYPES['dendrite']):
                continue
                if pi0 > i0:
                    num_true += 1
                else:
                    num_false += 1
            # for std statistics
            try:
                if shape_type == 'intensity':
                    stds.append(nss.ins_std_dict[idx])
                elif shape_type == 'radius':
                    stds.append(nss.rad_std_dict[idx])
            except KeyError:
                continue

        nprocessed += 1
        if nprocessed % 50 == 0:
            print(f'==> #processed: {nprocessed}')
            break

    print(f'--> #true ({num_true}) / #false ({num_false}). \n    Correct ratio are: {1.0*num_true/(num_true+num_false)}')

    stds = np.array(stds)
    print(stds.mean(), stds.max(), stds.min(), stds.std())
    with open('temp.pkl', 'wb') as fp:
        pickle.dump(stds, fp)
    """   

