import sys
sys.path.append(r'D:\code\support_lib\pylib')
import swc_handler
import numpy as np
import os
from pathlib import Path
from multiprocessing.pool import Pool
from math_utils import min_distances_between_two_sets
from morph_topo.morphology import Morphology


def inbox(xyz, shape, center):
    xyz = np.abs(np.array(xyz) - center) * 2
    return (xyz < np.array(shape)).all()
    

def main(in_ntb):
    outpath = Path(in_ntb).relative_to(ntb_dir)
    in_app2 = swc_dir / outpath.with_name(outpath.name.replace('_neutube', ''))
    outpath = outpath.with_name(outpath.name.replace('_neutube', ''))
    if outpath.exists():
        return
    if not in_app2.exists():
        return
    t1 = swc_handler.parse_swc(in_app2)
    t2 = swc_handler.parse_swc(in_ntb)
    if len(t1) == 0 or len(t2) == 0:
        return
    center = np.array(t1[0][2:5])
    pts1 = np.array([i[2:5] for i in t1])
    pts2 = swc_handler.tree_to_voxels(t2, block_size)
    dist = min_distances_between_two_sets(pts1, pts2, reciprocal=False)
    node = set(t1[i][0] for i in np.argwhere(dist.reshape(-1) > min_dist).reshape(-1) if np.linalg.norm(t1[i][2:5] - center) > soma_radius)
    t = swc_handler.prune(t1, node)
    try:
        m = Morphology(t)
        frag_len_dict = m.calc_frag_lengths()[1]
        seg_dict = m.convert_to_topology_tree()[1]
        seglen_dict = m.calc_seg_path_lengths(seg_dict, frag_len_dict)
        t = m.prune_by_seg_length(seg_dict, seglen_dict, seg_length_thresh)
    except:
        pass
    if len(t) < min_node:
        return
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    swc_handler.write_swc(t, outpath)


if __name__ == '__main__':
    #ntb_dir = '../../230k_all/stage2_ntb'
    ntb_dir = '../ntb'
    #swc_dir = '../stage3_seg_prune_weak'
    swc_dir = '../seg_prune'
    #os.system(f"find {ntb_dir} -name *swc > all.txt")
    min_dist = 5
    seg_length_thresh = 10
    soma_radius = 15
    block_size = [256, 512, 512]
    min_node = 20
    with open('all.txt') as f:
        ntbs = [i.rstrip() for i in f.readlines()]
    with Pool(40) as pt:
        pt.map(main, ntbs, 10)