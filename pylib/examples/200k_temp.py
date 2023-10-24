import sys
sys.path.append(sys.path[0] + "/..")
import os
from file_io import load_image, save_markers
import swc_handler
from neuron_quality.find_break_crossing import CrossingFinder
from morph_topo.morphology import Morphology
import numpy as np
from scipy.interpolate import CubicSpline
import math
from neuron_quality.find_break_crossing import find_point_by_distance


def anchor_angles(center: np.ndarray, p: np.ndarray, ch, spacing=(1, 1, 4)):
    """
    modified from Jingzhou's code
    :param center: coordinate for the multifurcation center
    :param p: coordinate for the parent anchor
    :param ch: coordinates for the furcation anchors
    :param spacing: scale factor for each dimension
    """
    vec_p = (p - center) * spacing
    vec_ch = [(coord - center) * spacing for coord in ch]
    cos_ch = [np.sum(vec_p * vec) / (np.linalg.norm(vec_p) * np.linalg.norm(vec)) for vec in vec_ch]
    return [*map(lambda x: math.acos(x) * 180 / math.pi, cos_ch)]


def get_anchors(morph: Morphology, ind: list, dist: float):
    """
    get anchors for a set of swc nodes to calculate angles, suppose they are one, their center is their mean coordinate,
    getting anchors requires removing redundant protrudes
    :param dist: path distance thr
    :param ind: array of coordinates
    """
    coords = [morph.pos_dict[i][2:5] for i in ind]
    center = np.mean(coords, axis=0)
    protrude = set()
    for i in ind:
        if morph.child_dict.get(i) is not None:
            protrude = protrude.union(morph.child_dict[i])
    com_line = None
    for i in ind:
        line = [i]
        while line[-1] != -1:
            line.append(morph.pos_dict[line[-1]][-1])
        line.reverse()
        protrude = protrude - set(line[:-1])
        if com_line is None:
            com_line = line
        else:
            for j, t in enumerate(com_line):
                if t != line[j]:
                    com_line = com_line[:j]
                    break
    # nearest common parent of all indices
    com_node = com_line[-1]
    protrude = list(protrude)
    anchor_p = find_point_by_distance(center, com_node, True, morph, dist, False, stop_by_branch=False)
    anchor_ch = [find_point_by_distance(center, i, False, morph, dist, False) for i in protrude]
    return center, anchor_p, anchor_ch, protrude, com_node


def my_find_point_by_distance(pt, anchor_idx, is_parent, morph, dist, return_center_point=True, epsilon=1e-7):
    """
    Find the point of exact `dist` to the start pt on tree structure. args are:
    - pt: the start point, [coordinate]
    - anchor_idx: the first node on swc tree to trace, first child or parent node
    - is_parent: whether the anchor_idx is the parent of `pt`, otherwise child.
                 if an furcation points encounted, then break
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
                if len(morph.child_dict[anchor_idx]) > 1:
                    break
            else:
                if (anchor_idx not in morph.child_dict) or (len(morph.child_dict[anchor_idx]) > 1):
                    break
                else:
                    anchor_idx = morph.child_dict[anchor_idx][0]

    # interpolate to find the exact point
    dd = d - dist
    if dd < 0:
        pt_a = cc
    else:
        dcur = np.linalg.norm(cc - ci)
        assert (dcur - dd >= 0)
        pt_a = ci + (cc - ci) * (dcur - dd) / (dcur + epsilon)
        pts.append(pt_a)

    if return_center_point:
        pt_a = np.mean(pts, axis=0)

    return pts


def prune(morph: Morphology, ind_set: set):
    print('to remove:', ind_set)
    tree = morph.tree.copy()
    for i in ind_set:
        q = []
        ind = morph.index_dict[i]
        if tree[ind] is None:
            continue
        tree[ind] = None
        if morph.child_dict.get(i) is not None:
            q.extend(morph.child_dict[i])
        while len(q) > 0:
            head = q.pop(0)
            ind = morph.index_dict[head]
            if tree[ind] is None:
                continue
            tree[ind] = None
            if morph.child_dict.get(head) is not None:
                q.extend(morph.child_dict[head])
    return filter(lambda x: x is not None, tree)


def spl_gray(morph, ind, parent, spl_dist):
    pts = find_point_by_distance(morph.pos_dict[ind][2:5], morph.pos_dict[ind][6], parent, morph, spl_dist, False)
    if len(pts) == 1:
        new_pts = np.array(pts)
    else:
        cs = CubicSpline(np.arange(len(pts)), np.array(pts))
        xs = np.arange(0, len(pts), len(pts) / spl_sample)
        new_pts = cs(xs)
    return np.median([img[0, p[2], p[1], p[0]] for p in new_pts])


if __name__ == '__main__':
    # test_dir = 'D:/PengLab/200k_testdata'
    test_dir = r'C:\Users\Zuohan Zhao\Desktop\test'
    # swc_dir = os.path.join(test_dir, 'Img_X_5922.13_Y_7785.15_Z_1800.97.swc')
    swc_dir = os.path.join(test_dir, 'Img_X_6840.64_Y_11892.9_Z_3605.37.swc')
    # img_dir = os.path.join(test_dir, 'Img_X_5922.13_Y_7785.15_Z_1800.97.v3draw')
    save_path = os.path.join(test_dir, 'out.swc')
    marker_path = os.path.join(test_dir, 'out.marker')
    tree = swc_handler.parse_swc(swc_dir)
    lower, upper = 10, 1000
    angle_thr = 90

    if len(tree) < lower or len(tree) > upper:
        exit()

    # img = load_image(img_dir)
    morph = Morphology(tree)

    # prune crossing first
    cf = CrossingFinder(morph, 15, 2)
    tri, db = cf.find_crossing_pairs()
    rm_ind = set()
    spline_dist = 15
    spl_sample = 10

    for i in tri:
        # angle
        center, anchor_p, anchor_ch, protrude, com_node = get_anchors(morph, [i], spline_dist)
        angles = anchor_angles(center, np.array(anchor_p), np.array(anchor_ch), spacing=(1, 1, 4))
        angles = np.array(angles)
        to = np.argmin(abs(angles - 180))
        rm_ind |= (set(protrude) - {protrude[to]})

        # # intensity
        # gray_par = spl_gray(morph, i, True, spline_dist)
        # gray_ch = []
        # for j in morph.child_dict[i]:
        #     gray_ch.append(spl_gray(morph, j, False, spline_dist))

    for down, up, dist in db:
        # angle
        center, anchor_p, anchor_ch, protrude, com_node = get_anchors(morph, [up, down], spline_dist)
        angles = anchor_angles(center, np.array(anchor_p), np.array(anchor_ch), spacing=(1, 1, 4))
        angles = np.array(angles)
        to = np.argmin(abs(angles - 180))
        j = protrude[to]
        jj = morph.pos_dict[j][6]
        while j != com_node:
            for k in morph.child_dict[jj]:
                if k != j:
                    rm_ind.add(k)
            j = jj
            jj = morph.pos_dict[j][6]

    awry_angle = set()
    for n, l in morph.pos_dict.items():
        if morph.child_dict.get(n) is not None and len(morph.child_dict[n]) > 1 and morph.pos_dict[n][6] != -1:
            center, anchor_p, anchor_ch, protrude, com_node = get_anchors(morph, [n], spline_dist)
            angles = anchor_angles(center, np.array(anchor_p), np.array(anchor_ch), spacing=(1, 1, 4))
            angles = np.array(angles)
            protrude = np.array(protrude)
            awry_angle |= set(protrude[angles < angle_thr])

    rm_ind |= awry_angle

    nt = prune(morph, rm_ind)
    swc_handler.write_swc(nt, save_path)
    markers = [morph.pos_dict[i][2:5] for i in tri]
    markers.extend(morph.pos_dict[up][2:5] for down, up, dist in db)
    markers.extend([morph.pos_dict[i][2:5] for i in awry_angle])
    save_markers(marker_path, markers)