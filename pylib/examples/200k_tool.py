import sys

sys.path.append(sys.path[0] + "/..")
import os
import glob
import math
from itertools import compress

import numpy as np
from scipy.interpolate import pchip_interpolate, interp1d
import fire
from concurrent.futures import ProcessPoolExecutor
from scipy.stats import ttest_ind
from sklearn.cluster import DBSCAN

from file_io import load_image
import swc_handler
from neuron_quality.find_break_crossing import CrossingFinder
from morph_topo.morphology import Morphology
from neuron_quality.find_break_crossing import find_point_by_distance


class HidePrint:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def anchor_angles(center: np.ndarray, p: np.ndarray, ch, spacing=(1, 1, 1)):
    """
    modified from Jingzhou's code
    :param center: coordinate for the multifurcation center
    :param p: coordinate for the parent anchor
    :param ch: coordinates for the furcation anchors
    :param spacing: scale factor for each dimension
    """
    vec_p = (p - center) * spacing
    vec_ch = [(coord - center) * spacing for coord in ch]
    cos_ch = [vec_p.dot(vec) / (np.linalg.norm(vec_p) * np.linalg.norm(vec)) for vec in vec_ch]
    out = [*map(lambda x: math.acos(max(min(x, 1), -1)) * 180 / math.pi, cos_ch)]
    return out


def get_anchors(morph: Morphology, ind: list, dist: float, spacing=(1, 1, 1), epsilon=1e-7):
    """
    get anchors for a set of swc nodes to calculate angles, suppose they are one, their center is their mean coordinate,
    getting anchors requires removing redundant protrudes
    :param dist: path distance thr
    :param ind: array of coordinates
    """
    center = np.mean([morph.pos_dict[i][2:5] for i in ind], axis=0)
    center_radius = np.mean([morph.pos_dict[i][5] for i in ind])
    protrude = set()
    for i in ind:
        if i in morph.child_dict:
            protrude |= set(morph.child_dict[i])
    com_line = None
    for i in ind:
        line = [i]
        while line[-1] != -1:
            line.append(morph.pos_dict[line[-1]][6])
        line.reverse()
        protrude -= set(line[1:])
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
    # com_node == center can cause problem for spline
    # for finding anchor_p, you must input sth different from the center to get the right pt list
    if np.linalg.norm((center - morph.pos_dict[com_node][2:5]) * spacing) <= epsilon:
        p = morph.pos_dict[com_node][6]
        # p can be -1 if the com_node is root
        # but when this happens, com_node can hardly == center
        # this case is dispelled when finding crossings
    else:
        p = com_node
    anchor_p, pts_p, rad_p = find_point_by_distance(center, p, True, morph, dist, False, spacing=spacing,
                                                    stop_by_branch=False, only_tgt_pt=False, radius=True,
                                                    pt_rad=center_radius)
    res = [find_point_by_distance(center, i, False, morph, dist, False, spacing=spacing,
                                  stop_by_branch=False, only_tgt_pt=False, radius=True, pt_rad=center_radius)
           for i in protrude]
    anchor_ch, pts_ch, rad_ch = [i[0] for i in res], [i[1] for i in res], [i[2] for i in res]
    return center, anchor_p, anchor_ch, protrude, com_node, pts_p, pts_ch, rad_p, rad_ch


def gray_sampling(pts: list, img: np.ndarray, sampling=10, pix_win_radius=1, spacing=(1, 1, 1)):
    """
    interpolate based on the coordinated point list given,
    aggregate the grayscale of the points in the image,
    can output either a median or the list of grayscale data
    there are at least 2 nodes in the input list
    sampling == 0 means no interpolation
    """
    pts = np.array(pts)
    if sampling > 0:
        dist = [np.linalg.norm((pts[i] - pts[i - 1]) * spacing) for i in range(1, len(pts))]
        dist.insert(0, 0)
        dist_cum = np.cumsum(dist)
        f = interp1d(dist_cum, np.array(pts), axis=0)
        pts = f(np.arange(dist_cum[-1] / sampling, dist_cum[-1], dist_cum[-1] / sampling).clip(0, dist_cum[-1]))
        # pts = pchip_interpolate(dist_cum, np.array(pts),
        #                         np.arange(dist_cum[-1] / sampling, dist_cum[-1], dist_cum[-1] / sampling))
    start = (pts.round() - pix_win_radius).astype(int).clip(0)
    end = (pts.round() + pix_win_radius + 1).astype(int).clip(None, np.array(img.shape[-1:-4:-1]) - 1)
    gray = [img[s[2]: e[2], s[1]: e[1], s[0]: e[0]].max()
            for s, e in zip(start, end) if (e - s > 0).all()]
    return gray


def radius_sampling(pts: list, rads: list, sampling=10, spacing=(1, 1, 1)):
    """
    sampling == 0 means returning the original radius list
    """
    if sampling > 0:
        dist = [np.linalg.norm((pts[i] - pts[i - 1]) * spacing) for i in range(1, len(pts))]
        dist.insert(0, 0)
        dist_cum = np.cumsum(dist)
        f = interp1d(dist_cum, rads)
        rads = f(np.arange(dist_cum[-1] / sampling, dist_cum[-1], dist_cum[-1] / sampling).clip(0, dist_cum[-1]))
        # rads = pchip_interpolate(dist_cum, rads,
        #                          np.arange(dist_cum[-1] / sampling, dist_cum[-1], dist_cum[-1] / sampling))
    return rads


def crossing_prune(args):
    """
    detect crossings and prune
    """
    swc_path, img_path, save_path, anchor_dist, soma_radius, dist_thr, spacing, sampling, pix_win_radius = args
    try:
        tree = swc_handler.parse_swc(swc_path)
        morph = Morphology(tree)
        cf = CrossingFinder(morph, soma_radius, dist_thr, spacing)
        img = load_image(img_path)
        if len(img.shape) == 4:
            img = img[0]
        with HidePrint():
            tri, db = cf.find_crossing_pairs()
        rm_ind = set()

        def scoring(center, anchor_p, anchor_ch, protrude, com_node, pts_p, pts_ch, rad_p: list, rad_ch: list):
            angles = anchor_angles(center, np.array(anchor_p), np.array(anchor_ch), spacing=spacing)
            angle_diff = [abs(a - 180) + 1 for a in angles]
            gray_p_med = np.median(gray_sampling(pts_p, img, sampling, pix_win_radius, spacing=spacing))
            gray_ch_med = [np.median(gray_sampling(pts, img, sampling, pix_win_radius, spacing=spacing)) for pts in
                           pts_ch]
            gray_diff = [abs(g - gray_p_med) + 1 for g in gray_ch_med]
            radius_p_med = np.median(radius_sampling(pts_p, rad_p, sampling, spacing=spacing))
            radius_ch_med = [np.median(radius_sampling(pts, rad, sampling, spacing=spacing))
                             for pts, rad in zip(pts_ch, rad_ch)]
            radius_diff = [abs(r - radius_p_med) + 1 for r in radius_ch_med]
            # to = np.argmin(abs(angles - 180))     # used to pick the smallest angle difference for them
            score = [a * g * r for a, g, r in zip(angle_diff, gray_diff, radius_diff)]
            return com_node, protrude, score

        for i in tri:
            # angle
            with HidePrint():
                com_node, protrude, score = scoring(*get_anchors(morph, [i], anchor_dist, spacing))
            rm_ind |= (set(protrude) - {protrude[np.argmin(score)]})
        for down, up, dist in db:
            # angle
            with HidePrint():
                com_node, protrude, score = scoring(*get_anchors(morph, [up, down], anchor_dist, spacing))
            j = protrude[np.argmin(score)]
            jj = morph.pos_dict[j][6]
            while j != com_node:
                for k in morph.child_dict[jj]:
                    if k != j:
                        rm_ind.add(k)
                j = jj
                jj = morph.pos_dict[j][6]
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        swc_handler.write_swc(swc_handler.prune(morph.tree, rm_ind), save_path)
        return True
    except Exception as e:
        print(str(e))
        print("The above error occurred in crossing prune. Proceed anyway: " + swc_path)
        return False


def branch_prune(args):
    """
    prune any awry child for every branch
    """
    swc_path, img_path, save_path, angle_thr, gray_pvalue, radius_pvalue, anchor_dist, \
    soma_radius, spacing, sampling, pix_win_radius = args
    try:
        tree = swc_handler.parse_swc(swc_path)
        morph = Morphology(tree)
        img = load_image(img_path)
        if len(img.shape) == 4:
            img = img[0]
        awry_node = set()
        cs = np.array(morph.pos_dict[morph.idx_soma][2:5])
        for n in morph.bifurcation | morph.multifurcation:
            if np.linalg.norm((morph.pos_dict[n][2:5] - cs) * spacing) <= soma_radius:
                continue
            # branch, no soma, away from soma
            with HidePrint():
                center, anchor_p, anchor_ch, protrude, com_node, pts_p, pts_ch, rad_p, rad_ch = \
                    get_anchors(morph, [n], anchor_dist, spacing)
            angles = anchor_angles(center, np.array(anchor_p), np.array(anchor_ch), spacing=spacing)
            gray_p = gray_sampling(pts_p, img, sampling, pix_win_radius, spacing=spacing)
            gray_ch = [gray_sampling(pts, img, sampling, pix_win_radius, spacing=spacing) for pts in pts_ch]
            radius_p = radius_sampling(pts_p, rad_p, sampling, spacing=spacing)
            radius_ch = [radius_sampling(pts, rad, sampling, spacing=spacing) for pts, rad in zip(pts_ch, rad_ch)]
            gray_pv = [ttest_ind(gray_p, ch, equal_var=False, random_state=1, alternative='less').pvalue
                       for ch in gray_ch]
            radius_pv = [ttest_ind(radius_p, ch, equal_var=False, random_state=1, alternative='less').pvalue
                         for ch in radius_ch]
            protrude = np.array(protrude)
            awry_node |= set(protrude[(np.array(angles) < angle_thr)
                                      | (np.array(gray_pv) < gray_pvalue)
                                      | (np.array(radius_pv) < radius_pvalue)])
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        swc_handler.write_swc(swc_handler.prune(morph.tree, awry_node), save_path)
        return True
    except Exception as e:
        print(str(e))
        print("The above error occurred in branch prune. Proceed anyway: " + swc_path)
        return False


def soma_limit_filter_radius(args):
    """
    limit the number of soma in an swc, based on radius of nodes
    """
    swc_path, soma_radius, max_count, min_radius, pass_rate, min_radius_keep, eps, spacing = args
    try:
        tree = swc_handler.parse_swc(swc_path)
        morph = Morphology(tree)
        dist = morph.get_distances_to_soma(spacing)
        soma_r = np.max([t[5] for t, d in zip(tree, dist) if d <= soma_radius])
        if soma_r < min_radius:
            return min_radius_keep
        pass_r = soma_r * pass_rate
        db = DBSCAN(eps=eps, min_samples=1)
        dt = [t[2:5] for t in tree if t[5] >= pass_r] * np.array(spacing)
        labs = db.fit_predict(dt)
        cluster = [dt[np.where(labs == l)].mean(axis=0) for l in np.unique(labs) if l != -1]
        ct = morph.pos_dict[morph.idx_soma][2:5] * np.array(spacing)
        outer_soma = [p for p in cluster if np.linalg.norm(p - ct) > soma_radius]
        return len(outer_soma) <= max_count - 1
    except Exception as e:
        print(str(e))
        print("The above error occurred in soma limit filter radius. Proceed anyway and leave this swc out: " + swc_path)
        return False


def soma_limit_filter_gray(args):
    """
    limit the number of soma in an swc
    """
    swc_path, img_path, soma_radius, min_radius, max_count, win_radius, pass_rate, eps, spacing = args
    try:
        tree = swc_handler.parse_swc(swc_path)
        img = load_image(img_path)
        if len(img.shape) == 4:
            img = img[0]
        morph = Morphology(tree)

        def stat(ct):
            ind = [np.clip([a - b, a + b + 1], 0, c - 1) for a, b, c in zip(ct, win_radius, img.shape[-1:-4:-1])]
            ind = np.array(ind, dtype=int)
            return img[ind[2][0]:ind[2][1], ind[1][0]:ind[1][1], ind[0][0]:ind[0][1]].flatten()

        ct = morph.pos_dict[morph.idx_soma][2:5]
        soma_stat = stat(ct)
        if soma_stat.size <= np.dot(win_radius, (1, 1, 1)):
            return False
        ct *= np.array(spacing)
        pts = morph.bifurcation | morph.multifurcation | set(t[0] for t in tree if t[5] >= min_radius)
        pts = np.array([morph.pos_dict[i][2:5] * np.array(spacing) for i in pts])
        if pts.size < 1:
            return False
        db = DBSCAN(eps=eps, min_samples=1)
        labs = db.fit_predict(pts)
        max_count -= 1
        cluster = [pts[np.where(labs == i)] for i in np.unique(labs) if i != -1]
        cluster.sort(key=len, reverse=True)
        for i in cluster:
            if max_count < 0:
                break
            c = i.mean(axis=0)
            if np.linalg.norm(c - ct) <= soma_radius:
                continue
            t_stat = stat(c / spacing)
            if t_stat.size <= np.dot(win_radius, (1, 1, 1)):
                continue
            if t_stat.mean() >= pass_rate * soma_stat.mean():
                max_count -= 1
        return max_count >= 0
    except Exception as e:
        print(str(e))
        print("The above error occurred in soma limit filter gray. Proceed anyway and leave this swc out: " + swc_path)
        return False


def node_limit_filter(args):
    path, downer, upper = args
    try:
        tree = swc_handler.parse_swc(path)
        return len(tree) >= downer and (upper < 0 or len(tree) <= upper)
    except Exception as e:
        print(str(e))
        print("The above error occurred in node limit filter. Proceed anyway and leave this swc out: " + path)
        return False


def save_path_gen(paths: list, in_dir: str, out_dir: str, suffix: str):
    out = [(f.split('.swc' if f.endswith('.swc') else '.eswc')[0]) + suffix + '.swc' for f in paths]
    if out_dir is None:
        if suffix == "":
            print('warning: replacing the input or the intermediate files.')
    else:
        out = [os.path.join(out_dir, os.path.relpath(f, in_dir)) for f in out]
    return out


def remove_disconnected(args):
    swc_path, save_path, center = args
    try:
        t = swc_handler.parse_swc(swc_path)
        roots = [t for t in t if t[6] == -1]
        min_r = None
        min_dist = None
        for r in roots:
            dist = np.linalg.norm(np.array(r[2:5]) - center)
            if min_r is None or dist < min_dist:
                min_dist = dist
                min_r = r
        t = swc_handler.rm_disconnected(t, min_r[0])
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        swc_handler.write_swc(t, save_path)
        return True
    except Exception as e:
        print(e)
        print("The above error occurred in remove disconnected. Proceed anyway: " + swc_path)
        return False


def soma_prune(args):
    swc_path, img_path, save_path, pass_rate, non_detect_radius, win_radius, \
        spacing, min_radius, anchor_dist, sampling, pix_win_radius, soma_radius = args
    try:
        img = load_image(img_path)
        if len(img.shape) == 4:
            img = img[0]
        tree = swc_handler.parse_swc(swc_path)
        morph = Morphology(tree)
        dist = morph.get_distances_to_soma(spacing)
        soma_r = np.max([t[5] for t, d in zip(tree, dist) if d <= soma_radius])
        pass_r = max(min_radius, soma_r) * pass_rate

        def stat(ct, win_radius):
            ind = [np.clip([a - b, a + b + 1], 0, c - 1) for a, b, c in zip(ct, win_radius, img.shape[-1:-4:-1])]
            ind = np.array(ind, dtype=int)
            return img[ind[2][0]:ind[2][1], ind[1][0]:ind[1][1], ind[0][0]:ind[0][1]].flatten()

        ct = morph.pos_dict[morph.idx_soma][2:5]
        soma_stat = stat(ct, win_radius)
        if soma_stat.size <= np.dot(win_radius, (1, 1, 1)):
            return False
        pass_gray = soma_stat.mean() * pass_rate
        awry_node = [t[0] for t, d in zip(tree, dist)
                     if d > non_detect_radius and (t[5] > pass_r or stat(t[2:5], win_radius).mean() > pass_gray)]
        rm_ind = set()
        cs = np.array(morph.pos_dict[morph.idx_soma][2:5])
        for i in awry_node:
            pts = [morph.pos_dict[i]]
            skip = False
            cs_fake = np.array([pts[0][2:5]])
            while pts[-1][6] != -1:
                if pts[-1][0] in rm_ind:
                    skip = True
                    break
                pts.append(morph.pos_dict[pts[-1][6]])
            if skip or len(pts) < 3:
                continue
            tmp_morph = Morphology(pts)
            gof1 = []
            len1 = []
            gof2 = []
            len2 = []
            gray2 = []
            for j, p in enumerate(pts):
                if j == 0 or j == len(pts) - 1 or np.linalg.norm((p[2:5] - cs) * spacing) <= soma_radius:
                    gof1.append(0)
                    len1.append(0)
                    gof2.append(0)
                    len2.append(0)
                    gray2.append(0)
                    continue
                with HidePrint():
                    center, anchor_p, anchor_ch, protrude, com_node, pts_p, pts_ch, rad_p, rad_ch = \
                        get_anchors(tmp_morph, [p[0]], anchor_dist, spacing)
                gof1.append(anchor_angles(center, cs, np.array(anchor_ch), spacing=spacing)[0])
                len1.append(np.linalg.norm((np.array(p[2:5]) - pts[j - 1][2:5]) * spacing))
                gof2.append(anchor_angles(center, anchor_p, cs_fake, spacing=spacing)[0])
                len2.append(np.linalg.norm((np.array(p[2:5]) - pts[j + 1][2:5]) * spacing))
                gray2.append(np.exp(10 * (1 - stat(p[2:5], [2, 2, 1]).mean() / 255) ** 2))
            cdf1 = [] # soma to fake
            cdf2 = [] # fake to soma
            gof1.reverse()
            len1.reverse()
            gray1 = gray2.copy()
            gray1.reverse()
            for a, b, c, d, e, f in zip(gof1, gof2, len1, len2, gray1, gray2):
                cdf1.append(a * c * e)
                if len(cdf1) > 1:
                    cdf1[-1] += cdf1[-2]
                cdf2.append(b * d * f)
                if len(cdf2) > 1:
                    cdf2[-1] += cdf2[-2]
            cdf1.reverse()
            k = 0
            for a, b in zip(cdf1, cdf2):
                if a < b:
                    break
                k += 1
            if len(morph.child_dict[pts[k][0]]) > 1:
                k -= 1
            rm_ind.add(pts[k][0])

        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        swc_handler.write_swc(swc_handler.prune(morph.tree, rm_ind), save_path)
        return True
    except Exception as e:
        print(str(e))
        print("The above error occurred in soma prune. Proceed anyway: " + swc_path)
        return False


class CLI200k:
    """
    BICCN 200K SWC post-processing CLI tools
    """

    def __init__(self, swc=None, img=None, ano=None, output_dir=None,
                 jobs=1, chunk_size=1000, spacing=(1, 1, 2), soma_radius=15, sampling=10, pix_win_radius=1):
        """
        The script allows 3 input manners for swc and img files, which can work in parallel.

        1. specify swc-files and img-files as 2 lists like [a,b,c], which should have the same length to map 1 by 1.

        2. specify swc-files ad img-files as 2 wildcards. This assumes swc and images share their number as well as
        their path pattern, meaning the file lists are generated and sorted by name.

        3. specify ano with same number of lines of SWCFILE and RAWIMG for everything, and they are matched by
        the order of input, or it can be a wildcard for a set of ano with the same requirement.

        Whichever way, the output paths are inferred based on the swc-files.

        soma-radius: nodes within soma-radius are not pruned
        spacing: swc and image scaling over the dimensions
        jobs: number of processes
        """
        assert type(swc) == type(img)
        self.pix_win_radius = pix_win_radius
        self.sampling = sampling
        self.spacing = spacing
        self.soma_radius = soma_radius
        self.swc_files = []
        self.img_files = []
        self.output_dir = output_dir
        if swc is not None:
            if type(swc) is str:
                self.swc_files.extend(sorted(glob.glob(swc, recursive=True)))
                self.img_files.extend(sorted(glob.glob(img, recursive=True)))
            elif type(swc) is list:
                self.swc_files.extend(swc)
                self.img_files.extend(img)
            else:
                raise "swc-files and img-files should be specified as either lists or strings"
            assert len(self.swc_files) == len(self.img_files) > 0
        if ano is not None:
            if type(ano) is str:
                ano_files = glob.glob(ano, recursive=True)
            elif type(ano) is list:
                ano_files = ano
            else:
                raise "ano should be specified as either a list or a string"
            for a in ano_files:
                with open(a) as f:
                    ano_dir = os.path.dirname(a)
                    lines = [i.rstrip().split('=') for i in f.readlines()]
                    swc = [i[1] if os.path.isabs(i[1]) else os.path.abspath(os.path.join(ano_dir, i[1]))
                           for i in lines if i[0] == "SWCFILE"]
                    img = [i[1] if os.path.isabs(i[1]) else os.path.abspath(os.path.join(ano_dir, i[1]))
                           for i in lines if i[0] == "RAWIMG"]
                    assert len(swc) == len(img) > 0
                    self.swc_files.extend(swc)
                    self.img_files.extend(img)
        assert len(self.swc_files) > 0
        self.root = os.path.commonpath([os.path.dirname(f) for f in self.swc_files])
        if self.output_dir is not None and self.root == "":
            print('warning: no common dir for all the swc, will save to the original folder.')
            self.output_dir = None
        self.jobs = jobs
        self.chunk = chunk_size
        # with ProcessPoolExecutor(max_workers=self.jobs) as pool:
        #     self.trees = list(pool.map(swc_handler.parse_swc, self.swc_files, chunksize=self.chunk))
        print('number of swc: {}'.format(len(self.swc_files)))

    def node_limit_filter(self, downer=1, upper=-1):
        """
        given a closed interval, return a list of swc paths that meet the range
        upper can be negative(default) to disable
        """
        with ProcessPoolExecutor(max_workers=self.jobs) as pool:
            bool_list = list(
                pool.map(node_limit_filter,
                         [(path, downer, upper) for path in self.swc_files],
                         chunksize=self.chunk)
            )
            self.swc_files = list(compress(self.swc_files, bool_list))
            self.img_files = list(compress(self.img_files, bool_list))
            print('swc node filtering finished.')
            print('remaining number of swc: {}'.format(len(self.swc_files)))
        return self

    def soma_limit_filter_radius(self, max_count=3, min_radius=5, pass_rate=0.8, min_radius_keep=False, eps=10, invert=False):
        """
        filter swc to have a limited number of soma, based on soma radius
        """
        with ProcessPoolExecutor(max_workers=self.jobs) as pool:
            bool_list = list(
                pool.map(soma_limit_filter_radius,
                         [(paths, self.soma_radius, max_count, min_radius, pass_rate,
                           min_radius_keep, eps, self.spacing) for paths in self.swc_files], chunksize=self.chunk)
            )
            if invert:
                bool_list = [not b for b in bool_list]
            self.swc_files = list(compress(self.swc_files, bool_list))
            self.img_files = list(compress(self.img_files, bool_list))
            print('swc soma number filtering finished.')
            print("remaining number of swc: {}".format(len(self.swc_files)))
        return self

    def soma_limit_filter_gray(self, min_radius=3, max_count=3, win_radius=(6, 6, 4), pass_rate=0.5, eps=10, invert=False):
        """
        filter swc to have a limited number of somas
        """
        with ProcessPoolExecutor(max_workers=self.jobs) as pool:
            bool_list = list(
                pool.map(soma_limit_filter_gray,
                         [(*paths, self.soma_radius, min_radius, max_count, win_radius, pass_rate, eps, self.spacing)
                          for paths in zip(self.swc_files, self.img_files)], chunksize=self.chunk)
            )
            if invert:
                bool_list = [not b for b in bool_list]
            self.swc_files = list(compress(self.swc_files, bool_list))
            self.img_files = list(compress(self.img_files, bool_list))
            print('swc soma number filtering finished.')
            print("remaining number of swc: {}".format(len(self.swc_files)))
        return self

    def remove_disconnected(self, center=(64,64,32), suffix="_rm_disc"):
        swc_saves = save_path_gen(self.swc_files, self.root, self.output_dir, suffix)
        with ProcessPoolExecutor(max_workers=self.jobs) as pool:
            bool_list = list(
                pool.map(remove_disconnected, [(*paths, center) for paths in zip(self.swc_files, swc_saves)])
            )
        self.swc_files = list(compress(swc_saves, bool_list))
        self.img_files = list(compress(self.img_files, bool_list))
        self.root = self.output_dir
        print('remove disconnected done.')
        return self

    def soma_prune(self, pass_rate=0.5, min_radius=3, anchor_dist=10, win_radius=(6, 6, 4),
                   non_detect_radius=25, suffix="_spruned"):
        swc_saves = save_path_gen(self.swc_files, self.root, self.output_dir, suffix)
        with ProcessPoolExecutor(max_workers=self.jobs) as pool:
            bool_list = list(
                pool.map(soma_prune,
                         [(*paths, pass_rate, non_detect_radius, win_radius, self.spacing, min_radius,
                           anchor_dist, self.sampling, self.pix_win_radius, self.soma_radius)
                          for paths in zip(self.swc_files, self.img_files, swc_saves)],
                         chunksize=self.chunk)
            )
        self.swc_files = list(compress(swc_saves, bool_list))
        self.img_files = list(compress(self.img_files, bool_list))
        self.root = self.output_dir
        print('soma prune done.')
        return self

    def crossing_prune(self, anchor_dist=10, dist_thr=5, suffix="_xpruned"):
        swc_saves = save_path_gen(self.swc_files, self.root, self.output_dir, suffix)
        with ProcessPoolExecutor(max_workers=self.jobs) as pool:
            bool_list = list(
                pool.map(crossing_prune,
                         [(*paths, anchor_dist, self.soma_radius, dist_thr,
                           self.spacing, self.sampling, self.pix_win_radius)
                          for paths in zip(self.swc_files, self.img_files, swc_saves)],
                         chunksize=self.chunk)
            )
        self.swc_files = list(compress(swc_saves, bool_list))
        self.img_files = list(compress(self.img_files, bool_list))
        self.root = self.output_dir
        print('crossing prune done.')
        return self

    def branch_prune(self, angle_thr=90, gray_pvalue=0.01, radius_pvalue=0.01, anchor_dist=10, suffix="_ypruned"):
        swc_saves = save_path_gen(self.swc_files, self.root, self.output_dir, suffix)
        with ProcessPoolExecutor(max_workers=self.jobs) as pool:
            bool_list = list(
                pool.map(branch_prune,
                         [(*paths, angle_thr, gray_pvalue, radius_pvalue, anchor_dist,
                           self.soma_radius, self.spacing, self.sampling, self.pix_win_radius)
                          for paths in zip(self.swc_files, self.img_files, swc_saves)],
                         chunksize=self.chunk)
            )
        self.swc_files = list(compress(swc_saves, bool_list))
        self.img_files = list(compress(self.img_files, bool_list))
        self.root = self.output_dir
        print('branch prune done.')
        return self

    def write_ano(self, ano=""):
        """
        output existing swc & img list to a linker file
        """
        assert ano
        os.makedirs(os.path.dirname(ano), exist_ok=True)
        with open(ano, 'w') as f:
            f.writelines(["SWCFILE=" + f + "\n" for f in self.swc_files])
            f.writelines(["RAWIMG=" + f + "\n" for f in self.img_files])
            print('ano writing finished')
        return self

    def end(self):
        print("Done.")


if __name__ == '__main__':
    fire.Fire(CLI200k)
