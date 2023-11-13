import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path
import json

import zmq

sys.path.append(sys.path[0] + "/..")
import os
import math

import swc_handler
from morph_topo.morphology import Morphology
import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import ttest_ind
from scipy.spatial import distance_matrix

from neuron_quality.find_break_crossing import CrossingFinder, find_point_by_distance
from file_io import load_image


class HidePrint:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


class Pruner:
    def __init__(self, tree: list, img: np.ndarray, spacing=(1, 1, 2), soma_radius=15, sampling=10, pix_win_radius=1,
                 anchor_dist=10, epsilon=1e-7):
        self.morph = Morphology(tree)
        self.img = img
        self.pix_win_radius = pix_win_radius
        self.sampling = sampling
        self.spacing = spacing
        self.soma_radius = soma_radius
        self.epsilon = epsilon
        self.anchor_dist = anchor_dist

    def anchor_angles(self, center: np.ndarray, p: np.ndarray, ch):
        """
        modified from Jingzhou's code
        :param center: coordinate for the multifurcation center
        :param p: coordinate for the parent anchor
        :param ch: coordinates for the furcation anchors
        :param spacing: scale factor for each dimension
        """
        vec_p = (p - center) * self.spacing
        vec_ch = [(coord - center) * self.spacing for coord in ch]
        cos_ch = [vec_p.dot(vec) / (np.linalg.norm(vec_p) * np.linalg.norm(vec)) for vec in vec_ch]
        out = [*map(lambda x: math.acos(max(min(x, 1), -1)) * 180 / math.pi, cos_ch)]
        return out

    def gray_sampling(self, pts: list):
        """
        interpolate based on the coordinated point list given,
        aggregate the grayscale of the points in the image,
        can output either a median or the list of grayscale data
        there are at least 2 nodes in the input list
        sampling == 0 means no interpolation
        """
        pts = np.array(pts)
        if self.sampling > 0:
            dist = [np.linalg.norm((pts[i] - pts[i - 1]) * self.spacing) for i in range(1, len(pts))]
            dist.insert(0, 0)
            dist_cum = np.cumsum(dist)
            f = interp1d(dist_cum, np.array(pts), axis=0)
            pts = f(np.arange(dist_cum[-1] / self.sampling, dist_cum[-1], dist_cum[-1] / self.sampling).clip(0,
                                                                                                             dist_cum[
                                                                                                                 -1]))
        start = (pts.round() - self.pix_win_radius).astype(int).clip(0)
        end = (pts.round() + self.pix_win_radius + 1).astype(int).clip(None, np.array(self.img.shape[-1:-4:-1]) - 1)
        gray = [self.img[s[2]: e[2], s[1]: e[1], s[0]: e[0]].max()
                for s, e in zip(start, end) if (e - s > 0).all()]
        return gray

    def radius_sampling(self, pts: list, rads: list):
        """
        sampling == 0 means returning the original radius list
        """
        if self.sampling > 0:
            dist = [np.linalg.norm((pts[i] - pts[i - 1]) * self.spacing) for i in range(1, len(pts))]
            dist.insert(0, 0)
            dist_cum = np.cumsum(dist)
            f = interp1d(dist_cum, rads)
            rads = f(np.arange(dist_cum[-1] / self.sampling, dist_cum[-1], dist_cum[-1] / self.sampling).clip(0,
                                                                                                              dist_cum[
                                                                                                                  -1]))
            # rads = pchip_interpolate(dist_cum, rads,
            #                          np.arange(dist_cum[-1] / sampling, dist_cum[-1], dist_cum[-1] / sampling))
        return rads

    def get_anchors(self, morph: Morphology, ind: list):
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
        if np.linalg.norm((center - morph.pos_dict[com_node][2:5]) * self.spacing) <= self.epsilon:
            p = morph.pos_dict[com_node][6]
            # p can be -1 if the com_node is root
            # but when this happens, com_node can hardly == center
            # this case is dispelled when finding crossings
        else:
            p = com_node
        anchor_p, pts_p, rad_p = find_point_by_distance(center, p, True, morph, self.anchor_dist, False,
                                                        spacing=self.spacing,
                                                        stop_by_branch=False, only_tgt_pt=False, radius=True,
                                                        pt_rad=center_radius)
        res = [find_point_by_distance(center, i, False, morph, self.anchor_dist, False, spacing=self.spacing,
                                      stop_by_branch=False, only_tgt_pt=False, radius=True, pt_rad=center_radius)
               for i in protrude]
        anchor_ch, pts_ch, rad_ch = [i[0] for i in res], [i[1] for i in res], [i[2] for i in res]
        return center, anchor_p, anchor_ch, protrude, com_node, pts_p, pts_ch, rad_p, rad_ch

    def scoring(self, center, anchor_p, anchor_ch, protrude, com_node, pts_p, pts_ch, rad_p: list, rad_ch: list):
        angles = self.anchor_angles(center, np.array(anchor_p), np.array(anchor_ch))
        angle_diff = [abs(a - 180) + 1 for a in angles]
        gray_p_med = np.median(self.gray_sampling(pts_p))
        gray_ch_med = [np.median(self.gray_sampling(pts)) for pts in
                       pts_ch]
        gray_diff = [abs(g - gray_p_med) + 1 for g in gray_ch_med]
        radius_p_med = np.median(self.radius_sampling(pts_p, rad_p))
        radius_ch_med = [np.median(self.radius_sampling(pts, rad)) for pts, rad in zip(pts_ch, rad_ch)]
        radius_diff = [abs(r - radius_p_med) + 1 for r in radius_ch_med]
        # to = np.argmin(abs(angles - 180))     # used to pick the smallest angle difference for them
        score = [a * g * r for a, g, r in zip(angle_diff, gray_diff, radius_diff)]
        return com_node, protrude, score

    def crossing_prune(self, dist_thr=5):
        """
        detect crossings and prune
        """
        cf = CrossingFinder(self.morph, self.soma_radius, dist_thr)
        with HidePrint():
            tri, db = cf.find_crossing_pairs()
        rm_ind = set()
        for i in tri:
            # angle
            with HidePrint():
                com_node, protrude, score = self.scoring(*self.get_anchors(self.morph, [i]))
            rm_ind |= (set(protrude) - {protrude[np.argmin(score)]})
        for down, up, dist in db:
            # angle
            with HidePrint():
                com_node, protrude, score = self.scoring(*self.get_anchors(self.morph, [up, down]))
            j = protrude[np.argmin(score)]
            jj = self.morph.pos_dict[j][6]
            while j != com_node:
                for k in self.morph.child_dict[jj]:
                    if k != j:
                        rm_ind.add(k)
                j = jj
                jj = self.morph.pos_dict[j][6]
        tree = swc_handler.prune(self.morph.tree, rm_ind)
        self.morph = Morphology(tree)

    def branch_prune(self, angle_thr=90, gray_pvalue=0.01, radius_pvalue=0.01):
        """
        prune any awry child for every branch
        """
        rm_ind = set()
        cs = np.array(self.morph.pos_dict[self.morph.idx_soma][2:5])
        for n in self.morph.bifurcation | self.morph.multifurcation:
            if np.linalg.norm((self.morph.pos_dict[n][2:5] - cs) * self.spacing) <= self.soma_radius:
                continue
            # branch, no soma, away from soma
            with HidePrint():
                center, anchor_p, anchor_ch, protrude, com_node, pts_p, pts_ch, rad_p, rad_ch = \
                    self.get_anchors(self.morph, [n])
            angles = self.anchor_angles(center, np.array(anchor_p), np.array(anchor_ch))
            gray_p = self.gray_sampling(pts_p)
            gray_ch = [self.gray_sampling(pts) for pts in pts_ch]
            radius_p = self.radius_sampling(pts_p, rad_p)
            radius_ch = [self.radius_sampling(pts, rad) for pts, rad in zip(pts_ch, rad_ch)]
            gray_pv = [ttest_ind(gray_p, ch, equal_var=False, random_state=1, alternative='less').pvalue
                       for ch in gray_ch]
            radius_pv = [ttest_ind(radius_p, ch, equal_var=False, random_state=1, alternative='less').pvalue
                         for ch in radius_ch]
            protrude = np.array(protrude)
            rm_ind |= set(protrude[(np.array(angles) < angle_thr)
                                   | (np.array(gray_pv) < gray_pvalue)
                                   | (np.array(radius_pv) < radius_pvalue)])
        self.morph = Morphology(swc_handler.prune(self.morph.tree, rm_ind))

    def sub_block(self, ct, win_radius):
        ind = [np.clip([a - b, a + b + 1], 0, c - 1) for a, b, c in zip(ct, win_radius, self.img.shape[-1:-4:-1])]
        ind = np.array(ind, dtype=int)
        return self.img[ind[2][0]:ind[2][1], ind[1][0]:ind[1][1], ind[0][0]:ind[0][1]].flatten()

    def soma_prune(self, pass_rate=0.5, non_detect_radius=25, win_radius=(6, 6, 4), min_radius=3, min_gray=0):
        dist = self.morph.get_distances_to_soma(self.spacing)
        soma_r = np.max([t[5] for t, d in zip(self.morph.tree, dist) if d <= self.soma_radius])
        pass_r = max(min_radius, soma_r * pass_rate)
        ct = self.morph.pos_dict[self.morph.idx_soma][2:5]
        soma_stat = self.sub_block(ct, win_radius)
        if soma_stat.size <= np.dot(win_radius, (1, 1, 1)):
            return False
        pass_gray = max(min_gray, soma_stat.mean() * pass_rate)
        awry_node = [t[0] for t, d in zip(self.morph.tree, dist)
                     if
                     d > non_detect_radius and (t[5] > pass_r or self.sub_block(t[2:5], win_radius).mean() > pass_gray)]
        rm_ind = set()
        cs = np.array(self.morph.pos_dict[self.morph.idx_soma][2:5])
        for i in awry_node:
            pts = [self.morph.pos_dict[i]]
            skip = False
            cs_fake = np.array([pts[0][2:5]])
            while pts[-1][6] != -1:
                if pts[-1][0] in rm_ind:
                    skip = True
                    break
                pts.append(self.morph.pos_dict[pts[-1][6]])
            if skip or len(pts) < 3:
                continue
            tmp_morph = Morphology(pts)
            gof1 = []
            len1 = []
            gof2 = []
            len2 = []
            gray2 = []
            for j, p in enumerate(pts):
                if j == 0 or j == len(pts) - 1 or np.linalg.norm((p[2:5] - cs) * self.spacing) <= self.soma_radius:
                    gof1.append(0)
                    len1.append(0)
                    gof2.append(0)
                    len2.append(0)
                    gray2.append(0)
                    continue
                with HidePrint():
                    center, anchor_p, anchor_ch, protrude, com_node, pts_p, pts_ch, rad_p, rad_ch = \
                        self.get_anchors(tmp_morph, [p[0]])
                gof1.append(self.anchor_angles(center, cs, np.array(anchor_ch))[0])
                len1.append(np.linalg.norm((np.array(p[2:5]) - pts[j - 1][2:5]) * self.spacing))
                gof2.append(self.anchor_angles(center, anchor_p, cs_fake)[0])
                len2.append(np.linalg.norm((np.array(p[2:5]) - pts[j + 1][2:5]) * self.spacing))
                gray2.append(np.exp(10 * (1 - self.sub_block(p[2:5], [2, 2, 1]).mean() / 255) ** 2))
            cdf1 = []  # soma to fake
            cdf2 = []  # fake to soma
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
            if len(self.morph.child_dict[pts[k][0]]) > 1:
                k -= 1
            rm_ind.add(pts[k][0])
        self.morph = Morphology(swc_handler.prune(self.morph.tree, rm_ind))
        return pass_r, pass_gray

    def ref_prune(self, tree: list, min_dist=2):
        r1 = [i[5] for i in self.morph.tree]
        pts1 = np.array([i[2:5] for i in self.morph.tree])
        pts2 = np.array([i[2:5] for i in tree])
        if len(pts1) == 0 or len(pts2) == 0:
            return
        dm = distance_matrix(pts1, pts2).min(axis=1)
        node = set(self.morph.tree[i][0] for i in np.argwhere((dm > r1) & (dm > min_dist)).reshape(-1))
        self.morph = Morphology(swc_handler.prune(self.morph.tree, node))


def main(args):
    client_addr, recved = args
    repl = {}
    try:
        tmp_dir = Path(recved['tmp'])
        raw_path = Path(recved['raw'])
        swc_path = Path(recved['swc'])
        dimX = recved['dimX']
        dimY = recved['dimY']
        dimZ = recved['dimZ']
        decoy_path = tmp_dir.joinpath("decoy.marker")
        enh_path = raw_path.with_stem(raw_path.stem + "_enh")
        soma_path = tmp_dir.joinpath("soma_stats.json")

        # image enhancement
        subprocess.check_output(
            f"{xvfb} {vaa3d_path} -x imPreProcess -f im_enhancement -i {raw_path} -o {enh_path}",
            timeout=timeout, shell=True
        )

        raw = load_image(str(raw_path))
        enh = load_image(str(enh_path))

        raw_thr = raw.mean() + raw.std() * 0.5
        enh_thr = enh.mean() + enh.std() * 0.5

        if soma_path.exists():
            with open(soma_path, 'r') as f:
                soma = json.load(f)
        else:
            soma = {"raw_thr": raw_thr}

        raw_thr = 0.5 * (raw_thr + soma['raw_thr'])

        # reconstructions
        with open(decoy_path, 'w') as f:
            f.write("##x,y,z,radius,shape,name,comment,color_r,color_g,color_b\n")
            f.write(f"{dimX / 2},{dimY / 2},{dimZ / 2},1,1,0,0,255,0,0")
        subprocess.check_output(
            f"{xvfb} {vaa3d_path} -x vn2 -f app2 -i {raw_path} -p {decoy_path} 0 {raw_thr} 0 1 0 0 5 0 0",
            timeout=timeout, shell=True
        )
        subprocess.check_output(
            f"{xvfb} {vaa3d_path} -x vn2 -f app2 -i {enh_path} -p {decoy_path} 0 {enh_thr} 0 1 0 0 5 0 0",
            timeout=timeout, shell=True
        )
        subprocess.check_output(
            f"{xvfb} {vaa3d_path} -x neuTube -f neutube_trace -i {raw_path} -p 1",
            timeout=timeout, shell=True
        )
        subprocess.check_output(
            f"{xvfb} {vaa3d_path} -x neuTube -f neutube_trace -i {enh_path} -p 1",
            timeout=timeout, shell=True
        )

        raw_app2 = raw_path.with_stem(raw_path.stem + "*_app2").with_suffix(".swc")
        enh_app2 = enh_path.with_stem(enh_path.stem + "*_app2").with_suffix(".swc")
        raw_neutube = raw_path.with_stem(raw_path.stem + "*_neutube").with_suffix(".swc")
        enh_neutube = enh_path.with_stem(enh_path.stem + "*_neutube").with_suffix(".swc")
        # consensus
        consensus_path = raw_path.with_stem(raw_path.stem + "_consensus").with_suffix(".swc")
        subprocess.check_output(
            f"{xvfb} {vaa3d_path} -x consensus_swc -f consensus_swc "
            f"-i {raw_app2} {enh_app2} {raw_neutube} {enh_neutube} -o {consensus_path} -p 2 3 0",
            timeout=timeout, shell=True
        )

        # post-pruning
        if len(raw.shape) == 4:
            raw = raw[0]
        p = Pruner(swc_handler.parse_swc(str(enh_app2)), raw, soma_radius=15 if soma_path.exists() else 5)
        p.crossing_prune()
        p.branch_prune()
        if not soma_path.exists():
            soma['min_radius'], soma['min_gray'] = p.soma_prune(min_radius=3, min_gray=0)
        else:
            p.soma_prune(non_detect_radius=5, min_radius=soma['min_radius'], min_gray=soma['min_gray'])
        p.ref_prune(swc_handler.parse_swc(str(consensus_path)))
        swc_handler.write_swc(p.morph.tree, swc_path)
        repl['state'] = 1
    except:
        repl['state'] = 0
    socket.send_multipart([client_addr, repl])


if __name__ == '__main__':
    context = zmq.Context()
    socket = context.socket(zmq.ROUTER)
    socket.bind("ipc:///tmp/feeds")
    poller = zmq.Poller()
    poller.register(socket, zmq.POLLIN)
    n_proc = 50
    timeout = 100
    vaa3d_path = "/PBshare/SEU-ALLEN/Users/zuohan/vaa3d_for_neutube/start_vaa3d.sh"
    xvfb = "xvfb-run -a"
    with Pool(n_proc) as pool:
        while True:
            socks = dict(poller.poll(1000))
            if socket in socks and socks[socket] == zmq.POLLIN:
                pool.apply_async(main, socket.recv_multipart())
