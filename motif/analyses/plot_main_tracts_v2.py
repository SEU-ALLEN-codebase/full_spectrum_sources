#!/usr/bin/env python

#================================================================
#   Copyright (C) 2022 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : valid.py
#   Author       : Yufeng Liu
#   Date         : 2022-08-06
#   Description  : 
#
#================================================================
import sys
sys.path.append(r'D:\code\support_lib\pylib')
import os
import glob
import cv2
import numpy as np
import random
from scipy.linalg import norm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec

from sklearn.decomposition import PCA
from sklearn.neighbors import KDTree
from sklearn.cluster import KMeans

from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

from swc_handler import parse_swc


def estimate_radius2d(pts, method=0):
    '''
    :params method:
                    0: mean distance as radius
                    1: median distance as radius
                    2: 75% percentile as radius
                    3: max distance as radius
    '''
    pca = PCA()
    tpts = pca.fit_transform(pts)
    pts2d = tpts[:,:2]
    center = pts2d.mean(axis=0)
    shift = pts2d - center
    lengths = np.linalg.norm(shift, axis=1)
    if method == 0:
        radius = lengths.mean()
    elif method == 1:
        radius = np.median(lengths)
    elif method == 2:
        radius = np.percentile(lengths, 75)
    elif method == 3:
        radius = np.max(lengths)

    return radius

def convert_to_proj_name(class_name):
    cn_split = class_name.split('-')
    ptype = cn_split[0]
    stype = '-'.join(cn_split[1:])
    proj_name = f'{stype} --> {ptype}'
    return proj_name


def truncated_cone(ax, p0, p1, R0, R1, color, alpha=0.5):
    """
    adapted from https://stackoverflow.com/questions/48703275/3d-truncated-cone-in-python
    """
    # vector in direction of axis
    v = p1 - p0
    # find magnitude of vector
    mag = norm(v)
    # unit vector in direction of axis
    v = v / mag
    # make some vector not in the same direction as v
    not_v = np.array([1, 1, 0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])
    # make vector perpendicular to v
    n1 = np.cross(v, not_v)
    # print n1,'\t',norm(n1)
    # normalize n1
    n1 /= norm(n1)
    # make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)
    # surface ranges over t from 0 to length of axis and 0 to 2*pi
    n = 80
    t = np.linspace(0, mag, n)
    theta = np.linspace(0, 2 * np.pi, n)
    # use meshgrid to make 2d arrays
    t, theta = np.meshgrid(t, theta)
    R = np.linspace(R0, R1, n)
    # generate coordinates for surface
    X, Y, Z = [p0[i] + v[i] * t + R *
               np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
    ax.plot_surface(X, Y, Z, color=color, linewidth=0, antialiased=False, alpha=alpha)


def calc_best_viewpoint(pts):
    pca = PCA()
    pca.fit(pts)
    x,y,z = pca.components_[2]
    elev = np.rad2deg(np.arcsin(z))
    azim = np.rad2deg(np.arctan(x/y))
    return elev, azim


if 0:
    def plot_main_tracts(class_name, key='tract.swc'):
        # directory containing the primary axonal tracts (in form of SWC file)
        mpath_dir = '../main_tracts_types'
        show_pts = 200  # Each tract will be sub-sampled by 200 evenly distributed points
        show_instances = 50 # Maximal number of neurons to be plotted
        figname = f'{class_name}_main_tract_vis.png'    # output figname
        scale = 1000    # transform to millimeter

        # intitalize the fig
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111, projection='3d')


        # load all the main paths
        pts = []
        termini = []
        somata = []
        files = list(glob.glob(os.path.join(mpath_dir, f'{class_name}*{key}')))
        for pathfile in files[:show_instances]:
            tree = parse_swc(pathfile)
            coords = np.array([node[2:5] for node in tree][::-1]) # soma to terminal
            coords /= scale

            nnodes = len(coords)
            step = max(nnodes // show_pts, 1)
            coords_sub = coords[::step]
            pts.extend(coords)
            termini.append(coords[-1])
            somata.append(coords[0])

            # plot
            ax.scatter(coords_sub[0,0], coords_sub[0,1], coords_sub[0,2], marker='o', color='r', s=150)
            ax.scatter(coords_sub[-1,0], coords_sub[-1,1], coords_sub[-1,2], marker='^', color='k', s=150)
            ax.plot(coords_sub[:,0], coords_sub[:,1], coords_sub[:,2], alpha=0.5)
        pts = np.array(pts)

        if False:
            # plot the projection cone
            vertice = np.array(somata).mean(axis=0)
            termini = np.array(termini)
            cone_b = termini.mean(axis=0)
            R = estimate_radius2d(termini, method=2)
            print(R)
            truncated_cone(ax, vertice, cone_b, 0, R, 'blue', alpha=0.2)
     
        fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
        ax.set_adjustable("box")
        xmin, ymin, zmin = pts.min(axis=0)
        xmax, ymax, zmax = pts.max(axis=0)
        offset = 0.5
        if class_name == 'TH_core-SMT':
            offset = 0.4
        elif class_name == 'CTX_IT-VISp':
            offset = 0.8

        if class_name.startswith('CTX_ET'):
            best_proj = False
        else:
            best_proj = False

        if best_proj:
            # get the best view orientation
            elev, azim = calc_best_viewpoint(pts)
            ax.view_init(elev, azim)

        # scalebar for better visualization when axis is of
        fontprops = fm.FontProperties(size=35)
        r = np.sqrt((xmax-xmin)**2 + (ymax-ymin)**2 + (zmax-zmin)**2)
        l = 11.05 / r * 0.02    # 11.05 for CTX_ET-SSp-ul, estimated value
        scalebar = AnchoredSizeBar(ax.transData,
                           l, f'1mm', 'lower left', 
                           pad=2,
                           color='k',
                           frameon=False,
                           size_vertical=0.001,
                           label_top=True,
                           fontproperties=fontprops)
        ax.add_artist(scalebar)

        print(r)
        ax.set_xlim(xmin+offset, xmax-offset)
        ax.set_ylim(ymin+offset, ymax-offset)
        ax.set_zlim(zmin+offset, zmax-offset)
        
        #title = convert_to_proj_name(class_name)
        if class_name[:2] == 'CP':
            title = class_name.split('-')[0]
        else:
            title = '-'.join(class_name.split('-')[1:])

        ax.set_title(title, fontsize=50, y=1.0)
        tick_label_size = 18
        #ax.set_xlabel(r'X-coord ({}$\mu$m)'.format(scale), fontsize=tick_label_size)
        #ax.set_ylabel(r'Y-coord ({}$\mu$m)'.format(scale), fontsize=tick_label_size)
        #ax.set_zlabel(r'Z-coord ({}$\mu$m)'.format(scale), fontsize=tick_label_size)
        #ax.tick_params(axis='both', which='major', labelsize=tick_label_size)

        ax.set_axis_off()
        #plt.legend()
        #plt.tight_layout()
        plt.savefig(figname, dpi=200)
        plt.close('all')

        # from rgb to rgba
        img = cv2.imread(figname)
        img_rgba = np.zeros((img.shape[0], img.shape[1], 4), dtype=img.dtype)
        mask = img.sum(axis=-1) != 255*3
        img_rgba[:,:, :3] = img
        img_rgba[mask, 3] = 255
        cv2.imwrite(figname, img_rgba)

        return len(files)

    ctypes = ['CP_GPe-CP', 'CP_SNr-CP', 'CTX_ET-MOp', 'CTX_ET-MOs', 'CTX_ET-RSPv', 'CTX_ET-SSp-bfd', 'CTX_ET-SSp-m', 'CTX_ET-SSp-n', 'CTX_ET-SSp-ul', 'CTX_ET-SSs', 'CTX_IT-MOp', 'CTX_IT-MOs', 'CTX_IT-SSp-bfd', 'CTX_IT-SSp-m', 'CTX_IT-SSp-n', 'CTX_IT-SSs', 'CTX_IT-VISp', 'TH_core-LGd', 'TH_core-MG', 'TH_core-SMT', 'TH_core-VPL', 'TH_core-VPLpc', 'TH_core-VPM', 'TH_matrix-LP', 'TH_matrix-VM']
    #ctypes = ['TH_matrix-LP']
    
    n = 0
    i = 0
    for class_name in ctypes:
        print(f'--> Plotting for {class_name}...')
        n += plot_main_tracts(class_name)
        i += 1
    print(n, i)


if 1:    # radius estimation along main tracts
    # params
    #mpath_dir = '../main_tracts_types'
    mpath_dir = '../CTX_ET-SSp-m-subclasses'
    show_pts = 200
    show_instances = 200

    np.set_printoptions(precision=2)


    # plot path evolution
    def calc_tract_radii1(class_name):
        pathfiles = list(glob.glob(os.path.join(mpath_dir, f'{class_name}*_tract.swc')))
        if show_instances < len(pathfiles):
            pathfiles = random.sample(pathfiles, show_instances)
        print(f'--> {class_name}: {len(pathfiles)}')

        pts = []
        min_npts = show_pts
        for pathfile in pathfiles:
            tree = parse_swc(pathfile)
            coords = np.array([node[2:5] for node in tree][::-1]) # soma to terminal
            
            nnodes = len(coords)
            step = max(nnodes // show_pts, 1)
            coords_sub = coords[::step]

            pts.append(coords_sub)
            if coords_sub.shape[0] < min_npts:
                min_npts = coords_sub.shape[0]
        
        pts = np.array(list([pt[:min_npts] for pt in pts]))
        radii = [estimate_radius2d(pts[:,i], method=2) for i in range(min_npts)]

        return radii

    def calc_tract_radii2(class_name):
        pathfiles = list(glob.glob(os.path.join(mpath_dir, f'{class_name}*_tract.swc')))
        if show_instances < len(pathfiles):
            pathfiles = random.sample(pathfiles, show_instances)
        print(f'--> {class_name}: {len(pathfiles)}')

        paths = []
        for pathfile in pathfiles:
            tree = parse_swc(pathfile)
            coords = np.array([node[2:5] for node in tree][::-1]) # soma to terminal
            paths.append(coords)

        path_select = 'median'
        if path_select == 'median':
            # find the path with median length
            termini = np.array([c[-1] for c in paths])
            tindices = np.array([len(c)-1 for c in paths])
            tm = termini.mean(axis=0)
            tdists = np.linalg.norm(termini - tm, axis=1)
            midx = np.argmin(tdists) #np.argsort(tdists)[(len(tdists) - 1) // 2]
            mpath = paths[midx]
        elif path_select == 'longest':
            lengths = np.array([len(c) for c in paths])
            tindices = np.array([len(c)-1 for c in paths])
            midx = np.argmax(lengths)
            mpath = paths[midx]

        
        # calculate the radius through nearest points matching
        dmins, imins = [], []
        for i, path in enumerate(paths):
            kdtree = KDTree(path, leaf_size=2)
            dmin, imin = kdtree.query(mpath, k=1)
            dmins.append(dmin[:,0])
            imins.append(imin[:,0])
        dmins = np.array(dmins).transpose()
        imins = np.array(imins).transpose()

        radii = []
        for pt, dmin, imin in zip(mpath, dmins, imins):
            tflags = imin != tindices
            if tflags.sum() < len(imin) / 2:
                break
            dists = dmin[tflags]
            dist = np.median(dists)
            radii.append(dist)
        #for dmin, imin in zip(dmins, imins):
        #    cur_pts = []
        #    for i, idx in enumerate(imin):
        #        pt = paths[i][idx]
        #        cur_pts.append(pt)
        #    radius = estimate_radius2d(cur_pts, method=2)
        #    radii.append(radius)

        print(f'{len(radii)} / {len(mpath)}')

        return radii

    def plot_tracts_radii(ctypes, figname, radius_type=2):
        scale = 1000
        fig = plt.figure(figsize=(8,8))

        cn_map = {
            '0_CTX_ET-SSp-m': 'C1',
            '1_CTX_ET-SSp-m': 'C2',
            '2_CTX_ET-SSp-m': 'C3',
        }
        for class_name in ctypes:
            print(class_name)
            if radius_type == 1:
                radii = np.array(calc_tract_radii1(class_name)) / scale
            else:
                radii = np.array(calc_tract_radii2(class_name)) / scale

            if figname == 'CTX_ET-SSp-m':
                label = cn_map[class_name]
            else:
                if class_name[:2] == 'CP':
                    label = class_name.split('-')[0]
                else:
                    label = '-'.join(class_name.split('-')[1:])
            plt.plot(np.linspace(0, 1, len(radii)), radii, label=label, linewidth=3)
        
        #title = convert_to_proj_name(figname)
        title = figname
        
        plt.title(title, fontsize=40, loc='center', y=1.0, pad=-40)
        plt.xticks([])
        
        axis_label_size = 35
        if radius_type == 1:
            plt.yticks([0,1,2,3], fontsize=axis_label_size)
            plt.ylim([0,3.5])
            plt.ylabel(r'Radius (mm)', fontsize=axis_label_size)
        else:
            plt.yticks([0,1,2], fontsize=axis_label_size)
            plt.ylim([0,2])
            plt.ylabel(r'Radius (mm)', fontsize=axis_label_size)
        plt.xlim([0,1])

        plt.xlabel('Normalized tract distance', fontsize=axis_label_size)
        if figname == 'CTX_ET-SSp-m':
            plt.legend(loc='center left', frameon=False, fontsize=25)
        else:
            plt.legend(loc=2, frameon=False, fontsize=25)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['left'].set_linewidth(3)
        plt.gca().spines['bottom'].set_linewidth(3)
        plt.gca().spines['left'].set_alpha(.7)
        plt.gca().spines['bottom'].set_alpha(.7)

        plt.tight_layout()
        plt.savefig(f'{figname}_radii.png', dpi=300)
        plt.close()

    def plot_tracts_radii_all(type_dict, figname, radius_type=2):
        scale = 1000
        fig = plt.figure(figsize=(8,12))
        gs = fig.add_gridspec(6, 4, 
                      left=0.1, right=0.98, bottom=0.05, top=0.98,
                      wspace=0.05, hspace=0.05,
                        )

        cls_list = ['CP', 'CTX_ET', 'CTX_IT', 'TH_core', 'TH_matrix']
        axis_label_size = 20
        ytick_size = 20
        xlabel = 'Main tract'
        ylabel = 'Radius (mm)'
        for cls in cls_list:
            if cls == 'CP':
                ax = fig.add_subplot(gs[:2,1:3])
                ax.set_ylabel(ylabel, fontsize=axis_label_size)
                ax.set_yticks([0,1,2,3])
                ax.tick_params(axis='y', which='major', labelsize=ytick_size)
            elif cls == 'CTX_ET':
                ax = fig.add_subplot(gs[2:4,:2])
                ax.set_ylabel(ylabel, fontsize=axis_label_size)
                ax.set_yticks([0,1,2,3])
                ax.tick_params(axis='y', which='major', labelsize=ytick_size)
            elif cls == 'CTX_IT':
                ax = fig.add_subplot(gs[2:4,2:4])
                ax.set_yticks([])
            elif cls == 'TH_core':
                ax = fig.add_subplot(gs[4:6,:2])
                ax.set_xlabel(xlabel, fontsize=axis_label_size)
                ax.set_ylabel(ylabel, fontsize=axis_label_size)
                ax.set_yticks([0,1,2,3])
                ax.tick_params(axis='y', which='major', labelsize=ytick_size)
            elif cls == 'TH_matrix':
                ax = fig.add_subplot(gs[4:6,2:4])
                ax.set_yticks([])
                ax.set_xlabel(xlabel, fontsize=axis_label_size)

            for class_name in type_dict[cls]:
                if radius_type == 1:
                    radii = np.array(calc_tract_radii1(class_name)) / scale
                else:
                    radii = np.array(calc_tract_radii2(class_name)) / scale

                cname = '-'.join(class_name.split('-')[1:])
                if cname == 'CP':
                    cname = class_name.split('-')[0]
                ax.plot(np.linspace(0, 1, len(radii)), radii, label=cname, linewidth=2)
                ax.legend(loc=2, frameon=False, fontsize=12)

            ax.set_xlim([0,1])
            ax.set_xticks([])
            ax.set_ylim([0,3.5])
            ax.spines.right.set_visible(False)
            ax.spines.top.set_visible(False)
            ax.spines['left'].set_linewidth(3)
            ax.spines['left'].set_alpha(0.3)
            ax.spines['bottom'].set_linewidth(3)
            ax.spines['bottom'].set_alpha(0.3)
            ax.set_title(cls, loc='center', y=1.0, pad=-25, fontsize=25)
        
        #plt.suptitle(figname, fontsize=30)

        plt.savefig(f'{figname}_radii.png', dpi=600)
        plt.close()
    

    # calculating
    '''
    type_dict = {
        'CP': ['CP_GPe-CP', 'CP_SNr-CP'],
        'CTX_ET': ['CTX_ET-MOp', 'CTX_ET-MOs', 'CTX_ET-RSPv', 'CTX_ET-SSp-bfd', 'CTX_ET-SSp-m', 'CTX_ET-SSp-n', 'CTX_ET-SSp-ul', 'CTX_ET-SSs'],
        'CTX_IT': ['CTX_IT-MOp', 'CTX_IT-MOs', 'CTX_IT-SSp-bfd', 'CTX_IT-SSp-m', 'CTX_IT-SSp-n', 'CTX_IT-SSs', 'CTX_IT-VISp'],
        'TH_core': ['TH_core-LGd', 'TH_core-MG', 'TH_core-SMT', 'TH_core-VPL', 'TH_core-VPLpc', 'TH_core-VPM'],
        'TH_matrix': ['TH_matrix-LP', 'TH_matrix-VM']
    }
    '''
    type_dict = {'CTX_ET-SSp-m': ['0_CTX_ET-SSp-m', '1_CTX_ET-SSp-m', '2_CTX_ET-SSp-m']}

    #plot_tracts_radii_all(type_dict, 'all', radius_type=1)
    for key, value in type_dict.items():
        print(key, value)
        plot_tracts_radii(value, key, radius_type=1)
 

if 0:    # clustering and divide neurons into sub-types

    # params
    mpath_dir = '../main_tracts_types'
    out_dir = '../CTX_ET-SSp-m-subclasses'
    show_instances = 200
    scale = 1000

    def plot_CTX_ET_neurons(class_name):
        pathfiles = list(glob.glob(os.path.join(mpath_dir, f'{class_name}*_tract.swc')))
        if show_instances < len(pathfiles):
            pathfiles = random.sample(pathfiles, show_instances)
        print(f'--> {class_name}: {len(pathfiles)}')

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111, projection='3d')
        fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
            hspace = 0, wspace = 0)

        pts = []
        for pathfile in pathfiles:
            tree = parse_swc(pathfile)
            pts.append(tree[0][2:5])
        pts = np.array(pts) / scale
        # standardize
        pts = (pts - pts.mean(axis=0)) / pts.std(axis=0)

        # clustering according to terminal points
        kmeans = KMeans(3).fit(pts)
        labels = kmeans.labels_

        # copy to file
        for pf, l in zip(pathfiles, labels):
            pfile = os.path.split(pf)[-1]
            outfile = os.path.join(out_dir, f'{l}_{pfile}')
            os.system(f'cp {pf} {outfile}')
        
        for i, c in zip(range(3), ['r', 'g', 'b']):
            pts_ = pts[labels == i]
            ax.scatter(pts_[:,0], pts_[:,1], pts_[:,2], marker='^', color=c, s=200)

        ax.set_adjustable("box")
        xmin, ymin, zmin = pts.min(axis=0)
        xmax, ymax, zmax = pts.max(axis=0)
        offset = 0.2
        ax.set_xlim(xmin+offset, xmax-offset)
        ax.set_ylim(ymin+offset, ymax-offset)
        ax.set_zlim(zmin+offset, zmax-offset)

        label_size = 25
        ax.set_xlabel('X (mm)', fontsize=label_size)
        ax.set_ylabel('Y (mm)', fontsize=label_size)
        ax.set_zlabel('Z (mm)', fontsize=label_size)

        ax.tick_params(axis='both', which='minor', labelsize=label_size)

        plt.savefig(f'{class_name}.png', dpi=200)
        plt.close()


    ctx_type = 'CTX_ET-SSp-m'
    plot_CTX_ET_neurons(ctx_type)
    



