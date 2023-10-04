import sys
sys.path.append("/PBshare/SEU-ALLEN/Users/zuohan/pylib")

from pathlib import Path
import os
import numpy as np
import pickle
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
from multiprocessing.pool import Pool


def get_matrix(pts):
    pts = np.array(pts)
    # make a grid, points in the grid and nearby grids are considered for distance
    grid_dims = (pts.max(axis=0) / grid_size).astype(int) + 1   # in case grid be smaller than needed, plus 1
    grid = [[[[] for z in range(grid_dims[2])] for y in range(grid_dims[1])] for x in range(grid_dims[0])]
    for k, p in enumerate(pts):
        xyz = (p / grid_size).astype(int)
        grid[xyz[0]][xyz[1]][xyz[2]].append(k)  # only index
    mtx = -np.ones([pts.shape[0]] * 2)
    for x, y, z in np.ndindex(*grid_dims):
        ind = []
        xyz = np.array([x, y, z])
        for xx, yy, zz in np.ndindex(3, 3, 3):
            if (xyz + (xx, yy, zz) > 0).all() and (xyz + (xx, yy, zz) <= grid_dims).all():
                ind.extend(grid[x+xx-1][y+yy-1][z+zz-1])
        if len(ind) == 0:
            continue
        sub = pts[ind]
        t1, t2 = np.mgrid[0:len(ind), 0:len(ind)]
        t1 = [ind[t] for t in t1.ravel()]
        t2 = [ind[t] for t in t2.ravel()]
        mtx[t1, t2] = distance_matrix(sub, sub, p=2).ravel()
    return mtx


if __name__ == '__main__':
    all = Path("/home/vkzohj/data/enhanced_soma_images/all.txt")
    coords_pik='coords.pickle'
    
    # count files
    try:
        with open(coords_pik, 'rb') as fp:
            brains = pickle.load(fp)
    except:
        print('generating coords pickle')
        brains = {}
        with open(all) as fp:
            for line in fp.readlines():
                path = Path(line.rstrip()).relative_to('..')
                path = Path('/home/vkzohj/data') / path
                brain = path.parent.parts[-1]
                if not path.exists():
                    continue
                if not brain in brains:
                    brains[brain] = []
                parts = str(path.stem).split('_')
                xyz = np.array([float(parts[2]), float(parts[4]), float(parts[6])])
                brains[brain].append(xyz)
            with open(coords_pik, 'wb') as fp:
                pickle.dump(brains, fp)
    
    # count distances
    mtx_pik = "mtx.pickle"
    try:
        with open(mtx_pik, 'rb') as fp:
            matrices = pickle.load(fp)
    except:
        print('generating matrix pickle')
        grid_size = (256, 256, 128)
        with Pool(20) as pt:
            matrices = pt.map(get_matrix, brains.values(), 5)
        matrices = dict(zip(brains.keys(), matrices))
        with open(mtx_pik, 'wb') as fp:
            pickle.dump(dict(matrices), fp)
    
    # plot histogram
    plt.figure(figsize=(8, 6))
    for b, mtx in matrices.items():
        up = mtx[np.triu_indices_from(mtx, k=1)]
        plt.hist(up[up > 0], bins=100, alpha=0.2, label=b)
    plt.title("Histograms of Soma-soma Distances for Each Brain")
    plt.xlabel("Distance/(pixel by 2nd res)", size=14)
    plt.ylabel("Count", size=14)
    plt.savefig('distribution.png')
