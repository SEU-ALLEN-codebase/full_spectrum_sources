import sys
sys.path.append("/PBshare/SEU-ALLEN/Users/zuohan/pylib")

from pathlib import Path
import os
import numpy as np
import pickle


if __name__ == '__main__':
    all = Path("../230k/stage2_app2/all.txt")
    
    coords_pik='coords.pickle'
    with open(coords_pik, 'rb') as fp:
        brains = pickle.load(fp)

    mtx_pik = "mtx.pickle"
    with open(mtx_pik, 'rb') as fp:
        matrices = pickle.load(fp)

    b150 = os.listdir('/home/vkzohj/data/soma_images/out/150k_0602')
    
    out='passed.txt'
    with open(out, 'w') as f:
        for b, mtx in matrices.items():
            bo = ((mtx > 0) & (mtx < 256)).sum(axis=0) < 5
            passed = np.array(brains[b])[bo]
            bb = f"150k_0602/{b}" if b in b150 else f"50k_0602/{b}"
            passed = [[str(i[0]).rstrip('0').rstrip('.'), 
                       str(i[1]).rstrip('0').rstrip('.'), 
                       str(i[2]).rstrip('0').rstrip('.')] for i in passed]
            f.writelines([f"out/{bb}/Img_X_{i[0]}_Y_{i[1]}_Z_{i[2]}.v3dpbd\n" for i in passed])