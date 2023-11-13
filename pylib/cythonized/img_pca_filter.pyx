import cython
cimport cython
import numpy as np
cimport numpy as np
from numpy.linalg import eigvalsh
from skimage.filters import gaussian


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef np.ndarray[np.float32_t, ndim=3] img_pca_test(np.ndarray[np.uint16_t, ndim=3] img, grid_radius=(2, 8, 8), stride=(2, 8, 8), float sigma=1.):
    cdef:
        int[3] grid
        int[3] gr = grid_radius
        int[3] st = stride
        int[3] dim = [img.shape[0], img.shape[1], img.shape[2]]
        int i, j, k
        int[3] ct = stride

    for i in range(3):
        grid[i] = dim[i] / st[i]
        assert grid[i] > 0

    cdef np.ndarray[np.float32_t, ndim=3] g = np.zeros(grid, dtype=np.float32)
    ct[0] = st[0]
    for i in range(grid[0]):
        ct[1] = st[1]
        for j in range(grid[1]):
            ct[2] = st[2]
            for k in range(grid[2]):
                g[i, j, k] = block_pca_test(img, ct, gr)
                ct[2] += st[2]
            ct[1] += st[1]
        ct[0] += st[0]

    g = gaussian(g, sigma)
    i, j, k = g.shape[:3]
    out = np.empty((i, st[0], j, st[1], k, st[2]), np.float32)
    out[...] = g[:, None, :, None, :, None]
    i, j, k = img.shape[:3]
    return out.reshape(i, j, k)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double block_pca_test(unsigned short[:, :, :] img, int[3] c, int[3] r):
    cdef:
        int i, j, k, count = 1
        int[3] start
        int[3] end
        double[3] r2

    for i in range(3):
        start[i] = max(c[i] - r[i], 0)
        end[i] = min(c[i] + r[i] + 1, img.shape[i])
        count *= end[i] - start[i]
        r2[i] = (r[i] + 1) * (r[i] + 1)

    # center of mass
    cdef:
        double[3] ct = [0, 0, 0]
        double mv, i2, j2, k2, tmpd, s = 0, sss = 0, w
    count = 0
    for i in range(start[0], end[0]):
        for j in range(start[1], end[1]):
            for k in range(start[2], end[2]):
                w = img[i, j, k]
                s += w
                w = w*w*w
                sss += w
                ct[0] += w * i
                ct[1] += w * j
                ct[2] += w * k
                count += 1
    if s > 0:
        ct[0] /= sss
        ct[1] /= sss
        ct[2] /= sss
        mv = s / count
    else:
        return 0

    # pca
    cdef:
        double cc11 = 0, cc12 = 0, cc13 = 0, cc22 = 0, cc23 = 0, cc33 = 0, ns = 0
        double[3] df = [0, 0, 0]
    for i in range(start[0], end[0]):
        df[0] = (i - ct[0]) / r[0]
        for j in range(start[1], end[1]):
            df[1] = (j - ct[1]) / r[1]
            for k in range(start[2], end[2]):
                df[2] = (k - ct[2]) / r[2]
                w = img[i, j, k]
                if w > 0:
                    w = w*w*w
                    cc11 += w * df[0] * df[0]
                    cc12 += w * df[0] * df[1]
                    cc13 += w * df[0] * df[2]
                    cc22 += w * df[1] * df[1]
                    cc23 += w * df[1] * df[2]
                    cc33 += w * df[2] * df[2]
                    ns += w
    cc11 /= ns
    cc12 /= ns
    cc13 /= ns
    cc22 /= ns
    cc23 /= ns
    cc33 /= ns
    cdef np.ndarray[np.float64_t, ndim = 2] cov = np.zeros([3, 3], dtype=np.float)
    cov[0, 0] = cc11
    cov[1, 1] = cc22
    cov[2, 2] = cc33
    cov[0, 1] = cc12
    cov[0, 2] = cc13
    cov[1, 2] = cc23
    try:
        val = eigvalsh(cov, 'U')
        return (val[2] - val[1]) / (val[2] + val[1]) * (val[2] - val[0]) / (val[2] + val[0])
    except:
        return 0
