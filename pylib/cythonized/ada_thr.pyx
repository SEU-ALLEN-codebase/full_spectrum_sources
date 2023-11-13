import cython
cimport cython
import numpy as np
cimport numpy as np


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef np.ndarray[np.uint16_t, ndim=3] adaptive_threshold(np.ndarray[np.uint16_t, ndim=3] img, int h=5, int d=3):
    cdef:
        int[3] dim = [img.shape[0], img.shape[1], img.shape[2]]
        int i, j, k, n, dist, count, m
        np.ndarray[np.uint16_t, ndim=3] out = np.empty(dim, dtype=np.uint16)
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                m = 0
                dist = count = 0
                for n in range(d):
                    dist += h
                    if i - dist >= 0:
                        m += img[i - dist, j, k]
                        count += 1
                    if j - dist >= 0:
                        m += img[i, j - dist, k]
                        count += 1
                    if k - dist >= 0:
                        m += img[i, j, k - dist]
                        count += 1
                    if i + dist < dim[0]:
                        m += img[i + dist, j, k]
                        count += 1
                    if j + dist < dim[1]:
                        m += img[i, j + dist, k]
                        count += 1
                    if k + dist < dim[2]:
                        m += img[i, j, k + dist]
                        count += 1
                out[i, j, k] = <int>max(0, img[i, j, k] - 1. * m / count)
    return out