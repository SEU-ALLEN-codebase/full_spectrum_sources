#!/usr/bin/env python

#================================================================
#   Copyright (C) 2021 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : image_utility.py
#   Author       : Yufeng Liu
#   Date         : 2021-07-03
#   Description  : 
#
#================================================================
import numpy as np
from file_io import load_image, save_image

def get_mip_image(img3d, axis=0, mode='MAX'):
    if mode == 'MAX':
        img2d = img3d.max(axis=axis)
    elif mode == 'MIN':
        img2d = img3d.min(axis=axis)
    else:
        raise ValueError

    return img2d

def image_histeq(image, number_bins=256):
    # from http://www.janeriksolem.net/histogram-equalization-with-python-and.html
    # get image histogram
    image_histogram, bins = np.histogram(image.flatten(), number_bins, density=True)
    cdf = image_histogram.cumsum() # cumulative distribution function
    cdf = (number_bins-1) * cdf / cdf[-1] # normalize

    # use linear interpolation of cdf to find new pixel values
    image_equalized = np.interp(image.flatten(), bins[:-1], cdf)

    return image_equalized.reshape(image.shape), cdf


class AbastractCropImage:
    def __init__(self):
        pass

