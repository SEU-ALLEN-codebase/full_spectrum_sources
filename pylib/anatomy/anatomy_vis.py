#!/usr/bin/env python

#================================================================
#   Copyright (C) 2023 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : anatomy_vis.py
#   Author       : Yufeng Liu
#   Date         : 2023-02-22
#   Description  : 
#
#================================================================

import numpy as np
from image_utils import get_mip_image

def detect_edges2d(img2d):
    img = img2d.astype(float)
    gx, gy = np.gradient(img)
    edges = (gx*gx + gy*gy) != 0
    return edges

def get_section_boundary(mask, axis=0, c=None, v=255):
    if c == None:
        c = mask.shape[axis] // 2
    if axis == 0:
        section = mask[c]
    elif axis == 1:
        section = mask[:,c]
    elif axis == 2:
        section = mask[:,:,c]
    else:
        raise NotImplementedError(f'Argument axis supports only 0/1/2, but got {axis}!')

    boundary = detect_edges2d(section)
    if v == 1:
        return boundary
    else:
        return boundary.astype(np.uint8) * v

def get_brain_outline2d(mask, axis=0, v=255):
    mask = mask > 0
    mask2d = get_mip_image(mask, axis=axis)
    outline = detect_edges2d(mask2d)
    if v == 1:
        return outline
    else:
        return outline.astype(np.uint8) * v

def get_brain_mask2d(mask, axis=0, v=255):
    mask = mask > 0
    mask2d = get_mip_image(mask, axis=axis)
    if v == 1:
        return mask2d
    else:
        return mask2d.astype(np.uint8) * v

def get_section_boundary_with_outline(mask, axis=0, sectionX=None, v=255, fuse=True):
    '''
    Args are:
    @mask:      3D CCF mask, with each region has unique value, uint8 3d array
    @axis:      section along which axis
    @sectionX:  plane/section coordiate along the axis
    @v:         mask value for edges
    @fuse:      whether to fuse the regional boundaries and brain outline
    '''

    boundary = get_section_boundary(mask, axis=axis, c=sectionX, v=v)
    outline = get_brain_outline2d(mask, axis=axis, v=v)

    if fuse:
        return np.maximum(boundary, outline)
    else:
        return boundary, outline


if __name__ == '__main__':
    from file_io import load_image
    import cv2
    
    mask_file = './resources/annotation_25.nrrd'
    axis = 0
    c = None
    v = 255
    
    mask = load_image(mask_file)
    bo = get_section_boundary_with_outline(mask, axis=axis, sectionX=c, v=v)
    cv2.imwrite(f'boundary_axis{axis}.png', bo)

