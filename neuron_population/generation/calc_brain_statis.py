#!/usr/bin/env python

#================================================================
#   Copyright (C) 2022 Yufeng Liu (Braintell, Southeast University). All rights reserved.
#   
#   Filename     : calc_brain_statis.py
#   Author       : Yufeng Liu
#   Date         : 2022-10-20
#   Description  : 
#
#================================================================

import os
import sys
sys.path.append(r'D:\code\support_lib\pylib')
import glob
import time
import csv
from collections import Counter
from skimage.filters import meijering
from cythonized.img_pca_filter import img_pca_test
from cythonized.ada_thr import adaptive_threshold

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import torch
import torch.nn.functional as F

from file_io import load_image, save_image, get_tera_res_path
from image_utils import get_mip_image



def get_zeng_threshs(thresh_file):
    df = pd.read_csv(thresh_file, skiprows=1)
    brains = df['brain_id'].to_numpy()
    threshs = df[['cur_thresh', 'supposed']].to_numpy()
    thresh_dict = dict(zip(brains, threshs))
    return thresh_dict

fMOST_Zeng_THRESH = get_zeng_threshs(thresh_file='../brain_statistics/statis_out_old/statis_out_adaThr_all/fMOST-Zeng/threshold_merged_20230424_fmost1.txt')


def get_filesize(tera_dir, res_id=-3, outdir=None):
    np.random.seed(1024)

    res_path = get_tera_res_path(tera_dir, res_id, False)
    fsizes = []
    counter = 0
    t0 = time.time()
    for imgfile in glob.glob(os.path.join(res_path, '*/*/[0-9]*[0-9].tif')):
        fs = os.path.getsize(imgfile) / 1000. / 1000.
        fsizes.append(fs)

        # randomly save 2D MIP images for checking and debugging
        output_min_size = 0.4
        output_max_size = 10.0
        prob = 0.2
        if (fs < output_max_size) and (fs > output_min_size):
            if np.random.random() < prob:
                img = load_image(imgfile)
                img2d = get_mip_image(img)
                max_pv = img2d.max()
                min_pv = img2d.min()
                img2d = ((img2d - min_pv) / (max_pv - min_pv + 1e-10) * 255).astype(np.uint8)
                
                fname = os.path.split(imgfile)[-1]
                prefix = os.path.splitext(fname)[0]
                outfile = os.path.join(outdir, f'{prefix}_fs{fs:.3f}_vmax{max_pv}.png')
                save_image(outfile, img2d)

        counter += 1
        if counter % 100 == 0:
            print(f'--> Processed {counter} files in {time.time() - t0:.4f} seconds')

    fsizes = np.array(fsizes)
    sns.histplot(fsizes, bins=200, binrange=(fsizes.min(),4))
    plt.ylim([0,300])
    plt.xlabel('TeraFly block image size')
    plt.savefig(f'filesize_distr.png', dpi=300)

def get_block_counts(zdim, ydim, xdim, h=5, d=3):
    dh = (d+1) * h
    yv, zv, xv = np.meshgrid(np.arange(ydim), np.arange(zdim), np.arange(xdim))
    yc = np.maximum(dh - yv - 1, 0)//h + np.maximum(dh - (ydim - yv), 0)//h
    zc = np.maximum(dh - zv - 1, 0)//h + np.maximum(dh - (zdim - zv), 0)//h
    xc = np.maximum(dh - xv - 1, 0)//h + np.maximum(dh - (xdim - xv), 0)//h
    bc = 2 * d * 3 - xc - yc - zc
    #print(bc.max(), bc.min(), bc.mean(), bc[:11,:11,:11], yc[:11,:11,:11], zc[:11,:11,:11], xc[:11,:11,:11])
    bc = torch.from_numpy(np.expand_dims(np.expand_dims(bc, 0), 0).astype(np.float32))
    return bc

def ada_thresholding(img, block_counts, fg_thresh, h=5, d=3, cuda=True):
    if img.ndim == 4:
        imgt = torch.from_numpy(img.astype(np.float32)).unsqueeze(0)
    elif img.ndim == 3:
        imgt = torch.from_numpy(img.astype(np.float32)).unsqueeze(0).unsqueeze(0)
    k = 2 * d + 1
    weight = torch.zeros((1,1,k,k,k))
    weight[0,0,:,d,d] = 1
    weight[0,0,d,:,d] = 1
    weight[0,0,d,d,:] = 1
    weight[0,0,d,d,d] = 0
    
    if cuda:
        imgt = imgt.cuda()
        weight = weight.cuda()
    conved = F.conv3d(imgt, weight, dilation=h, padding=d*h)
    if (conved.shape[2] != block_counts.shape[2]) or \
        (conved.shape[3] != block_counts.shape[3]) or \
        (conved.shape[4] != block_counts.shape[4]):
        block_counts = get_block_counts(*img.shape[-3:])
        if cuda:
            block_counts = block_counts.cuda()
    conved = conved / block_counts
    diff = imgt - conved
    diff[diff < fg_thresh] = 0
    diff[diff >= fg_thresh] = 1
    # remove isolated digits
    #diff_conv = F.conv3d(diff, weight, padding=d)
    #diff[diff_conv < 3] = 0

    diff = diff.bool()
    
    if cuda:
        diff = diff.cpu()
    diff = diff.numpy()[0][0]

    return diff

class CalcBrainStatis(object):
    def __init__(self, tera_dir, res_id_statis=-3, mip_dir='', cuda=True, fmt='tif'):
        self.tera_dir = tera_dir
        self.res_id_statis = res_id_statis
        assert (res_id_statis < 0)  # -1: the highest resolution, -2: second-highest, -3,...,
        self.multiplier = np.power(2, np.fabs(res_id_statis)-1)
        self.mip_dir = mip_dir
        self.cuda = cuda
        self.fmt = fmt

    def set_region_mask(self, region_mask, max_res_dims, mask_dims):
        """
        max_res_dims & mask_dims: np.ndarray, in order [y, x, z]. Do not get it wrong!
        """
        self.region_mask = region_mask
        self.mask_dims = mask_dims
        # self.coord_factor = max_res_dims / mask_dims / self.multiplier
        self.coord_factor = max_res_dims / mask_dims 

    def set_strides_through_block(self, block):
        """
        dimension of each block, in order [y, x, z]. WARNING: The order is unconventional! 
        """
        zs, ys, xs = block.shape[-3:]  # in case multi-channel images
        self.strides = np.array([ys, xs, zs])   # to numpy array batch-computing convenience
        # initialize the local coordinates of current block
        #zz, yy, xx = np.meshgrid(np.arange(ys), np.arange(zs), np.arange(xs))

    def get_image_range(self, idx=0):
        """
        Pre-estimate the overall statistics of the brain, using the lowest resolution to speed up
        """
        res_path = get_tera_res_path(self.tera_dir, idx, False)
        vmin = 2**16
        vmax = 0
        vmean = []
        vstd = []
        for imgfile in glob.glob(os.path.join(res_path, f'*/*/[0-9]*[0-9].{self.fmt}')):
            img = load_image(imgfile)
            #print(img.dtype)
            vmin_, vmax_ = img.min(), img.max()
            vmean.append(img.mean())
            vstd.append(img.std())
            vmax = max(vmax_, vmax)
            vmin = min(vmin_, vmin)

        vmean = np.median(vmean)
        vstd = np.median(vstd)
        return vmax, vmin, vmean, vstd


    def block_statis(self, block_path, img_bin, start_x=0, start_y=0, start_z=0):
        """
        """
        # infer the coordinate from the block_path
        fname = os.path.splitext(os.path.split(block_path)[-1])[0]
        start_hres = np.array(list(map(int, fname.split('_')))) # y-x-z order
        start_res = start_hres
        # assert int(start_hres[0] / self.multiplier) == start_hres[0] / self.multiplier
        # start_res = start_hres // self.multiplier
        # print(start_hres)
        #img = load_image(block_path)    # Reminder: TeraFly use uint16 by default
        assert img_bin.ndim == 3, "Only 3 dimensional image is supported!"
        # now you should get brain regions from the mask image
        fg_pos = np.nonzero(img_bin)

        pos_mapped_z = np.round((fg_pos[0] * self.multiplier  + start_res[2]/10 + start_z) / self.coord_factor[2]).astype(np.int32)

        pos_mapped_z = np.clip(pos_mapped_z, 0, self.mask_dims[2]-1)
        pos_mapped_y = np.round((fg_pos[1] * self.multiplier  + start_res[0]/10 + start_y) / self.coord_factor[0]).astype(np.int32)
        pos_mapped_y = np.clip(pos_mapped_y, 0, self.mask_dims[0]-1)
        pos_mapped_x = np.round((fg_pos[2] * self.multiplier  + start_res[1]/10 + start_x) / self.coord_factor[1]).astype(np.int32)
        pos_mapped_x = np.clip(pos_mapped_x, 0, self.mask_dims[1]-1)

        regions = self.region_mask[0,pos_mapped_z, pos_mapped_y, pos_mapped_x]
        # then you can use the regions of foreground voxels
        region_counter = Counter(regions)
        return region_counter

    # @staticmethod
    # def detect_anisotropy(image):
        

    def brain_statis(self, filesize_thresh=1.7, vmax_thresh=300, fg_thresh=300, save_mip=True, start_x=0, start_y=0, start_z=0):
        if start_x > 0 and start_y > 0:
            # This is LSFM-Osten data, which miss the sec-highest resolution images
            # This if is very ugly!!!
            res_path = get_tera_res_path(self.tera_dir, res_ids=self.res_id_statis+1, bracket_escape=False)
        else:
            res_path = get_tera_res_path(self.tera_dir, res_ids=self.res_id_statis, bracket_escape=False)

        n_processed = 0
        idx_file = 0
        n_small_block = 0
        n_lowQ_block = 0
        n_highQ_block = 0
        display_freq = 10
        t0 = time.time()
        brain_counter = Counter()
        for block_file in glob.glob(os.path.join(res_path, f'*/*/[0-9]*[0-9].{self.fmt}')):
            n_processed += 1
            if n_processed % display_freq == 0:
                print(f'--> Procssed: {n_processed} in {time.time() - t0:.4f}s. \t\tStatis: small: {n_small_block}, low-quality: {n_lowQ_block}, high-quality: {n_highQ_block}')

            fs = os.path.getsize(block_file) / 1000. / 1000.
            if fs < filesize_thresh:
                n_small_block += 1
                continue

            # filter with vmax
            img = load_image(block_file)
            if img.ndim == 4: img = img[0]
            if idx_file == 0:   # get the dimension of each block, only once!
                self.set_strides_through_block(img)
                # self.bc = get_block_counts(self.strides[2], self.strides[0], self.strides[1])
                # if self.cuda:
                    # self.bc = self.bc.cuda()
                idx_file += 1

            vmax = img.max()
            if vmax < vmax_thresh:
                n_lowQ_block += 1
                continue

            n_highQ_block += 1
            # do pca filter
            # pca block by block, default window radius = 2, 8, 8 (zyx), stride = 2, 8, 8
            # it comes out as a downsampled image, so resample back
            # for better speed, the raw size must be proportional to the stride, so do a padding beforehand
            
            # do ada_thresholding
            # img_a = ada_thresholding(img, self.bc, fg_thresh=fg_thresh, cuda=self.cuda)
            img_a = adaptive_threshold(img.astype(np.uint16))
            img_adtThr = img_a >= fg_thresh
            
            # dim modulo stride, and get remaining distance to the bigger nearest multiple of stride
            piz = 2 - (img.shape[0] & 1) & 1
            piy = 8 - (img.shape[1] & 7) & 7
            pix = 8 - (img.shape[2] & 7) & 7
            arr = np.empty([img.shape[0] + piz, img.shape[1] + piy, img.shape[2] + pix], dtype=img.dtype)
            arr[:img.shape[0], :img.shape[1], :img.shape[2]] = img[...]
            if piz > 0:
                arr[:-piz-1:-1, :, :] = arr[img.shape[0]-piz:img.shape[0]]
            if piy > 0:
                arr[:, :-piy-1:-1, :] = arr[:, img.shape[1]-piy:img.shape[1]]
            if pix > 0:
                arr[..., :-pix-1:-1] = arr[..., img.shape[2]-pix:img.shape[2]]
            
            # do pca on the padded array, only allow unsigned short
            pca = img_pca_test(arr.astype(np.uint16))[:img.shape[0], :img.shape[1], :img.shape[2]]
            # do a rescale based on pca on the array, rescale back and save it to another image
            img_a = (img_a * pca).astype(img_a.dtype)
            thr = 0.1 * fg_thresh
            
            # # get a mask to estimate a new thershold, because pca can cause an overall drop down of intensity
            # # the mask uses a similar approach to the original method, but the threshold is doubled
            # # to make sure not too much background is included
            # mask = adaptive_threshold(img.astype(np.uint16)) >= fg_thresh * 2
            # # the image can be mostly blank, which needs a different thresholding
            # if mask.sum() > mask.size * 0.01:
            #     # usually we threshold based on a mean drop down on the foreground
            #     thr = (img2 / (img + 1))[mask].mean() * fg_thresh
            # else:
            #     # otherwise we directly turn down the threshold, 0.5 here may still be too big
            #     thr = fg_thresh * 0.5
            # img = img2
            
            img_a[img_a < thr] = 0
            img_a[img_a >= thr] = 1
            img_a = img_a.astype(bool)
            cur_counter = self.block_statis(block_file, img_a, start_x=start_x, start_y=start_y, start_z=start_z)
            brain_counter = brain_counter + cur_counter
            # do pca filtering
            # import ipdb; ipdb.set_trace()

            if save_mip and np.random.random() < .1:
                # save mip for inspection
                img2d = get_mip_image(img)
                max_pv = img2d.max()
                min_pv = img2d.min()
                img_thr = get_mip_image(img_a.astype(np.uint8) * 255)
                img_adtThr2d = get_mip_image(img_adtThr).astype(np.uint8) * 255

                # normalize for better visualization
                imga2d = ((img2d - min_pv) / (max_pv - min_pv + 1e-10) * 255).astype(np.uint8)

                # imgc = get_mip_image(imgc); imgc = ((imgc - imgc.min()) / (imgc.max() - imgc.min() + 1e-10) * 255).astype(np.uint8)
                img2d_merge = np.hstack((imga2d, img_thr))
                #img2d_merge = np.hstack((imga2d, img_adtThr2d, img_thr))

                fname = os.path.split(block_file)[-1]
                prefix = os.path.splitext(fname)[0]
                outfile = os.path.join(self.mip_dir, f'{prefix}_fs{fs:.3f}_vmax{max_pv}_thresh{vmax_thresh:.1f}.png')
                save_image(outfile, img2d_merge)

        print(f'<-Summary of current brain: -> small: {n_small_block}, low-quality: {n_lowQ_block}, high-quality: {n_highQ_block}')

        return brain_counter

def brain_statis_wrapper(tera_dir, mask_file_dir, out_dir, max_res_dims, mask_dims, filesize_thresh, brain_id, 
                         res_ids=-3, source='fMOST-Zeng', cuda=True):
    csv_out = os.path.join(out_dir, f'{brain_id}.csv')
    if os.path.exists(csv_out):
        return 

    print(f'===> Processing {brain_id} / {source}')
    mask_file = os.path.join(mask_file_dir, f'{brain_id}.v3draw')
    if source == 'LSFM-Osten':
        if brain_id == 'A2ds221_sitich_terafly_ver4.0':
            start_x, start_y, start_z = 6802, 1715, 0
            mask_file = os.path.join(mask_file_dir, 'AVP-IHC-A2_16um.v3draw')
        elif brain_id == 'A3_Pos_p5um_nv6_dsk_correction_result_9.29ver_terafly':
            start_x, start_y, start_z = 7886, 3561, 0
            mask_file = os.path.join(mask_file_dir, 'A3_16um_final.v3draw')
        
    else:
        start_x, start_y, start_z = 0, 0, 0
    
    print(mask_file)
    mask = load_image(mask_file)

    mip_dir = os.path.join(out_dir, f'mip2d_{brain_id}')
    if not os.path.exists(mip_dir):
        os.mkdir(mip_dir)

    # statistics
    if source == 'LSFM-Wu':
        fmt = 'raw'
    else:
        fmt = 'tif'

    cbs = CalcBrainStatis(tera_dir, mip_dir=mip_dir, cuda=cuda, res_id_statis=res_ids, fmt=fmt)
    cbs.set_region_mask(mask, max_res_dims, mask_dims)
    if source == 'fMOST-Zeng':
        _, _, vmean, vstd = cbs.get_image_range()
        vmax_thresh = min(max(vmean + 1.5 * vstd, 400), 1000)
        fg_thresh = vmax_thresh * 0.9
        
        # use manual assigned threshold according to empirical values
        if brain_id not in fMOST_Zeng_THRESH:
            return
        
        ths = fMOST_Zeng_THRESH[brain_id]
        vmax_thresh = ths[1]
        fg_thresh = vmax_thresh * 0.9


    elif source == 'fMOST-Huang':
        vmax_thresh = 800
        fg_thresh = vmax_thresh * 0.9
    elif source == 'LSFM-Wu':
        vmax_thresh = 1000
        fg_thresh = 600
    elif source == 'LSFM-Dong':
        vmax_thresh = 70
        fg_thresh = 70
    elif source == 'LSFM-Osten':
        vmax_thresh = 800
        fg_thresh = 250

    brain_counter = cbs.brain_statis(filesize_thresh=filesize_thresh, vmax_thresh=vmax_thresh, fg_thresh=fg_thresh, 
                                     start_x=start_x, start_y=start_y, start_z=start_z)
    #print(f'{brain_id}: {cbs.get_image_range()}')
    
    if len(brain_counter) > 0:
        with open(csv_out, 'w', newline='') as f:
            header = ['region','count'] #csv列名
            writer = csv.writer(f)
            writer.writerow(header)
            for key,value in brain_counter.items():
                writer.writerow([key,value])
        



if __name__ == '__main__':
    from multiprocessing.pool import Pool

    tera_downsize_file = '../brain_statistics/ccf_info/TeraDownsampleSize.csv'
    mask_file_dir = '/PBshare/SEU-ALLEN/Users/ZhixiYun/data/registration/Inverse'
    for source in ['LSFM-Wu','LSFM-Osten']:
        out_dir = f'./unzipped/{source}'
        res_ids = -3
        filesize_thresh = 1.7
        nproc = 2
 
        if source == 'fMOST-Zeng':
            match_str = 'mouse[0-9]*'
            tera_path = '/PBshare/TeraconvertedBrain'
            #tera_path = '/PBshare/BrainRaw/Unzipped_Brains'
        elif source == 'fMOST-Huang':
            match_str = 'mouse*[0-9]'
            tera_path = '/PBshare/Huang_Brains'
        elif source == 'STPT-Huang':
            match_str = '[1-9]*processed'
            tera_path = '/PBshare/Huang_Brains'
        elif source == 'LSFM-Wu':
            tera_path = '/PBshare/Zhuhao_Wu'
            match_str = 'WHOLE_mouse_B*'
            res_ids = -1    # use -1 as it have relative low resolution
        elif source == 'LSFM-Osten':
            tera_path = '/PBshare/SEU-ALLEN/Projects/A2_A3_MouseBrain'
            match_str = 'A[2-3]*terafly*' 
            res_ids = -3    # 
        elif source == 'LSFM-Dong':
            tera_path = '/PBshare/DongHW_Brains/20220315_SW220203_03_LS_6x_1000z'
            match_str = '*TeraFly'
            res_ids = -1
        
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
       
        dim_f = pd.read_csv(tera_downsize_file, index_col='ID', sep=',')
        args_list = []
        for tera_dir in glob.glob(os.path.join(tera_path, match_str)):
            if source == 'fMOST-Zeng':
                brain_id = int(os.path.split(tera_dir)[-1].split('_')[0][5:])
            elif source == 'fMOST-Huang':
                if os.path.exists(os.path.join(tera_dir, 'terafly')):
                    tera_dir = os.path.join(tera_dir, 'terafly')
                    brain_folder = os.path.split(os.path.split(tera_dir)[0])[-1]
                else:
                    brain_folder = os.path.split(tera_dir)[-1]
                brain_id = int(brain_folder[-6:])
            elif source == 'LSFM-Wu':
                tera_dir = '/PBshare/Zhuhao_Wu/WHOLE_mouse_Brain_AZ10_SR3B_6_A_647'
                brain_id = os.path.split(tera_dir)[-1]
            elif source == 'LSFM-Dong':
                tera_dir = '/PBshare/DongHW_Brains/20220315_SW220203_03_LS_6x_1000z/Ex_642_Em_680_TeraFly'
                brain_id = 'Ex_642_Em_680'
            elif source == 'LSFM-Osten':
                brain_id = os.path.split(tera_dir)[-1]
            
            
            if source == 'LSFM-Osten':
                if os.path.split(tera_dir)[-1][:2] == 'A2':
                    max_res_dims = np.array([25727, 26235, 19710])
                    mask_dims = np.array([828,824,514])
                    # start_y = 1715, start_x = 6802
                    #Z04-45Y08-20 of Z01-45Y01-27
                else:
                    max_res_dims = np.array([27897, 22781, 18598])
                    mask_dims = np.array([766,797,588])
                    # start_y = 3561, start_x = 7886
                    # Z07-42Y10-19 of Z01-47Y01-26
            else:
                try:
                    max_res_dims = np.array([dim_f.loc[str(brain_id)][0],dim_f.loc[str(brain_id)][1],dim_f.loc[str(brain_id)][2]])
                    mask_dims = np.array([dim_f.loc[str(brain_id)][3],dim_f.loc[str(brain_id)][4],dim_f.loc[str(brain_id)][5]])
                except KeyError:
                    continue
  
            print(brain_id)           
            # if brain_id != 17298: continue
            args = tera_dir, mask_file_dir, out_dir, max_res_dims, mask_dims, filesize_thresh, brain_id, res_ids, source, False
            # brain_statis_wrapper(*args)
            args_list.append(args)
 
        print('\n\n')
        print(args_list)
        print('\n\n')
  
        # sys.exit()
        print(f'Number of brains to process: {len(args_list)}')
        with Pool(nproc) as pt:
            pt.starmap(brain_statis_wrapper, args_list)



