from .virtual_volume import VirtualVolume
from pathlib import Path
from .config import *
import struct
from file_io import load_image
from .block import Block
import numpy as np


class TiledVolume(VirtualVolume):
    def __init__(self, root_dir):
        super(TiledVolume, self).__init__(root_dir)
        self.VXL_1 = self.VXL_2 = self.VXL_3 = 0
        self.N_ROWS = self.N_COLS = 0
        self.BLOCKS = None
        self.reference_system_first = self.reference_system_second = self.reference_system_thrid = None

        mdata_filepath = Path(root_dir) / MDATA_BIN_FILE_NAME
        if mdata_filepath.is_file():
            self.load(mdata_filepath)
            self.init_channels()
        else:
            raise f"in TiledVolume: unable to find metadata file at {mdata_filepath}"

    def load(self, mdata_filepath):
        with open(mdata_filepath, 'rb') as f:
            mdata_version_read = struct.unpack('f', f.read(4))[0]
            if mdata_version_read != MDATA_BIN_FILE_VERSION:
                f.seek(0)
                str_size = struct.unpack('H', f.read(2))[0]
                f.read(str_size)
            self.reference_system_first, self.reference_system_second, self.reference_system_thrid, \
            self.VXL_1, self.VXL_2, self.VXL_3, self.VXL_V, self.VXL_H, self.VXL_D, \
            self.ORG_V, self.ORG_H, self.ORG_D, self.DIM_V, self.DIM_H, self.DIM_D, \
            self.N_ROWS, self.N_COLS = struct.unpack('iiifffffffffIIIHH', f.read(64))
            self.BLOCKS = [[Block(self, i, j, f) for j in range(self.N_COLS)] for i in range(self.N_ROWS)]

    def init_channels(self):
        self.DIM_C = self.BLOCKS[0][0].N_CHANS
        self.BYTESxCHAN = int(self.BLOCKS[0][0].N_BYTESxCHAN)
        self.n_active = self.DIM_C
        self.active = list(range(self.n_active))

    def load_subvolume(self, v0=-1, v1=-1, h0=-1, h1=-1, d0=-1, d1=-1):
        v0, h0, d0 = max(0, v0), max(0, h0), max(0, d0)
        v1 = v1 if 0 <= v1 <= self.DIM_V else self.DIM_V
        h1 = h1 if 0 <= h1 <= self.DIM_H else self.DIM_H
        d1 = d1 if 0 <= d1 <= self.DIM_D else self.DIM_D
        assert v1 > v0 and h1 > h0 and d1 > d0
        sbv_height = v1 - v0
        sbv_width = h1 - h0
        sbv_depth = d1 - d0
        subvol_area = Rect(h0, v0, h1, v1)
        first_time = True
        intersect_segm = self.BLOCKS[0][0].intersects_segm(d0, d1)
        if intersect_segm is not None:
            for row in range(self.N_ROWS):
                for col in range(self.N_COLS):
                    block = self.BLOCKS[row][col]
                    intersect_area = block.intersects_rect(subvol_area)
                    if intersect_area is not None:
                        for k in range(intersect_segm.ind0, intersect_segm.ind1 + 1):
                            if first_time:
                                first_time = False
                                sbv_channels = self.DIM_C
                                sbv_bytes_chan = self.BYTESxCHAN
                                assert sbv_channels == 1  # multi channel not supported here
                                if sbv_bytes_chan == 1:
                                    dt = np.uint8
                                elif sbv_bytes_chan == 2:
                                    dt = np.uint16
                                elif sbv_bytes_chan == 4:
                                    dt = np.float32
                                else:
                                    raise "unsupported datatype"
                                buffer = np.zeros((sbv_depth, sbv_height, sbv_width), dtype=dt)
                            slice_fullpath = Path(self.root_dir) / block.DIR_NAME / block.FILENAMES[k]
                            # vertices of file block
                            sv0 = 0 if v0 < intersect_area.V0 else v0 - block.ABS_V
                            sv1 = block.HEIGHT if v1 > intersect_area.V1 else v1 - block.ABS_V
                            sh0 = 0 if h0 < intersect_area.H0 else h0 - block.ABS_H
                            sh1 = block.WIDTH if h1 > intersect_area.H1 else h1 - block.ABS_H
                            sd0 = 0 if d0 < block.BLOCK_ABS_D[k] else d0 - block.BLOCK_ABS_D[k]
                            sd1 = block.BLOCK_SIZE[k] if d1 > block.BLOCK_ABS_D[k] + block.BLOCK_SIZE[k] \
                                else d1 - block.BLOCK_ABS_D[k]

                            # vertices of buffer block
                            bv0 = 0 if v0 > intersect_area.V0 else intersect_area.V0 - v0
                            bh0 = 0 if h0 > intersect_area.H0 else intersect_area.H0 - h0
                            bd0 = 0 if d0 > block.BLOCK_ABS_D[k] else block.BLOCK_ABS_D[k] - d0

                            if "NULL.tif" in str(slice_fullpath):
                                continue

                            img = load_image(slice_fullpath, flip_tif=False)
                            if len(img.shape) == 4:
                                img = img[0]
                            buffer[bd0:bd0 + sd1 - sd0, bv0:bv0 + sv1 - sv0, bh0:bh0 + sh1 - sh0] = \
                                img[sd0:sd1, sv0:sv1, sh0:sh1]
        else:
            raise "TiledVolume: depth interval out of range"
        return buffer
