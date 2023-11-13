from .virtual_volume import VirtualVolume
import struct
from .config import *


class Segm:
    def __init__(self, d0, d1, ind0, ind1):
        self.D0 = d0
        self.D1 = d1
        self.ind0 = ind0
        self.ind1 = ind1


class Block(VirtualVolume):
    def __init__(self, container, row_index, col_index, file):
        super(Block, self).__init__()
        self.CONTAINER = container
        self.ROW_INDEX = row_index
        self.COL_INDEX = col_index

        self.DIR_NAME = None
        self.FILENAMES = None
        self.HEIGHT = self.WIDTH = self.DEPTH = 0
        self.N_BLOCKS = 0
        self.N_CHANS = 0
        self.N_BYTESxCHAN = 0
        self.ABS_V = self.ABS_H = 0
        self.BLOCK_SIZE = 0
        self.BLOCK_ABS_D = 0

        self.unbinarize_from(file)

    def unbinarize_from(self, file):
        self.HEIGHT, self.WIDTH, self.DEPTH, self.N_BLOCKS, self.N_CHANS, self.ABS_V, self.ABS_H, str_size = \
            struct.unpack('IIIIIiiH', file.read(30))
        self.DIR_NAME = file.read(str_size).decode('utf-8').rstrip("\x00")
        self.FILENAMES = []
        self.BLOCK_SIZE = []
        self.BLOCK_ABS_D = []
        for i in range(self.N_BLOCKS):
            str_size = struct.unpack('H', file.read(2))[0]
            self.FILENAMES.append(file.read(str_size).decode('utf-8').rstrip("\x00"))
            self.BLOCK_SIZE.append(struct.unpack('I', file.read(4))[0])
            self.BLOCK_ABS_D.append(struct.unpack('i', file.read(4))[0])
        self.N_BYTESxCHAN = struct.unpack('I', file.read(4))[0]

    def intersects_segm(self, d0: int, d1: int):
        if d0 >= self.BLOCK_ABS_D[self.N_BLOCKS - 1] + self.BLOCK_SIZE[self.N_BLOCKS - 1] or d1 <= 0:
            return None # no intersection

        i0 = 0
        while i0 < self.N_BLOCKS - 1:
            if d0 < self.BLOCK_ABS_D[i0 + 1]:
                break
            i0 += 1

        i1 = self.N_BLOCKS - 1
        while i1 > 0:
            if d1 > self.BLOCK_ABS_D[i1]:
                break
            i1 -= 1

        return Segm(max(d0, 0), min(d1, self.DEPTH), i0, i1)

    def intersects_rect(self, area: Rect):
        if area.H0 < self.ABS_H + self.WIDTH and area.H1 > self.ABS_H and \
                area.V0 < self.ABS_V + self.HEIGHT and area.V1 > self.ABS_V:
            return Rect(
                max(self.ABS_H, area.H0),
                max(self.ABS_V, area.V0),
                min(self.ABS_H + self.WIDTH, area.H1),
                min(self.ABS_V + self.HEIGHT, area.V1)
            )
        else:
            return None
