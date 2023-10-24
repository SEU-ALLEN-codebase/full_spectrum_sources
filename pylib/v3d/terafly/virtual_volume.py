class VirtualVolume:
    def __init__(self, root_dir=None, vxl_1=0, vxl_2=0, vxl_3=0):
        self.root_dir = root_dir
        self.VXL_V = vxl_1
        self.VXL_H = vxl_2
        self.VXL_D = vxl_3
        self.ORG_V = self.ORG_H = self.ORG_D = 0.0
        self.DIM_V = self.DIM_H = self.DIM_D = self.DIM_C = 0
        self.BYTESxCHAN = 0
        self.active = None
        self.n_active = 0
        self.t0 = self.t1 = 0
        self.DIM_T = 1
