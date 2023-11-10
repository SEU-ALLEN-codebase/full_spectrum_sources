from .tiled_volume import TiledVolume
from pathlib import Path
from .config import *


class Interface:
    loaded = {}

    @staticmethod
    def instance(path):
        """
        now it only supports TiledVolume (Tiled Tiff/V3DRaw 3D)
        :param path: teraconverted brain / resolution
        :return: an instance of the image manager from the path
        """
        path = Path(path)
        volume = None
        if path.is_dir():
            format = "unknown"
            if (path / FORMAT_MDATA_FILE_NAME).is_file():
                try:
                    with open(path / FORMAT_MDATA_FILE_NAME) as f:
                        format = f.readline().rstrip()
                    if format == TILED_MC_FORMAT:
                        pass
                    elif format == STACKED_FORMAT:
                        pass
                    elif format == TILED_FORMAT:
                        volume = TiledVolume(path)
                    elif format == SIMPLE_FORMAT:
                        pass
                    elif format == SIMPLE_RAW_FORMAT:
                        pass
                    elif format == TIME_SERIES:
                        pass
                    else:
                        raise f"cannot recognize format {format}"
                except:
                    raise f"cannot import {format} at {path}"
            if volume is None:
                # try:
                #     volume = TiledMCVolume(path)
                # except:
                #     print(f"cannot import TiledMCVolume at {path}")
                #     try:
                #         volume = StackedVolume(path)
                #     except:
                #         print(f"cannot import StackedVolume at {path}")
                try:
                    volume = TiledVolume(path)
                except Exception as e:
                    print(e)
                    print(f"cannot import TiledVolume at {path}")
                    # try:
                    #     volume = SimpleVolume(path)
                    # except:
                    #     print(f"cannot import SimpleVolume at {path}")
                    #     try:
                    #         volume = SimpleVolumeRaw(path)
                    #     except:
                    #         print(f"cannot import SimpleVolumeRaw at {path}")
                    #         try:
                    #             volume = TimeSeries(path)
                    #         except:
                    #             print(f"cannot import TimeSeries at {path}")
                    raise f"generic error occurred when importing volume at {path}"
        elif path.is_file():
            if path.suffix in ['.raw', '.v3draw', '.RAW', '.V3DRAW', '.tif', '.tiff', '.TIF', '.TIFF']:
                pass
            elif path.suffix in ['.xml', '.XML']:
                pass
            else:
                raise f"Unsupported file extensions for {path}"
        else:
            raise f"path {path} does not exist"
        return volume

    @classmethod
    def get_dim(cls, path):
        """
        get teraconvert dim range.
        it uses static member `loaded` to store metadata (the image manager instance), no need to load it each time.
        if you need clearing, refer to get_subvolume.loaded.
        you can also manage the image manager at your own will.
        :param path: teraconverted brain / resolution
        :return: (max x, max y, max z)
        """
        if path not in cls.loaded:
            cls.loaded[path] = cls.instance(path)
        vol = cls.loaded[path]
        return vol.DIM_H, vol.DIM_V, vol.DIM_D, vol.DIM_C

    @classmethod
    def get_subvolume(cls, path, x0, x1, y0, y1, z0, z1):
        """
        different from vaa3d, it returns the image of its original pixel size.
        it uses static member `loaded` to store metadata (the image manager instance), no need to load it each time.
        if you need clearing, refer to get_subvolume.loaded.
        you can also manage the image manager at your own will.
        the indexing is pixel-wise, meaning you have to use different coordinates for different resolutions
        if the loaded image is tif, it's already flipped
        :param path: teraconverted brain / resolution
        :param x0: starting x
        :param x1: ending x
        :param y0: starting y
        :param y1: ending y
        :param z0: starting z
        :param z1: ending z
        :return: the cropped image
        """
        if path not in cls.loaded:
            cls.loaded[path] = cls.instance(path)
        return cls.loaded[path].load_subvolume(y0, y1, x0, x1, z0, z1)




