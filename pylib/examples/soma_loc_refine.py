import sys
sys.path.append(sys.path[0] + "/..")
from file_io import load_image
import numpy as np
import SimpleITK as sitk
import fire
import glob


class SomaLocRefine(object):

    def __init__(self, win_radius=(20, 20, 20), kernel_radius=(3, 3, 3), sigma_xy=10, sigma_z=10, iterations=30):
        """
        :param win-size: The radius of the window for mean-shift, as [x,y,z].
        :param kernel_radius: the 3D kernel size for morphological opening, as [x,y,z]
        :param sigma_xy: The Gaussian sigma for xy-plane mean-shift.
        :param sigma_z: The Gaussian sigma for z-axis mean-shift.
        :param iterations: Number of mean-shift iterations.
        """
        self._win_radius = win_radius
        self._sigma_xy = sigma_xy
        self._sigma_z = sigma_z
        self._kernel_radius = kernel_radius
        self._iterations = iterations

    def find(self, path: str, start="AUTO"):
        """
        Refine the soma location for a single 3D image.
        :param path: The path to your input image.
        :param start: The 3D coordinate(list) of the initial soma location. Default as the average of the DT maxima.
        """
        img = load_image(path)
        if len(img.shape) == 4:
            img = img[0]
        # bin & fill
        a = sitk.GetImageFromArray(img)
        a = sitk.BinaryThreshold(a, img.mean() + 0.5 * img.std())
        a = sitk.BinaryFillhole(a)
        a = sitk.BinaryMorphologicalOpening(a, kernelRadius=self._kernel_radius[-1:-4:-1])  # [x,y,z] -> [z,y,x]
        a = sitk.GetArrayFromImage(a)
        img = (a > 0) * img
        t = img.max(axis=0)
        a = sitk.GetImageFromArray(t)
        a = sitk.BinaryThreshold(a, t.mean() + 0.5 * t.std())
        a = sitk.BinaryFillhole(a)
        a = sitk.BinaryMorphologicalOpening(a, kernelRadius=self._kernel_radius[-2:-4:-1])    # [x,y,z] -> [y,x]
        a = sitk.ApproximateSignedDistanceMap(a)
        a = sitk.GetArrayFromImage(a)
        a = -a
        a[a < 0] = 0
        xy = a ** 3
        if start == "AUTO":
            pos = np.argwhere(xy == xy.max()).mean(axis=0)
        elif start == "CENTER":
            pos = np.array(img.shape[1:]) // 2
        else:
            pos = start[-2:-4:-1]   # [x,y,z] -> [y,x]
        rr = self._win_radius[-2:-4:-1]     # [x,y,z] -> [y,x]
        for k in range(self._iterations):
            ind = np.array(
                [np.clip(np.round(pos - rr), 0, None), np.clip(np.round(pos + rr + 1), None, img.shape[1:])],
                dtype=int)
            yy, xx = np.meshgrid(*[range(*ind[:, j]) for j in range(2)], indexing='ij')
            coord = np.array([yy, xx]).reshape(2, -1)
            weight = coord - np.repeat([pos], coord.shape[1], axis=0).T
            weight = np.linalg.norm(weight, axis=0) / self._sigma_xy
            weight = xy[yy, xx].reshape(-1) * np.exp(-weight ** 2)
            pos = coord.dot(weight) / weight.sum()
        t = img[:, int(pos[0]), int(pos[1])]
        a = sitk.GetImageFromArray([t])
        a = sitk.BinaryThreshold(a, t.mean() + 0.5 * t.std())
        a = sitk.BinaryMorphologicalOpening(a, kernelRadius=(1, self._kernel_radius[2]))
        a = sitk.ApproximateSignedDistanceMap(a)
        a = sitk.GetArrayFromImage(a)
        a = -a
        a[a < 0] = 0
        z = a[0] ** 3
        h = np.argwhere(z == z.max()).mean()
        for k in range(self._iterations):
            zz, = np.meshgrid(range(np.clip(np.round(h - self._win_radius[2]), 0, None).astype(int),
                                    np.clip(np.round(h + self._win_radius[2] + 1), None, img.shape[0]).astype(int)))
            w = zz - [h] * len(zz)
            w /= self._sigma_z
            w = z[zz] * np.exp(-w ** 2)
            h = zz.dot(w) / w.sum()
        return "{2},{1},{0}".format(h, *pos)

    def batch(self, pattern=".", start="CENTER"):
        """
        :param pattern: filepath pattern for all input images, used by glob
        :param start:
        """
        for fp in glob.iglob(pattern):
            try:
                print(fp + "\t" + self.find(fp, start=start))
            except:
                print(fp + "\tfailure")


    def wrapper(self, v3d, apo, attenuationArgs=()):
        pass


if __name__ == '__main__':
    fire.Fire(SomaLocRefine)
