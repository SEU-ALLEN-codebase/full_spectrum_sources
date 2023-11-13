from file_io import load_v3draw
import os
import numpy as np
import SimpleITK as sitk


def main(dir, r, sigma, sigma2):
    """
    r: z, y, x
    """
    for b in os.listdir(dir):
        for i in os.listdir(os.path.join(dir, b)):
            if not i.endswith('v3draw'):
                continue
            # if not i.startswith('11874'):
            #     continue
            img = load_v3draw(os.path.join(dir, b, i)).sum(axis=0)
            # bin & fill
            a = sitk.GetImageFromArray(img)
            a = sitk.BinaryThreshold(a, img.mean() + 0.5 * img.std())
            a = sitk.BinaryFillhole(a)
            a = sitk.BinaryMorphologicalOpening(a, kernelRadius=(3, 3, 3))
            a = sitk.GetArrayFromImage(a)
            img = (a > 0) * img
            t = img.max(axis=0)
            a = sitk.GetImageFromArray(t)
            a = sitk.BinaryThreshold(a, t.mean() + 0.5 * t.std())
            a = sitk.BinaryFillhole(a)
            a = sitk.BinaryMorphologicalOpening(a, kernelRadius=(3, 3))
            a = sitk.ApproximateSignedDistanceMap(a)
            a = sitk.GetArrayFromImage(a)
            a = -a
            a[a < 0] = 0
            # a = (a*255/a.max()).astype('uint8')
            # sitk.WriteImage(sitk.GetImageFromArray(a), dir + "test.tiff")
            xy = a ** 3
            # pos = np.array(np.unravel_index(xy.argmax(), xy.shape))
            pos = np.argwhere(xy == xy.max()).mean(axis=0)
            rr = r[1:]
            for k in range(30):
                # print(pos)
                ind = np.array(
                    [np.clip(np.round(pos - rr), 0, None), np.clip(np.round(pos + rr + 1), None, img.shape[1:])],
                    dtype=int)
                yy, xx = np.meshgrid(*[range(*ind[:, j]) for j in range(2)], indexing='ij')
                coord = np.array([yy, xx]).reshape(2, -1)
                weight = coord - np.repeat([pos], coord.shape[1], axis=0).T
                weight = np.linalg.norm(weight, axis=0) / sigma
                weight = xy[yy, xx].reshape(-1) * np.exp(-weight ** 2)
                if weight.sum() == 0:
                    print('sth wrong with', i)
                    break
                pos = coord.dot(weight) / weight.sum()
            t = img[:, int(pos[0]), int(pos[1])]
            a = sitk.GetImageFromArray([t])
            a = sitk.BinaryThreshold(a, t.mean() + 0.5 * t.std())
            a = sitk.BinaryMorphologicalOpening(a, kernelRadius=(1, 3))
            a = sitk.ApproximateSignedDistanceMap(a)
            a = sitk.GetArrayFromImage(a)
            a = -a
            a[a < 0] = 0
            z = a[0] ** 3
            h = np.argwhere(z == z.max()).mean()
            for k in range(30):
                zz, = np.meshgrid(range(np.clip(np.round(h - r[0]), 0, None).astype(int),
                                        np.clip(np.round(h + r[0] + 1), None, img.shape[0]).astype(int)))
                w = zz - [h] * len(zz)
                w /= sigma2
                w = z[zz] * np.exp(-w ** 2)
                if w.sum() == 0:
                    print('sth wrong with', i)
                    break
                h = zz.dot(w) / w.sum()
            with open(os.path.join(dir, b, i).removesuffix('.v3draw') + '.marker', 'w') as f:
                print(b, i)
                f.write(
                    "##x,y,z,radius,shape,name,comment,color_r,color_g,color_b\n{2},{1},{0},1,1,0,0,255,0,0\n".format(h,
                                                                                                                      *pos))


if __name__ == '__main__':
    dir = "C:/Users/Zuohan Zhao/Desktop/"
    main(dir + 'raw', [20, 20, 20], 10, 10)
    # main(dir + 'z2', [20, 20, 20], 10, 10)
