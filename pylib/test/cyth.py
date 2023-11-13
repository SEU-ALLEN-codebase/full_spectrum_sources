import unittest
from file_io import load_image, save_image
from cythonized.ada_thr import adaptive_threshold
from cythonized.img_pca_filter import img_pca_filter
import time


class MyTestCase(unittest.TestCase):
    def test_pca(self):
        img = load_image(r"D:\GitHub\Segment-Splicing-Search\data\gt\512_128_tiff\14976_13296_5475.tiff")
        t1 = time.time()
        out = img_pca_filter(img)
        t2 = time.time()
        print(f'used {t2-t1}s.')
        save_image('d:/test.tiff', out)

    def test_ada(self):
        img = load_image(r"D:\GitHub\Segment-Splicing-Search\data\gt\512_128_tiff\14976_13296_5475.tiff")
        t1 = time.time()
        out = adaptive_threshold(img)
        t2 = time.time()
        print(f'used {t2-t1}s.')
        save_image('d:/test.tiff', out)


if __name__ == '__main__':
    unittest.main()
