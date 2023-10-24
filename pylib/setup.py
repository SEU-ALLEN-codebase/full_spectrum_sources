from setuptools import setup, Extension
import numpy

setup(
    name='pylib',
    version='',
    packages=['v3d', 'v3d.terafly', 'graph', 'anatomy', 'morph_topo', 'neuron_quality'],
    url='',
    license='MIT',
    author='Yufeng Liu',
    author_email='',
    description='',
    ext_modules=[
        Extension(
            'cythonized.img_pca_filter',
            sources=['cythonized/img_pca_filter.pyx'],
            language='c++',
            include_dirs=[numpy.get_include()],
            library_dirs=[],
            libraries=[],
            extra_compile_args=[],
            extra_link_args=[]
            ),
        Extension(
            'cythonized.ada_thr',
            sources=['cythonized/ada_thr.pyx'],
            language='c++',
            include_dirs=[numpy.get_include()],
            library_dirs=[],
            libraries=[],
            extra_compile_args=[],
            extra_link_args=[]
            )
    ]
)
