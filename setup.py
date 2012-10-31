from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os

ext_modules = [
    Extension(name="pysamstats", 
              sources=["pysamstats.pyx"],
              include_dirs=['pysam', 'samtools']
              )

]

setup(
  name = 'pysamstats',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
