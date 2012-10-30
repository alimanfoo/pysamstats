from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("pysamstats", ["pysamstats.pyx"])]

setup(
  name = 'pysamstats',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
