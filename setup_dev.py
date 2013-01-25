from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import pysam

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("pysamstats", 
                             ["pysamstats.pyx"], 
#                             include_dirs=pysam.get_include(),
                             include_dirs=['pysam', 'samtools'],
                             define_macros=pysam.get_defines())],
)
