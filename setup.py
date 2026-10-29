import pysam
from setuptools import Extension, setup


setup(
    ext_modules=[
        Extension(
            "pysamstats.opt",
            sources=["src/pysamstats/opt.pyx"],
            include_dirs=pysam.get_include(),
            define_macros=pysam.get_defines(),
        )
    ]
)
