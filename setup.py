try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from ast import literal_eval
from distutils.extension import Extension
from Cython.Build import cythonize
import pysam


def get_version(source='pysamstats.pyx'):
    with open(source) as f:
        for line in f:
            if line.startswith('__version__'):
                return literal_eval(line.partition('=')[2].lstrip())
    raise ValueError("__version__ not found")


extensions = [Extension('pysamstats',
                        sources=['pysamstats.pyx'],
                        include_dirs=pysam.get_include(),
                        define_macros=pysam.get_defines())]


setup(
    name='pysamstats',
    version=get_version(),
    author='Alistair Miles',
    author_email='alimanfoo@googlemail.com',
    url='https://github.com/alimanfoo/pysamstats',
    license='MIT Licenses',
    description='A Python utility for calculating statistics against genome '
                'position based on sequence alignments from a SAM or BAM file.',
    scripts=['scripts/pysamstats'],
    classifiers=[
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    ext_modules=cythonize(extensions),
    install_requires=['pysam>=0.8.1'],
)
