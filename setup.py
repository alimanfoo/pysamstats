from distutils.core import setup
from distutils.extension import Extension
from ast import literal_eval


try:
    import pysam
    from Cython.Distutils import build_ext # Cython should be installed via pysam
except ImportError:
    raise Exception('please install pysam first, e.g.: pip install --upgrade pysam')


try:
    import numpy as np
except ImportError:
    raise Exception('please install numpy first')


def get_version(source='pysamstats.pyx'):
    with open(source) as f:
        for line in f:
            if line.startswith('__version__'):
                return literal_eval(line.partition('=')[2].lstrip())
    raise ValueError("__version__ not found")


setup(
    name='pysamstats',
    version=get_version(),
    author='Alistair Miles',
    author_email='alimanfoo@googlemail.com',
    url='https://github.com/alimanfoo/pysamstats',
    license='MIT Licenses',
    description='A Python utility for calculating statistics against genome position based on sequence alignments from a SAM or BAM file.',
    scripts=['pysamstats'],
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension('pysamstats',
                           sources=['pysamstats.pyx'],
                           include_dirs=[np.get_include()] + pysam.get_include(),
                           define_macros=pysam.get_defines()),
                 ],
    classifiers=['Intended Audience :: Developers',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python',
                 'Topic :: Software Development :: Libraries :: Python Modules'
                 ],
)


