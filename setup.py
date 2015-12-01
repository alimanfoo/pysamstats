from distutils.core import setup
from distutils.extension import Extension


# require pysam is pre-installed
try:
    import pysam
except ImportError:
    raise Exception('pysam not found; please install pysam first')
from distutils.version import StrictVersion
required_pysam_version = '0.8.4'
if StrictVersion(pysam.__version__) < StrictVersion(required_pysam_version):
    raise Exception('pysam version >= %s is required; found %s' %
                    (required_pysam_version, pysam.__version__))


def get_version():
    """Extract version number from source file."""
    from ast import literal_eval
    with open('pysamstats.pyx') as f:
        for line in f:
            if line.startswith('__version__'):
                return literal_eval(line.partition('=')[2].lstrip())
    raise ValueError("__version__ not found")


extensions = [Extension('pysamstats',
                        sources=['pysamstats.c'],
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
                'position based on sequence alignments from a SAM, '
                'BAM or CRAM file.',
    scripts=['scripts/pysamstats'],
    classifiers=[
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    ext_modules=extensions,
)
