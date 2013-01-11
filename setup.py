from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name='pysamstats',
    version='0.3',
    author='Alistair Miles',
    author_email='alimanfoo@googlemail.com',
    url='https://github.com/alimanfoo/pysamstats',
    license='MIT Licenses',
    description='A small Python utility for calculating statistics per genome position based on pileups from a SAM or BAM file.',
    scripts=['pysamstats'],
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("pysamstats", ["pysamstats.pyx"], include_dirs=["pysam"])],
    classifiers=['Intended Audience :: Developers',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python',
                 'Topic :: Software Development :: Libraries :: Python Modules'
                 ],
    requires=['pysam (>=0.7)'],
)
