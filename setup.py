from distutils.core import setup
from distutils.extension import Extension

setup(
    name='pysamstats',
    version='0.4.4',
    author='Alistair Miles',
    author_email='alimanfoo@googlemail.com',
    url='https://github.com/alimanfoo/pysamstats',
    license='MIT Licenses',
    description='A small Python utility for calculating statistics per genome position based on pileups from a SAM or BAM file.',
    scripts=['pysamstats'],
    ext_modules = [Extension("pysamstats", ["pysamstats.c"], include_dirs=["pysam"])],
    classifiers=['Intended Audience :: Developers',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python',
                 'Topic :: Software Development :: Libraries :: Python Modules'
                 ],
    install_requires=['pysam>=0.7', 'numpy>=1.6'],
)
