import sys
import time
import pstats
import cProfile
import timeit


from pysam import Samfile


def test():
    samfile = Samfile('fixture/test.bam')
    count = 0
    for _ in samfile.pileup('Pf3D7_01_v3', start=0, end=1000):
        count += 1

cProfile.runctx('test()', globals(), locals(), 'profile/test.prof')
s = pstats.Stats('profile/test.prof')
s.strip_dirs().sort_stats('time').print_stats()
print timeit.repeat('test()', number=20, repeat=3, setup='from __main__ import test')


