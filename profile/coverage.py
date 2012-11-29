import sys
import time
import pstats
import cProfile
import timeit


from pysam import Samfile
import pyximport; pyximport.install()


sys.path.append('.')
import pysamstats


print 'raw pileup'

def test1():
    samfile = Samfile('fixture/test.bam')
    count = 0
    for col in samfile.pileup('Pf3D7_01_v3'):
        reads = col.pileups
        count += 1
    print count

print timeit.repeat('test1()', number=1, repeat=3, setup='from __main__ import test1')
cProfile.runctx('test1()', globals(), locals(), 'test1.prof')
s = pstats.Stats('test1.prof')
s.strip_dirs().sort_stats('time').print_stats()


print 'stat_coverage'

def test2():
    samfile = Samfile('fixture/test.bam')
    count = 0
    for rec in pysamstats.stat_coverage(samfile, 'Pf3D7_01_v3'):
        count += 1
    print count

print timeit.repeat('test2()', number=1, repeat=3, setup='from __main__ import test2')
cProfile.runctx('test2()', globals(), locals(), 'test2.prof')
s = pstats.Stats('test2.prof')
s.strip_dirs().sort_stats('time').print_stats()


print 'stat_coverage_strand'

def test3():
    samfile = Samfile('fixture/test.bam')
    count = 0
    for rec in pysamstats.stat_coverage_strand(samfile, 'Pf3D7_01_v3'):
        count += 1
    print count

print timeit.repeat('test3()', number=1, repeat=3, setup='from __main__ import test3')
cProfile.runctx('test3()', globals(), locals(), 'test3.prof')
s = pstats.Stats('test3.prof')
s.strip_dirs().sort_stats('time').print_stats()


