#!/usr/bin/env python

import sys
import time
import pstats
import cProfile
import timeit


from pysam import Samfile, Fastafile


sys.path.append('.')
import pysamstats


def test(fun):
    samfile = Samfile('fixture/test.bam')
    count = 0
    f = getattr(pysamstats, fun)
    for _ in f(samfile, chrom='Pf3D7_01_v3', start=0, end=1000):
        count += 1

def test_withrefseq(fun):
    samfile = Samfile('fixture/test.bam')
    fafile = Fastafile('fixture/ref.fa')
    count = 0
    f = getattr(pysamstats, fun)
    for _ in f(samfile, fafile, chrom='Pf3D7_01_v3', start=0, end=1000):
        count += 1

stats_types_requiring_fasta = ('variation', 'variation_strand', 'baseq_ext', 'baseq_ext_strand')

fun = sys.argv[1]
if fun in stats_types_requiring_fasta:
    cmd = 'test_withrefseq("stat_%s")' % fun
else:
    cmd = 'test("stat_%s")' % fun
    
cProfile.runctx(cmd, globals(), locals(), 'test.prof')
s = pstats.Stats('test.prof')
s.strip_dirs().sort_stats('time').print_stats()
print timeit.repeat(cmd, number=20, repeat=3, setup='from __main__ import test, test_withrefseq')


