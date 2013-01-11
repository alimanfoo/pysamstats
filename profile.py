#!/usr/bin/env python

import sys
import time
import pstats
import cProfile
import timeit


from pysam import Samfile, Fastafile


sys.path.append('.')
import pysamstats


def profile(fun):
    samfile = Samfile('fixture/test.bam')
    count = 0
    f = getattr(pysamstats, fun)
    for _ in f(samfile, chrom='Pf3D7_01_v3', start=0, end=1000):
        count += 1

def profile_withrefseq(fun):
    samfile = Samfile('fixture/test.bam')
    fafile = Fastafile('fixture/ref.fa')
    count = 0
    f = getattr(pysamstats, fun)
    for _ in f(samfile, fafile, chrom='Pf3D7_01_v3', start=0, end=1000):
        count += 1

stats_types_requiring_fasta = ('variation', 
                               'variation_strand', 
                               'baseq_ext', 
                               'baseq_ext_strand', 
                               'coverage_gc',
                               'coverage_normed_gc')

fun = sys.argv[1]
if fun in stats_types_requiring_fasta:
    cmd = 'profile_withrefseq("stat_%s")' % fun
else:
    cmd = 'profile("stat_%s")' % fun
    
prof_fn = '%s.prof' % fun
cProfile.runctx(cmd, globals(), locals(), prof_fn)
s = pstats.Stats(prof_fn)
s.strip_dirs().sort_stats('time').print_stats()
print timeit.repeat(cmd, number=20, repeat=3, setup='from __main__ import profile, profile_withrefseq')


