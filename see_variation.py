#!/usr/bin/env python

from pysamstats import VariationStatsTable
t = VariationStatsTable('fixture/test.bam', 'fixture/ref.fa',
                       'Pf3D7_01_v3', 0, 10000)
from petl import *
print see(t, 5)
