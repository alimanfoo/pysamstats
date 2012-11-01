#!/usr/bin/env python

from pysamstats import NormedCoverageStatsTable
t = NormedCoverageStatsTable('fixture/test.bam',
                             'Pf3D7_01_v3', 5000, 10000)
from petl import *
print see(t, 10)
