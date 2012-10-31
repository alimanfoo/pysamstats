#!/usr/bin/env python

from pysamstats import TlenStatsTable
t = TlenStatsTable('fixture/test.bam', 
                    'Pf3D7_01_v3', 0, 10000)
from petl import *
print see(t, 5)
