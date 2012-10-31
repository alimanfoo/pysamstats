#!/usr/bin/env python

from pysamstats import MapqStatsTable
t = MapqStatsTable('fixture/test.bam',
                    'Pf3D7_01_v3', 0, 10000)

from petl import *
nrows(progress(t, 1000))
