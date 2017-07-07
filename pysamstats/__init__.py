# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division


__version__ = '1.0.0.dev0'


from .pileup import *
from .binned import *

stats_types_noref = ('coverage',
                     'coverage_strand',
                     'coverage_ext',
                     'coverage_ext_strand',
                     'tlen',
                     'tlen_strand',
                     'mapq',
                     'mapq_strand',
                     'baseq',
                     'baseq_strand',
                     'mapq_binned',
                     'alignment_binned',
                     'tlen_binned')

stats_types_withref = ('variation',
                       'variation_strand',
                       'baseq_ext',
                       'baseq_ext_strand',
                       'coverage_gc',
                       'coverage_binned',
                       'coverage_ext_binned')

stats_types = sorted(stats_types_noref + stats_types_withref)
