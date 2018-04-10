# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import logging
import sys


from nose.tools import eq_
from pysam import Samfile, Fastafile
import numpy as np


# from .util import normalise_coords, fwd, rev, pp, mean, std, rms, vmax, compare_iterators
# import pysamstats
from pysamstats.io import write_hdf5
import tables

logger = logging.getLogger(__name__)
debug = logger.debug

def test_pileup_csize():

    import tempfile
    label = "AS2_scf7180000696055"
    tmp = tempfile.mktemp(suffix=".h5")
    print("xxx")
    write_hdf5("coverage", tmp, Samfile("fixture/z2.bam"), chrom=label)
    print(tmp)
    eq_(1, 1)
    #tables.File(tmp)

