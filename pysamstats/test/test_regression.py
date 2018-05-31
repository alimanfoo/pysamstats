# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import logging
import sys

from pysam import Samfile
import numpy as np
from os.path import isfile

from pysamstats.config import stats_types, stats_types_withref
import pysamstats

logger = logging.getLogger(__name__)
debug = logger.debug

# PY2/3 compatibility
PY2 = sys.version_info[0] == 2


# no test_prefix so not run
def generate_fixtures():

    bampath = "fixture/test.bam"
    fastapath = "fixture/ref.fa"
    archive = "fixture/regression.npz"

    # simple stats
    dat = {}
    for q in stats_types:
        if q in stats_types_withref:
            dat[q] = getattr(pysamstats, "load_" + q)(Samfile(bampath), fafile=fastapath)
        else:
            dat[q] = getattr(pysamstats, "load_" + q)(Samfile(bampath))

    assert not isfile(archive), "Attempting to regenerate static test archive."
    np.savez_compressed(archive, **dat)


def test_against_fixtures():

    # load fixtures from numpy array.
    bampath = "fixture/test.bam"
    fastapath = "fixture/ref.fa"
    archive = "fixture/regression.npz"

    testset = np.load(archive)

    for q in stats_types:
        if q in stats_types_withref:
            x = getattr(pysamstats, "load_" + q)(Samfile(bampath), fafile=fastapath)
        else:
            x = getattr(pysamstats, "load_" + q)(Samfile(bampath))

        # Loop through all fields.
        for key in testset[q].dtype.fields.keys():
            np.testing.assert_array_equal(testset[q][key], x[key], err_msg=key)
