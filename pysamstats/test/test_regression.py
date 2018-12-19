# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import logging
import sys

from pysam import Samfile
import numpy as np
from os.path import isfile


import pysamstats
from pysamstats.config import stats_types, stats_types_withref


# no test_prefix so not run during unit tests
def generate_fixtures():

    bampath = "fixture/test.bam"
    fastapath = "fixture/ref.fa"
    archive = "fixture/regression.npz"
    assert not isfile(archive)

    # simple stats
    dat = {}
    for q in stats_types:
        if q in stats_types_withref:
            dat[q] = getattr(pysamstats, "load_" + q)(Samfile(bampath), fafile=fastapath)
        else:
            dat[q] = getattr(pysamstats, "load_" + q)(Samfile(bampath))

    np.savez_compressed(archive, **dat)


def test_against_fixtures():

    # load fixtures from numpy array
    bampath = "fixture/test.bam"
    fastapath = "fixture/ref.fa"
    archive = "fixture/regression.npz"

    testset = np.load(archive)

    for q in stats_types:
        if q in stats_types_withref:
            x = getattr(pysamstats, "load_" + q)(Samfile(bampath), fafile=fastapath)
        else:
            x = getattr(pysamstats, "load_" + q)(Samfile(bampath))

        # loop through all fields
        for key in testset[q].dtype.names:
            expect = testset[q][key]
            actual = x[key]
            try:
                np.testing.assert_array_equal(expect, actual, err_msg=key)
            except AssertionError:
                print(expect[expect != actual])
                print(actual[expect != actual])
                raise
