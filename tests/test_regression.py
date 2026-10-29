# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division

from pysam import Samfile
import numpy as np


import pysamstats
from pysamstats.config import stats_types, stats_types_withref
from .util import fixture_path


def test_against_fixtures():

    # load fixtures from numpy array
    bampath = fixture_path("test.bam")
    fastapath = fixture_path("ref.fa")
    archive = fixture_path("regression.npz")

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
