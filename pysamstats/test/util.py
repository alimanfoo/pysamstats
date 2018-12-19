# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import logging
import sys
from math import sqrt


import numpy as np
from numpy import around as round
from nose.tools import eq_, assert_almost_equal


from pysam import Samfile, Fastafile


logger = logging.getLogger(__name__)
debug = logger.debug


# PY2/3 compatibility
PY2 = sys.version_info[0] == 2
if PY2:
    # noinspection PyUnresolvedReferences
    from itertools import izip_longest
else:
    from itertools import zip_longest as izip_longest


def compare_iterators(expected, actual):
    for e, a in izip_longest(expected, actual, fillvalue=None):
        assert e is not None, ('expected value is None', e, a)
        assert a is not None, ('actual value is None', e, a)
        for k, v in e.items():
            try:
                if isinstance(v, float):
                    assert_almost_equal(v, a[k])
                else:
                    eq_(v, a[k])
            except:
                debug('mismatch %r, expected %r, found %r' % (k, v, a[k]))
                debug('expected: %r' % sorted(e.items()))
                debug('actual: %r' % sorted(a.items()))
                raise
        for k in a:  # check no unexpected fields
            try:
                assert k in e
            except:
                debug('missing %r' % k)
                debug('expected: %r' % sorted(e.items()))
                debug('actual: %r' % sorted(a.items()))
                raise


def normalise_coords(one_based, start, end):
    """Normalise start and end coordinates.

    Parameters
    ----------
    one_based : bool
    start : int
    end : int

    Returns
    -------
    start : int
    end : int

    """
    if one_based:
        start = start - 1 if start is not None else None
        end = end - 1 if end is not None else None
    return start, end


def fwd(reads):
    return [read for read in reads if not read.alignment.is_reverse]


def rev(reads):
    return [read for read in reads if read.alignment.is_reverse]


def pp(reads):
    return [read for read in reads if read.alignment.is_proper_pair]


def rms(a):
    if a:
        return int(round(sqrt(np.mean(np.power(a, 2)))))
    else:
        return 0


def mean(a):
    if a:
        return int(round(np.mean(a)))
    else:
        return 0


def std(a):
    if len(a) >= 2:
        std = np.std(a, ddof=1)
        if np.isnan(std):
            return 0
        return int(round(std))
    else:
        return 0


def vmax(a):
    if a:
        return max(a)
    else:
        return 0


def rootmean(sqsum, count):
    if count > 0:
        return int(round(sqrt(sqsum / count)))
    else:
        return 0
