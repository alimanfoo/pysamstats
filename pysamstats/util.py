# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import sys
import csv
import time
import itertools
from operator import itemgetter


def flatten(recs, *fields):
    """Convert a record (dict) iterator to a row (tuple) iterator.

    Parameters
    ----------

    recs : iterator of dicts
        records generator
    fields : list of strings
        names of fields to select

    Returns
    -------

    rows : iterator of tuples
        rows generator

    """

    getter = itemgetter(*fields)
    it = (getter(rec) for rec in recs)
    return it


def tabulate(stat, *args, **kwargs):
    """Tabulate statistics.

    Parameters
    ----------

    stat : string
        statistics type
    *args
        passed through to statistics function
    fields : list of strings
        names of fields to select
    **args
        passed through to statistics function

    Returns
    -------

    table : row container

    """

    return _StatsTable(stat, *args, **kwargs)


class _StatsTable(object):

    def __init__(self, stats_type, *args, **kwargs):
        try:
            self.stats_function = globals()['stat_' + stats_type]
            self.fields = kwargs.pop('fields', None)
            if self.fields is None:
                self.fields = globals()['fields_' + stats_type]
            self.args = args
            self.kwargs = kwargs
        except KeyError:
            raise Exception('statistics type not found: %r' % stats_type)

    def __iter__(self):
        recs = self.stats_function(*self.args, **self.kwargs)
        fields = tuple(self.fields)
        rows = flatten(recs, *fields)
        yield fields
        for row in rows:
            yield row


def load_stats(statfun, default_dtype, user_dtype, user_fields, **kwargs):

    import numpy as np

    # determine fields to load
    default_fields = [t[0] for t in default_dtype]
    if user_fields is None:
        fields = default_fields
    else:
        fields = user_fields
        if any([f not in default_fields for f in fields]):
            raise ValueError('invalid fields: %r' % fields)

    # determine dtype
    dtype = dict(default_dtype)
    if user_dtype is not None:
        dtype.update(dict(user_dtype))

    # handle single field requested
    if len(fields) == 1:
        dtype = dtype[fields[0]]
    else:
        dtype = [(f, dtype[f]) for f in fields]

    # setup record generator
    recs = statfun(**kwargs)

    # flatten records
    it = flatten(recs, *fields)

    # load into a Numpy array
    a = np.fromiter(it, dtype=dtype)

    # view as recarray for convenience
    if len(fields) > 1:
        a = a.view(np.recarray)

    return a
