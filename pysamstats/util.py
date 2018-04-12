# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
from operator import itemgetter
from pysam import AlignmentFile


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

    # check if contig label dtype is appropriate length
    max_seqid_len = determine_max_seqid(kwargs["alignmentfile"])
    dtype["chrom"] = "a{0}".format(max_seqid_len)

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


def determine_max_seqid(alignmentfile):

    if isinstance(alignmentfile, str):
        alignmentfile = AlignmentFile(alignmentfile)

    return max([len(x) for x in alignmentfile.references])
