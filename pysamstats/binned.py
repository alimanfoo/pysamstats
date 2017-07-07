# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import functools


import pysamstats.opt as opt
import pysamstats.util as util


# TODO factor out docstring parameters


# noinspection PyShadowingBuiltins
def stat_binned(type,
                alignmentfile,
                fafile,
                chrom=None,
                start=None,
                end=None,
                one_based=False,
                window_size=300,
                window_offset=None):
    """Generate statistics per genome window, based on all reads whose alignment starts within
    the window.

    Parameters
    ----------
    type : string
        Statistics type. One of "coverage", "coverage_ext", "mapq", "alignment", "tlen".
    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path.
    fafile : pysam.FastaFile or string
        FASTA file or file path.
    chrom : string
        Chromosome/contig.
    start : int
        Start position.
    end : int
        End position.
    one_based : bool
        Coordinate system, False if zero-based (default), True if one-based.
    window_size : int
        Window size to use.
    window_offset : int
        Distance from window start to record position.

    Returns
    -------
    recs : iterator
        An iterator yielding dict objects, where each dict holds data for a single window.

    """

    try:
        stat = statsobj_binned[type]
    except KeyError:
        raise ValueError('unsupported statistics type: %r' % type)

    return opt.iter_binned(stat, alignmentfile=alignmentfile, fafile=fafile, chrom=chrom,
                           start=start, end=end, one_based=one_based, window_size=window_size,
                           window_offset=window_offset)


# noinspection PyShadowingBuiltins
def load_binned(type,
                alignmentfile,
                fafile,
                chrom=None,
                start=None,
                end=None,
                one_based=False,
                window_size=300,
                window_offset=None):
    """Load statistics per genome window, based on all reads whose alignment starts within
    the window.

    Parameters
    ----------
    type : string
        Statistics type. One of "coverage", "coverage_ext", "mapq", "alignment", "tlen".
    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path.
    fafile : pysam.FastaFile or string
        FASTA file or file path.
    chrom : string
        Chromosome/contig.
    start : int
        Start position.
    end : int
        End position.
    one_based : bool
        Coordinate system, False if zero-based (default), True if one-based.
    window_size : int
        Window size to use.
    window_offset : int
        Distance from window start to record position.

    Returns
    -------
    ra : numpy structured array
        A structured array.

    """

    stat = functools.partial(stat_binned, type)
    try:
        dtype = dtypes_binned[type]
    except KeyError:
        raise ValueError('unsupported statistics type: %r' % type)

    return util.load_stats(stat, dtype, alignmentfile=alignmentfile, fafile=fafile, chrom=chrom,
                           start=start, end=end, one_based=one_based, window_size=window_size,
                           window_offset=window_offset)


statsobj_binned = {
    'coverage': opt.CoverageBinned(),
}


dtypes_binned = {
    'coverage': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('gc', 'u1'),
        ('reads_all', 'i4'),
        ('reads_pp', 'i4')
    ],
}


# backwards compatibility
#########################


_stat_doc_lines = stat_binned.__doc__.split('\n')
_load_doc_lines = load_binned.__doc__.split('\n')
# strip "type" parameter
_stat_doc = '\n'.join(_stat_doc_lines[:4] + _stat_doc_lines[6:])
_load_doc = '\n'.join(_load_doc_lines[:4] + _stat_doc_lines[6:])
# named functions
stat_coverage_binned = functools.partial(stat_binned, 'coverage')
stat_coverage_binned.__doc__ = _stat_doc
load_coverage_binned = functools.partial(load_binned, 'coverage')
load_coverage_binned.__doc__ = _load_doc
