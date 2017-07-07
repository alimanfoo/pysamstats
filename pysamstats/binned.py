# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import functools


import pysamstats.opt as opt
import pysamstats.util as util


_doc_params = """
    Parameters
    ----------
    type : string
        Statistics type. One of "coverage", "coverage_ext", "mapq", "alignment", "tlen".
    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path, not required for all statistics types.
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
"""


# noinspection PyShadowingBuiltins
def stat_binned(type,
                alignmentfile,
                fafile=None,
                chrom=None,
                start=None,
                end=None,
                one_based=False,
                window_size=300,
                window_offset=None):
    """Generate statistics per genome window, based on all reads whose alignment starts within
    the window.
    {params}
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


stat_binned.__doc__ = stat_binned.__doc__.format(params=_doc_params)


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
    {params}
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


load_binned.__doc__ = load_binned.__doc__.format(params=_doc_params)


statsobj_binned = {
    'coverage': opt.CoverageBinned(),
    'coverage_ext': opt.CoverageExtBinned(),
    'mapq': opt.MapqBinned(),
    'alignment': opt.AlignmentBinned(),
    'tlen': opt.TlenBinned(),
}


dtypes_binned = {
    'coverage': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('gc', 'u1'),
        ('reads_all', 'i4'),
        ('reads_pp', 'i4')
    ],
    'coverage_ext': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('gc', 'u1'),
        ('reads_all', 'i4'),
        ('reads_pp', 'i4'),
        ('reads_mate_unmapped', 'i4'),
        ('reads_mate_other_chr', 'i4'),
        ('reads_mate_same_strand', 'i4'),
        ('reads_faceaway', 'i4'),
        ('reads_softclipped', 'i4'),
        ('reads_duplicate', 'i4')
    ],
    'mapq': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_mapq0', 'i4'),
        ('rms_mapq', 'i4'),
    ],
    'alignment': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('bases_all', 'i4'),
        ('M', 'i4'),
        ('I', 'i4'),
        ('D', 'i4'),
        ('N', 'i4'),
        ('S', 'i4'),
        ('H', 'i4'),
        ('P', 'i4'),
        ('=', 'i4'),
        ('X', 'i4')
    ],
    'tlen': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_pp', 'i4'),
        ('mean_tlen', 'i4'),
        ('mean_tlen_pp', 'i4'),
        ('rms_tlen', 'i4'),
        ('rms_tlen_pp', 'i4'),
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
stat_coverage_ext_binned = functools.partial(stat_binned, 'coverage_ext')
stat_coverage_ext_binned.__doc__ = _stat_doc
load_coverage_ext_binned = functools.partial(load_binned, 'coverage_ext')
load_coverage_ext_binned.__doc__ = _load_doc
stat_mapq_binned = functools.partial(stat_binned, 'mapq')
stat_mapq_binned.__doc__ = _stat_doc
load_mapq_binned = functools.partial(load_binned, 'mapq')
load_mapq_binned.__doc__ = _load_doc
stat_alignment_binned = functools.partial(stat_binned, 'alignment')
stat_alignment_binned.__doc__ = _stat_doc
load_alignment_binned = functools.partial(load_binned, 'alignment')
load_alignment_binned.__doc__ = _load_doc
stat_tlen_binned = functools.partial(stat_binned, 'tlen')
stat_tlen_binned.__doc__ = _stat_doc
load_tlen_binned = functools.partial(load_binned, 'tlen')
load_tlen_binned.__doc__ = _load_doc
