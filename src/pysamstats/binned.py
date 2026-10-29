# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import functools


import pysamstats.opt as opt
import pysamstats.util as util
import pysamstats.config as config


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
    min_mapq : int, optional
        Only reads with mapping quality equal to or greater than this value will be counted (0
        by default).
    no_dup : bool, optional
        If True, don't count reads flagged as duplicate."""



# noinspection PyShadowingBuiltins
def stat_binned(type,
                alignmentfile,
                fafile=None,
                chrom=None,
                start=None,
                end=None,
                one_based=False,
                window_size=300,
                window_offset=None,
                min_mapq=0,
                no_dup=False):
    """Generate statistics per genome window, based on all reads whose alignment starts within
    the window.
    {params}

    Returns
    -------
    recs : iterator
        An iterator yielding dict objects, where each dict holds data for a single window.

    """

    try:
        stat = stats_classes_binned[type]()
    except KeyError:
        raise ValueError('unsupported statistics type: %r' % type)

    return opt.iter_binned(stat, alignmentfile=alignmentfile, fafile=fafile, chrom=chrom,
                           start=start, end=end, one_based=one_based, window_size=window_size,
                           window_offset=window_offset, min_mapq=min_mapq, no_dup=no_dup)


stat_binned.__doc__ = stat_binned.__doc__.format(params=_doc_params)


# noinspection PyShadowingBuiltins
def load_binned(type,
                alignmentfile,
                fafile=None,
                chrom=None,
                start=None,
                end=None,
                one_based=False,
                window_size=300,
                window_offset=None,
                min_mapq=0,
                no_dup=False,
                dtype=None,
                fields=None):
    """Load statistics per genome window, based on all reads whose alignment starts within
    the window.
    {params}
    dtype : dtype
        Override default dtype.
    fields : string or list of strings
        Select a subset of fields to load.

    Returns
    -------
    ra : numpy structured array
        A structured array.

    """

    statfun = functools.partial(stat_binned, type)
    try:
        default_dtype = getattr(config, 'dtype_' + type + '_binned')
    except KeyError:
        raise ValueError('unsupported statistics type: %r' % type)

    return util.load_stats(statfun, user_dtype=dtype, default_dtype=default_dtype,
                           user_fields=fields, alignmentfile=alignmentfile, fafile=fafile,
                           chrom=chrom, start=start, end=end, one_based=one_based,
                           window_size=window_size, window_offset=window_offset,
                           min_mapq=min_mapq, no_dup=no_dup)


load_binned.__doc__ = load_binned.__doc__.format(params=_doc_params)


stats_classes_binned = {
    'coverage': opt.CoverageBinned,
    'coverage_ext': opt.CoverageExtBinned,
    'mapq': opt.MapqBinned,
    'alignment': opt.AlignmentBinned,
    'tlen': opt.TlenBinned,
}


# backwards compatibility
#########################


_stat_doc_lines = stat_binned.__doc__.split('\n')
_load_doc_lines = load_binned.__doc__.split('\n')
# strip "type" parameter
_stat_doc = '\n'.join(_stat_doc_lines[:5] + _stat_doc_lines[7:])
_load_doc = '\n'.join(_load_doc_lines[:5] + _load_doc_lines[7:])


def _specialize(type):
    statfun = functools.partial(stat_binned, type)
    statfun.__doc__ = _stat_doc
    statfun.__name__ = 'stat_' + type
    loadfun = functools.partial(load_binned, type)
    loadfun.__doc__ = _load_doc
    loadfun.__name__ = 'load_' + type
    return statfun, loadfun


# named functions
stat_coverage_binned, load_coverage_binned = _specialize('coverage')
stat_coverage_ext_binned, load_coverage_ext_binned = _specialize('coverage_ext')
stat_mapq_binned, load_mapq_binned = _specialize('mapq')
stat_alignment_binned, load_alignment_binned = _specialize('alignment')
stat_tlen_binned, load_tlen_binned = _specialize('tlen')
