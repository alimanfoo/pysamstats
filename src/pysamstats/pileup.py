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
        Statistics type. One of "coverage", "coverage_strand", "coverage_ext",
        "coverage_ext_strand", "variation", "variation_strand", "tlen", "tlen_strand", "mapq",
        "mapq_strand", "baseq", "baseq_strand", "baseq_ext", "baseq_ext_strand", "coverage_gc".
    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path.
    fafile : pysam.FastaFile or string, optional
        FASTA file or file path, only required for some statistics types.
    chrom : string, optional
        Chromosome/contig.
    start : int, optional
        Start position.
    end : int, optional
        End position.
    one_based : bool, optional
        Coordinate system, False if zero-based (default), True if one-based.
    truncate : bool, optional
        If True, truncate output to selected region.
    pad : bool, optional
        If True, emit records for every position, even if no reads are aligned.
    max_depth : int, optional
        Maximum depth to allow in pileup column.
    window_size : int, optional
        Window size to use for percent GC calculation (only applies to coverage_gc).
    window_offset : int, optional
        Distance from window start to record position (only applies to coverage_gc).
    min_mapq : int, optional
        Only reads with mapping quality equal to or greater than this value will be counted (0
        by default).
    min_baseq : int, optional
        Only reads with base quality equal to or greater than this value will be counted (0 by
        default).
    no_del : bool, optional
        If True, don't count reads aligned with a deletion at the current position.
    no_dup : bool, optional
        If True, don't count reads flagged as duplicate."""


# noinspection PyShadowingBuiltins
def stat_pileup(type,
                alignmentfile,
                fafile=None,
                chrom=None,
                start=None,
                end=None,
                one_based=False,
                truncate=False,
                stepper="all",
                pad=False,
                max_depth=8000,
                window_size=300,
                window_offset=None,
                min_mapq=0,
                min_baseq=0,
                no_del=False,
                no_dup=False):
    """Generate statistics per genome position, based on read pileups.
    {params}

    Returns
    -------
    recs : iterator
        An iterator yielding dict objects, where each dict holds data for a single genome position.

    """

    if type in config.stats_types_withref and fafile is None:
        raise ValueError('reference sequence is required; please provide fafile argument')

    try:
        if type == 'coverage_gc':
            stat = stats_classes_pileup[type](window_size=window_size, window_offset=window_offset)
        else:
            stat = stats_classes_pileup[type]()
    except KeyError:
        raise ValueError('unsupported statistics type: %r' % type)

    return opt.iter_pileup(stat, alignmentfile=alignmentfile, fafile=fafile, chrom=chrom,
                           start=start, end=end, one_based=one_based, truncate=truncate, stepper=stepper, pad=pad,
                           max_depth=max_depth, min_mapq=min_mapq, min_baseq=min_baseq,
                           no_del=no_del, no_dup=no_dup)


stat_pileup.__doc__ = stat_pileup.__doc__.format(params=_doc_params)


# noinspection PyShadowingBuiltins
def load_pileup(type,
                alignmentfile,
                fafile=None,
                chrom=None,
                start=None,
                end=None,
                one_based=False,
                truncate=False,
                stepper="all",
                pad=False,
                max_depth=8000,
                window_size=300,
                window_offset=None,
                min_mapq=0,
                min_baseq=0,
                no_del=False,
                no_dup=False,
                dtype=None,
                fields=None):
    """Load statistics per genome position, based on read pileups.
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

    statfun = functools.partial(stat_pileup, type)
    try:
        default_dtype = getattr(config, 'dtype_' + type)
    except AttributeError:
        raise ValueError('unsupported statistics type: %r' % type)

    return util.load_stats(statfun, user_dtype=dtype, default_dtype=default_dtype,
                           user_fields=fields, alignmentfile=alignmentfile, fafile=fafile,
                           chrom=chrom, start=start, end=end, one_based=one_based,
                           truncate=truncate, stepper=stepper, pad=pad, max_depth=max_depth, window_size=window_size,
                           window_offset=window_offset, min_mapq=min_mapq, min_baseq=min_baseq, no_del=no_del, no_dup=no_dup)


load_pileup.__doc__ = load_pileup.__doc__.format(params=_doc_params)


stats_classes_pileup = {
    'coverage': opt.Coverage,
    'coverage_strand': opt.CoverageStrand,
    'coverage_ext': opt.CoverageExt,
    'coverage_ext_strand': opt.CoverageExtStrand,
    'variation': opt.Variation,
    'variation_strand': opt.VariationStrand,
    'tlen': opt.Tlen,
    'tlen_strand': opt.TlenStrand,
    'mapq': opt.Mapq,
    'mapq_strand': opt.MapqStrand,
    'baseq': opt.Baseq,
    'baseq_strand': opt.BaseqStrand,
    'baseq_ext': opt.BaseqExt,
    'baseq_ext_strand': opt.BaseqExtStrand,
    'coverage_gc': opt.CoverageGC,
}


# backwards compatibility
#########################


_stat_doc_lines = stat_pileup.__doc__.split('\n')
_load_doc_lines = load_pileup.__doc__.split('\n')
# strip "type" parameter
_stat_doc = '\n'.join(_stat_doc_lines[:4] + _stat_doc_lines[8:])
_load_doc = '\n'.join(_load_doc_lines[:4] + _load_doc_lines[8:])


def _specialize(type):
    stat = functools.partial(stat_pileup, type)
    stat.__doc__ = _stat_doc
    stat.__name__ = 'stat_' + type
    load = functools.partial(load_pileup, type)
    load.__doc__ = _load_doc
    load.__name__ = 'load_' + type
    return stat, load


# named functions
stat_coverage, load_coverage = _specialize('coverage')
stat_coverage_strand, load_coverage_strand = _specialize('coverage_strand')
stat_coverage_ext, load_coverage_ext = _specialize('coverage_ext')
stat_coverage_ext_strand, load_coverage_ext_strand = _specialize('coverage_ext_strand')
stat_variation, load_variation = _specialize('variation')
stat_variation_strand, load_variation_strand = _specialize('variation_strand')
stat_tlen, load_tlen = _specialize('tlen')
stat_tlen_strand, load_tlen_strand = _specialize('tlen_strand')
stat_mapq, load_mapq = _specialize('mapq')
stat_mapq_strand, load_mapq_strand = _specialize('mapq_strand')
stat_baseq, load_baseq = _specialize('baseq')
stat_baseq_strand, load_baseq_strand = _specialize('baseq_strand')
stat_baseq_ext, load_baseq_ext = _specialize('baseq_ext')
stat_baseq_ext_strand, load_baseq_ext_strand = _specialize('baseq_ext_strand')
stat_coverage_gc, load_coverage_gc = _specialize('coverage_gc')
