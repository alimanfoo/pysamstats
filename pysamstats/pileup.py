# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import functools


import pysamstats.opt as opt
import pysamstats.util as util


# noinspection PyShadowingBuiltins
def stat_pileup(type,
                alignmentfile,
                fafile=None,
                chrom=None,
                start=None,
                end=None,
                one_based=False,
                truncate=False,
                pad=False,
                max_depth=8000,
                window_size=300,
                window_offset=None):
    """Generate statistics per genome position, based on read pileups.

    Parameters
    ----------
    type : string
        Statistics type. One of "coverage", "coverage_strand", "coverage_ext",
        "coverage_ext_strand", "variation", "variation_strand", "tlen", "tlen_strand", "mapq",
        "mapq_strand", "baseq", "baseq_strand", "baseq_ext", "baseq_ext_strand", @@TODO
    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path.
    fafile : pysam.FastaFile or string
        FASTA file or file path, only required for some statistics types.
    chrom : string
        Chromosome/contig.
    start : int
        Start position.
    end : int
        End position.
    one_based : bool
        Coordinate system, False if zero-based (default), True if one-based.
    truncate : bool
        If True, truncate output to selected region.
    pad : bool
        If True, emit records for every position, even if no reads are aligned.
    max_depth : int
        Maximum depth to allow in pileup column.
    window_size : int
        Window size to use for percent GC calculation (only applies to coverage_gc).
    window_offset : int
        Distance from window start to record position (only applies to coverage_gc).

    Returns
    -------
    recs : iterator
        An iterator yielding dict objects, where each dict holds data for a single genome position.

    """

    try:
        if type == 'coverage_gc':
            # special case needed to handle window parameters
            rec, rec_pad = opt.frecs_coverage_gc(window_size=window_size,
                                                 window_offset=window_offset)
        else:
            rec, rec_pad = frecs[type]
    except KeyError:
        raise ValueError('unsupported statistics type: %r' % type)

    return opt.iter_pileup(rec, rec_pad, alignmentfile=alignmentfile, fafile=fafile, chrom=chrom,
                           start=start, end=end, one_based=one_based, truncate=truncate, pad=pad,
                           max_depth=max_depth)


# noinspection PyShadowingBuiltins
def load_pileup(type,
                alignmentfile,
                fafile=None,
                chrom=None,
                start=None,
                end=None,
                one_based=False,
                truncate=False,
                pad=False,
                max_depth=8000,
                window_size=300,
                window_offset=None):
    """Load statistics per genome position, based on read pileups.

    Parameters
    ----------
    type : string
        Statistics type. One of "coverage", "coverage_strand", "coverage_ext",
        "coverage_ext_strand", @@TODO
    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path.
    fafile : pysam.FastaFile or string
        FASTA file or file path, only required for some statistics types.
    chrom : string
        Chromosome/contig.
    start : int
        Start position.
    end : int
        End position.
    one_based : bool
        Coordinate system, False if zero-based (default), True if one-based.
    truncate : bool
        If True, truncate output to selected region.
    pad : bool
        If True, emit records for every position, even if no reads are aligned.
    max_depth : int
        Maximum depth to allow in pileup column.
    window_size : int
        Window size to use for percent GC calculation (only applies to coverage_gc).
    window_offset : int
        Distance from window start to record position (only applies to coverage_gc).

    Returns
    -------
    ra : numpy structured array
        A structured array.

    """

    stat = functools.partial(stat_pileup, type)
    try:
        dtype = dtypes[type]
    except KeyError:
        raise ValueError('unsupported statistics type: %r' % type)

    return util.load_stats(stat, dtype, alignmentfile=alignmentfile, fafile=fafile, chrom=chrom,
                           start=start, end=end, one_based=one_based, truncate=truncate, pad=pad,
                           max_depth=max_depth, window_size=window_size,
                           window_offset=window_offset)


frecs = {
    'coverage': (opt.rec_coverage, opt.rec_coverage_pad),
    'coverage_strand': (opt.rec_coverage_strand, opt.rec_coverage_strand_pad),
    'coverage_ext': (opt.rec_coverage_ext, opt.rec_coverage_ext_pad),
    'coverage_ext_strand': (opt.rec_coverage_ext_strand, opt.rec_coverage_ext_strand_pad),
    'variation': (opt.rec_variation, opt.rec_variation_pad),
    'variation_strand': (opt.rec_variation_strand, opt.rec_variation_strand_pad),
    'tlen': (opt.rec_tlen, opt.rec_tlen_pad),
    'tlen_strand': (opt.rec_tlen_strand, opt.rec_tlen_strand_pad),
    'mapq': (opt.rec_mapq, opt.rec_mapq_pad),
    'mapq_strand': (opt.rec_mapq_strand, opt.rec_mapq_strand_pad),
    'baseq': (opt.rec_baseq, opt.rec_baseq_pad),
    'baseq_strand': (opt.rec_baseq_strand, opt.rec_baseq_strand_pad),
    'baseq_ext': (opt.rec_baseq_ext, opt.rec_baseq_ext_pad),
    'baseq_ext_strand': (opt.rec_baseq_ext_strand, opt.rec_baseq_ext_strand_pad),
}


dtypes = {
    'coverage': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_pp', 'i4')
    ],
    'coverage_strand': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_fwd', 'i4'),
        ('reads_rev', 'i4'),
        ('reads_pp', 'i4'),
        ('reads_pp_fwd', 'i4'),
        ('reads_pp_rev', 'i4'),
    ],
    'coverage_ext': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_pp', 'i4'),
        ('reads_mate_unmapped', 'i4'),
        ('reads_mate_other_chr', 'i4'),
        ('reads_mate_same_strand', 'i4'),
        ('reads_faceaway', 'i4'),
        ('reads_softclipped', 'i4'),
        ('reads_duplicate', 'i4')
    ],
    'coverage_ext_strand': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_fwd', 'i4'),
        ('reads_rev', 'i4'),
        ('reads_pp', 'i4'),
        ('reads_pp_fwd', 'i4'),
        ('reads_pp_rev', 'i4'),
        ('reads_mate_unmapped', 'i4'),
        ('reads_mate_unmapped_fwd', 'i4'),
        ('reads_mate_unmapped_rev', 'i4'),
        ('reads_mate_other_chr', 'i4'),
        ('reads_mate_other_chr_fwd', 'i4'),
        ('reads_mate_other_chr_rev', 'i4'),
        ('reads_mate_same_strand', 'i4'),
        ('reads_mate_same_strand_fwd', 'i4'),
        ('reads_mate_same_strand_rev', 'i4'),
        ('reads_faceaway', 'i4'),
        ('reads_faceaway_fwd', 'i4'),
        ('reads_faceaway_rev', 'i4'),
        ('reads_softclipped', 'i4'),
        ('reads_softclipped_fwd', 'i4'),
        ('reads_softclipped_rev', 'i4'),
        ('reads_duplicate', 'i4'),
        ('reads_duplicate_fwd', 'i4'),
        ('reads_duplicate_rev', 'i4'),
    ],
    'variation': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('ref', 'a1'),
        ('reads_all', 'i4'),
        ('reads_pp', 'i4'),
        ('matches', 'i4'),
        ('matches_pp', 'i4'),
        ('mismatches', 'i4'),
        ('mismatches_pp', 'i4'),
        ('deletions', 'i4'),
        ('deletions_pp', 'i4'),
        ('insertions', 'i4'),
        ('insertions_pp', 'i4'),
        ('A', 'i4'),
        ('A_pp', 'i4'),
        ('C', 'i4'),
        ('C_pp', 'i4'),
        ('T', 'i4'),
        ('T_pp', 'i4'),
        ('G', 'i4'),
        ('G_pp', 'i4'),
        ('N', 'i4'),
        ('N_pp', 'i4')
    ],
    'variation_strand': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('ref', 'a1'),
        ('reads_all', 'i4'),
        ('reads_fwd', 'i4'),
        ('reads_rev', 'i4'),
        ('reads_pp', 'i4'),
        ('reads_pp_fwd', 'i4'),
        ('reads_pp_rev', 'i4'),
        ('matches', 'i4'),
        ('matches_fwd', 'i4'),
        ('matches_rev', 'i4'),
        ('matches_pp', 'i4'),
        ('matches_pp_fwd', 'i4'),
        ('matches_pp_rev', 'i4'),
        ('mismatches', 'i4'),
        ('mismatches_fwd', 'i4'),
        ('mismatches_rev', 'i4'),
        ('mismatches_pp', 'i4'),
        ('mismatches_pp_fwd', 'i4'),
        ('mismatches_pp_rev', 'i4'),
        ('deletions', 'i4'),
        ('deletions_fwd', 'i4'),
        ('deletions_rev', 'i4'),
        ('deletions_pp', 'i4'),
        ('deletions_pp_fwd', 'i4'),
        ('deletions_pp_rev', 'i4'),
        ('insertions', 'i4'),
        ('insertions_fwd', 'i4'),
        ('insertions_rev', 'i4'),
        ('insertions_pp', 'i4'),
        ('insertions_pp_fwd', 'i4'),
        ('insertions_pp_rev', 'i4'),
        ('A', 'i4'), ('A_fwd', 'i4'), ('A_rev', 'i4'),
        ('A_pp', 'i4'), ('A_pp_fwd', 'i4'), ('A_pp_rev', 'i4'),
        ('C', 'i4'), ('C_fwd', 'i4'), ('C_rev', 'i4'),
        ('C_pp', 'i4'), ('C_pp_fwd', 'i4'), ('C_pp_rev', 'i4'),
        ('T', 'i4'), ('T_fwd', 'i4'), ('T_rev', 'i4'),
        ('T_pp', 'i4'), ('T_pp_fwd', 'i4'), ('T_pp_rev', 'i4'),
        ('G', 'i4'), ('G_fwd', 'i4'), ('G_rev', 'i4'),
        ('G_pp', 'i4'), ('G_pp_fwd', 'i4'), ('G_pp_rev', 'i4'),
        ('N', 'i4'), ('N_fwd', 'i4'), ('N_rev', 'i4'),
        ('N_pp', 'i4'), ('N_pp_fwd', 'i4'), ('N_pp_rev', 'i4')
    ],
    'tlen': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_paired', 'i4'),
        ('reads_pp', 'i4'),
        ('mean_tlen', 'i4'),
        ('mean_tlen_pp', 'i4'),
        ('rms_tlen', 'i4'),
        ('rms_tlen_pp', 'i4'),
        ('std_tlen', 'i4'),
        ('std_tlen_pp', 'i4')
    ],
    'tlen_strand': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_fwd', 'i4'),
        ('reads_rev', 'i4'),
        ('reads_paired', 'i4'),
        ('reads_paired_fwd', 'i4'),
        ('reads_paired_rev', 'i4'),
        ('reads_pp', 'i4'),
        ('reads_pp_fwd', 'i4'),
        ('reads_pp_rev', 'i4'),
        ('mean_tlen', 'i4'),
        ('mean_tlen_fwd', 'i4'),
        ('mean_tlen_rev', 'i4'),
        ('mean_tlen_pp', 'i4'),
        ('mean_tlen_pp_fwd', 'i4'),
        ('mean_tlen_pp_rev', 'i4'),
        ('rms_tlen', 'i4'),
        ('rms_tlen_fwd', 'i4'),
        ('rms_tlen_rev', 'i4'),
        ('rms_tlen_pp', 'i4'),
        ('rms_tlen_pp_fwd', 'i4'),
        ('rms_tlen_pp_rev', 'i4'),
        ('std_tlen', 'i4'),
        ('std_tlen_fwd', 'i4'),
        ('std_tlen_rev', 'i4'),
        ('std_tlen_pp', 'i4'),
        ('std_tlen_pp_fwd', 'i4'),
        ('std_tlen_pp_rev', 'i4')
    ],
    'mapq': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_pp', 'i4'),
        ('reads_mapq0', 'i4'),
        ('reads_mapq0_pp', 'i4'),
        ('rms_mapq', 'i4'),
        ('rms_mapq_pp', 'i4'),
        ('max_mapq', 'i4'),
        ('max_mapq_pp', 'i4')
    ],
    'mapq_strand': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_fwd', 'i4'),
        ('reads_rev', 'i4'),
        ('reads_pp', 'i4'),
        ('reads_pp_fwd', 'i4'),
        ('reads_pp_rev', 'i4'),
        ('reads_mapq0', 'i4'),
        ('reads_mapq0_fwd', 'i4'),
        ('reads_mapq0_rev', 'i4'),
        ('reads_mapq0_pp', 'i4'),
        ('reads_mapq0_pp_fwd', 'i4'),
        ('reads_mapq0_pp_rev', 'i4'),
        ('rms_mapq', 'i4'),
        ('rms_mapq_fwd', 'i4'),
        ('rms_mapq_rev', 'i4'),
        ('rms_mapq_pp', 'i4'),
        ('rms_mapq_pp_fwd', 'i4'),
        ('rms_mapq_pp_rev', 'i4'),
        ('max_mapq', 'i4'),
        ('max_mapq_fwd', 'i4'),
        ('max_mapq_rev', 'i4'),
        ('max_mapq_pp', 'i4'),
        ('max_mapq_pp_fwd', 'i4'),
        ('max_mapq_pp_rev', 'i4'),
    ],
    'baseq': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_pp', 'i4'),
        ('rms_baseq', 'i4'),
        ('rms_baseq_pp', 'i4'),
    ],
    'baseq_strand': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('reads_all', 'i4'),
        ('reads_fwd', 'i4'),
        ('reads_rev', 'i4'),
        ('reads_pp', 'i4'),
        ('reads_pp_fwd', 'i4'),
        ('reads_pp_rev', 'i4'),
        ('rms_baseq', 'i4'),
        ('rms_baseq_fwd', 'i4'),
        ('rms_baseq_rev', 'i4'),
        ('rms_baseq_pp', 'i4'),
        ('rms_baseq_pp_fwd', 'i4'),
        ('rms_baseq_pp_rev', 'i4'),
    ],
    'baseq_ext': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('ref', 'a1'),
        ('reads_all', 'i4'),
        ('reads_pp', 'i4'),
        ('matches', 'i4'),
        ('matches_pp', 'i4'),
        ('mismatches', 'i4'),
        ('mismatches_pp', 'i4'),
        ('rms_baseq', 'i4'),
        ('rms_baseq_pp', 'i4'),
        ('rms_baseq_matches', 'i4'),
        ('rms_baseq_matches_pp', 'i4'),
        ('rms_baseq_mismatches', 'i4'),
        ('rms_baseq_mismatches_pp', 'i4'),
    ],
    'baseq_ext_strand': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('ref', 'a1'),
        ('reads_all', 'i4'),
        ('reads_fwd', 'i4'),
        ('reads_rev', 'i4'),
        ('reads_pp', 'i4'),
        ('reads_pp_fwd', 'i4'),
        ('reads_pp_rev', 'i4'),
        ('matches', 'i4'),
        ('matches_fwd', 'i4'),
        ('matches_rev', 'i4'),
        ('matches_pp', 'i4'),
        ('matches_pp_fwd', 'i4'),
        ('matches_pp_rev', 'i4'),
        ('mismatches', 'i4'),
        ('mismatches_fwd', 'i4'),
        ('mismatches_rev', 'i4'),
        ('mismatches_pp', 'i4'),
        ('mismatches_pp_fwd', 'i4'),
        ('mismatches_pp_rev', 'i4'),
        ('rms_baseq', 'i4'),
        ('rms_baseq_fwd', 'i4'),
        ('rms_baseq_rev', 'i4'),
        ('rms_baseq_pp', 'i4'),
        ('rms_baseq_pp_fwd', 'i4'),
        ('rms_baseq_pp_rev', 'i4'),
        ('rms_baseq_matches', 'i4'),
        ('rms_baseq_matches_fwd', 'i4'),
        ('rms_baseq_matches_rev', 'i4'),
        ('rms_baseq_matches_pp', 'i4'),
        ('rms_baseq_matches_pp_fwd', 'i4'),
        ('rms_baseq_matches_pp_rev', 'i4'),
        ('rms_baseq_mismatches', 'i4'),
        ('rms_baseq_mismatches_fwd', 'i4'),
        ('rms_baseq_mismatches_rev', 'i4'),
        ('rms_baseq_mismatches_pp', 'i4'),
        ('rms_baseq_mismatches_pp_fwd', 'i4'),
        ('rms_baseq_mismatches_pp_rev', 'i4')
    ],
    'coverage_gc': [
        ('chrom', 'a12'),
        ('pos', 'i4'),
        ('gc', 'u1'),
        ('reads_all', 'i4'),
        ('reads_pp', 'i4')
    ],
}


# backwards compatibility
#########################


_stat_doc_lines = stat_pileup.__doc__.split('\n')
_load_doc_lines = load_pileup.__doc__.split('\n')
# strip "type" parameter
_stat_doc = '\n'.join(_stat_doc_lines[:4] + _stat_doc_lines[8:])
_load_doc = '\n'.join(_load_doc_lines[:4] + _stat_doc_lines[8:])
# named functions
stat_coverage = functools.partial(stat_pileup, 'coverage')
stat_coverage.__doc__ = _stat_doc
load_coverage = functools.partial(load_pileup, 'coverage')
load_coverage.__doc__ = _load_doc
stat_coverage_strand = functools.partial(stat_pileup, 'coverage_strand')
stat_coverage_strand.__doc__ = _stat_doc
load_coverage_strand = functools.partial(load_pileup, 'coverage_strand')
load_coverage_strand.__doc__ = _load_doc
stat_coverage_ext = functools.partial(stat_pileup, 'coverage_ext')
stat_coverage_ext.__doc__ = _stat_doc
load_coverage_ext = functools.partial(load_pileup, 'coverage_ext')
load_coverage_ext.__doc__ = _load_doc
stat_coverage_ext_strand = functools.partial(stat_pileup, 'coverage_ext_strand')
stat_coverage_ext_strand.__doc__ = _stat_doc
load_coverage_ext_strand = functools.partial(load_pileup, 'coverage_ext_strand')
load_coverage_ext_strand.__doc__ = _load_doc
stat_variation = functools.partial(stat_pileup, 'variation')
stat_variation.__doc__ = _stat_doc
load_variation = functools.partial(load_pileup, 'variation')
load_variation.__doc__ = _load_doc
stat_variation_strand = functools.partial(stat_pileup, 'variation_strand')
stat_variation_strand.__doc__ = _stat_doc
load_variation_strand = functools.partial(load_pileup, 'variation_strand')
load_variation_strand.__doc__ = _load_doc
stat_tlen = functools.partial(stat_pileup, 'tlen')
stat_tlen.__doc__ = _stat_doc
load_tlen = functools.partial(load_pileup, 'tlen')
load_tlen.__doc__ = _load_doc
stat_tlen_strand = functools.partial(stat_pileup, 'tlen_strand')
stat_tlen_strand.__doc__ = _stat_doc
load_tlen_strand = functools.partial(load_pileup, 'tlen_strand')
load_tlen_strand.__doc__ = _load_doc
stat_mapq = functools.partial(stat_pileup, 'mapq')
stat_mapq.__doc__ = _stat_doc
load_mapq = functools.partial(load_pileup, 'mapq')
load_mapq.__doc__ = _load_doc
stat_mapq_strand = functools.partial(stat_pileup, 'mapq_strand')
stat_mapq_strand.__doc__ = _stat_doc
load_mapq_strand = functools.partial(load_pileup, 'mapq_strand')
load_mapq_strand.__doc__ = _load_doc
stat_baseq = functools.partial(stat_pileup, 'baseq')
stat_baseq.__doc__ = _stat_doc
load_baseq = functools.partial(load_pileup, 'baseq')
load_baseq.__doc__ = _load_doc
stat_baseq_strand = functools.partial(stat_pileup, 'baseq_strand')
stat_baseq_strand.__doc__ = _stat_doc
load_baseq_strand = functools.partial(load_pileup, 'baseq_strand')
load_baseq_strand.__doc__ = _load_doc
stat_baseq_ext = functools.partial(stat_pileup, 'baseq_ext')
stat_baseq_ext.__doc__ = _stat_doc
load_baseq_ext = functools.partial(load_pileup, 'baseq_ext')
load_baseq_ext.__doc__ = _load_doc
stat_baseq_ext_strand = functools.partial(stat_pileup, 'baseq_ext_strand')
stat_baseq_ext_strand.__doc__ = _stat_doc
load_baseq_ext_strand = functools.partial(load_pileup, 'baseq_ext_strand')
load_baseq_ext_strand.__doc__ = _load_doc
stat_coverage_gc = functools.partial(stat_pileup, 'coverage_gc')
stat_coverage_gc.__doc__ = _stat_doc
load_coverage_gc = functools.partial(load_pileup, 'coverage_gc')
load_coverage_gc.__doc__ = _load_doc





