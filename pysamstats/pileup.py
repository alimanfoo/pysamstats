# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division


import pysamstats.opt as opt
import pysamstats.util as util


_doc_stat_params_noref = """
    Parameters
    ----------
    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path.
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
"""


_doc_stat_returns = """
    Returns
    -------
    recs : iterator
        An iterator yielding dict objects, where each dict holds data for a single genome position.
"""


_doc_load_returns = """
    Returns
    -------
    ra : NumPy structured array
        A structured array.
"""


# coverage
##########


dtype_coverage = [('chrom', 'a12'),
                  ('pos', 'i4'),
                  ('reads_all', 'i4'),
                  ('reads_pp', 'i4')]


fields_coverage = [t[0] for t in dtype_coverage]


_doc_coverage = """basic coverage statistics"""


def stat_coverage(*args, **kwargs):
    """Generate {info} per genome position.
    {params}{returns}
    """
    return opt.iter_pileup(opt.rec_coverage, opt.rec_coverage_pad, *args, **kwargs)


stat_coverage.__doc__ = stat_coverage.__doc__.format(
    info=_doc_coverage,
    params=_doc_stat_params_noref,
    returns=_doc_stat_returns
)


def load_coverage(*args, **kwargs):
    """Load {info} per genome position.
    {params}{returns}
    """
    return util.load_stats(stat_coverage, dtype_coverage, *args, **kwargs)


load_coverage.__doc__ = load_coverage.__doc__.format(
    info=_doc_coverage,
    params=_doc_stat_params_noref,
    returns=_doc_load_returns
)


# coverage_strand
#################


dtype_coverage_strand = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_fwd', 'i4'),
    ('reads_rev', 'i4'),
    ('reads_pp', 'i4'),
    ('reads_pp_fwd', 'i4'),
    ('reads_pp_rev', 'i4'),
]


fields_coverage_strand = [t[0] for t in dtype_coverage_strand]


_doc_coverage_strand = """basic coverage statistics by strand"""


def stat_coverage_strand(*args, **kwargs):
    """Generate {info} per genome position.
    {params}{returns}
    """
    return opt.iter_pileup(opt.rec_coverage_strand, opt.rec_coverage_strand_pad, *args, **kwargs)


stat_coverage_strand.__doc__ = stat_coverage_strand.__doc__.format(
    info=_doc_coverage_strand,
    params=_doc_stat_params_noref,
    returns=_doc_stat_returns
)


def load_coverage_strand(*args, **kwargs):
    """Load {info} per genome position.
    {params}{returns}
    """
    return util.load_stats(stat_coverage_strand, dtype_coverage_strand, *args, **kwargs)


load_coverage_strand.__doc__ = load_coverage_strand.__doc__.format(
    info=_doc_coverage_strand,
    params=_doc_stat_params_noref,
    returns=_doc_load_returns
)


# coverage_ext
##############


dtype_coverage_ext = [
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
]


fields_coverage_ext = [t[0] for t in dtype_coverage_ext]


_doc_coverage_ext = """extended coverage statistics"""


def stat_coverage_ext(*args, **kwargs):
    """Generate {info} per genome position.
    {params}{returns}
    """
    return opt.iter_pileup(opt.rec_coverage_ext, opt.rec_coverage_ext_pad, *args, **kwargs)


stat_coverage_ext.__doc__ = stat_coverage_ext.__doc__.format(
    info=_doc_coverage_ext,
    params=_doc_stat_params_noref,
    returns=_doc_stat_returns
)


def load_coverage_ext(*args, **kwargs):
    """Load {info} per genome position.
    {params}{returns}
    """
    return util.load_stats(stat_coverage_ext, dtype_coverage_ext, *args, **kwargs)


load_coverage_ext.__doc__ = load_coverage_ext.__doc__.format(
    info=_doc_coverage_ext,
    params=_doc_stat_params_noref,
    returns=_doc_load_returns
)


# coverage_ext_strand
#####################


dtype_coverage_ext_strand = [
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
]


fields_coverage_ext_strand = [t[0] for t in dtype_coverage_ext_strand]


_doc_coverage_ext_strand = """extended coverage statistics by strand"""


def stat_coverage_ext_strand(*args, **kwargs):
    """Generate {info} per genome position.
    {params}{returns}
    """
    return opt.iter_pileup(opt.rec_coverage_ext_strand, opt.rec_coverage_ext_strand_pad, *args,
                           **kwargs)


stat_coverage_ext_strand.__doc__ = stat_coverage_ext_strand.__doc__.format(
    info=_doc_coverage_ext_strand,
    params=_doc_stat_params_noref,
    returns=_doc_stat_returns
)


def load_coverage_ext_strand(*args, **kwargs):
    """Load {info} per genome position.
    {params}{returns}
    """
    return util.load_stats(stat_coverage_ext_strand, dtype_coverage_ext_strand, *args, **kwargs)


load_coverage_ext_strand.__doc__ = load_coverage_ext_strand.__doc__.format(
    info=_doc_coverage_ext_strand,
    params=_doc_stat_params_noref,
    returns=_doc_load_returns
)


# variation
###########


dtype_variation = [
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
]


fields_variation = [t[0] for t in dtype_variation]


_doc_variation = """variation statistics"""


def stat_variation(*args, **kwargs):
    """Generate {info} per genome position.
    {params}{returns}
    """
    return opt.iter_pileup(opt.rec_variation, opt.rec_variation_pad, *args, **kwargs)


stat_variation.__doc__ = stat_variation.__doc__.format(
    info=_doc_variation,
    params=_doc_stat_params_noref,
    returns=_doc_stat_returns
)


def load_variation(*args, **kwargs):
    """Load {info} per genome position.
    {params}{returns}
    """
    return util.load_stats(stat_variation, dtype_variation, *args, **kwargs)


load_variation.__doc__ = load_variation.__doc__.format(
    info=_doc_variation,
    params=_doc_stat_params_noref,
    returns=_doc_load_returns
)


# variation_strand
##################


dtype_variation_strand = [
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
]


fields_variation_strand = [t[0] for t in dtype_variation_strand]


_doc_variation_strand = """variation statistics by strand"""


def stat_variation_strand(*args, **kwargs):
    """Generate {info} per genome position.
    {params}{returns}
    """
    return opt.iter_pileup(opt.rec_variation_strand, opt.rec_variation_strand_pad, *args, **kwargs)


stat_variation_strand.__doc__ = stat_variation_strand.__doc__.format(
    info=_doc_variation_strand,
    params=_doc_stat_params_noref,
    returns=_doc_stat_returns
)


def load_variation_strand(*args, **kwargs):
    """Load {info} per genome position.
    {params}{returns}
    """
    return util.load_stats(stat_variation_strand, dtype_variation_strand, *args, **kwargs)


load_variation_strand.__doc__ = load_variation_strand.__doc__.format(
    info=_doc_variation_strand,
    params=_doc_stat_params_noref,
    returns=_doc_load_returns
)
