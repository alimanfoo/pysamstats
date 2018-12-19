# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
from itertools import chain
from collections import Counter
import logging


from pysam import Samfile, Fastafile
from nose.tools import eq_
from numpy import around as round


import pysamstats
from .util import normalise_coords, mean, rms, rootmean, compare_iterators


logger = logging.getLogger(__name__)
debug = logger.debug


def compare_stats(impl, refimpl):
    # no read filters
    kwargs = {'chrom': 'Pf3D7_01_v3',
              'start': 0,
              'end': 2000,
              'one_based': False}
    expected = refimpl(Samfile('fixture/test.bam'), **kwargs)
    actual = impl(Samfile('fixture/test.bam'), **kwargs)
    compare_iterators(expected, actual)
    # read filters
    kwargs['min_mapq'] = 1
    kwargs['no_dup'] = True
    expected = refimpl(Samfile('fixture/test.bam'), **kwargs)
    actual = impl(Samfile('fixture/test.bam'), **kwargs)
    compare_iterators(expected, actual)


def compare_stats_withref(impl, refimpl, bam_fn='fixture/test.bam',
                          fasta_fn='fixture/ref.fa'):
    # no read filters
    kwargs = {'chrom': 'Pf3D7_01_v3',
              'start': 0,
              'end': 2000,
              'one_based': False}
    expected = refimpl(Samfile(bam_fn), Fastafile(fasta_fn), **kwargs)
    actual = impl(Samfile(bam_fn), Fastafile(fasta_fn), **kwargs)
    compare_iterators(expected, actual)
    # read filters
    kwargs['min_mapq'] = 1
    kwargs['no_dup'] = True
    expected = refimpl(Samfile(bam_fn), Fastafile(fasta_fn), **kwargs)
    actual = impl(Samfile(bam_fn), Fastafile(fasta_fn), **kwargs)
    compare_iterators(expected, actual)


def filter_alignments(alignments, min_mapq, no_dup):
    if min_mapq > 0:
        alignments = (a for a in alignments if a.mapq >= min_mapq)
    if no_dup:
        alignments = (a for a in alignments if not a.is_duplicate)
    return alignments


def stat_coverage_binned_refimpl(samfile, fastafile, chrom=None, start=None,
                                 end=None, one_based=False, window_size=300,
                                 window_offset=150, min_mapq=0, no_dup=False):
    if chrom is None:
        # noinspection PyTypeChecker
        it = chain(*[
            iter_coverage_binned(samfile, fastafile, chrom, None, None, one_based, window_size,
                                 window_offset, min_mapq, no_dup)
            for chrom in samfile.references
        ])
    else:
        it = iter_coverage_binned(samfile, fastafile, chrom, start, end, one_based, window_size,
                                  window_offset, min_mapq, no_dup)
    return it


def gc_content(fastafile, chrom, start, end):
    seq = fastafile.fetch(chrom, start, end).lower()
    nc = Counter(seq)
    gc = int(round((nc['g'] + nc['c']) * 100. / (end-start)))
    return gc


def iter_coverage_binned(samfile, fastafile, chrom, start, end, one_based,
                         window_size, window_offset, min_mapq, no_dup):
    assert chrom is not None
    start, end = normalise_coords(one_based, start, end)
    chrlen = samfile.lengths[samfile.references.index(chrom)]
    if start is None:
        start = 0
    if end is None:
        end = chrlen
    if end > chrlen:
        end = chrlen
    # setup first bin
    bin_start = start
    bin_end = bin_start + window_size
    reads_all = reads_pp = 0

    # iterate over reads
    alignments = samfile.fetch(chrom, start, end)
    alignments = filter_alignments(alignments, min_mapq, no_dup)
    for aln in alignments:
        while aln.pos > bin_end:  # end of bin
            gc = gc_content(fastafile, chrom, bin_start, bin_end)
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'gc': gc, 'reads_all': reads_all,
                   'reads_pp': reads_pp}
            yield rec
            reads_all = reads_pp = 0
            bin_start = bin_end
            bin_end = bin_start + window_size
        if not aln.is_unmapped:
            reads_all += 1
            if aln.is_proper_pair:
                reads_pp += 1

    # deal with last non-empty bin
    gc = gc_content(fastafile, chrom, bin_start, bin_end)
    pos = bin_start + window_offset
    if one_based:
        pos += 1
    rec = {'chrom': chrom, 'pos': pos,
           'gc': gc, 'reads_all': reads_all, 'reads_pp': reads_pp}
    yield rec

    # deal with empty bins up to explicit end
    if end is not None:
        while bin_end < end:
            reads_all = reads_pp = 0
            bin_start = bin_end
            bin_end = bin_start + window_size
            gc = gc_content(fastafile, chrom, bin_start, bin_end)
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'gc': gc, 'reads_all': reads_all,
                   'reads_pp': reads_pp}
            yield rec


def test_stat_coverage_binned():
    compare_stats_withref(pysamstats.stat_coverage_binned,
                          stat_coverage_binned_refimpl)


def test_stat_coverage_binned_uppercase_fasta():
    compare_stats_withref(pysamstats.stat_coverage_binned,
                          stat_coverage_binned_refimpl,
                          fasta_fn='fixture/ref.upper.fa')


def stat_coverage_ext_binned_refimpl(samfile, fastafile, chrom=None, start=None, end=None,
                                     one_based=False, window_size=300, window_offset=150,
                                     min_mapq=0, no_dup=False):
    if chrom is None:
        # noinspection PyTypeChecker
        it = chain(*[
            iter_coverage_ext_binned(samfile, fastafile, chrom, None, None, one_based,
                                     window_size, window_offset, min_mapq, no_dup)
            for chrom in samfile.references
        ])
    else:
        it = iter_coverage_ext_binned(samfile, fastafile, chrom, start, end, one_based,
                                      window_size, window_offset, min_mapq, no_dup)
    return it


def iter_coverage_ext_binned(samfile, fastafile, chrom, start, end, one_based, window_size,
                             window_offset, min_mapq, no_dup):
    assert chrom is not None
    start, end = normalise_coords(one_based, start, end)
    chrlen = samfile.lengths[samfile.references.index(chrom)]
    if start is None:
        start = 0
    if end is None:
        end = chrlen
    if end > chrlen:
        end = chrlen
    # setup first bin
    bin_start = start
    bin_end = bin_start + window_size
    reads_all = reads_pp = reads_mate_unmapped = reads_mate_other_chr = \
        reads_mate_same_strand = reads_faceaway = reads_softclipped = \
        reads_duplicate = 0

    # iterate over reads
    alignments = samfile.fetch(chrom, start, end)
    alignments = filter_alignments(alignments, min_mapq, no_dup)
    for aln in alignments:
        while aln.pos > bin_end:  # end of bin
            gc = gc_content(fastafile, chrom, bin_start, bin_end)
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'gc': gc,
                   'reads_all': reads_all,
                   'reads_pp': reads_pp,
                   'reads_mate_unmapped': reads_mate_unmapped,
                   'reads_mate_other_chr': reads_mate_other_chr,
                   'reads_mate_same_strand': reads_mate_same_strand,
                   'reads_faceaway': reads_faceaway,
                   'reads_softclipped': reads_softclipped,
                   'reads_duplicate': reads_duplicate}
            yield rec
            reads_all = reads_pp = reads_mate_unmapped = reads_mate_other_chr\
                = reads_mate_same_strand = reads_faceaway = reads_softclipped\
                = reads_duplicate = 0
            bin_start = bin_end
            bin_end = bin_start + window_size
#        debug(reads, reads.cigar, repr(reads.cigarstring))
        if not aln.is_unmapped:
            reads_all += 1
            if aln.is_proper_pair:
                reads_pp += 1
            if aln.is_duplicate:
                reads_duplicate += 1
            if aln.cigar is not None and any((op[0] == 4) for op in aln.cigar):
                reads_softclipped += 1
            # should be mutually exclusive
            if aln.mate_is_unmapped:
                reads_mate_unmapped += 1
            elif aln.tid != aln.rnext:
                reads_mate_other_chr += 1
            elif aln.is_reverse == aln.mate_is_reverse:
                reads_mate_same_strand += 1
            elif (
                # mapped to reverse strand but leftmost
                (aln.is_reverse and aln.tlen > 0)
                # mapped to fwd strand but rightmost
                or (not aln.is_reverse and aln.tlen < 0)
            ):
                reads_faceaway += 1

    # deal with last non-empty bin
    gc = gc_content(fastafile, chrom, bin_start, bin_end)
    pos = bin_start + window_offset
    if one_based:
        pos += 1
    rec = {'chrom': chrom, 'pos': pos,
           'gc': gc,
           'reads_all': reads_all,
           'reads_pp': reads_pp,
           'reads_mate_unmapped': reads_mate_unmapped,
           'reads_mate_other_chr': reads_mate_other_chr,
           'reads_mate_same_strand': reads_mate_same_strand,
           'reads_faceaway': reads_faceaway,
           'reads_softclipped': reads_softclipped,
           'reads_duplicate': reads_duplicate}
    yield rec

    # deal with empty bins up to explicit end
    if end is not None:
        while bin_end < end:
            reads_all = reads_pp = reads_mate_unmapped = reads_mate_other_chr\
                = reads_mate_same_strand = reads_faceaway = reads_softclipped\
                = reads_duplicate = 0
            bin_start = bin_end
            bin_end = bin_start + window_size
            gc = gc_content(fastafile, chrom, bin_start, bin_end)
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'gc': gc,
                   'reads_all': reads_all,
                   'reads_pp': reads_pp,
                   'reads_mate_unmapped': reads_mate_unmapped,
                   'reads_mate_other_chr': reads_mate_other_chr,
                   'reads_mate_same_strand': reads_mate_same_strand,
                   'reads_faceaway': reads_faceaway,
                   'reads_softclipped': reads_softclipped,
                   'reads_duplicate': reads_duplicate}
            yield rec


def test_stat_coverage_ext_binned():
    compare_stats_withref(pysamstats.stat_coverage_ext_binned,
                          stat_coverage_ext_binned_refimpl)


def test_stat_coverage_ext_binned_uppercase_fasta():
    compare_stats_withref(pysamstats.stat_coverage_ext_binned,
                          stat_coverage_ext_binned_refimpl,
                          fasta_fn='fixture/ref.upper.fa')


def stat_mapq_binned_refimpl(samfile, chrom=None, start=None, end=None, one_based=False,
                             window_size=300, window_offset=150, min_mapq=0, no_dup=False):
    if chrom is None:
        # noinspection PyTypeChecker
        it = chain(*[iter_mapq_binned(samfile, chrom, None, None, one_based, window_size,
                                      window_offset, min_mapq, no_dup)
                     for chrom in samfile.references])
    else:
        it = iter_mapq_binned(samfile, chrom, start, end, one_based, window_size, window_offset,
                              min_mapq, no_dup)
    return it


def iter_mapq_binned(samfile, chrom, start, end, one_based, window_size, window_offset, min_mapq,
                     no_dup):
    assert chrom is not None
    start, end = normalise_coords(one_based, start, end)
    chrlen = samfile.lengths[samfile.references.index(chrom)]
    if start is None:
        start = 0
    if end is None:
        end = chrlen
    if end > chrlen:
        end = chrlen
    # setup first bin
    bin_start = start
    bin_end = bin_start + window_size
    reads_all = reads_mapq0 = mapq_square_sum = 0

    # iterate over reads
    alignments = samfile.fetch(chrom, start, end)
    alignments = filter_alignments(alignments, min_mapq, no_dup)
    for aln in alignments:
        while aln.pos > bin_end:  # end of bin
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'reads_all': reads_all,
                   'reads_mapq0': reads_mapq0,
                   'rms_mapq': rootmean(mapq_square_sum, reads_all)}
            yield rec
            reads_all = reads_mapq0 = mapq_square_sum = 0
            bin_start = bin_end
            bin_end = bin_start + window_size
        if not aln.is_unmapped:
            reads_all += 1
            mapq_square_sum += aln.mapq**2
            if aln.mapq == 0:
                reads_mapq0 += 1

    # deal with last non-empty bin
    pos = bin_start + window_offset
    if one_based:
        pos += 1
    rec = {'chrom': chrom, 'pos': pos,
           'reads_all': reads_all,
           'reads_mapq0': reads_mapq0,
           'rms_mapq': rootmean(mapq_square_sum, reads_all)}
    yield rec

    # deal with empty bins up to explicit end
    if end is not None:
        while bin_end < end:
            reads_all = reads_mapq0 = mapq_square_sum = 0
            bin_start = bin_end
            bin_end = bin_start + window_size
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'reads_all': reads_all,
                   'reads_mapq0': reads_mapq0,
                   'rms_mapq': rootmean(mapq_square_sum, reads_all)}
            yield rec


def test_stat_mapq_binned():
    compare_stats(pysamstats.stat_mapq_binned, stat_mapq_binned_refimpl)


def stat_alignment_binned_refimpl(samfile, chrom=None, start=None, end=None, one_based=False,
                                  window_size=300, window_offset=150, min_mapq=0, no_dup=False):
    if chrom is None:
        # noinspection PyTypeChecker
        it = chain(*[
            iter_alignment_binned(samfile, chrom, None, None, one_based, window_size,
                                  window_offset, min_mapq, no_dup)
            for chrom in samfile.references]
        )
    else:
        it = iter_alignment_binned(samfile, chrom, start, end, one_based, window_size,
                                   window_offset, min_mapq, no_dup)
    return it


CIGAR = 'MIDNSHP=X'


def iter_alignment_binned(samfile, chrom, start, end, one_based, window_size, window_offset,
                          min_mapq, no_dup):
    assert chrom is not None
    start, end = normalise_coords(one_based, start, end)
    chrlen = samfile.lengths[samfile.references.index(chrom)]
    if start is None:
        start = 0
    if end is None:
        end = chrlen
    if end > chrlen:
        end = chrlen
    # setup first bin
    bin_start = start
    bin_end = bin_start + window_size
    c = Counter()
    reads_all = 0

    # iterate over reads
    alignments = samfile.fetch(chrom, start, end)
    alignments = filter_alignments(alignments, min_mapq, no_dup)
    for aln in alignments:
        while aln.pos > bin_end:  # end of bin
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos, 'reads_all': reads_all}
            for i in range(len(CIGAR)):
                rec[CIGAR[i]] = c[i]
#            rec['NM'] = c['NM']
            rec['bases_all'] = c[0] + c[1] + c[4] + c[7] + c[8]
            yield rec
            c = Counter()
            reads_all = 0
            bin_start = bin_end
            bin_end = bin_start + window_size
#        debug(aln.cigar)
        if not aln.is_unmapped:
            reads_all += 1
            if aln.cigar is not None:
                for op, l in aln.cigar:
                    c[op] += l
            # add edit distance
    #        tags = dict(aln.tags)
    #        if 'NM' in tags:
    #            c['NM'] += tags['NM']

    # deal with last non-empty bin
    pos = bin_start + window_offset
    if one_based:
        pos += 1
    rec = {'chrom': chrom, 'pos': pos, 'reads_all': reads_all}
    for i in range(len(CIGAR)):
        rec[CIGAR[i]] = c[i]
#            rec['NM'] = c['NM']
    rec['bases_all'] = c[0] + c[1] + c[4] + c[7] + c[8]
    yield rec

    # deal with empty bins up to explicit end
    if end is not None:
        while bin_end < end:
            c = Counter()
            reads_all = 0
            bin_start = bin_end
            bin_end = bin_start + window_size
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos, 'reads_all': reads_all}
            for i in range(len(CIGAR)):
                rec[CIGAR[i]] = c[i]
        #            rec['NM'] = c['NM']
            rec['bases_all'] = c[0] + c[1] + c[4] + c[7] + c[8]
            yield rec


def test_stat_alignment_binned():
    compare_stats(pysamstats.stat_alignment_binned, stat_alignment_binned_refimpl)


def stat_tlen_binned_refimpl(samfile, chrom=None, start=None, end=None, one_based=False,
                             window_size=300, window_offset=150, min_mapq=0, no_dup=False):
    if chrom is None:
        # noinspection PyTypeChecker
        it = chain(*[iter_tlen_binned(samfile, chrom, None, None, one_based, window_size,
                                      window_offset, min_mapq, no_dup)
                     for chrom in samfile.references])
    else:
        it = iter_tlen_binned(samfile, chrom, start, end, one_based, window_size, window_offset,
                              min_mapq, no_dup)
    return it


def iter_tlen_binned(samfile, chrom, start, end, one_based, window_size, window_offset, min_mapq,
                     no_dup):
    assert chrom is not None
    start, end = normalise_coords(one_based, start, end)
    chrlen = samfile.lengths[samfile.references.index(chrom)]
    if start is None:
        start = 0
    if end is None:
        end = chrlen
    if end > chrlen:
        end = chrlen
    # setup first bin
    bin_start = start
    bin_end = bin_start + window_size
    reads_all = reads_pp = 0
    tlens = []
    tlens_pp = []

    # iterate over reads
    alignments = samfile.fetch(chrom, start, end)
    alignments = filter_alignments(alignments, min_mapq, no_dup)
    for aln in alignments:
        while aln.pos > bin_end:  # end of bin
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'reads_all': reads_all,
                   'reads_pp': reads_pp,
                   'mean_tlen': mean(tlens),
                   'rms_tlen': rms(tlens),
                   'mean_tlen_pp': mean(tlens_pp),
                   'rms_tlen_pp': rms(tlens_pp),
                   }
            yield rec
            reads_all = reads_pp = 0
            tlens = []
            tlens_pp = []
            bin_start = bin_end
            bin_end = bin_start + window_size
        if not aln.is_unmapped:
            reads_all += 1
            tlens.append(aln.tlen)
            if aln.is_proper_pair:
                reads_pp += 1
                tlens_pp.append(aln.tlen)

    # deal with last non-empty bin
    pos = bin_start + window_offset
    if one_based:
        pos += 1
    rec = {'chrom': chrom, 'pos': pos,
           'reads_all': reads_all,
           'reads_pp': reads_pp,
           'mean_tlen': mean(tlens),
           'rms_tlen': rms(tlens),
           'mean_tlen_pp': mean(tlens_pp),
           'rms_tlen_pp': rms(tlens_pp),
           }
    yield rec

    # deal with empty bins up to explicit end
    if end is not None:
        while bin_end < end:
            reads_all = reads_pp = 0
            tlens = []
            tlens_pp = []
            bin_start = bin_end
            bin_end = bin_start + window_size
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'reads_all': reads_all,
                   'reads_pp': reads_pp,
                   'mean_tlen': mean(tlens),
                   'rms_tlen': rms(tlens),
                   'mean_tlen_pp': mean(tlens_pp),
                   'rms_tlen_pp': rms(tlens_pp),
                   }
            yield rec


def test_stat_tlen_binned():
    compare_stats(pysamstats.stat_tlen_binned, stat_tlen_binned_refimpl)


binned_functions = [
    (pysamstats.load_coverage_binned, 1),
    (pysamstats.load_coverage_ext_binned, 1),
    (pysamstats.load_mapq_binned, 0),
    (pysamstats.load_alignment_binned, 0),
    (pysamstats.load_tlen_binned, 0),
]


def test_binned_pad_region():
    kwargs = {'chrom': 'Pf3D7_01_v3',
              'start': 1000,
              'end': 20000,
              'one_based': False,
              'window_size': 200,
              'window_offset': 100}
    for f, needs_ref in binned_functions:
        debug(f.__name__)
        if needs_ref:
            a = f(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'),
                  **kwargs)
        else:
            a = f(Samfile('fixture/test.bam'), **kwargs)
        assert set(a['chrom']) == {b'Pf3D7_01_v3'}
        eq_(1100, a['pos'][0])
        eq_(19900, a['pos'][-1])


def test_binned_pad_wg():
    expected = stat_coverage_binned_refimpl(
        Samfile('fixture/test.bam'),
        Fastafile('fixture/ref.fa'))

    actual = pysamstats.stat_coverage_binned(Samfile('fixture/test.bam'),
                                             Fastafile('fixture/ref.fa'))
    compare_iterators(expected, actual)
    kwargs = {'window_size': 200,
              'window_offset': 100}
    for f, needs_ref in binned_functions:
        debug(f.__name__)
        if needs_ref:
            a = f(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'),
                  **kwargs)
        else:
            a = f(Samfile('fixture/test.bam'), **kwargs)
        assert sorted(set(a['chrom'])) == [b'Pf3D7_01_v3', b'Pf3D7_02_v3',
                                           b'Pf3D7_03_v3']
        eq_(100, a[a['chrom'] == b'Pf3D7_01_v3']['pos'][0])
        eq_(50100, a[a['chrom'] == b'Pf3D7_01_v3']['pos'][-1])
        eq_(100, a[a['chrom'] == b'Pf3D7_02_v3']['pos'][0])
        eq_(60100, a[a['chrom'] == b'Pf3D7_02_v3']['pos'][-1])
        eq_(100, a[a['chrom'] == b'Pf3D7_03_v3']['pos'][0])
        eq_(70100, a[a['chrom'] == b'Pf3D7_03_v3']['pos'][-1])
