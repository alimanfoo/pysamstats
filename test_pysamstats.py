"""
Tests for the pysamstats module.

The strategy here is to compare the outputs of the functions under test with
unoptimised, pure-python reference implementations of the same functions, over
an example dataset. 

"""


from pysam import Samfile
from nose.tools import eq_
import pysamstats


def _test(impl, refimpl):
    kwargs = {'chrom': 'Pf3D7_01_v3',
              'start': 0,
              'end': 10000,
              'one_based': False}
    expected = refimpl(Samfile('fixture/test.bam'), **kwargs)
    actual = impl(Samfile('fixture/test.bam'), **kwargs)
    for e, a in zip(expected, actual):
        for k, v in e.items():
            try:
                eq_(v, a[k])
            except:
                print k
                print e, a
                raise


def stat_coverage_refimpl(samfile, chrom=None, start=None, end=None, one_based=False):
    if one_based:
        start = start - 1 if start is not None else None
        end = end - 1 if end is not None else None
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads_all = col.n
        reads_pp = sum(1 for read in col.pileups if read.alignment.is_proper_pair)
        yield {'chr': chrom, 'pos': pos, 'reads_all': reads_all, 'reads_pp': reads_pp}
        

def test_stat_coverage():
    _test(pysamstats.stat_coverage, stat_coverage_refimpl)

        
def stat_coverage_strand_refimpl(samfile, chrom=None, start=None, end=None, one_based=False):
    if one_based:
        start = start - 1 if start is not None else None
        end = end - 1 if end is not None else None
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads_all = col.n
        reads = col.pileups
        reads_fwd = sum(1 for read in reads if not read.alignment.is_reverse)
        reads_rev = sum(1 for read in reads if read.alignment.is_reverse)
        reads_pp = sum(1 for read in reads if read.alignment.is_proper_pair)
        reads_pp_fwd = sum(1 for read in reads if read.alignment.is_proper_pair and not read.alignment.is_reverse)
        reads_pp_rev = sum(1 for read in reads if read.alignment.is_proper_pair and read.alignment.is_reverse)
        yield {'chr': chrom, 'pos': pos, 
               'reads_all': reads_all, 'reads_fwd': reads_fwd, 'reads_rev': reads_rev,
               'reads_pp': reads_pp, 'reads_pp_fwd': reads_pp_fwd, 'reads_pp_rev': reads_pp_rev}
        

def test_stat_coverage_strand():
    _test(pysamstats.stat_coverage_strand, stat_coverage_strand_refimpl)


def stat_coverage_ext_refimpl(samfile, chrom=None, start=None, end=None, one_based=False):
    if one_based:
        start = start - 1 if start is not None else None
        end = end - 1 if end is not None else None
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads_all = col.n
        reads = col.pileups
        reads_pp = sum(1 for read in reads if read.alignment.is_proper_pair)
        reads_mate_unmapped = sum(1 for read in reads if read.alignment.mate_is_unmapped)
        reads_mate_other_chr = sum(1 for read in reads
                                   if not read.alignment.mate_is_unmapped 
                                   and col.tid != read.alignment.rnext)
        reads_mate_same_strand = sum(1 for read in reads
                                     if not read.alignment.mate_is_unmapped
                                     and col.tid == read.alignment.rnext
                                     and (read.alignment.is_reverse == read.alignment.mate_is_reverse))
        reads_faceaway = sum(1 for read in reads
                             if not read.alignment.mate_is_unmapped
                             and ((read.alignment.is_reverse and read.alignment.tlen > 0) # mapped to reverse strand but leftmost
                                  or (not read.alignment.is_reverse and read.alignment.tlen < 0)) # mapped to fwd strand but rightmost
                             )
        reads_softclipped = sum(1 for read in reads
                                if any((op[0] == 4) for op in read.alignment.cigar))
        yield {'chr': chrom, 'pos': pos, 
               'reads_all': reads_all, 
               'reads_pp': reads_pp,
               'reads_mate_unmapped': reads_mate_unmapped,
               'reads_mate_other_chr': reads_mate_other_chr,
               'reads_mate_same_strand': reads_mate_same_strand,
               'reads_faceaway': reads_faceaway,
               'reads_softclipped': reads_softclipped}
        

def test_stat_coverage_ext():
    _test(pysamstats.stat_coverage_ext, stat_coverage_ext_refimpl)
