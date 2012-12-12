"""
Tests for the pysamstats module.

The strategy here is to compare the outputs of the functions under test with
unoptimised, pure-python reference implementations of the same functions, over
an example dataset. 

"""


# TODO simplify reference implementations by building lists and using len()

from pysam import Samfile, Fastafile
from nose.tools import eq_, assert_almost_equal
import numpy as np
from math import sqrt
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
                if isinstance(v, float):
                    assert_almost_equal(v, a[k])
                else:
                    eq_(v, a[k])
            except:
                print k
                print e
                print a
                raise


def _test_withrefseq(impl, refimpl):
    kwargs = {'chrom': 'Pf3D7_01_v3',
              'start': 0,
              'end': 10000,
              'one_based': False}
    expected = refimpl(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'), **kwargs)
    actual = impl(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'), **kwargs)
    for e, a in zip(expected, actual):
        for k, v in e.items():
            try:
                eq_(v, a[k])
            except:
                print k
                print e, a
                raise


def normalise_coords(one_based, start, end):
    if one_based:
        start = start - 1 if start is not None else None
        end = end - 1 if end is not None else None
    return start, end


def stat_coverage_refimpl(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(one_based, start, end)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads = col.pileups
        reads_pp = [read for read in reads if read.alignment.is_proper_pair]
        yield {'chr': chrom, 'pos': pos, 'reads_all': len(reads), 'reads_pp': len(reads_pp)}
        

def test_stat_coverage():
    _test(pysamstats.stat_coverage, stat_coverage_refimpl)


def fwd(reads):
    return [read for read in reads if not read.alignment.is_reverse]


def rev(reads):
    return [read for read in reads if read.alignment.is_reverse]

        
def stat_coverage_strand_refimpl(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(one_based, start, end)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads = col.pileups
        reads_fwd = fwd(reads)
        reads_rev = rev(reads)
        reads_pp = [read for read in reads if read.alignment.is_proper_pair]
        reads_pp_fwd = fwd(reads_pp)
        reads_pp_rev = rev(reads_pp)
        yield {'chr': chrom, 'pos': pos, 
               'reads_all': len(reads), 'reads_fwd': len(reads_fwd), 'reads_rev': len(reads_rev),
               'reads_pp': len(reads_pp), 'reads_pp_fwd': len(reads_pp_fwd), 'reads_pp_rev': len(reads_pp_rev)}
        

def test_stat_coverage_strand():
    _test(pysamstats.stat_coverage_strand, stat_coverage_strand_refimpl)


def stat_coverage_ext_refimpl(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(one_based, start, end)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads = col.pileups
        reads_pp = [read for read in reads if read.alignment.is_proper_pair]
        reads_mate_unmapped = [read for read in reads if read.alignment.mate_is_unmapped]
        reads_mate_mapped = [read for read in reads if not read.alignment.mate_is_unmapped]
        reads_mate_other_chr = [read for read in reads_mate_mapped
                                if col.tid != read.alignment.rnext]
        reads_mate_same_strand = [read for read in reads_mate_mapped
                                  if col.tid == read.alignment.rnext
                                  and (read.alignment.is_reverse == read.alignment.mate_is_reverse)]
        reads_faceaway = [read for read in reads_mate_mapped
                          if read.alignment.is_reverse != read.alignment.mate_is_reverse
                          and ((read.alignment.is_reverse and read.alignment.tlen > 0) # mapped to reverse strand but leftmost
                               or (not read.alignment.is_reverse and read.alignment.tlen < 0)) # mapped to fwd strand but rightmost
                          ]
        reads_softclipped = [read for read in reads
                             if any((op[0] == 4) for op in read.alignment.cigar)]
        yield {'chr': chrom, 'pos': pos, 
               'reads_all': len(reads), 
               'reads_pp': len(reads_pp),
               'reads_mate_unmapped': len(reads_mate_unmapped),
               'reads_mate_other_chr': len(reads_mate_other_chr),
               'reads_mate_same_strand': len(reads_mate_same_strand),
               'reads_faceaway': len(reads_faceaway),
               'reads_softclipped': len(reads_softclipped)}
        

def test_stat_coverage_ext():
    _test(pysamstats.stat_coverage_ext, stat_coverage_ext_refimpl)


def stat_coverage_ext_strand_refimpl(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(one_based, start, end)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads = col.pileups
        reads_pp = [read for read in reads if read.alignment.is_proper_pair]
        reads_mate_unmapped = [read for read in reads if read.alignment.mate_is_unmapped]
        reads_mate_mapped = [read for read in reads if not read.alignment.mate_is_unmapped]
        reads_mate_other_chr = [read for read in reads_mate_mapped
                                if col.tid != read.alignment.rnext]
        reads_mate_same_strand = [read for read in reads_mate_mapped
                                  if col.tid == read.alignment.rnext
                                  and (read.alignment.is_reverse == read.alignment.mate_is_reverse)]
        reads_faceaway = [read for read in reads_mate_mapped
                          if read.alignment.is_reverse != read.alignment.mate_is_reverse
                          and ((read.alignment.is_reverse and read.alignment.tlen > 0) # mapped to reverse strand but leftmost
                               or (not read.alignment.is_reverse and read.alignment.tlen < 0)) # mapped to fwd strand but rightmost
                          ]
        reads_softclipped = [read for read in reads
                             if any((op[0] == 4) for op in read.alignment.cigar)]
        reads_fwd =fwd(reads)
        reads_rev = rev(reads)
        reads_pp_fwd = fwd(reads_pp)
        reads_pp_rev = rev(reads_pp)
        reads_mate_unmapped_fwd = fwd(reads_mate_unmapped)
        reads_mate_unmapped_rev = rev(reads_mate_unmapped)
        reads_mate_other_chr_fwd = fwd(reads_mate_other_chr)
        reads_mate_other_chr_rev = rev(reads_mate_other_chr)
        reads_mate_same_strand_fwd = fwd(reads_mate_same_strand)
        reads_mate_same_strand_rev = rev(reads_mate_same_strand)
        reads_faceaway_fwd = fwd(reads_faceaway)
        reads_faceaway_rev = rev(reads_faceaway)
        reads_softclipped_fwd = fwd(reads_softclipped)
        reads_softclipped_rev = rev(reads_softclipped)
        yield {'chr': chrom, 'pos': pos, 
               'reads_all': len(reads), 
               'reads_fwd': len(reads_fwd),
               'reads_rev': len(reads_rev),
               'reads_pp': len(reads_pp),
               'reads_pp_fwd': len(reads_pp_fwd),
               'reads_pp_rev': len(reads_pp_rev),
               'reads_mate_unmapped': len(reads_mate_unmapped),
               'reads_mate_unmapped_fwd': len(reads_mate_unmapped_fwd),
               'reads_mate_unmapped_rev': len(reads_mate_unmapped_rev),
               'reads_mate_other_chr': len(reads_mate_other_chr),
               'reads_mate_other_chr_fwd': len(reads_mate_other_chr_fwd),
               'reads_mate_other_chr_rev': len(reads_mate_other_chr_rev),
               'reads_mate_same_strand': len(reads_mate_same_strand),
               'reads_mate_same_strand_fwd': len(reads_mate_same_strand_fwd),
               'reads_mate_same_strand_rev': len(reads_mate_same_strand_rev),
               'reads_faceaway': len(reads_faceaway),
               'reads_faceaway_fwd': len(reads_faceaway_fwd),
               'reads_faceaway_rev': len(reads_faceaway_rev),
               'reads_softclipped': len(reads_softclipped),
               'reads_softclipped_fwd': len(reads_softclipped_fwd),
               'reads_softclipped_rev': len(reads_softclipped_rev)}
        

def test_stat_coverage_ext_strand():
    _test(pysamstats.stat_coverage_ext_strand, stat_coverage_ext_strand_refimpl)


def stat_variation_refimpl(samfile, fafile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(one_based, start, end)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads = col.pileups
        reads_nodel = [read for read in reads if not read.is_del]
        reads_pp = [read for read in reads if read.alignment.is_proper_pair]
        reads_pp_nodel = [read for read in reads if read.alignment.is_proper_pair and not read.is_del]
        ref = fafile.fetch(chrom, col.pos, col.pos+1).upper()
        matches = [read for read in reads_nodel
                      if read.alignment.seq[read.qpos] == ref]
        matches_pp = [read for read in reads_pp_nodel
                         if read.alignment.seq[read.qpos] == ref]
        mismatches = [read for read in reads_nodel
                         if read.alignment.seq[read.qpos] != ref]
        mismatches_pp = [read for read in reads_pp_nodel
                            if read.alignment.seq[read.qpos] != ref]
        deletions = [read for read in reads
                        if read.is_del]
        deletions_pp = [read for read in reads_pp
                           if read.is_del]
        insertions = [read for read in reads
                         if read.indel > 0]
        insertions_pp = [read for read in reads_pp
                            if read.indel > 0]
        A = [read for read in reads_nodel
                if read.alignment.seq[read.qpos] == 'A']
        A_pp = [read for read in reads_pp_nodel
                   if read.alignment.seq[read.qpos] == 'A']
        C = [read for read in reads_nodel
                if read.alignment.seq[read.qpos] == 'C']
        C_pp = [read for read in reads_pp_nodel
                   if read.alignment.seq[read.qpos] == 'C']
        T = [read for read in reads_nodel
                if read.alignment.seq[read.qpos] == 'T']
        T_pp = [read for read in reads_pp_nodel
                   if read.alignment.seq[read.qpos] == 'T']
        G = [read for read in reads_nodel
                if read.alignment.seq[read.qpos] == 'G']
        G_pp = [read for read in reads_pp_nodel
                   if read.alignment.seq[read.qpos] == 'G']
        N = [read for read in reads_nodel
                if read.alignment.seq[read.qpos] == 'N']
        N_pp = [read for read in reads_pp_nodel
                   if read.alignment.seq[read.qpos] == 'N']
        yield {'chr': chrom, 'pos': pos, 'ref': ref,
               'reads_all': len(reads), 
               'reads_pp': len(reads_pp),
               'matches': len(matches),
               'matches_pp': len(matches_pp),
               'mismatches': len(mismatches),
               'mismatches_pp': len(mismatches_pp),
               'deletions': len(deletions),
               'deletions_pp': len(deletions_pp),
               'insertions': len(insertions),
               'insertions_pp': len(insertions_pp),
               'A': len(A), 'A_pp': len(A_pp),
               'C': len(C), 'C_pp': len(C_pp),
               'T': len(T), 'T_pp': len(T_pp),
               'G': len(G), 'G_pp': len(G_pp),
               'N': len(N), 'N_pp': len(N_pp)}
        

def test_stat_variation():
    _test_withrefseq(pysamstats.stat_variation, stat_variation_refimpl)

        
def stat_tlen_refimpl(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(one_based, start, end)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads = [read for read in col.pileups if not read.alignment.mate_is_unmapped and read.alignment.rnext == col.tid]
        tlen = [read.alignment.tlen for read in reads]
        rms_tlen = sqrt(np.mean(np.power(tlen, 2)))
        std_tlen = np.std(tlen)
        reads_pp = [read for read in reads if read.alignment.is_proper_pair]
        tlen_pp = [read.alignment.tlen for read in reads_pp]
        rms_tlen_pp = sqrt(np.mean(np.power(tlen_pp, 2))) 
        std_tlen_pp = np.std(tlen_pp)
        yield {'chr': chrom, 'pos': pos, 
               'reads_all': col.n, 
               'reads_pp': len(reads_pp),
               'rms_tlen': int(round(rms_tlen)),
               'rms_tlen_pp': int(round(rms_tlen_pp)),
               'std_tlen': int(round(std_tlen)),
               'std_tlen_pp': int(round(std_tlen_pp))}
        

def test_stat_tlen():
    _test(pysamstats.stat_tlen, stat_tlen_refimpl)

        
