"""
Tests for the pysamstats module.

The strategy here is to compare the outputs of the functions under test with
unoptimised, pure-python reference implementations of the same functions, over
an example dataset. 

"""


from pprint import pprint 
from pysam import Samfile, Fastafile
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
                pprint(k)
                pprint(e)
                pprint(a)
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
                             and read.alignment.is_reverse != read.alignment.mate_is_reverse
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


def stat_coverage_ext_strand_refimpl(samfile, chrom=None, start=None, end=None, one_based=False):
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
        reads_pp_rev = sum(1 for read in reads 
                           if read.alignment.is_reverse and read.alignment.is_proper_pair)
        reads_pp_fwd = sum(1 for read in reads 
                           if not read.alignment.is_reverse and read.alignment.is_proper_pair)
        reads_mate_unmapped = sum(1 for read in reads if read.alignment.mate_is_unmapped)
        reads_mate_unmapped_rev = sum(1 for read in reads 
                                      if read.alignment.is_reverse and read.alignment.mate_is_unmapped)
        reads_mate_unmapped_fwd = sum(1 for read in reads 
                                      if not read.alignment.is_reverse and read.alignment.mate_is_unmapped)
        reads_mate_other_chr = sum(1 for read in reads
                                   if not read.alignment.mate_is_unmapped 
                                   and col.tid != read.alignment.rnext)
        reads_mate_other_chr_rev = sum(1 for read in reads
                                       if read.alignment.is_reverse
                                       and not read.alignment.mate_is_unmapped 
                                       and col.tid != read.alignment.rnext)
        reads_mate_other_chr_fwd = sum(1 for read in reads
                                       if not read.alignment.is_reverse
                                       and not read.alignment.mate_is_unmapped 
                                       and col.tid != read.alignment.rnext)
        reads_mate_same_strand = sum(1 for read in reads
                                     if not read.alignment.mate_is_unmapped
                                     and col.tid == read.alignment.rnext
                                     and (read.alignment.is_reverse == read.alignment.mate_is_reverse))
        reads_mate_same_strand_rev = sum(1 for read in reads
                                         if read.alignment.is_reverse
                                         and not read.alignment.mate_is_unmapped
                                         and col.tid == read.alignment.rnext
                                         and (read.alignment.is_reverse == read.alignment.mate_is_reverse))
        reads_mate_same_strand_fwd = sum(1 for read in reads
                                         if not read.alignment.is_reverse
                                         and not read.alignment.mate_is_unmapped
                                         and col.tid == read.alignment.rnext
                                         and (read.alignment.is_reverse == read.alignment.mate_is_reverse))
        reads_faceaway = sum(1 for read in reads
                             if not read.alignment.mate_is_unmapped
                             and read.alignment.is_reverse != read.alignment.mate_is_reverse
                             and ((read.alignment.is_reverse and read.alignment.tlen > 0) # mapped to reverse strand but leftmost
                                  or (not read.alignment.is_reverse and read.alignment.tlen < 0)) # mapped to fwd strand but rightmost
                             )
        reads_faceaway_rev = sum(1 for read in reads
                                 if read.alignment.is_reverse
                                 and not read.alignment.mate_is_unmapped
                                 and read.alignment.is_reverse != read.alignment.mate_is_reverse
                                 and ((read.alignment.is_reverse and read.alignment.tlen > 0) # mapped to reverse strand but leftmost
                                      or (not read.alignment.is_reverse and read.alignment.tlen < 0)) # mapped to fwd strand but rightmost
                                 )
        reads_faceaway_fwd = sum(1 for read in reads
                                 if not read.alignment.is_reverse
                                 and not read.alignment.mate_is_unmapped
                                 and read.alignment.is_reverse != read.alignment.mate_is_reverse
                                 and ((read.alignment.is_reverse and read.alignment.tlen > 0) # mapped to reverse strand but leftmost
                                      or (not read.alignment.is_reverse and read.alignment.tlen < 0)) # mapped to fwd strand but rightmost
                                 )
        reads_softclipped = sum(1 for read in reads
                                if any((op[0] == 4) for op in read.alignment.cigar))
        reads_softclipped_rev = sum(1 for read in reads
                                if read.alignment.is_reverse
                                and any((op[0] == 4) for op in read.alignment.cigar))
        reads_softclipped_fwd = sum(1 for read in reads
                                if not read.alignment.is_reverse
                                and any((op[0] == 4) for op in read.alignment.cigar))
        yield {'chr': chrom, 'pos': pos, 
               'reads_all': reads_all, 
               'reads_fwd': reads_fwd,
               'reads_rev': reads_rev,
               'reads_pp': reads_pp,
               'reads_pp_fwd': reads_pp_fwd,
               'reads_pp_rev': reads_pp_rev,
               'reads_mate_unmapped': reads_mate_unmapped,
               'reads_mate_unmapped_fwd': reads_mate_unmapped_fwd,
               'reads_mate_unmapped_rev': reads_mate_unmapped_rev,
               'reads_mate_other_chr': reads_mate_other_chr,
               'reads_mate_other_chr_fwd': reads_mate_other_chr_fwd,
               'reads_mate_other_chr_rev': reads_mate_other_chr_rev,
               'reads_mate_same_strand': reads_mate_same_strand,
               'reads_mate_same_strand_fwd': reads_mate_same_strand_fwd,
               'reads_mate_same_strand_rev': reads_mate_same_strand_rev,
               'reads_faceaway': reads_faceaway,
               'reads_faceaway_fwd': reads_faceaway_fwd,
               'reads_faceaway_rev': reads_faceaway_rev,
               'reads_softclipped': reads_softclipped,
               'reads_softclipped_fwd': reads_softclipped_fwd,
               'reads_softclipped_rev': reads_softclipped_rev}
        

def test_stat_coverage_ext_strand():
    _test(pysamstats.stat_coverage_ext_strand, stat_coverage_ext_strand_refimpl)


def stat_variation_refimpl(samfile, fafile, chrom=None, start=None, end=None, one_based=False):
    if one_based:
        start = start - 1 if start is not None else None
        end = end - 1 if end is not None else None
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        reads_all = col.n
        reads = col.pileups
        reads_pp = sum(1 for read in reads if read.alignment.is_proper_pair)
        ref = fafile.fetch(chrom, col.pos, col.pos+1).upper()
        matches = sum(1 for read in reads
                      if not read.is_del
#                      and read.alignment.seq[read.qpos] in 'ACTG'
                      and read.alignment.seq[read.qpos] == ref)
        matches_pp = sum(1 for read in reads
                         if not read.is_del
                         and read.alignment.is_proper_pair
 #                        and read.alignment.seq[read.qpos] in 'ACTG'
                         and read.alignment.seq[read.qpos] == ref)
        mismatches = sum(1 for read in reads
                         if not read.is_del
#                         and read.alignment.seq[read.qpos] in 'ACTG'
                         and read.alignment.seq[read.qpos] != ref)
        mismatches_pp = sum(1 for read in reads
                            if not read.is_del
                            and read.alignment.is_proper_pair
#                            and read.alignment.seq[read.qpos] in 'ACTG'
                            and read.alignment.seq[read.qpos] != ref)
        deletions = sum(1 for read in reads
                        if read.is_del)
        deletions_pp = sum(1 for read in reads
                           if read.is_del
                           and read.alignment.is_proper_pair)
        insertions = sum(1 for read in reads
                         if read.indel > 0)
        insertions_pp = sum(1 for read in reads
                            if read.indel > 0
                            and read.alignment.is_proper_pair)
        A = sum(1 for read in reads
                if not read.is_del
                and read.alignment.seq[read.qpos] == 'A')
        A_pp = sum(1 for read in reads
                   if not read.is_del
                   and read.alignment.is_proper_pair
                   and read.alignment.seq[read.qpos] == 'A')
        C = sum(1 for read in reads
                if not read.is_del
                and read.alignment.seq[read.qpos] == 'C')
        C_pp = sum(1 for read in reads
                   if not read.is_del
                   and read.alignment.is_proper_pair
                   and read.alignment.seq[read.qpos] == 'C')
        T = sum(1 for read in reads
                if not read.is_del
                and read.alignment.seq[read.qpos] == 'T')
        T_pp = sum(1 for read in reads
                   if not read.is_del
                   and read.alignment.is_proper_pair
                   and read.alignment.seq[read.qpos] == 'T')
        G = sum(1 for read in reads
                if not read.is_del
                and read.alignment.seq[read.qpos] == 'G')
        G_pp = sum(1 for read in reads
                   if not read.is_del
                   and read.alignment.is_proper_pair
                   and read.alignment.seq[read.qpos] == 'G')
        N = sum(1 for read in reads
                if not read.is_del
                and read.alignment.seq[read.qpos] == 'N')
        N_pp = sum(1 for read in reads
                   if not read.is_del
                   and read.alignment.is_proper_pair
                   and read.alignment.seq[read.qpos] == 'N')
        yield {'chr': chrom, 'pos': pos, 'ref': ref,
               'reads_all': reads_all, 
               'reads_pp': reads_pp,
               'matches': matches,
               'matches_pp': matches_pp,
               'mismatches': mismatches,
               'mismatches_pp': mismatches_pp,
               'deletions': deletions,
               'deletions_pp': deletions_pp,
               'insertions': insertions,
               'insertions_pp': insertions_pp,
               'A': A, 'A_pp': A_pp,
               'C': C, 'C_pp': C_pp,
               'T': T, 'T_pp': T_pp,
               'G': G, 'G_pp': G_pp,
               'N': N, 'N_pp': N_pp}
        

def test_stat_variation():
    _test_withrefseq(pysamstats.stat_variation, stat_variation_refimpl)

        
