from mock import Mock
from nose.tools import eq_


import pyximport; pyximport.install()
from pysamstats import construct_rec_coverage, construct_rec_coverage_strand


def test_construct_rec_coverage():
    
    # mock samfile
    samfile = Mock()
    samfile.getrname = Mock(return_value='chr1')

    # mock reads
    read1 = Mock()
    read1.alignment.is_proper_pair = True
    read2 = Mock()
    read2.alignment.is_proper_pair = False
    reads = [read1, read2]

    # mock pileup proxy (a.k.a. pileup column)
    col = Mock()
    col.tid = 0
    col.pos = 0
    col.n = len(reads)
    col.pileups = reads

    # call function under test
    chrom, pos, reads_all, reads_pp = construct_rec_coverage(samfile, col)

    # assertions
    eq_('chr1', chrom)
    eq_(0, pos)
    eq_(2, reads_all)
    eq_(1, reads_pp)


def test_construct_rec_coverage_strand():
    
    # mock samfile
    samfile = Mock()
    samfile.getrname = Mock(return_value='chr1')

    # mock reads
    read1 = Mock()
    read1.alignment.is_reverse = False
    read1.alignment.is_proper_pair = True
    read2 = Mock()
    read2.alignment.is_reverse = True
    read2.alignment.is_proper_pair = True
    read3 = Mock()
    read3.alignment.is_reverse = False
    read3.alignment.is_proper_pair = False
    read4 = Mock()
    read4.alignment.is_reverse = True
    read4.alignment.is_proper_pair = False
    reads = [read1, read2, read3, read4]

    # mock pileup proxy (a.k.a. pileup column)
    col = Mock()
    col.tid = 0
    col.pos = 0
    col.n = len(reads)
    col.pileups = reads

    # call function under test
    chrom, pos, reads_all, reads_fwd, reads_rev, reads_pp, reads_pp_fwd, reads_pp_rev = construct_rec_coverage_strand(samfile, col)

    # assertions
    eq_('chr1', chrom)
    eq_(0, pos)
    eq_(4, reads_all)
    eq_(2, reads_fwd)
    eq_(2, reads_rev)
    eq_(2, reads_pp)
    eq_(1, reads_pp_fwd)
    eq_(1, reads_pp_rev)


            
    
