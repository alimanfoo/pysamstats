from mock import Mock
from nose.tools import eq_


import pyximport; pyximport.install()
from pysamstats import construct_rec_coverage, construct_rec_coverage_strand, \
    construct_rec_coverage_ext


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
    rec = construct_rec_coverage(samfile, col)

    # assertions
    eq_('chr1', rec['chr'])
    eq_(0, rec['pos'])
    eq_(2, rec['reads_all'])
    eq_(1, rec['reads_pp'])


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
    rec = construct_rec_coverage_strand(samfile, col)

    # assertions
    eq_('chr1', rec['chr'])
    eq_(0, rec['pos'])
    eq_(4, rec['reads_all'])
    eq_(2, rec['reads_fwd'])
    eq_(2, rec['reads_rev'])
    eq_(2, rec['reads_pp'])
    eq_(1, rec['reads_pp_fwd'])
    eq_(1, rec['reads_pp_rev'])


def test_construct_rec_coverage_ext():
    
    # mock samfile
    samfile = Mock()
    samfile.getrname = Mock(return_value='chr1')

    # mock reads

    # good read
    read1 = Mock() 
    read1.alignment.is_proper_pair = True
    read1.alignment.is_reverse = False
    read1.alignment.mate_is_unmapped = False
    read1.alignment.tid = 0
    read1.alignment.rnext = 0 # mate same chromosome
    read1.alignment.mate_is_reverse = True
    read1.alignment.tlen = 170 # leftmost
    read1.alignment.opt = Mock(side_effect=lambda tag: 2 if tag == 'NM' else None)
    read1.alignment.cigar = ((0, 10),) # no soft clipping

    # not properly paired
    read2 = Mock()
    read2.alignment.is_proper_pair = False
    read2.alignment.is_reverse = False
    read2.alignment.mate_is_unmapped = False
    read2.alignment.tid = 0
    read2.alignment.rnext = 0 # mate same chromosome
    read2.alignment.mate_is_reverse = True
    read2.alignment.tlen = 17000 # leftmost
    read2.alignment.opt = Mock(side_effect=lambda tag: 2 if tag == 'NM' else None)
    read2.alignment.cigar = ((0, 10),) # no soft clipping
    
    # mate is unmapped
    read3 = Mock()
    read3.alignment.is_proper_pair = False
    read3.alignment.is_reverse = False
    read3.alignment.mate_is_unmapped = True
    read3.alignment.opt = Mock(side_effect=lambda tag: 2 if tag == 'NM' else None)
    read3.alignment.cigar = ((0, 10),) # no soft clipping

    # mate other chromosome
    read4 = Mock()
    read4.alignment.is_proper_pair = False
    read4.alignment.is_reverse = False
    read4.alignment.mate_is_unmapped = False
    read4.alignment.tid = 0
    read4.alignment.rnext = 1 # mate other chromosome
    read4.alignment.opt = Mock(side_effect=lambda tag: 2 if tag == 'NM' else None)
    read4.alignment.cigar = ((0, 10),) # no soft clipping

    # mate same strand    
    read5 = Mock() 
    read5.alignment.is_proper_pair = False
    read5.alignment.is_reverse = False
    read5.alignment.mate_is_unmapped = False
    read5.alignment.tid = 0
    read5.alignment.rnext = 0 # mate same chromosome
    read5.alignment.mate_is_reverse = False
    read5.alignment.tlen = 170 # leftmost
    read5.alignment.opt = Mock(side_effect=lambda tag: 2 if tag == 'NM' else None)
    read5.alignment.cigar = ((0, 10),) # no soft clipping
    
    # faceaway
    read6 = Mock() 
    read6.alignment.is_proper_pair = False
    read6.alignment.is_reverse = False
    read6.alignment.mate_is_unmapped = False
    read6.alignment.tid = 0
    read6.alignment.rnext = 0 # mate same chromosome
    read6.alignment.mate_is_reverse = True
    read6.alignment.tlen = -170 # rightmost
    read6.alignment.opt = Mock(side_effect=lambda tag: 2 if tag == 'NM' else None)
    read6.alignment.cigar = ((0, 10),) # no soft clipping

    # edit distance 0
    read7 = Mock() 
    read7.alignment.is_proper_pair = True
    read7.alignment.is_reverse = False
    read7.alignment.mate_is_unmapped = False
    read7.alignment.tid = 0
    read7.alignment.rnext = 0 # mate same chromosome
    read7.alignment.mate_is_reverse = True
    read7.alignment.tlen = 170 # leftmost
    read7.alignment.opt = Mock(side_effect=lambda tag: 0 if tag == 'NM' else None)
    read7.alignment.cigar = ((0, 10),) # no soft clipping

    # soft clipped
    read8 = Mock() 
    read8.alignment.is_proper_pair = True
    read8.alignment.is_reverse = False
    read8.alignment.mate_is_unmapped = False
    read8.alignment.tid = 0
    read8.alignment.rnext = 0 # mate same chromosome
    read8.alignment.mate_is_reverse = True
    read8.alignment.tlen = 170 # leftmost
    read8.alignment.opt = Mock(side_effect=lambda tag: 2 if tag == 'NM' else None)
    read8.alignment.cigar = ((0, 5), (4, 2)) # soft clipping

    reads = [read1, read2, read3, read4, read5, read6, read7, read8]

    # mock pileup proxy (a.k.a. pileup column)
    col = Mock()
    col.tid = 0
    col.pos = 0
    col.n = len(reads)
    col.pileups = reads

    # call function under test
    rec = construct_rec_coverage_ext(samfile, col)

    # assertions
    eq_('chr1', rec['chr'])
    eq_(0, rec['pos'])
    eq_(8, rec['reads_all'])
    eq_(3, rec['reads_pp'])
    eq_(1, rec['reads_mate_unmapped'])
    eq_(1, rec['reads_mate_other_chr'])
    



            
    
