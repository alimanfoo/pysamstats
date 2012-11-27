from mock import Mock
from nose.tools import eq_


import pyximport; pyximport.install()
from pysamstats import construct_rec_coverage


def test_construct_rec_coverage():
    
    # mock samfile
    samfile = Mock()
    samfile.getrname = Mock(return_value='chr1')

    # mock pileup proxy (a.k.a. pileup column)
    col = Mock()
    col.tid = 0
    col.pos = 0
    col.n = 4

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
    col.pileups = Mock(return_value=[read1, read2, read3, read4])

    # call function under test
    rec = construct_rec_coverage(samfile, col)

    # assertions
    eq_('chr1', rec['chr'])
    eq_(0, rec['pos'])
    eq_(4, rec['reads']['all'])
    eq_(2, rec['reads']['fwd'])
    eq_(2, rec['reads']['rev'])
    eq_(2, rec['reads']['pp'])
    eq_(1, rec['reads']['pp_fwd'])
    eq_(1, rec['reads']['pp_rev'])


            
    
