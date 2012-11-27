# cython: profile=True


import numpy as np
cimport numpy as np


cdef struct CountStrand:
    unsigned short all
    unsigned short fwd
    unsigned short rev


cdef struct CountPpStrand:
    unsigned short all
    unsigned short fwd
    unsigned short rev
    unsigned short pp
    unsigned short pp_fwd
    unsigned short pp_rev


cdef struct RecCoverage:
    char* chr
    unsigned int pos
    CountPpStrand reads


cdef inline void init_pp_strand(CountPpStrand *c):
    c.all = 0
    c.rev = 0
    c.fwd = 0
    c.pp = 0
    c.pp_fwd = 0
    c.pp_rev = 0


cdef inline void incr_pp_strand(CountPpStrand *c, bint is_reverse, bint is_proper_pair):
    c.all += 1
    if is_reverse:
        c.rev += 1
        if is_proper_pair:
            c.pp += 1
            c.pp_rev += 1
    else:
        c.fwd += 1
        if is_proper_pair:
            c.pp += 1
            c.pp_fwd += 1


cpdef RecCoverage construct_rec_coverage(samfile, col):

    # statically typed variables
    cdef RecCoverage rec # return value
    cdef Py_ssize_t i # loop index
    cdef Py_ssize_t n # total number of reads in column
    cdef bint is_reverse # whether read is mapped to reverse strand
    cdef bint is_proper_pair # whether the read is mapped in a proper pair
    init_pp_strand(&rec.reads) # initialise counts to zero

    # set chromosome name and position
    chr = samfile.getrname(col.tid)
    rec.chr = chr
    pos = col.pos
    rec.pos = pos
    
    # loop over reads
    n = col.n
    reads = col.pileups
    for i in range(n):
        read = reads[i]
        aln = read.alignment
        is_reverse = aln.is_reverse
        is_proper_pair = aln.is_proper_pair
        incr_pp_strand(&rec.reads, is_reverse, is_proper_pair)

    return rec


def stat_coverage(samfile, chr=None, start=None, end=None):
    for col in samfile.pileup(chr, start, end):
        yield construct_rec_coverage(samfile, col)


cpdef RecCoverage construct_rec_coverage2(samfile, col) except *:

    # statically typed variables
    cdef RecCoverage rec # return value
    cdef Py_ssize_t i # loop index
    cdef Py_ssize_t n # total number of reads in column
    cdef bint b_is_reverse, b_is_proper_pair
    cdef np.ndarray[np.uint8_t, ndim=1] is_reverse # whether read is mapped to reverse strand
    cdef np.ndarray[np.uint8_t, ndim=1] is_proper_pair # whether the read is mapped in a proper pair

    # initialise
    n = col.n
    is_reverse = np.zeros((n,), dtype=np.uint8)
    is_proper_pair = np.zeros((n,), dtype=np.uint8)
    init_pp_strand(&rec.reads) # initialise counts to zero

    # set chromosome name and position
    chr = samfile.getrname(col.tid)
    rec.chr = chr
    pos = col.pos
    rec.pos = pos
    
    # loop over reads
    reads = col.pileups
    for i in range(n):
        read = reads[i]
        aln = read.alignment
        b_is_reverse = aln.is_reverse
        is_reverse[i] = b_is_reverse
        b_is_proper_pair = aln.is_proper_pair
        is_proper_pair[i] = b_is_proper_pair

    # determine counts
    is_reverse.dtype = np.bool
    is_proper_pair.dtype = np.bool
    is_forward = ~is_reverse
    is_pp_fwd = is_forward & is_proper_pair
    is_pp_rev = is_reverse & is_proper_pair
    rec.reads.all = n
    rec.reads.fwd = np.count_nonzero(is_forward)
    rec.reads.rev = np.count_nonzero(is_reverse)
    rec.reads.pp = np.count_nonzero(is_proper_pair)
    rec.reads.pp_fwd = np.count_nonzero(is_pp_fwd)
    rec.reads.pp_rev = np.count_nonzero(is_pp_rev)

    return rec


def stat_coverage2(samfile, chr=None, start=None, end=None):
    for col in samfile.pileup(chr, start, end):
        yield construct_rec_coverage2(samfile, col)

