

cdef struct CountStrand:
    unsigned int all
    unsigned int fwd
    unsigned int rev


cdef struct CountPpStrand:
    unsigned int all
    unsigned int fwd
    unsigned int rev
    unsigned int pp
    unsigned int pp_fwd
    unsigned int pp_rev


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
    reads = col.pileups()
    for i in range(n):
        read = reads[i]
        is_reverse = read.alignment.is_reverse
        is_proper_pair = read.alignment.is_proper_pair
        incr_pp_strand(&rec.reads, is_reverse, is_proper_pair)

    return rec


def stat_coverage(samfile, chr=None, start=None, end=None):
    for col in samfile.pileup(chr, start, end):
        yield construct_rec_coverage(samfile, col)

