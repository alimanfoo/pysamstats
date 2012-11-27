

cdef struct CountStrand:
    int all
    int fwd
    int rev


cdef struct CountPpStrand:
    int all
    int fwd
    int rev
    int pp
    int pp_fwd
    int pp_rev


cdef struct RecCoverage:
    char* chr
    int pos
    CountPpStrand reads


cdef void init_pp_strand(CountPpStrand *c):
    c.all = 0
    c.rev = 0
    c.fwd = 0
    c.pp = 0
    c.pp_fwd = 0
    c.pp_rev = 0


cdef void incr_pp_strand(CountPpStrand *c, bint is_reverse, bint is_proper_pair):
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
    cdef int i # loop index
    cdef int n # total number of reads in column
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
