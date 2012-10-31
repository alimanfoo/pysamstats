# cython: profile=True

"""
TODO doc me

"""

# standard library imports


# 3rd party imports
from csamtools cimport Samfile, PileupRead, AlignedRead, PileupProxy


cdef class CountAggregator:
    cdef int n
    cdef int m
    def __cinit__(self, n):
        self.n = n
        self.m = 0
    cdef get_result(self):
        return self.m


cdef class AggReadsFwd(CountAggregator):
    cdef add_read(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse):
        if not is_reverse:
            self.m += 1
    

cdef class AggReadsRev(CountAggregator):
    cdef add_read(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse):
        if is_reverse:
            self.m += 1
    

cdef class AggReadsProperPair(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair:
            self.m += 1
    

cdef class AggReadsProperPairFwd(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair and not is_reverse:
            self.m += 1
    

cdef class AggReadsProperPairRev(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair and is_reverse:
            self.m += 1
    
    
cdef class AggReadsMateUnmapped(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if aln.mate_is_unmapped:
            self.m += 1
    
    
cdef class AggReadsMateUnmappedFwd(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not is_reverse:
            if aln.mate_is_unmapped:
                self.m += 1
    
    
cdef class AggReadsMateUnmappedRev(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_reverse:
            if aln.mate_is_unmapped:
                self.m += 1
    
    
cdef class AggReadsMateOtherChr(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped and aln.rnext != aln.tid: 
            self.m += 1
        
    
cdef class AggReadsMateOtherChrFwd(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not is_reverse:
            if not aln.mate_is_unmapped and aln.rnext != aln.tid: 
                self.m += 1
    
    
cdef class AggReadsMateOtherChrRev(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_reverse:
            if not aln.mate_is_unmapped and aln.rnext != aln.tid:
                self.m += 1
    
    
cdef class AggReadsMateSameStrand(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if is_reverse:
                if aln.mate_is_reverse: 
                    self.m += 1
            elif not aln.mate_is_reverse:
                    self.m += 1
    
    
cdef class AggReadsMateSameStrandFwd(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if not is_reverse:
                if not aln.mate_is_reverse:
                    self.m += 1
    
    
cdef class AggReadsMateSameStrandRev(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if is_reverse:
                if aln.mate_is_reverse:
                    self.m += 1
    
    
cdef class AggReadsFaceaway(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if is_reverse:
                if aln.tlen > 0: # mapped to reverse strand but leftmost
                    self.m += 1
            elif aln.tlen < 0: # mapped to fwd strand but rightmost
                    self.m += 1

    
cdef class AggReadsFaceawayFwd(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if not is_reverse:
                if aln.tlen < 0: # mapped to fwd strand but rightmost
                    self.m += 1
    
    
cdef class AggReadsFaceawayRev(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if is_reverse:
                if aln.tlen > 0: # mapped to rev strand but leftmost
                    self.m += 1

    
cdef class AggReadsMapq0(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if aln.mapq == 0:
            self.m += 1

    
cdef class AggReadsMapq0Fwd(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not is_reverse:
            if aln.mapq == 0:
                self.m += 1

    
cdef class AggReadsMapq0Rev(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_reverse:
            if aln.mapq == 0:
                self.m += 1

    
cdef class AggReadsMapq0ProperPair(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair:
            if aln.mapq == 0:
                self.m += 1

    
cdef class AggReadsMapq0ProperPairFwd(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair and not is_reverse:
            if aln.mapq == 0:
                self.m += 1

    
cdef class AggReadsMapq0ProperPairRev(CountAggregator):
    cdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair and is_reverse:
            if aln.mapq == 0:
                self.m += 1


cpdef build_coverage_stats(PileupProxy col):
    cdef int n = col.n
    cdef int ri
    cdef bint is_proper_pair
    cdef bint is_reverse
    cdef PileupRead read
    cdef AlignedRead aln

    reads_fwd = AggReadsFwd(n)
    reads_rev = AggReadsRev(n)
    
    reads = col.pileups

    # iterate over reads in the column
    for ri in range(n):
        read = reads[ri]
        aln = read._alignment
        is_proper_pair = aln.is_proper_pair
        is_reverse = aln.is_reverse
        # pass reads to aggregators
        reads_fwd.add_read(read, aln, is_proper_pair, is_reverse)
        reads_rev.add_read(read, aln, is_proper_pair, is_reverse)
            
    # construct output row
    data = {
            'reads_fwd': reads_fwd.get_result(),
            'pc_reads_fwd': reads_fwd.get_result() * 100. / n,
            'reads_rev': reads_rev.get_result(),
            'pc_reads_rev': reads_rev.get_result() * 100. / n,
            }
    return data

    
class CoverageStatsTable(object):
    
    def __init__(self, samfn, chr=None, start=None, end=None):
        self.samfn = samfn
        self.chr = chr
        self.start = start
        self.end = end
        
    def __iter__(self):

        # define header
        fixed_variables = ['chr', 'pos', 'reads']
        computed_variables = ['reads_fwd', 'reads_rev', 'pc_reads_fwd', 'pc_reads_rev']
        header = fixed_variables + computed_variables
        yield header
        
        # open sam file
        sam = Samfile(self.samfn)
        
        # run pileup
        for col in sam.pileup(self.chr, self.start, self.end):
            
            # fixed variables            
            chr = sam.getrname(col.tid)
            pos = col.pos + 1 # 1-based
            row = [chr, pos, col.n] # start with raw read depth
            
            # computed variables
            data = build_coverage_stats(col)
            row.extend(data[v] for v in computed_variables) 
            yield row


