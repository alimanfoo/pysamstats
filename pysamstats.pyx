"""
TODO doc me

"""

# standard library imports
from collections import OrderedDict


# 3rd party imports
import pysam


cdef class CountAggregator(object):
    cdef int n
    cdef int m
    def __init__(self, n):
        self.n = n
        self.m = 0
    cpdef get_stats(self):
        cdef float p
        p = self.m * 100. / self.n
        return self.m, p


cdef class AggReadsFwd(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not is_reverse:
            self.m += 1
    

cdef class AggReadsRev(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_reverse:
            self.m += 1
    

cdef class AggReadsProperPair(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair:
            self.m += 1
    

cdef class AggReadsProperPairFwd(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair and not is_reverse:
            self.m += 1
    

cdef class AggReadsProperPairRev(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair and is_reverse:
            self.m += 1
    
    
cdef class AggReadsMateUnmapped(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if aln.mate_is_unmapped:
            self.m += 1
    
    
cdef class AggReadsMateUnmappedFwd(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not is_reverse:
            if aln.mate_is_unmapped:
                self.m += 1
    
    
cdef class AggReadsMateUnmappedRev(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_reverse:
            if aln.mate_is_unmapped:
                self.m += 1
    
    
cdef class AggReadsMateOtherChr(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped and aln.rnext != aln.tid: 
            self.m += 1
        
    
cdef class AggReadsMateOtherChrFwd(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not is_reverse:
            if not aln.mate_is_unmapped and aln.rnext != aln.tid: 
                self.m += 1
    
    
cdef class AggReadsMateOtherChrRev(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_reverse:
            if not aln.mate_is_unmapped and aln.rnext != aln.tid:
                self.m += 1
    
    
cdef class AggReadsMateSameStrand(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if is_reverse:
                if aln.mate_is_reverse: 
                    self.m += 1
            elif not aln.mate_is_reverse:
                    self.m += 1
    
    
cdef class AggReadsMateSameStrandFwd(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if not is_reverse:
                if not aln.mate_is_reverse:
                    self.m += 1
    
    
cdef class AggReadsMateSameStrandRev(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if is_reverse:
                if aln.mate_is_reverse:
                    self.m += 1
    
    
cdef class AggReadsFaceaway(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if is_reverse:
                if aln.tlen > 0: # mapped to reverse strand but leftmost
                    self.m += 1
            elif aln.tlen < 0: # mapped to fwd strand but rightmost
                    self.m += 1

    
cdef class AggReadsFaceawayFwd(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if not is_reverse:
                if aln.tlen < 0: # mapped to fwd strand but rightmost
                    self.m += 1
    
    
cdef class AggReadsFaceawayRev(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not aln.mate_is_unmapped:
            if is_reverse:
                if aln.tlen > 0: # mapped to rev strand but leftmost
                    self.m += 1

    
cdef class AggReadsMapq0(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if aln.mapq == 0:
            self.m += 1

    
cdef class AggReadsMapq0Fwd(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if not is_reverse:
            if aln.mapq == 0:
                self.m += 1

    
cdef class AggReadsMapq0Rev(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_reverse:
            if aln.mapq == 0:
                self.m += 1

    
cdef class AggReadsMapq0ProperPair(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair:
            if aln.mapq == 0:
                self.m += 1

    
cdef class AggReadsMapq0ProperPairFwd(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair and not is_reverse:
            if aln.mapq == 0:
                self.m += 1

    
cdef class AggReadsMapq0ProperPairRev(CountAggregator):
    cpdef add_read(self, read, aln, bint is_proper_pair, bint is_reverse):
        if is_proper_pair and is_reverse:
            if aln.mapq == 0:
                self.m += 1

    
class CoverageStatsTable(object):
    
    def __init__(self, samfn, chr=None, start=None, end=None):
        self.samfn = samfn
        self.chr = chr
        self.start = start
        self.end = end
        
    def __iter__(self):
        cdef:
            int alen
            int n
            int ri
            int ai
            bint is_proper_pair
            bint is_reverse
        
        # set up aggregators
        a = OrderedDict()
        a['reads_fwd'] = AggReadsFwd
        a['reads_rev'] = AggReadsRev
        a['reads_pp'] = AggReadsProperPair
        a['reads_pp_fwd'] = AggReadsProperPairFwd 
        a['reads_pp_rev'] = AggReadsProperPairRev 
        a['reads_mate_unmapped'] = AggReadsMateUnmapped
        a['reads_mate_unmapped_fwd'] = AggReadsMateUnmappedFwd
        a['reads_mate_unmapped_rev'] = AggReadsMateUnmappedRev
        a['reads_mate_other_chr'] = AggReadsMateOtherChr
        a['reads_mate_other_chr_fwd'] = AggReadsMateOtherChrFwd
        a['reads_mate_other_chr_rev'] = AggReadsMateOtherChrRev
        a['reads_mate_same_strand'] = AggReadsMateSameStrand
        a['reads_mate_same_strand_fwd'] = AggReadsMateSameStrandFwd
        a['reads_mate_same_strand_rev'] = AggReadsMateSameStrandRev
        a['reads_faceaway'] = AggReadsFaceaway
        a['reads_faceaway_fwd'] = AggReadsFaceawayFwd
        a['reads_faceaway_rev'] = AggReadsFaceawayRev
        a['reads_mapq0'] = AggReadsMapq0
        a['reads_mapq0_fwd'] = AggReadsMapq0Fwd
        a['reads_mapq0_rev'] = AggReadsMapq0Rev
        a['reads_mapq0_pp'] = AggReadsMapq0ProperPair
        a['reads_mapq0_pp_fwd'] = AggReadsMapq0ProperPairFwd
        a['reads_mapq0_pp_rev'] = AggReadsMapq0ProperPairRev
        alen = len(a)
        
        # define header
        header = ['reads'] 
        header.extend(a.keys()) # count fields
        header.extend('pc_%s' % f for f in a.keys()) # percentage fields
        yield header
        
        # open sam file
        sam = pysam.Samfile(self.samfn)
        
        # run pileup
        for col in sam.pileup(self.chr, self.start, self.end):
            
            n = col.n
            reads = col.pileups
            
            # instantiate aggregators
            ags = [c(n) for c in a.values()]
            
            # iterate over reads in the column
            for ri in range(n):
                read = reads[ri]
                aln = read.alignment
                is_proper_pair = aln.is_proper_pair
                is_reverse = aln.is_reverse
                # iterate over aggregators
                for ai in range(alen):
                    ag = ags[ai]
                    ag.add_read(read, aln, is_proper_pair, is_reverse)
                    
            # construct output row
            row = [n] # start with raw read depth
            data = [ag.get_stats() for ag in ags]
            row.extend(d[0] for d in data) # add counts
            row.extend(d[1] for d in data) # add percentages
            yield tuple(row)
                    
    
        
        
    
        
