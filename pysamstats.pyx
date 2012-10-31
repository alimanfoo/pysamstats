# cython: profile=True

"""
TODO doc me

"""

# standard library imports


# 3rd party imports
from csamtools cimport Samfile, PileupRead, AlignedRead, PileupProxy


cdef class CountAggregator:
    cdef int n
    cdef int all
    cdef int fwd
    cdef int rev
    
    def __cinit__(self, n):
        self.n = n
        self.all = 0
        self.fwd = 0
        self.rev = 0
    
    # to be overridden in subclasses    
    cdef bint test(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse, bint mate_is_unmapped):
        return 0

    cdef add(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse, bint mate_is_unmapped):
        if self.test(read, aln, is_proper_pair, is_reverse, mate_is_unmapped):
            self.all += 1
            if is_reverse:
                self.rev += 1
            else:
                self.fwd += 1
        

cdef class ProperPairCountAggregator(CountAggregator):
    cdef int pp
    cdef int pp_fwd
    cdef int pp_rev
    
    def __cinit__(self, n):
        self.n = n
        self.all = 0
        self.fwd = 0
        self.rev = 0
        self.pp = 0
        self.pp_fwd = 0
        self.pp_rev = 0
    
    cdef add(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse, bint mate_is_unmapped):
        if self.test(read, aln, is_proper_pair, is_reverse, mate_is_unmapped):
            self.all += 1
            if is_reverse:
                self.rev += 1
            else:
                self.fwd += 1
            if is_proper_pair:
                self.pp += 1
                if is_reverse:
                    self.pp_rev += 1
                else:
                    self.pp_fwd += 1
        

cdef class AggReads(ProperPairCountAggregator):
    cdef bint test(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse, bint mate_is_unmapped):
        return 1


cdef class AggReadsMateUnmapped(CountAggregator):
    cdef bint test(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse, bint mate_is_unmapped):
        return mate_is_unmapped
    
    
cdef class AggReadsMateOtherChr(CountAggregator):
    cdef bint test(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse, bint mate_is_unmapped):
        return not mate_is_unmapped and aln.rnext != aln.tid 
        
    
cdef class AggReadsMateSameStrand(CountAggregator):
    cdef bint test(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse, bint mate_is_unmapped):
        cdef bint mate_is_reverse
        if not mate_is_unmapped:
            mate_is_reverse = aln.mate_is_reverse
            return (is_reverse and mate_is_reverse) or (not is_reverse and not mate_is_reverse)
        else:
            return 0
    
    
cdef class AggReadsFaceaway(CountAggregator):
    cdef bint test(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse, bint mate_is_unmapped):
        cdef int tlen
        if not mate_is_unmapped:
            tlen = aln.tlen
            return ((is_reverse and tlen > 0) # mapped to reverse strand but leftmost
                    or (not is_reverse and tlen < 0)) # mapped to fwd strand but rightmost
        else:
            return 0
    
    
cdef class AggReadsEdit0(ProperPairCountAggregator):
    cdef bint test(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse, bint mate_is_unmapped):
        return aln.opt('NM') == 0
    
    
cdef class AggReadsSoftClipped(ProperPairCountAggregator):
    cdef bint test(self, PileupRead read, AlignedRead aln, bint is_proper_pair, bint is_reverse, bint mate_is_unmapped):
        cigar = aln.cigar
        cdef int i = 0
        for i in range(len(cigar)):
            op = cigar[i]
            if op[0] == 4: # softclip code
                return 1
        return 0
    
    
cpdef build_coverage_stats(PileupProxy col):
    cdef int n = col.n
    cdef int ri
    cdef bint is_proper_pair
    cdef bint is_reverse
    cdef bint mate_is_unmapped
    cdef PileupRead read
    cdef AlignedRead aln

    # create aggregators
    agg_reads = AggReads(n)
    agg_reads_mate_unmapped = AggReadsMateUnmapped(n)
    agg_reads_mate_other_chr = AggReadsMateOtherChr(n)
    agg_reads_mate_same_strand = AggReadsMateSameStrand(n)
    agg_reads_faceaway = AggReadsFaceaway(n)
    agg_reads_edit0 = AggReadsEdit0(n)
    agg_reads_softclipped = AggReadsSoftClipped(n)
    
    # access reads
    reads = col.pileups

    # iterate over reads in the column
    for ri in range(n):
        read = reads[ri]
        aln = read.alignment
        
        # optimisation - access these now so done only once
        is_proper_pair = aln.is_proper_pair
        is_reverse = aln.is_reverse
        mate_is_unmapped = aln.mate_is_unmapped
        
        # pass reads to aggregators
        agg_reads.add(read, aln, is_proper_pair, is_reverse, mate_is_unmapped)
        agg_reads_mate_unmapped.add(read, aln, is_proper_pair, is_reverse, mate_is_unmapped)
        agg_reads_mate_other_chr.add(read, aln, is_proper_pair, is_reverse, mate_is_unmapped)
        agg_reads_mate_same_strand.add(read, aln, is_proper_pair, is_reverse, mate_is_unmapped)
        agg_reads_faceaway.add(read, aln, is_proper_pair, is_reverse, mate_is_unmapped)
        agg_reads_edit0.add(read, aln, is_proper_pair, is_reverse, mate_is_unmapped)
        agg_reads_softclipped.add(read, aln, is_proper_pair, is_reverse, mate_is_unmapped)
            
    # construct output row
    data = {
            'reads_fwd': agg_reads.fwd,
            'reads_rev': agg_reads.rev,
            'reads_pp': agg_reads.pp,
            'reads_pp_fwd': agg_reads.pp_fwd,
            'reads_pp_rev': agg_reads.pp_rev,
            'reads_mate_unmapped': agg_reads_mate_unmapped.all,
            'reads_mate_unmapped_fwd': agg_reads_mate_unmapped.fwd,
            'reads_mate_unmapped_rev': agg_reads_mate_unmapped.rev,
            'reads_mate_other_chr': agg_reads_mate_other_chr.all,
            'reads_mate_other_chr_fwd': agg_reads_mate_other_chr.fwd,
            'reads_mate_other_chr_rev': agg_reads_mate_other_chr.rev,
            'reads_mate_same_strand': agg_reads_mate_same_strand.all,
            'reads_mate_same_strand_fwd': agg_reads_mate_same_strand.fwd,
            'reads_mate_same_strand_rev': agg_reads_mate_same_strand.rev,
            'reads_faceaway': agg_reads_faceaway.all,
            'reads_faceaway_fwd': agg_reads_faceaway.fwd,
            'reads_faceaway_rev': agg_reads_faceaway.rev,
            'reads_edit0': agg_reads_edit0.all,                      
            'reads_edit0_fwd': agg_reads_edit0.fwd,
            'reads_edit0_rev': agg_reads_edit0.rev,                            
            'reads_edit0_pp': agg_reads_edit0.pp,     
            'reads_edit0_pp_fwd': agg_reads_edit0.pp_fwd,
            'reads_edit0_pp_rev': agg_reads_edit0.pp_rev,      
            'reads_softclipped': agg_reads_softclipped.all,                      
            'reads_softclipped_fwd': agg_reads_softclipped.fwd,
            'reads_softclipped_rev': agg_reads_softclipped.rev,                            
            'reads_softclipped_pp': agg_reads_softclipped.pp,     
            'reads_softclipped_pp_fwd': agg_reads_softclipped.pp_fwd,
            'reads_softclipped_pp_rev': agg_reads_softclipped.pp_rev,      
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
        computed_variables = [
                              'reads_fwd', 
                              'reads_rev',
                              'reads_pp',
                              'reads_pp_fwd',
                              'reads_pp_rev',
                              'reads_mate_unmapped',
                              'reads_mate_unmapped_fwd',
                              'reads_mate_unmapped_rev',
                              'reads_mate_other_chr',
                              'reads_mate_other_chr_fwd',
                              'reads_mate_other_chr_rev',
                              'reads_mate_same_strand',
                              'reads_mate_same_strand_fwd',
                              'reads_mate_same_strand_rev',
                              'reads_faceaway',
                              'reads_faceaway_fwd',
                              'reads_faceaway_rev',
                              'reads_edit0',                            
                              'reads_edit0_fwd',
                              'reads_edit0_rev',                            
                              'reads_edit0_pp',                            
                              'reads_edit0_pp_fwd',
                              'reads_edit0_pp_rev',
                              'reads_softclipped',                         
                              'reads_softclipped_fwd',                         
                              'reads_softclipped_rev',                         
                              'reads_softclipped_pp',
                              'reads_softclipped_pp_fwd',                         
                              'reads_softclipped_pp_rev',                         
                              ]
        header = fixed_variables + computed_variables
        yield header
        
        # open sam file
        sam = Samfile(self.samfn)
        
        # run pileup
        for col in sam.pileup(self.chr, self.start, self.end):
            
            # fixed variables            
            chr = sam.getrname(col.tid)
            pos = col.pos + 1 # 1-based
            row = [chr, pos, col.n] 
            
            # computed variables
            data = build_coverage_stats(col)
            row.extend(data[v] for v in computed_variables) 
            yield row


