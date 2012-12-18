# cython: profile=True


import sys
#import numpy as np
#import numpy.ma as ma
#cimport numpy as np
import time
import csv
from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t
from libc.math cimport sqrt
from cpython cimport PyBytes_FromStringAndSize
from csamtools cimport Samfile, Fastafile, PileupProxy, bam1_t, bam_pileup1_t, bam1_cigar, bam1_seq, bam1_qual


## These are bits set in the flag.
## have to put these definitions here, in csamtools.pxd they got ignored
## @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
DEF BAM_FPAIRED       =1
## @abstract the read is mapped in a proper pair */
DEF BAM_FPROPER_PAIR  =2
## @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
DEF BAM_FUNMAP        =4
## @abstract the mate is unmapped */
DEF BAM_FMUNMAP       =8
## @abstract the read is mapped to the reverse strand */
DEF BAM_FREVERSE      =16
## @abstract the mate is mapped to the reverse strand */
DEF BAM_FMREVERSE     =32
## @abstract this is read1 */
DEF BAM_FREAD1        =64
## @abstract this is read2 */
DEF BAM_FREAD2       =128
## @abstract not primary alignment */
DEF BAM_FSECONDARY   =256
## @abstract QC failure */
DEF BAM_FQCFAIL      =512
## @abstract optical or PCR duplicate */
DEF BAM_FDUP        =1024

# CIGAR operations
DEF BAM_CIGAR_SHIFT=4
DEF BAM_CIGAR_MASK=((1 << BAM_CIGAR_SHIFT) - 1)

DEF BAM_CMATCH     = 0
DEF BAM_CINS       = 1
DEF BAM_CDEL       = 2
DEF BAM_CREF_SKIP  = 3
DEF BAM_CSOFT_CLIP = 4
DEF BAM_CHARD_CLIP = 5
DEF BAM_CPAD       = 6
DEF BAM_CEQUAL     = 7
DEF BAM_CDIFF      = 8

cdef char * bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"


#############################
# BASIC COVERAGE STATISTICS #
#############################


cpdef object construct_rec_coverage(Samfile samfile, PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef unsigned int reads_pp = 0

    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        if is_proper_pair:
            reads_pp += 1

    return {'chr': chrom, 
            'pos': pos, 
            'reads_all': n, 
            'reads_pp': reads_pp}


def stat_coverage(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_coverage(samfile, col, one_based)
        
        
def write_coverage(*args, **kwargs):
    fieldnames = ('chr', 'pos', 'reads_all', 'reads_pp')
    write_stats(stat_coverage, fieldnames, *args, **kwargs)
    
    
################################
# STRANDED COVERAGE STATISTICS #
################################


cpdef object construct_rec_coverage_strand(Samfile samfile, PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef uint32_t flag
    cdef bint is_reverse 
    cdef bint is_proper_pair 
    cdef unsigned int reads_fwd = 0
    cdef unsigned int reads_rev = 0
    cdef unsigned int reads_pp = 0
    cdef unsigned int reads_pp_fwd = 0
    cdef unsigned int reads_pp_rev = 0
    
    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_reverse = <bint>(flag & BAM_FREVERSE)
        if is_reverse:
            reads_rev += 1
        else:
            reads_fwd += 1
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        if is_proper_pair:
            reads_pp += 1
            if is_reverse:
                reads_pp_rev += 1
            else:
                reads_pp_fwd += 1

    return {'chr': chrom, 
            'pos': pos, 
            'reads_all': n, 
            'reads_fwd': reads_fwd, 
            'reads_rev': reads_rev, 
            'reads_pp': reads_pp,
            'reads_pp_fwd': reads_pp_fwd,
            'reads_pp_rev': reads_pp_rev}


def stat_coverage_strand(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(chrom, start, end):
        yield construct_rec_coverage_strand(samfile, col, one_based)
        
        
def write_coverage_strand(*args, **kwargs):
    fieldnames = ('chr', 'pos', 
                  'reads_all', 'reads_fwd', 'reads_rev', 
                  'reads_pp', 'reads_pp_fwd', 'reads_pp_rev')
    write_stats(stat_coverage_strand, fieldnames, *args, **kwargs)
    
    
################################
# EXTENDED COVERAGE STATISTICS #
################################


cpdef object construct_rec_coverage_ext(Samfile samfile, PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef bint is_reverse 
    cdef bint is_proper_pair 
    cdef bint mate_is_unmappped 
    cdef bint mate_is_reverse
    cdef int tlen
    # counting variables 
    cdef unsigned int reads_pp = 0
    cdef unsigned int reads_mate_unmapped = 0
    cdef unsigned int reads_mate_other_chr = 0
    cdef unsigned int reads_mate_same_strand = 0
    cdef unsigned int reads_faceaway = 0
    cdef unsigned int reads_softclipped = 0

    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    tid = col.tid
    chrom = samfile.getrname(tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_reverse = <bint>(flag & BAM_FREVERSE)
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
        mate_is_reverse = <bint>(flag & BAM_FMREVERSE)
        tlen = aln.core.isize
        if is_proper_pair:
            reads_pp += 1
        if mate_is_unmapped:
            reads_mate_unmapped += 1
        elif tid != aln.core.mtid:
            reads_mate_other_chr += 1
        elif (is_reverse and mate_is_reverse) or (not is_reverse and not mate_is_reverse):
            reads_mate_same_strand += 1
        elif (is_reverse and tlen > 0) or (not is_reverse and tlen < 0):
            reads_faceaway += 1
        if is_softclipped(aln):
            reads_softclipped += 1
            
    return {'chr': chrom, 'pos': pos, 
               'reads_all': n, 
               'reads_pp': reads_pp,
               'reads_mate_unmapped': reads_mate_unmapped,
               'reads_mate_other_chr': reads_mate_other_chr,
               'reads_mate_same_strand': reads_mate_same_strand,
               'reads_faceaway': reads_faceaway,
               'reads_softclipped': reads_softclipped}


def stat_coverage_ext(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(chrom, start, end):
        yield construct_rec_coverage_ext(samfile, col, one_based)
        
        
def write_coverage_ext(*args, **kwargs):
    fieldnames = ('chr', 'pos', 
                  'reads_all', 
                  'reads_pp', 
                  'reads_mate_unmapped', 
                  'reads_mate_other_chr',
                  'reads_mate_same_strand',
                  'reads_faceaway', 
                  'reads_softclipped')
    write_stats(stat_coverage_ext, fieldnames, *args, **kwargs)
    
    
##########################################
# EXTENDED COVERAGE STATISTICS BY STRAND #
##########################################


cpdef object construct_rec_coverage_ext_strand(Samfile samfile, PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef bint is_reverse 
    cdef bint is_proper_pair 
    cdef bint mate_is_unmappped 
    cdef bint mate_is_reverse
    cdef int tlen
    # counting variables 
    cdef unsigned int reads_rev = 0
    cdef unsigned int reads_fwd = 0
    cdef unsigned int reads_pp = 0
    cdef unsigned int reads_pp_rev = 0
    cdef unsigned int reads_pp_fwd = 0
    cdef unsigned int reads_mate_unmapped = 0
    cdef unsigned int reads_mate_unmapped_rev = 0
    cdef unsigned int reads_mate_unmapped_fwd = 0
    cdef unsigned int reads_mate_other_chr = 0
    cdef unsigned int reads_mate_other_chr_rev = 0
    cdef unsigned int reads_mate_other_chr_fwd = 0
    cdef unsigned int reads_mate_same_strand = 0
    cdef unsigned int reads_mate_same_strand_rev = 0
    cdef unsigned int reads_mate_same_strand_fwd = 0
    cdef unsigned int reads_faceaway = 0
    cdef unsigned int reads_faceaway_rev = 0
    cdef unsigned int reads_faceaway_fwd = 0
    cdef unsigned int reads_softclipped = 0
    cdef unsigned int reads_softclipped_rev = 0
    cdef unsigned int reads_softclipped_fwd = 0

    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    tid = col.tid
    chrom = samfile.getrname(tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_reverse = <bint>(flag & BAM_FREVERSE)
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
        mate_is_reverse = <bint>(flag & BAM_FMREVERSE)
        tlen = aln.core.isize
        if is_reverse:
            reads_rev += 1
        else:
            reads_fwd += 1
        if is_proper_pair:
            reads_pp += 1
            if is_reverse:
                reads_pp_rev += 1
            else:
                reads_pp_fwd += 1
        if mate_is_unmapped:
            reads_mate_unmapped += 1
            if is_reverse:
                reads_mate_unmapped_rev += 1
            else:
                reads_mate_unmapped_fwd += 1
        elif tid != aln.core.mtid:
            reads_mate_other_chr += 1
            if is_reverse:
                reads_mate_other_chr_rev += 1
            else:
                reads_mate_other_chr_fwd += 1
        elif is_reverse and mate_is_reverse:
            reads_mate_same_strand += 1
            reads_mate_same_strand_rev += 1
        elif not is_reverse and not mate_is_reverse:
            reads_mate_same_strand += 1
            reads_mate_same_strand_fwd += 1
        elif (is_reverse and tlen > 0) or (not is_reverse and tlen < 0):
            reads_faceaway += 1
            if is_reverse:
                reads_faceaway_rev += 1
            else:
                reads_faceaway_fwd += 1
        if is_softclipped(aln):
            reads_softclipped += 1
            if is_reverse:
                reads_softclipped_rev += 1
            else:
                reads_softclipped_fwd += 1
            
    return {'chr': chrom, 'pos': pos, 
           'reads_all': n, 
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


def stat_coverage_ext_strand(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(chrom, start, end):
        yield construct_rec_coverage_ext_strand(samfile, col, one_based)
        
        
def write_coverage_ext_strand(*args, **kwargs):
    fieldnames = ('chr', 'pos', 
                  'reads_all', 
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
                  'reads_softclipped',
                  'reads_softclipped_fwd',
                  'reads_softclipped_rev')
    write_stats(stat_coverage_ext_strand, fieldnames, *args, **kwargs)


########################
# VARIATION STATISTICS #
########################


cpdef object construct_rec_variation(Samfile samfile, Fastafile fafile, 
                                     PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    # counting variables
    cdef unsigned int reads_pp = 0
    cdef unsigned int matches = 0
    cdef unsigned int matches_pp = 0
    cdef unsigned int mismatches = 0
    cdef unsigned int mismatches_pp = 0
    cdef unsigned int deletions = 0
    cdef unsigned int deletions_pp = 0
    cdef unsigned int insertions = 0
    cdef unsigned int insertions_pp = 0
    cdef unsigned int A = 0
    cdef unsigned int A_pp = 0
    cdef unsigned int C = 0
    cdef unsigned int C_pp = 0
    cdef unsigned int T = 0
    cdef unsigned int T_pp = 0
    cdef unsigned int G = 0
    cdef unsigned int G_pp = 0
    cdef unsigned int N = 0
    cdef unsigned int N_pp = 0
    
    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # reference base
    refbase = fafile.fetch(reference=chrom, start=col.pos, end=col.pos+1).upper()
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        # read.qpos
        # read.is_del
        # read.indel
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        if is_proper_pair:
            reads_pp += 1
        if read.is_del:
            deletions += 1
            if is_proper_pair:
                deletions_pp += 1
        else:
#            alnbase = get_seq_range(aln, 0, aln.core.l_qseq)[read.qpos]
            alnbase = get_seq_base(aln, read.qpos)
#            print refbase, alnbase
            if alnbase == 'A':
                A += 1
                if is_proper_pair:
                    A_pp += 1
            elif alnbase == 'T':
                T += 1
                if is_proper_pair:
                    T_pp += 1
            elif alnbase == 'C':
                C += 1
                if is_proper_pair:
                    C_pp += 1
            elif alnbase == 'G':
                G += 1
                if is_proper_pair:
                    G_pp += 1
            elif alnbase == 'N':
                N += 1
                if is_proper_pair:
                    N_pp += 1
            if read.indel > 0:
                insertions += 1
                if is_proper_pair:
                    insertions_pp += 1
            if alnbase == refbase:
                matches += 1
                if is_proper_pair:
                    matches_pp += 1
            else:
                mismatches += 1
                if is_proper_pair:
                    mismatches_pp += 1

    return {'chr': chrom, 'pos': pos, 'ref': refbase,
            'reads_all': n, 'reads_pp': reads_pp,
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


def stat_variation(samfile, fafile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_variation(samfile, fafile, col, one_based)
        
        
def write_variation(*args, **kwargs):
    fieldnames = ('chr', 'pos', 'ref', 
                  'reads_all', 'reads_pp',
                  'matches', 'matches_pp',
                  'mismatches', 'mismatches_pp',
                  'deletions', 'deletions_pp',
                  'insertions', 'insertions_pp',
                  'A', 'A_pp', 'C', 'C_pp', 'T', 'T_pp', 'G', 'G_pp', 'N', 'N_pp')
    write_stats(stat_variation, fieldnames, *args, **kwargs)
    
    
#################################
# STRANDED VARIATION STATISTICS #
#################################


cdef struct CountPpStrand:
    unsigned int all, pp, fwd, rev, pp_fwd, pp_rev


cdef inline init_pp_strand(CountPpStrand* c):
    c.all = c.fwd = c.rev = c.pp = c.pp_fwd = c.pp_rev = 0    


cdef inline incr_pp_strand(CountPpStrand* c, bint is_reverse, bint is_proper_pair):
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
                

cpdef object construct_rec_variation_strand(Samfile samfile, Fastafile fafile, 
                                            PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair, is_reverse
    # counting variables
    cdef CountPpStrand reads, matches, mismatches, deletions, insertions, A, C, T, G, N
    
    # initialise variables
    n = col.n
    plp = col.plp
    init_pp_strand(&reads)
    init_pp_strand(&matches)
    init_pp_strand(&mismatches)
    init_pp_strand(&deletions)
    init_pp_strand(&insertions)
    init_pp_strand(&A)
    init_pp_strand(&T)
    init_pp_strand(&C)
    init_pp_strand(&G)
    init_pp_strand(&N)

    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # reference base
    refbase = fafile.fetch(reference=chrom, start=col.pos, end=col.pos+1).upper()
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        # read.qpos
        # read.is_del
        # read.indel
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)
        incr_pp_strand(&reads, is_reverse, is_proper_pair)
        if read.is_del:
            incr_pp_strand(&deletions, is_reverse, is_proper_pair)
        else:
            alnbase = get_seq_base(aln, read.qpos)
            if alnbase == 'A':
                incr_pp_strand(&A, is_reverse, is_proper_pair)
            elif alnbase == 'T':
                incr_pp_strand(&T, is_reverse, is_proper_pair)
            elif alnbase == 'C':
                incr_pp_strand(&C, is_reverse, is_proper_pair)
            elif alnbase == 'G':
                incr_pp_strand(&G, is_reverse, is_proper_pair)
            elif alnbase == 'N':
                incr_pp_strand(&N, is_reverse, is_proper_pair)
            if read.indel > 0:
                incr_pp_strand(&insertions, is_reverse, is_proper_pair)
            if alnbase == refbase:
                incr_pp_strand(&matches, is_reverse, is_proper_pair)
            else:
                incr_pp_strand(&mismatches, is_reverse, is_proper_pair)

    return {'chr': chrom, 'pos': pos, 'ref': refbase,
            'reads_all': n, 'reads_fwd': reads.fwd, 'reads_rev': reads.rev, 
            'reads_pp': reads.pp, 'reads_pp_fwd': reads.pp_fwd, 'reads_pp_rev': reads.pp_rev,
            'matches': matches.all, 'matches_fwd': matches.fwd, 'matches_rev': matches.rev, 
            'matches_pp': matches.pp, 'matches_pp_fwd': matches.pp_fwd, 'matches_pp_rev': matches.pp_rev,
            'mismatches': mismatches.all, 'mismatches_fwd': mismatches.fwd, 'mismatches_rev': mismatches.rev, 
            'mismatches_pp': mismatches.pp, 'mismatches_pp_fwd': mismatches.pp_fwd, 'mismatches_pp_rev': mismatches.pp_rev,
            'deletions': deletions.all, 'deletions_fwd': deletions.fwd, 'deletions_rev': deletions.rev, 
            'deletions_pp': deletions.pp, 'deletions_pp_fwd': deletions.pp_fwd, 'deletions_pp_rev': deletions.pp_rev,
            'insertions': insertions.all, 'insertions_fwd': insertions.fwd, 'insertions_rev': insertions.rev, 
            'insertions_pp': insertions.pp, 'insertions_pp_fwd': insertions.pp_fwd, 'insertions_pp_rev': insertions.pp_rev,
            'A': A.all, 'A_fwd': A.fwd, 'A_rev': A.rev, 'A_pp': A.pp, 'A_pp_fwd': A.pp_fwd, 'A_pp_rev': A.pp_rev,
            'C': C.all, 'C_fwd': C.fwd, 'C_rev': C.rev, 'C_pp': C.pp, 'C_pp_fwd': C.pp_fwd, 'C_pp_rev': C.pp_rev,
            'T': T.all, 'T_fwd': T.fwd, 'T_rev': T.rev, 'T_pp': T.pp, 'T_pp_fwd': T.pp_fwd, 'T_pp_rev': T.pp_rev,
            'G': G.all, 'G_fwd': G.fwd, 'G_rev': G.rev, 'G_pp': G.pp, 'G_pp_fwd': G.pp_fwd, 'G_pp_rev': G.pp_rev,
            'N': N.all, 'N_fwd': N.fwd, 'N_rev': N.rev, 'N_pp': N.pp, 'N_pp_fwd': N.pp_fwd, 'N_pp_rev': N.pp_rev,
            }


def stat_variation_strand(samfile, fafile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_variation_strand(samfile, fafile, col, one_based)
        
        
def write_variation_strand(*args, **kwargs):
    fieldnames = ('chr', 'pos', 'ref', 
                'reads_all', 'reads_fwd', 'reads_rev', 'reads_pp', 'reads_pp_fwd', 'reads_pp_rev',
                'matches', 'matches_fwd', 'matches_rev', 'matches_pp', 'matches_pp_fwd', 'matches_pp_rev',
                'mismatches', 'mismatches_fwd', 'mismatches_rev', 'mismatches_pp', 'mismatches_pp_fwd', 'mismatches_pp_rev',
                'deletions', 'deletions_fwd', 'deletions_rev', 'deletions_pp', 'deletions_pp_fwd', 'deletions_pp_rev',
                'insertions', 'insertions_fwd', 'insertions_rev', 'insertions_pp', 'insertions_pp_fwd', 'insertions_pp_rev',
                'A', 'A_fwd', 'A_rev', 'A_pp', 'A_pp_fwd', 'A_pp_rev',
                'C', 'C_fwd', 'C_rev', 'C_pp', 'C_pp_fwd', 'C_pp_rev',
                'T', 'T_fwd', 'T_rev', 'T_pp', 'T_pp_fwd', 'T_pp_rev',
                'G', 'G_fwd', 'G_rev', 'G_pp', 'G_pp_fwd', 'G_pp_rev',
                'N', 'N_fwd', 'N_rev', 'N_pp', 'N_pp_fwd', 'N_pp_rev')
    write_stats(stat_variation_strand, fieldnames, *args, **kwargs)
    
    
##########################
# INSERT SIZE STATISTICS #
##########################


cpdef object construct_rec_tlen(Samfile samfile, PileupProxy col, 
                                bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef bint mate_is_unmappped 
    cdef bint mate_other_chr
    cdef unsigned int reads_p = 0 # reads "paired", i.e., mate is mapped to same chromosome, so tlen is meaningful
    cdef unsigned int reads_pp = 0 # reads "properly paired", as defined by aligner
    cdef int64_t tlen
    cdef int64_t tlen_squared
    cdef int64_t tlen_p_sum = 0
    cdef double tlen_p_mean = 0
    cdef double tlen_p_dev_squared
    cdef double tlen_p_dev_squared_sum = 0
    cdef int64_t tlen_p_squared_sum = 0
    cdef int64_t tlen_pp_sum = 0
    cdef double tlen_pp_mean = 0
    cdef double tlen_pp_dev_squared
    cdef double tlen_pp_dev_squared_sum = 0
    cdef int64_t tlen_pp_squared_sum = 0
    
    # initialise variables
    n = col.n
    plp = col.plp
    
    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # loop over reads
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag

        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
        mate_other_chr = <bint>(aln.core.tid != aln.core.mtid)

        # N.B., pysam exposes this property as 'tlen' rather than 'isize' so we 
        # follow their naming convention
        tlen = aln.core.isize 
        tlen_squared = tlen**2

        # N.B. insert size is only meaningful if mate is mapped to same chromosome
        if not mate_is_unmapped and not mate_other_chr:
            reads_p += 1
            tlen_p_sum += tlen
            tlen_p_squared_sum += tlen_squared
            if is_proper_pair:
                reads_pp += 1
                tlen_pp_sum += tlen
                tlen_pp_squared_sum += tlen_squared

    # calculate intermediate variables
    if reads_p > 0:
        tlen_p_mean = tlen_p_sum * 1. / reads_p
    if reads_pp > 0:
        tlen_pp_mean = tlen_pp_sum * 1. / reads_pp
        
    # loop over reads again to calculate variance (and hence std)
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
        mate_other_chr = <bint>(aln.core.tid != aln.core.mtid)
        tlen = aln.core.isize
        # N.B. insert size is only meaningful if mate is mapped to same chromosome
        if not mate_is_unmapped and not mate_other_chr:
            tlen_p_dev_squared = (tlen - tlen_p_mean)**2
            tlen_p_dev_squared_sum += tlen_p_dev_squared
            if is_proper_pair:
                tlen_pp_dev_squared = (tlen - tlen_pp_mean)**2
                tlen_pp_dev_squared_sum += tlen_pp_dev_squared

    # calculate output variables
    # N.B. round values to nearest integer, any finer precision is probably not
    # interesting    
    if reads_p > 0:
        mean_tlen = int(round(tlen_p_mean))
        rms_tlen = int(round(sqrt(tlen_p_squared_sum * 1. / reads_p)))
        variance_tlen = tlen_p_dev_squared_sum * 1. / reads_p
        std_tlen = int(round(sqrt(variance_tlen)))
    else:
        rms_tlen = std_tlen = mean_tlen = 'NA'
    if reads_pp > 0:
        mean_tlen_pp = int(round(tlen_pp_mean))
        rms_tlen_pp = int(round(sqrt(tlen_pp_squared_sum * 1. / reads_pp)))
        variance_tlen_pp = tlen_pp_dev_squared_sum * 1. / reads_pp
        std_tlen_pp = int(round(sqrt(variance_tlen_pp)))
    else:
        rms_tlen_pp = std_tlen_pp = mean_tlen_pp = 'NA'

    return {'chr': chrom, 
            'pos': pos, 
#            'reads_all': n, 
#            'reads_paired': reads_p,
#            'reads_pp': reads_pp,
            'mean_tlen': mean_tlen,
            'mean_tlen_pp': mean_tlen_pp,
            'rms_tlen': rms_tlen,
            'rms_tlen_pp': rms_tlen_pp,
            'std_tlen': std_tlen,
            'std_tlen_pp': std_tlen_pp}


def stat_tlen(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_tlen(samfile, col, one_based)
        
        
def write_tlen(*args, **kwargs):
    fieldnames = ('chr', 'pos', 
 #                 'reads_all', 'reads_paired', 'reads_pp', 
                  'mean_tlen', 'mean_tlen_pp',
                  'rms_tlen', 'rms_tlen_pp',
                  'std_tlen', 'std_tlen_pp')
    write_stats(stat_tlen, fieldnames, *args, **kwargs)
    
    
####################################
# INSERT SIZE STATISTICS BY STRAND #
####################################


cpdef object construct_rec_tlen_strand(Samfile samfile, PileupProxy col, 
                                       bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef bint mate_is_unmappped 
    cdef bint mate_other_chr
    
    # counting variables
    cdef unsigned int reads_fwd = 0
    cdef unsigned int reads_rev = 0
    cdef unsigned int reads_p = 0 # reads "paired", i.e., mate is mapped to same chromosome, so tlen is meaningful
    cdef unsigned int reads_p_fwd = 0
    cdef unsigned int reads_p_rev = 0
    cdef unsigned int reads_pp = 0 # reads "properly paired", as defined by aligner
    cdef unsigned int reads_pp_fwd = 0
    cdef unsigned int reads_pp_rev = 0
    
    cdef int64_t tlen
    cdef int64_t tlen_squared
    
    cdef int64_t tlen_p_sum = 0
    cdef double tlen_p_mean = 0
    cdef double tlen_p_dev_squared
    cdef double tlen_p_dev_squared_sum = 0
    cdef int64_t tlen_p_squared_sum = 0
    cdef int64_t tlen_p_fwd_sum = 0
    cdef double tlen_p_fwd_mean = 0
    cdef double tlen_p_fwd_dev_squared
    cdef double tlen_p_fwd_dev_squared_sum = 0
    cdef int64_t tlen_p_fwd_squared_sum = 0
    cdef int64_t tlen_p_rev_sum = 0
    cdef double tlen_p_rev_mean = 0
    cdef double tlen_p_rev_dev_squared
    cdef double tlen_p_rev_dev_squared_sum = 0
    cdef int64_t tlen_p_rev_squared_sum = 0

    cdef int64_t tlen_pp_sum = 0
    cdef double tlen_pp_mean = 0
    cdef double tlen_pp_dev_squared
    cdef double tlen_pp_dev_squared_sum = 0
    cdef int64_t tlen_pp_squared_sum = 0
    cdef int64_t tlen_pp_fwd_sum = 0
    cdef double tlen_pp_fwd_mean = 0
    cdef double tlen_pp_fwd_dev_squared
    cdef double tlen_pp_fwd_dev_squared_sum = 0
    cdef int64_t tlen_pp_fwd_squared_sum = 0
    cdef int64_t tlen_pp_rev_sum = 0
    cdef double tlen_pp_rev_mean = 0
    cdef double tlen_pp_rev_dev_squared
    cdef double tlen_pp_rev_dev_squared_sum = 0
    cdef int64_t tlen_pp_rev_squared_sum = 0
    
    # initialise variables
    n = col.n
    plp = col.plp
    
    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # loop over reads
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag

        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)
        mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
        mate_other_chr = <bint>(aln.core.tid != aln.core.mtid)
        
        # not sure these are really needed
        if is_reverse:
            reads_rev += 1
        else:
            reads_fwd += 1

        # N.B., pysam exposes this property as 'tlen' rather than 'isize' so we 
        # follow their naming convention
        tlen = aln.core.isize 
        tlen_squared = tlen**2

        # N.B. insert size is only meaningful if mate is mapped to same chromosome
        if not mate_is_unmapped and not mate_other_chr:
            reads_p += 1
            tlen_p_sum += tlen
            tlen_p_squared_sum += tlen_squared
            if is_reverse:
                reads_p_rev += 1
                tlen_p_rev_sum += tlen
                tlen_p_rev_squared_sum += tlen_squared
            else:
                reads_p_fwd += 1
                tlen_p_fwd_sum += tlen
                tlen_p_fwd_squared_sum += tlen_squared
                
            if is_proper_pair:
                reads_pp += 1
                tlen_pp_sum += tlen
                tlen_pp_squared_sum += tlen_squared
                if is_reverse:
                    reads_pp_rev += 1
                    tlen_pp_rev_sum += tlen
                    tlen_pp_rev_squared_sum += tlen_squared
                else:
                    reads_pp_fwd += 1
                    tlen_pp_fwd_sum += tlen
                    tlen_pp_fwd_squared_sum += tlen_squared

    # calculate intermediate variables
    if reads_p > 0:
        tlen_p_mean = tlen_p_sum * 1. / reads_p
        if reads_p_rev > 0:
            tlen_p_rev_mean = tlen_p_rev_sum * 1. / reads_p_rev
        if reads_p_fwd > 0:
            tlen_p_fwd_mean = tlen_p_fwd_sum * 1. / reads_p_fwd
    if reads_pp > 0:
        tlen_pp_mean = tlen_pp_sum * 1. / reads_pp
        if reads_pp_rev > 0:
            tlen_pp_rev_mean = tlen_pp_rev_sum * 1. / reads_pp_rev
        if reads_pp_fwd > 0:
            tlen_pp_fwd_mean = tlen_pp_fwd_sum * 1. / reads_pp_fwd
        
    # loop over reads again to calculate variance (and hence std)
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)
        mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
        mate_other_chr = <bint>(aln.core.tid != aln.core.mtid)
        tlen = aln.core.isize
        # N.B. insert size is only meaningful if mate is mapped to same chromosome
        if not mate_is_unmapped and not mate_other_chr:
            tlen_p_dev_squared = (tlen - tlen_p_mean)**2
            tlen_p_dev_squared_sum += tlen_p_dev_squared
            if is_reverse:
                tlen_p_rev_dev_squared = (tlen - tlen_p_rev_mean)**2
                tlen_p_rev_dev_squared_sum += tlen_p_rev_dev_squared
            else:
                tlen_p_fwd_dev_squared = (tlen - tlen_p_fwd_mean)**2
                tlen_p_fwd_dev_squared_sum += tlen_p_fwd_dev_squared
            if is_proper_pair:
                tlen_pp_dev_squared = (tlen - tlen_pp_mean)**2
                tlen_pp_dev_squared_sum += tlen_pp_dev_squared
                if is_reverse:
                    tlen_pp_rev_dev_squared = (tlen - tlen_pp_rev_mean)**2
                    tlen_pp_rev_dev_squared_sum += tlen_pp_rev_dev_squared
                else:
                    tlen_pp_fwd_dev_squared = (tlen - tlen_pp_fwd_mean)**2
                    tlen_pp_fwd_dev_squared_sum += tlen_pp_fwd_dev_squared
                    
    # calculate output variables
    # N.B. round values to nearest integer, any finer precision is probably not
    # interesting    
    if reads_p > 0:
        mean_tlen = int(round(tlen_p_mean))
        rms_tlen = int(round(sqrt(tlen_p_squared_sum * 1. / reads_p)))
        variance_tlen = tlen_p_dev_squared_sum * 1. / reads_p
        std_tlen = int(round(sqrt(variance_tlen)))
    else:
        rms_tlen = std_tlen = mean_tlen = 'NA'
    if reads_p_rev > 0:
        mean_tlen_rev = int(round(tlen_p_rev_mean))
        rms_tlen_rev = int(round(sqrt(tlen_p_rev_squared_sum * 1. / reads_p_rev)))
        variance_tlen_rev = tlen_p_rev_dev_squared_sum * 1. / reads_p_rev
        std_tlen_rev = int(round(sqrt(variance_tlen_rev)))
    else:
        rms_tlen_rev = std_tlen_rev = mean_tlen_rev = 'NA'
    if reads_p_fwd > 0:
        mean_tlen_fwd = int(round(tlen_p_fwd_mean))
        rms_tlen_fwd = int(round(sqrt(tlen_p_fwd_squared_sum * 1. / reads_p_fwd)))
        variance_tlen_fwd = tlen_p_fwd_dev_squared_sum * 1. / reads_p_fwd
        std_tlen_fwd = int(round(sqrt(variance_tlen_fwd)))
    else:
        rms_tlen_fwd = std_tlen_fwd = mean_tlen_fwd = 'NA'
    if reads_pp > 0:
        mean_tlen_pp = int(round(tlen_pp_mean))
        rms_tlen_pp = int(round(sqrt(tlen_pp_squared_sum * 1. / reads_pp)))
        variance_tlen_pp = tlen_pp_dev_squared_sum * 1. / reads_pp
        std_tlen_pp = int(round(sqrt(variance_tlen_pp)))
    else:
        rms_tlen_pp = std_tlen_pp = mean_tlen_pp = 'NA'
    if reads_pp_rev > 0:
        mean_tlen_pp_rev = int(round(tlen_pp_rev_mean))
        rms_tlen_pp_rev = int(round(sqrt(tlen_pp_rev_squared_sum * 1. / reads_pp_rev)))
        variance_tlen_pp_rev = tlen_pp_rev_dev_squared_sum * 1. / reads_pp_rev
        std_tlen_pp_rev = int(round(sqrt(variance_tlen_pp_rev)))
    else:
        rms_tlen_pp_rev = std_tlen_pp_rev = mean_tlen_pp_rev = 'NA'
    if reads_pp_fwd > 0:
        mean_tlen_pp_fwd = int(round(tlen_pp_fwd_mean))
        rms_tlen_pp_fwd = int(round(sqrt(tlen_pp_fwd_squared_sum * 1. / reads_pp_fwd)))
        variance_tlen_pp_fwd = tlen_pp_fwd_dev_squared_sum * 1. / reads_pp_fwd
        std_tlen_pp_fwd = int(round(sqrt(variance_tlen_pp_fwd)))
    else:
        rms_tlen_pp_fwd = std_tlen_pp_fwd = mean_tlen_pp_fwd = 'NA'

    return {'chr': chrom, 
            'pos': pos, 
#            'reads_all': n, 
#            'reads_fwd': reads_fwd, 
#            'reads_rev': reads_rev,
#            'reads_paired': reads_p,
#            'reads_paired_fwd': reads_p_fwd,
#            'reads_paired_rev': reads_p_rev,
#            'reads_pp': reads_pp,
#            'reads_pp_fwd': reads_pp_fwd,
#            'reads_pp_rev': reads_pp_rev,
            'mean_tlen': mean_tlen,
            'mean_tlen_fwd': mean_tlen_fwd,
            'mean_tlen_rev': mean_tlen_rev,
            'mean_tlen_pp': mean_tlen_pp,
            'mean_tlen_pp_fwd': mean_tlen_pp_fwd,
            'mean_tlen_pp_rev': mean_tlen_pp_rev,
            'rms_tlen': rms_tlen,
            'rms_tlen_fwd': rms_tlen_fwd,
            'rms_tlen_rev': rms_tlen_rev,
            'rms_tlen_pp': rms_tlen_pp,
            'rms_tlen_pp_fwd': rms_tlen_pp_fwd,
            'rms_tlen_pp_rev': rms_tlen_pp_rev,
            'std_tlen': std_tlen,
            'std_tlen_fwd': std_tlen_fwd,
            'std_tlen_rev': std_tlen_rev,
            'std_tlen_pp': std_tlen_pp,
            'std_tlen_pp_fwd': std_tlen_pp_fwd,
            'std_tlen_pp_rev': std_tlen_pp_rev}


def stat_tlen_strand(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_tlen_strand(samfile, col, one_based)
        
        
def write_tlen_strand(*args, **kwargs):
    fieldnames = ('chr', 'pos', 
#                  'reads_all', 'reads_fwd', 'reads_rev', 
#                  'reads_paired', 'reads_paired_fwd', 'reads_paired_rev', 
#                  'reads_pp', 'reads_pp_fwd', 'reads_pp_rev', 
                  'mean_tlen', 'mean_tlen_fwd', 'mean_tlen_rev', 
                  'mean_tlen_pp', 'mean_tlen_pp_fwd', 'mean_tlen_pp_rev',
                  'rms_tlen', 'rms_tlen_fwd', 'rms_tlen_rev', 
                  'rms_tlen_pp', 'rms_tlen_pp_fwd', 'rms_tlen_pp_rev',
                  'std_tlen', 'std_tlen_fwd', 'std_tlen_rev', 
                  'std_tlen_pp', 'std_tlen_pp_fwd', 'std_tlen_pp_rev')
    write_stats(stat_tlen_strand, fieldnames, *args, **kwargs)
    
    
##############################
# MAPPING QUALITY STATISTICS #
##############################


cpdef object construct_rec_mapq(Samfile samfile, PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef uint32_t flag
    cdef uint64_t mapq
    cdef uint64_t mapq_max = 0
    cdef uint64_t mapq_pp_max = 0
    cdef uint64_t mapq_squared
    cdef uint64_t mapq_squared_sum = 0
    cdef uint64_t mapq_pp_squared_sum = 0
    cdef bint is_proper_pair
    cdef unsigned int reads_pp = 0
    cdef unsigned int reads_mapq0 = 0
    cdef unsigned int reads_mapq0_pp = 0

    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        mapq = aln.core.qual
        mapq_squared = mapq**2
        mapq_squared_sum += mapq_squared
        if mapq == 0:
            reads_mapq0 += 1
        if mapq > mapq_max:
            mapq_max = mapq
        if is_proper_pair:
            reads_pp += 1
            mapq_pp_squared_sum += mapq_squared
            if mapq > mapq_pp_max:
                mapq_pp_max = mapq
            if mapq == 0:
                reads_mapq0_pp += 1

    # construct output variables
    rms_mapq = int(round(sqrt(mapq_squared_sum * 1. / n)))
    max_mapq = mapq_max
    if reads_pp > 0:
        rms_mapq_pp = int(round(sqrt(mapq_pp_squared_sum * 1. / reads_pp)))
        max_mapq_pp = mapq_pp_max
    else:
        rms_mapq_pp = max_mapq_pp = 'NA'
        
    return {'chr': chrom, 
            'pos': pos, 
            'reads_all': n, 
            'reads_pp': reads_pp,
            'reads_mapq0': reads_mapq0,
            'reads_mapq0_pp': reads_mapq0_pp,
            'rms_mapq': rms_mapq,
            'rms_mapq_pp': rms_mapq_pp,
            'max_mapq': max_mapq,
            'max_mapq_pp': max_mapq_pp}


def stat_mapq(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_mapq(samfile, col, one_based)
        
        
def write_mapq(*args, **kwargs):
    fieldnames = ('chr', 'pos', 
                  'reads_all', 'reads_pp',
                  'reads_mapq0', 'reads_mapq0_pp',
                  'rms_mapq', 'rms_mapq_pp',
                  'max_mapq', 'max_mapq_pp')
    write_stats(stat_mapq, fieldnames, *args, **kwargs)
    
    
########################################
# MAPPING QUALITY STATISTICS BY STRAND #
########################################


cpdef object construct_rec_mapq_strand(Samfile samfile, PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column

    cdef uint32_t flag
    cdef uint64_t mapq
    cdef uint64_t mapq_squared
    cdef bint is_proper_pair
    cdef bint is_reverse

    cdef uint64_t mapq_max = 0
    cdef uint64_t mapq_rev_max = 0
    cdef uint64_t mapq_fwd_max = 0
    cdef uint64_t mapq_pp_max = 0
    cdef uint64_t mapq_pp_rev_max = 0
    cdef uint64_t mapq_pp_fwd_max = 0
    cdef uint64_t mapq_squared_sum = 0
    cdef uint64_t mapq_fwd_squared_sum = 0
    cdef uint64_t mapq_rev_squared_sum = 0
    cdef uint64_t mapq_pp_squared_sum = 0
    cdef uint64_t mapq_pp_fwd_squared_sum = 0
    cdef uint64_t mapq_pp_rev_squared_sum = 0

    cdef unsigned int reads_rev = 0
    cdef unsigned int reads_fwd = 0
    cdef unsigned int reads_pp = 0
    cdef unsigned int reads_pp_rev = 0
    cdef unsigned int reads_pp_fwd = 0
    cdef unsigned int reads_mapq0 = 0
    cdef unsigned int reads_mapq0_fwd = 0
    cdef unsigned int reads_mapq0_rev = 0
    cdef unsigned int reads_mapq0_pp = 0
    cdef unsigned int reads_mapq0_pp_fwd = 0
    cdef unsigned int reads_mapq0_pp_rev = 0

    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)
        mapq = aln.core.qual
        mapq_squared = mapq**2

        mapq_squared_sum += mapq_squared
        if mapq > mapq_max:
            mapq_max = mapq
        if is_reverse:
            reads_rev += 1
            mapq_rev_squared_sum += mapq_squared
            if mapq > mapq_rev_max:
                mapq_rev_max = mapq
        else:
            reads_fwd += 1
            mapq_fwd_squared_sum += mapq_squared
            if mapq > mapq_fwd_max:
                mapq_fwd_max = mapq

        if mapq == 0:
            reads_mapq0 += 1
            if is_reverse:
                reads_mapq0_rev += 1
            else:
                reads_mapq0_fwd += 1

        if is_proper_pair:
            reads_pp += 1
            mapq_pp_squared_sum += mapq_squared
            if mapq > mapq_pp_max:
                mapq_pp_max = mapq
            if is_reverse:
                reads_pp_rev += 1
                mapq_pp_rev_squared_sum += mapq_squared
                if mapq > mapq_pp_rev_max:
                    mapq_pp_rev_max = mapq
            else:
                reads_pp_fwd += 1
                mapq_pp_fwd_squared_sum += mapq_squared
                if mapq > mapq_pp_fwd_max:
                    mapq_pp_fwd_max = mapq
            if mapq == 0:
                reads_mapq0_pp += 1
                if is_reverse:
                    reads_mapq0_pp_rev += 1
                else:
                    reads_mapq0_pp_fwd += 1

    # construct output variables
    rms_mapq = int(round(sqrt(mapq_squared_sum * 1. / n)))
    max_mapq = mapq_max
    if reads_rev > 0:
        rms_mapq_rev = int(round(sqrt(mapq_rev_squared_sum * 1. / reads_rev)))
        max_mapq_rev = mapq_rev_max
    else:
        rms_mapq_rev = max_mapq_rev = 'NA'
    if reads_fwd > 0:
        rms_mapq_fwd = int(round(sqrt(mapq_fwd_squared_sum * 1. / reads_fwd)))
        max_mapq_fwd = mapq_fwd_max
    else:
        rms_mapq_fwd = max_mapq_fwd = 'NA'
    if reads_pp > 0:
        rms_mapq_pp = int(round(sqrt(mapq_pp_squared_sum * 1. / reads_pp)))
        max_mapq_pp = mapq_pp_max
    else:
        rms_mapq_pp = max_mapq_pp = 'NA'
    if reads_pp_fwd > 0:
        rms_mapq_pp_fwd = int(round(sqrt(mapq_pp_fwd_squared_sum * 1. / reads_pp_fwd)))
        max_mapq_pp_fwd = mapq_pp_fwd_max
    else:
        rms_mapq_pp_fwd = max_mapq_pp_fwd = 'NA'
    if reads_pp_rev > 0:
        rms_mapq_pp_rev = int(round(sqrt(mapq_pp_rev_squared_sum * 1. / reads_pp_rev)))
        max_mapq_pp_rev = mapq_pp_rev_max
    else:
        rms_mapq_pp_rev = max_mapq_pp_rev = 'NA'
        
    return {'chr': chrom, 
            'pos': pos, 
            'reads_all': n,
            'reads_fwd': reads_fwd, 
            'reads_rev': reads_rev, 
            'reads_pp': reads_pp,
            'reads_pp_fwd': reads_pp_fwd,
            'reads_pp_rev': reads_pp_rev,
            'reads_mapq0': reads_mapq0,
            'reads_mapq0_fwd': reads_mapq0_fwd, 
            'reads_mapq0_rev': reads_mapq0_rev, 
            'reads_mapq0_pp': reads_mapq0_pp,
            'reads_mapq0_pp_fwd': reads_mapq0_pp_fwd,
            'reads_mapq0_pp_rev': reads_mapq0_pp_rev,
            'rms_mapq': rms_mapq,
            'rms_mapq_fwd': rms_mapq_fwd,
            'rms_mapq_rev': rms_mapq_rev,
            'rms_mapq_pp': rms_mapq_pp,
            'rms_mapq_pp_fwd': rms_mapq_pp_fwd,
            'rms_mapq_pp_rev': rms_mapq_pp_rev,
            'max_mapq': max_mapq,
            'max_mapq_fwd': max_mapq_fwd,
            'max_mapq_rev': max_mapq_rev,
            'max_mapq_pp': max_mapq_pp,
            'max_mapq_pp_fwd': max_mapq_pp_fwd,
            'max_mapq_pp_rev': max_mapq_pp_rev,
            }


def stat_mapq_strand(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_mapq_strand(samfile, col, one_based)
        
        
def write_mapq_strand(*args, **kwargs):
    fieldnames = ('chr', 'pos', 
                  'reads_all', 'reads_fwd', 'reads_rev', 
                  'reads_pp', 'reads_pp_fwd', 'reads_pp_rev', 
                  'reads_mapq0', 'reads_mapq0_fwd', 'reads_mapq0_rev', 
                  'reads_mapq0_pp', 'reads_mapq0_pp_fwd', 'reads_mapq0_pp_rev', 
                  'rms_mapq', 'rms_mapq_fwd', 'rms_mapq_rev', 
                  'rms_mapq_pp', 'rms_mapq_pp_fwd', 'rms_mapq_pp_rev', 
                  'max_mapq', 'max_mapq_fwd', 'max_mapq_rev', 
                  'max_mapq_pp', 'max_mapq_pp_fwd', 'max_mapq_pp_rev', 
                  )
    write_stats(stat_mapq_strand, fieldnames, *args, **kwargs)
    
    
###########################
# BASE QUALITY STATISTICS #
###########################


cpdef object construct_rec_baseq(Samfile samfile, PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef unsigned int reads_nodel = 0
#    cdef unsigned int reads_pp = 0
    cdef unsigned int reads_pp_nodel = 0
    cdef uint64_t baseq, baseq_squared
    cdef uint64_t baseq_squared_sum = 0
    cdef uint64_t baseq_pp_squared_sum = 0

    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
 #       if is_proper_pair:
 #           reads_pp += 1
        # N.B., base quality only makes sense if the aligned read is not a deletion
        if not read.is_del:
            reads_nodel += 1
            baseq = bam1_qual(aln)[read.qpos]
            baseq_squared = baseq**2
            baseq_squared_sum += baseq_squared
            if is_proper_pair:
                reads_pp_nodel += 1
                baseq_pp_squared_sum += baseq_squared

    # output variables
    rms_baseq = rootmean(baseq_squared_sum, reads_nodel)
    rms_baseq_pp = rootmean(baseq_pp_squared_sum, reads_pp_nodel)

    return {'chr': chrom, 
            'pos': pos, 
 #           'reads_all': n, 
 #           'reads_pp': reads_pp,
            'rms_baseq': rms_baseq,
            'rms_baseq_pp': rms_baseq_pp}


def stat_baseq(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_baseq(samfile, col, one_based)
        
        
def write_baseq(*args, **kwargs):
    fieldnames = ('chr', 'pos', 
#                  'reads_all', 
#                  'reads_pp', 
                  'rms_baseq', 
                  'rms_baseq_pp',
                  )
    write_stats(stat_baseq, fieldnames, *args, **kwargs)
    
    
#####################################
# BASE QUALITY STATISTICS BY STRAND #
#####################################


cpdef object construct_rec_baseq_strand(Samfile samfile, PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column

    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef bint is_reverse
    cdef uint64_t baseq
    cdef uint64_t baseq_squared

    cdef uint64_t baseq_squared_sum = 0
    cdef uint64_t baseq_fwd_squared_sum = 0
    cdef uint64_t baseq_rev_squared_sum = 0
    cdef uint64_t baseq_pp_squared_sum = 0
    cdef uint64_t baseq_pp_fwd_squared_sum = 0
    cdef uint64_t baseq_pp_rev_squared_sum = 0

#    cdef unsigned int reads_rev = 0
#    cdef unsigned int reads_fwd = 0
#    cdef unsigned int reads_pp = 0
#    cdef unsigned int reads_pp_rev = 0
#    cdef unsigned int reads_pp_fwd = 0
    cdef unsigned int reads_nodel = 0
    cdef unsigned int reads_rev_nodel = 0
    cdef unsigned int reads_fwd_nodel = 0
    cdef unsigned int reads_pp_nodel = 0
    cdef unsigned int reads_pp_rev_nodel = 0
    cdef unsigned int reads_pp_fwd_nodel = 0

    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)

#        if is_reverse:
#            reads_rev += 1
#        else:
#            reads_fwd += 1
#        if is_proper_pair:
#            reads_pp += 1
#            if is_reverse:
#                reads_pp_rev += 1
#            else:
#                reads_pp_fwd += 1

        # N.B., baseq only makes sense if not a deletion
        if not read.is_del:
            reads_nodel += 1
            baseq = bam1_qual(aln)[read.qpos]
            baseq_squared = baseq**2
            baseq_squared_sum += baseq_squared
            if is_reverse:
                reads_rev_nodel += 1
                baseq_rev_squared_sum += baseq_squared
            else:
                reads_fwd_nodel += 1
                baseq_fwd_squared_sum += baseq_squared
            if is_proper_pair:
                reads_pp_nodel += 1
                baseq_pp_squared_sum += baseq_squared
                if is_reverse:
                    reads_pp_rev_nodel += 1
                    baseq_pp_rev_squared_sum += baseq_squared
                else:
                    reads_pp_fwd_nodel += 1
                    baseq_pp_fwd_squared_sum += baseq_squared

    # construct output variables
    rms_baseq = rootmean(baseq_squared_sum, reads_nodel)
    rms_baseq_rev = rootmean(baseq_rev_squared_sum, reads_rev_nodel)
    rms_baseq_fwd = rootmean(baseq_fwd_squared_sum, reads_fwd_nodel)
    rms_baseq_pp = rootmean(baseq_pp_squared_sum, reads_pp_nodel)
    rms_baseq_pp_fwd = rootmean(baseq_pp_fwd_squared_sum, reads_pp_fwd_nodel)
    rms_baseq_pp_rev = rootmean(baseq_pp_rev_squared_sum, reads_pp_rev_nodel)
        
    return {'chr': chrom, 
            'pos': pos, 
#            'reads_all': n,
#            'reads_fwd': reads_fwd, 
#            'reads_rev': reads_rev, 
#            'reads_pp': reads_pp,
#            'reads_pp_fwd': reads_pp_fwd,
#            'reads_pp_rev': reads_pp_rev,
            'rms_baseq': rms_baseq,
            'rms_baseq_fwd': rms_baseq_fwd,
            'rms_baseq_rev': rms_baseq_rev,
            'rms_baseq_pp': rms_baseq_pp,
            'rms_baseq_pp_fwd': rms_baseq_pp_fwd,
            'rms_baseq_pp_rev': rms_baseq_pp_rev,
            }


def stat_baseq_strand(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_baseq_strand(samfile, col, one_based)
        
        
def write_baseq_strand(*args, **kwargs):
    fieldnames = ('chr', 'pos', 
#                  'reads_all', 'reads_fwd', 'reads_rev', 
#                  'reads_pp', 'reads_pp_fwd', 'reads_pp_rev', 
                  'rms_baseq', 'rms_baseq_fwd', 'rms_baseq_rev', 
                  'rms_baseq_pp', 'rms_baseq_pp_fwd', 'rms_baseq_pp_rev', 
                  )
    write_stats(stat_baseq_strand, fieldnames, *args, **kwargs)
    
    
####################################
# EXTENDED BASE QUALITY STATISTICS #
####################################


cpdef object construct_rec_baseq_ext(Samfile samfile, Fastafile fafile, 
                                     PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    # counting variables
    cdef unsigned int reads_nodel = 0
#    cdef unsigned int reads_pp = 0
    cdef unsigned int reads_pp_nodel = 0
    cdef unsigned int matches = 0
    cdef unsigned int matches_pp = 0
    cdef unsigned int mismatches = 0
    cdef unsigned int mismatches_pp = 0

    cdef uint64_t baseq
    cdef uint64_t baseq_squared

    cdef uint64_t baseq_squared_sum = 0
    cdef uint64_t baseq_pp_squared_sum = 0
    cdef uint64_t baseq_matches_squared_sum = 0
    cdef uint64_t baseq_matches_pp_squared_sum = 0
    cdef uint64_t baseq_mismatches_squared_sum = 0
    cdef uint64_t baseq_mismatches_pp_squared_sum = 0
    
    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # reference base
    refbase = fafile.fetch(reference=chrom, start=col.pos, end=col.pos+1).upper()
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#        if is_proper_pair:
#            reads_pp += 1
        if not read.is_del:
            reads_nodel += 1
            baseq = bam1_qual(aln)[read.qpos]
            baseq_squared = baseq**2
            baseq_squared_sum += baseq_squared
            if is_proper_pair:
                reads_pp_nodel += 1
                baseq_pp_squared_sum += baseq_squared
            alnbase = get_seq_base(aln, read.qpos)
            if alnbase == refbase:
                matches += 1
                baseq_matches_squared_sum += baseq_squared
                if is_proper_pair:
                    matches_pp += 1
                    baseq_matches_pp_squared_sum += baseq_squared
            else:
                mismatches += 1
                baseq_mismatches_squared_sum += baseq_squared
                if is_proper_pair:
                    mismatches_pp += 1
                    baseq_mismatches_pp_squared_sum += baseq_squared

    # construct output variables
    rms_baseq = rootmean(baseq_squared_sum, reads_nodel)
    rms_baseq_pp = rootmean(baseq_pp_squared_sum, reads_pp_nodel)
    rms_baseq_matches = rootmean(baseq_matches_squared_sum, matches)
    rms_baseq_matches_pp = rootmean(baseq_matches_pp_squared_sum, matches_pp)
    rms_baseq_mismatches = rootmean(baseq_mismatches_squared_sum, mismatches)
    rms_baseq_mismatches_pp = rootmean(baseq_mismatches_pp_squared_sum, mismatches_pp)

    return {'chr': chrom, 'pos': pos, 'ref': refbase,
#            'reads_all': n, 'reads_pp': reads_pp,
#            'matches': matches,
#            'matches_pp': matches_pp,
#            'mismatches': mismatches,
#            'mismatches_pp': mismatches_pp,
            'rms_baseq': rms_baseq,
            'rms_baseq_pp': rms_baseq_pp,
            'rms_baseq_matches': rms_baseq_matches,
            'rms_baseq_matches_pp': rms_baseq_matches_pp,
            'rms_baseq_mismatches': rms_baseq_mismatches,
            'rms_baseq_mismatches_pp': rms_baseq_mismatches_pp,
            }


def stat_baseq_ext(samfile, fafile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_baseq_ext(samfile, fafile, col, one_based)
        
        
def write_baseq_ext(*args, **kwargs):
    fieldnames = ('chr', 'pos', 'ref', 
#                  'reads_all', 'reads_pp',
#                  'matches', 'matches_pp',
#                  'mismatches', 'mismatches_pp',
                  'rms_baseq', 'rms_baseq_pp',
                  'rms_baseq_matches', 'rms_baseq_matches_pp',
                  'rms_baseq_mismatches', 'rms_baseq_mismatches_pp',
                  )
    write_stats(stat_baseq_ext, fieldnames, *args, **kwargs)
    
    
##############################################
# EXTENDED BASE QUALITY STATISTICS BY STRAND #
##############################################


cpdef object construct_rec_baseq_ext_strand(Samfile samfile, Fastafile fafile, 
                                            PileupProxy col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int n # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef bint is_reverse
    
    # counting variables
#    cdef unsigned int reads_fwd = 0
#    cdef unsigned int reads_rev = 0
    cdef unsigned int reads_nodel = 0
    cdef unsigned int reads_fwd_nodel = 0
    cdef unsigned int reads_rev_nodel = 0
#    cdef unsigned int reads_pp = 0
#    cdef unsigned int reads_pp_fwd = 0
#    cdef unsigned int reads_pp_rev = 0
    cdef unsigned int reads_pp_nodel = 0
    cdef unsigned int reads_pp_fwd_nodel = 0
    cdef unsigned int reads_pp_rev_nodel = 0
    cdef unsigned int matches = 0
    cdef unsigned int matches_fwd = 0
    cdef unsigned int matches_rev = 0
    cdef unsigned int matches_pp = 0
    cdef unsigned int matches_pp_fwd = 0
    cdef unsigned int matches_pp_rev = 0
    cdef unsigned int mismatches = 0
    cdef unsigned int mismatches_fwd = 0
    cdef unsigned int mismatches_rev = 0
    cdef unsigned int mismatches_pp = 0
    cdef unsigned int mismatches_pp_fwd = 0
    cdef unsigned int mismatches_pp_rev = 0

    cdef uint64_t baseq
    cdef uint64_t baseq_squared

    cdef uint64_t baseq_squared_sum = 0
    cdef uint64_t baseq_fwd_squared_sum = 0
    cdef uint64_t baseq_rev_squared_sum = 0
    cdef uint64_t baseq_pp_squared_sum = 0
    cdef uint64_t baseq_pp_fwd_squared_sum = 0
    cdef uint64_t baseq_pp_rev_squared_sum = 0
    cdef uint64_t baseq_matches_squared_sum = 0
    cdef uint64_t baseq_matches_fwd_squared_sum = 0
    cdef uint64_t baseq_matches_rev_squared_sum = 0
    cdef uint64_t baseq_matches_pp_squared_sum = 0
    cdef uint64_t baseq_matches_pp_fwd_squared_sum = 0
    cdef uint64_t baseq_matches_pp_rev_squared_sum = 0
    cdef uint64_t baseq_mismatches_squared_sum = 0
    cdef uint64_t baseq_mismatches_fwd_squared_sum = 0
    cdef uint64_t baseq_mismatches_rev_squared_sum = 0
    cdef uint64_t baseq_mismatches_pp_squared_sum = 0
    cdef uint64_t baseq_mismatches_pp_fwd_squared_sum = 0
    cdef uint64_t baseq_mismatches_pp_rev_squared_sum = 0
    
    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    chrom = samfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # reference base
    refbase = fafile.fetch(reference=chrom, start=col.pos, end=col.pos+1).upper()
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)

#        if is_reverse:
#            reads_rev += 1
#        else:
#            reads_fwd += 1
#        if is_proper_pair:
#            reads_pp += 1
#            if is_reverse:
#                reads_pp_rev += 1
#            else:
#                reads_pp_fwd += 1

        if not read.is_del:
            reads_nodel += 1
            baseq = bam1_qual(aln)[read.qpos]
            baseq_squared = baseq**2
            baseq_squared_sum += baseq_squared
            if is_reverse:
                reads_rev_nodel += 1
                baseq_rev_squared_sum += baseq_squared
            else:
                reads_fwd_nodel += 1
                baseq_fwd_squared_sum += baseq_squared
            if is_proper_pair:
                reads_pp_nodel += 1
                baseq_pp_squared_sum += baseq_squared
                if is_reverse:
                    reads_pp_rev_nodel += 1
                    baseq_pp_rev_squared_sum += baseq_squared
                else:
                    reads_pp_fwd_nodel += 1
                    baseq_pp_fwd_squared_sum += baseq_squared
            alnbase = get_seq_base(aln, read.qpos)
            if alnbase == refbase:
                matches += 1
                baseq_matches_squared_sum += baseq_squared
                if is_reverse:
                    matches_rev += 1
                    baseq_matches_rev_squared_sum += baseq_squared
                else:
                    matches_fwd += 1
                    baseq_matches_fwd_squared_sum += baseq_squared
                    
                if is_proper_pair:
                    matches_pp += 1
                    baseq_matches_pp_squared_sum += baseq_squared
                    if is_reverse:
                        matches_pp_rev += 1
                        baseq_matches_pp_rev_squared_sum += baseq_squared
                    else:
                        matches_pp_fwd += 1
                        baseq_matches_pp_fwd_squared_sum += baseq_squared
            else:
                mismatches += 1
                baseq_mismatches_squared_sum += baseq_squared
                if is_reverse:
                    mismatches_rev += 1
                    baseq_mismatches_rev_squared_sum += baseq_squared
                else:
                    mismatches_fwd += 1
                    baseq_mismatches_fwd_squared_sum += baseq_squared
                    
                if is_proper_pair:
                    mismatches_pp += 1
                    baseq_mismatches_pp_squared_sum += baseq_squared
                    if is_reverse:
                        mismatches_pp_rev += 1
                        baseq_mismatches_pp_rev_squared_sum += baseq_squared
                    else:
                        mismatches_pp_fwd += 1
                        baseq_mismatches_pp_fwd_squared_sum += baseq_squared

    # construct output variables
    rms_baseq = rootmean(baseq_squared_sum, reads_nodel)
    rms_baseq_fwd = rootmean(baseq_fwd_squared_sum, reads_fwd_nodel)
    rms_baseq_rev = rootmean(baseq_rev_squared_sum, reads_rev_nodel)
    rms_baseq_pp = rootmean(baseq_pp_squared_sum, reads_pp_nodel)
    rms_baseq_pp_fwd = rootmean(baseq_pp_fwd_squared_sum, reads_pp_fwd_nodel)
    rms_baseq_pp_rev = rootmean(baseq_pp_rev_squared_sum, reads_pp_rev_nodel)
    rms_baseq_matches = rootmean(baseq_matches_squared_sum, matches)
    rms_baseq_matches_fwd = rootmean(baseq_matches_fwd_squared_sum, matches_fwd)
    rms_baseq_matches_rev = rootmean(baseq_matches_rev_squared_sum, matches_rev)
    rms_baseq_matches_pp = rootmean(baseq_matches_pp_squared_sum, matches_pp)
    rms_baseq_matches_pp_fwd = rootmean(baseq_matches_pp_fwd_squared_sum, matches_pp_fwd)
    rms_baseq_matches_pp_rev = rootmean(baseq_matches_pp_rev_squared_sum, matches_pp_rev)
    rms_baseq_mismatches = rootmean(baseq_mismatches_squared_sum, mismatches)
    rms_baseq_mismatches_fwd = rootmean(baseq_mismatches_fwd_squared_sum, mismatches_fwd)
    rms_baseq_mismatches_rev = rootmean(baseq_mismatches_rev_squared_sum, mismatches_rev)
    rms_baseq_mismatches_pp = rootmean(baseq_mismatches_pp_squared_sum, mismatches_pp)
    rms_baseq_mismatches_pp_fwd = rootmean(baseq_mismatches_pp_fwd_squared_sum, mismatches_pp_fwd)
    rms_baseq_mismatches_pp_rev = rootmean(baseq_mismatches_pp_rev_squared_sum, mismatches_pp_rev)

    return {'chr': chrom, 'pos': pos, 'ref': refbase,
#            'reads_all': n, 'reads_fwd': reads_fwd, 'reads_rev': reads_rev, 
#            'reads_pp': reads_pp, 'reads_pp_fwd': reads_pp_fwd, 'reads_pp_rev': reads_pp_rev,
#            'matches': matches, 'matches_fwd': matches_fwd, 'matches_rev': matches_rev, 
#            'matches_pp': matches_pp, 'matches_pp_fwd': matches_pp_fwd, 'matches_pp_rev': matches_pp_rev,
#            'mismatches': mismatches, 'mismatches_fwd': mismatches_fwd, 'mismatches_rev': mismatches_rev, 
#            'mismatches_pp': mismatches_pp, 'mismatches_pp_fwd': mismatches_pp_fwd, 'mismatches_pp_rev': mismatches_pp_rev,
            'rms_baseq': rms_baseq, 'rms_baseq_fwd': rms_baseq_fwd, 'rms_baseq_rev': rms_baseq_rev, 
            'rms_baseq_pp': rms_baseq_pp, 'rms_baseq_pp_fwd': rms_baseq_pp_fwd, 'rms_baseq_pp_rev': rms_baseq_pp_rev,
            'rms_baseq_matches': rms_baseq_matches, 'rms_baseq_matches_fwd': rms_baseq_matches_fwd, 'rms_baseq_matches_rev': rms_baseq_matches_rev, 
            'rms_baseq_matches_pp': rms_baseq_matches_pp, 'rms_baseq_matches_pp_fwd': rms_baseq_matches_pp_fwd, 'rms_baseq_matches_pp_rev': rms_baseq_matches_pp_rev,
            'rms_baseq_mismatches': rms_baseq_mismatches, 'rms_baseq_mismatches_fwd': rms_baseq_mismatches_fwd, 'rms_baseq_mismatches_rev': rms_baseq_mismatches_rev, 
            'rms_baseq_mismatches_pp': rms_baseq_mismatches_pp, 'rms_baseq_mismatches_pp_fwd': rms_baseq_mismatches_pp_fwd, 'rms_baseq_mismatches_pp_rev': rms_baseq_mismatches_pp_rev,
            }


def stat_baseq_ext_strand(samfile, fafile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_baseq_ext_strand(samfile, fafile, col, one_based)
        
        
def write_baseq_ext_strand(*args, **kwargs):
    fieldnames = ('chr', 'pos', 'ref', 
#                  'reads_all', 'reads_fwd', 'reads_rev', 
#                  'reads_pp', 'reads_pp_fwd', 'reads_pp_rev',
#                  'matches', 'matches_fwd', 'matches_rev', 
#                  'matches_pp', 'matches_pp_fwd', 'matches_pp_rev',
#                  'mismatches', 'mismatches_fwd', 'mismatches_rev', 
#                  'mismatches_pp', 'mismatches_pp_fwd', 'mismatches_pp_rev', 
                  'rms_baseq', 'rms_baseq_fwd', 'rms_baseq_rev', 
                  'rms_baseq_pp', 'rms_baseq_pp_fwd', 'rms_baseq_pp_rev',
                  'rms_baseq_matches', 'rms_baseq_matches_fwd', 'rms_baseq_matches_rev', 
                  'rms_baseq_matches_pp', 'rms_baseq_matches_pp_fwd', 'rms_baseq_matches_pp_rev',
                  'rms_baseq_mismatches', 'rms_baseq_mismatches_fwd', 'rms_baseq_mismatches_rev', 
                  'rms_baseq_mismatches_pp', 'rms_baseq_mismatches_pp_fwd', 'rms_baseq_mismatches_pp_rev',
                  )
    write_stats(stat_baseq_ext_strand, fieldnames, *args, **kwargs)
    
    
# TODO normed coverage
# TODO check tlen & mapq stats and anything else has NA where it should


#####################
# UTILITY FUNCTIONS #
#####################  


def normalise_coords(start, end, one_based):
    if one_based:
        start = start - 1 if start is not None else None
        end = end - 1 if end is not None else None
    return start, end

    
def write_stats(statfun, fieldnames, outfile, samfile, fafile=None,
                dialect=csv.excel_tab, write_header=True, 
                chrom=None, start=None, end=None, 
                one_based=False, progress=None):
    cdef long long counter = 0
    cdef long long modulus
    
    writer = csv.DictWriter(outfile, fieldnames, dialect=dialect)
    
    if write_header:
        writer.writeheader()

    if fafile is None:
        recs = statfun(samfile, chrom=chrom, start=start, end=end, one_based=one_based)
    else:
        recs = statfun(samfile, fafile, chrom=chrom, start=start, end=end, one_based=one_based)

    if progress is None:
        writer.writerows(recs)

    else:
        modulus = progress
        before = time.time()
        before_all = before
        for rec in recs:
            counter += 1
            writer.writerow(rec)
            if counter % modulus == 0:
                after = time.time()
                elapsed = after - before_all
                batch_elapsed = after - before
                print >>sys.stderr, '%s rows in %.2fs (%d rows/s); batch in %.2fs (%d rows/s)' % (counter, elapsed, counter/elapsed, batch_elapsed, progress/batch_elapsed)
                before = after
        after_all = time.time()
        elapsed_all = after_all - before_all
        print >>sys.stderr, '%s rows in %.2fs (%d rows/s)' % (counter, elapsed_all, counter/elapsed_all)
    
    
cdef inline bint is_softclipped(bam1_t * aln):
    cigar_p = bam1_cigar(aln);
    for k in range(aln.core.n_cigar):
        op = cigar_p[k] & BAM_CIGAR_MASK
        if op == BAM_CSOFT_CLIP:
            return 1
    return 0


cdef inline object get_seq_base(bam1_t *src, uint32_t k):
    cdef uint8_t * p
    cdef char * s

    if not src.core.l_qseq:
        return None

    seq = PyBytes_FromStringAndSize(NULL, 1)
    s   = <char*>seq
    p   = bam1_seq(src)

    # equivalent to bam_nt16_rev_table[bam1_seqi(s, i)] (see bam.c)
    # note: do not use string literal as it will be a python string
    s[0] = bam_nt16_rev_table[p[k/2] >> 4 * (1 - k%2) & 0xf]

    return seq


cdef inline object rootmean(uint64_t sqsum, unsigned int count):
    if count > 0:
        return int(round(sqrt(sqsum * 1. / count)))
    else:
        return 'NA'
    
    