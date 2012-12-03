# cython: profile=True


import sys
import numpy as np
cimport numpy as np
import time
import csv
from libc.stdint cimport uint32_t
from csamtools cimport Samfile, PileupProxy, bam1_t, bam_pileup1_t, bam1_cigar


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


def normalise_coords(start, end, one_based):
    if one_based:
        start = start - 1 if start is not None else None
        end = end - 1 if end is not None else None
    return start, end


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
        
        
def write_stats(statfun, outfile, fieldnames, samfile,
                dialect=csv.excel_tab, write_header=True, 
                chrom=None, start=None, end=None, 
                one_based=False, progress=None):
    cdef long long counter = 0
    cdef long long modulus
    
    writer = csv.DictWriter(outfile, fieldnames, dialect=dialect)
    
    if write_header:
        writer.writeheader()
    
    if progress is None:
        recs = statfun(samfile, chrom=chrom, start=start, end=end, 
                       one_based=one_based)
        writer.writerows(recs)

    else:
        modulus = progress
        before = time.time()
        before_all = before
        for rec in statfun(samfile, chrom=chrom, start=start, end=end, 
                           one_based=one_based):
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
    
    
def write_coverage(outfile, samfile, dialect=csv.excel_tab, write_header=True,
                   chrom=None, start=None, end=None, 
                   one_based=False, progress=None):
    fieldnames = ('chr', 'pos', 'reads_all', 'reads_pp')
    write_stats(stat_coverage, outfile, fieldnames, samfile, 
                dialect=dialect, write_header=write_header,
                chrom=chrom, start=start, end=end, 
                one_based=one_based, progress=progress)
    
    
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
        
        
def write_coverage_strand(outfile, samfile, dialect=csv.excel_tab, write_header=True,
                          chrom=None, start=None, end=None, 
                          one_based=False, progress=None):
    fieldnames = ('chr', 'pos', 
                  'reads_all', 'reads_fwd', 'reads_rev', 
                  'reads_pp', 'reads_pp_fwd', 'reads_pp_rev')
    write_stats(stat_coverage_strand, outfile, fieldnames, samfile, 
                dialect=dialect, write_header=write_header,
                chrom=chrom, start=start, end=end, 
                one_based=one_based, progress=progress)
    
    
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


cdef bint is_softclipped(bam1_t * aln):
    cigar_p = bam1_cigar(aln);
    for k in range(aln.core.n_cigar):
        op = cigar_p[k] & BAM_CIGAR_MASK
        if op == BAM_CSOFT_CLIP:
            return 1
    return 0


def stat_coverage_ext(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(chrom, start, end):
        yield construct_rec_coverage_ext(samfile, col, one_based)
        
        
def write_coverage_ext(outfile, samfile, dialect=csv.excel_tab, write_header=True,
                       chrom=None, start=None, end=None, 
                       one_based=False, progress=None):
    fieldnames = ('chr', 'pos', 
                  'reads_all', 'reads_pp', 
                  'reads_mate_unmapped', 
                  'reads_mate_other_chr',
                  'reads_mate_same_strand',
                  'reads_faceaway', 
                  'reads_softclipped')
    write_stats(stat_coverage_ext, outfile, fieldnames, samfile, 
                dialect=dialect, write_header=write_header,
                chrom=chrom, start=start, end=end, 
                one_based=one_based, progress=progress)
    
    
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
        
        
def write_coverage_ext_strand(outfile, samfile, dialect=csv.excel_tab, write_header=True,
                              chrom=None, start=None, end=None, 
                              one_based=False, progress=None):
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
                  'reads_softclipped_rev',
                  )
    write_stats(stat_coverage_ext_strand, outfile, fieldnames, samfile, 
                dialect=dialect, write_header=write_header,
                chrom=chrom, start=start, end=end, 
                one_based=one_based, progress=progress)
    
    
