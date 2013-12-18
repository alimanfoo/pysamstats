# cython: profile=False

# 0.9
# 0.9.1
__version__ = '0.11.1'


import sys
import numpy as np
cimport numpy as np
import time
import csv
from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t
from libc.math cimport sqrt
from cpython cimport PyBytes_FromStringAndSize
from pysam.csamtools cimport Samfile, Fastafile, PileupProxy, bam1_t, bam_pileup1_t, bam1_cigar, bam1_seq, bam1_qual, IteratorRowRegion


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

cdef char* CODE2CIGAR= "MIDNSHP=X"

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
    cdef int reads_pp = 0

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

    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': n, 
            'reads_pp': reads_pp}


cpdef object construct_rec_coverage_pad(chrom, pos, bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom,
            'pos': pos,
            'reads_all': 0,
            'reads_pp': 0}


def stat_pileup(frec, fpad, Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False, **kwargs):
    cdef PileupProxy col
    cdef int curpos
    start, end = normalise_coords(start, end, one_based)
    it = samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate)

    if pad:
        curpos = start
        curchr = chrom
        for col in it:
            curchr = samfile.getrname(col.tid)
            while curpos < col.pos:
                yield fpad(curchr, curpos, one_based)
                curpos += 1
            yield frec(samfile, col, one_based)
            curpos = col.pos + 1
        if chrom is not None and end is not None:
            while curpos < end:
                yield fpad(chrom, curpos, one_based)
                curpos += 1

    else:
        for col in it:
            yield frec(samfile, col, one_based)


def stat_coverage(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup(construct_rec_coverage, construct_rec_coverage_pad, samfile,
                       chrom=chrom, start=start, end=end, one_based=one_based,
                       truncate=truncate, pad=pad, **kwargs)


def write_coverage(*args, **kwargs):
    try:
        fields = kwargs['fields']
    except:
        fields = ('chrom', 'pos', 'reads_all', 'reads_pp')
    write_stats(stat_coverage, fields, *args, **kwargs)
    
    
def load_coverage(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'),
                     ('reads_pp', 'i4')]
    return load_stats(stat_coverage, default_dtype, *args, **kwargs)
    
    
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
    cdef int reads_fwd = 0
    cdef int reads_rev = 0
    cdef int reads_pp = 0
    cdef int reads_pp_fwd = 0
    cdef int reads_pp_rev = 0
    
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

    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': n, 
            'reads_fwd': reads_fwd, 
            'reads_rev': reads_rev, 
            'reads_pp': reads_pp,
            'reads_pp_fwd': reads_pp_fwd,
            'reads_pp_rev': reads_pp_rev}


cpdef object construct_rec_coverage_strand_pad(chrom, pos, bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom,
            'pos': pos,
            'reads_all': 0,
            'reads_fwd': 0,
            'reads_rev': 0,
            'reads_pp': 0,
            'reads_pp_fwd': 0,
            'reads_pp_rev': 0}


def stat_coverage_strand(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup(construct_rec_coverage_strand, construct_rec_coverage_strand_pad, samfile,
                       chrom=chrom, start=start, end=end, one_based=one_based,
                       truncate=truncate, pad=pad, **kwargs)


def write_coverage_strand(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 
                  'reads_all', 'reads_fwd', 'reads_rev', 
                  'reads_pp', 'reads_pp_fwd', 'reads_pp_rev')
    write_stats(stat_coverage_strand, fieldnames, *args, **kwargs)


def load_coverage_strand(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'),
                     ('reads_fwd', 'i4'),
                     ('reads_rev', 'i4'),
                     ('reads_pp', 'i4'),
                     ('reads_pp_fwd', 'i4'),
                     ('reads_pp_rev', 'i4'),
                     ]
    return load_stats(stat_coverage_strand, default_dtype, *args, **kwargs)
    
    
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
    cdef bint is_duplicate
    cdef bint mate_is_unmappped 
    cdef bint mate_is_reverse
    cdef int tlen
    # counting variables 
    cdef int reads_pp = 0
    cdef int reads_mate_unmapped = 0
    cdef int reads_mate_other_chr = 0
    cdef int reads_mate_same_strand = 0
    cdef int reads_faceaway = 0
    cdef int reads_softclipped = 0
    cdef int reads_duplicate = 0

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
        is_duplicate = <bint>(flag & BAM_FDUP)
        mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
        mate_is_reverse = <bint>(flag & BAM_FMREVERSE)
        tlen = aln.core.isize
        if is_duplicate:
            reads_duplicate += 1
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
            
    return {'chrom': chrom, 'pos': pos, 
               'reads_all': n, 
               'reads_pp': reads_pp,
               'reads_mate_unmapped': reads_mate_unmapped,
               'reads_mate_other_chr': reads_mate_other_chr,
               'reads_mate_same_strand': reads_mate_same_strand,
               'reads_faceaway': reads_faceaway,
               'reads_softclipped': reads_softclipped,
               'reads_duplicate': reads_duplicate}


# def stat_coverage_ext(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_coverage_ext(samfile, col, one_based)
        
        
cpdef object construct_rec_coverage_ext_pad(chrom, pos, bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 'pos': pos,
            'reads_all': 0,
            'reads_pp': 0,
            'reads_mate_unmapped': 0,
            'reads_mate_other_chr': 0,
            'reads_mate_same_strand': 0,
            'reads_faceaway': 0,
            'reads_softclipped': 0,
            'reads_duplicate': 0}


def stat_coverage_ext(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup(construct_rec_coverage_ext, construct_rec_coverage_ext_pad, samfile,
                       chrom=chrom, start=start, end=end, one_based=one_based,
                       truncate=truncate, pad=pad, **kwargs)


def write_coverage_ext(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 
                  'reads_all', 
                  'reads_pp', 
                  'reads_mate_unmapped', 
                  'reads_mate_other_chr',
                  'reads_mate_same_strand',
                  'reads_faceaway', 
                  'reads_softclipped',
                  'reads_duplicate')
    write_stats(stat_coverage_ext, fieldnames, *args, **kwargs)


def load_coverage_ext(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'), 
                     ('reads_pp', 'i4'), 
                     ('reads_mate_unmapped', 'i4'), 
                     ('reads_mate_other_chr', 'i4'),
                     ('reads_mate_same_strand', 'i4'),
                     ('reads_faceaway', 'i4'), 
                     ('reads_softclipped', 'i4'),
                     ('reads_duplicate', 'i4')
                     ]
    return load_stats(stat_coverage_ext, default_dtype, *args, **kwargs)
    
    
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
    cdef bint is_duplicate
    cdef bint mate_is_unmappped 
    cdef bint mate_is_reverse
    cdef int tlen
    # counting variables 
    cdef int reads_rev = 0
    cdef int reads_fwd = 0
    cdef int reads_pp = 0
    cdef int reads_pp_rev = 0
    cdef int reads_pp_fwd = 0
    cdef int reads_mate_unmapped = 0
    cdef int reads_mate_unmapped_rev = 0
    cdef int reads_mate_unmapped_fwd = 0
    cdef int reads_mate_other_chr = 0
    cdef int reads_mate_other_chr_rev = 0
    cdef int reads_mate_other_chr_fwd = 0
    cdef int reads_mate_same_strand = 0
    cdef int reads_mate_same_strand_rev = 0
    cdef int reads_mate_same_strand_fwd = 0
    cdef int reads_faceaway = 0
    cdef int reads_faceaway_rev = 0
    cdef int reads_faceaway_fwd = 0
    cdef int reads_softclipped = 0
    cdef int reads_softclipped_rev = 0
    cdef int reads_softclipped_fwd = 0
    cdef int reads_duplicate = 0
    cdef int reads_duplicate_rev = 0
    cdef int reads_duplicate_fwd = 0

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
        is_duplicate = <bint>(flag & BAM_FDUP)
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
        if is_duplicate:
            reads_duplicate += 1
            if is_reverse:
                reads_duplicate_rev += 1
            else:
                reads_duplicate_fwd += 1
            
    return {'chrom': chrom, 'pos': pos, 
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
           'reads_softclipped_rev': reads_softclipped_rev,
           'reads_duplicate': reads_duplicate,
           'reads_duplicate_fwd': reads_duplicate_fwd,
           'reads_duplicate_rev': reads_duplicate_rev,
           }


# def stat_coverage_ext_strand(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_coverage_ext_strand(samfile, col, one_based)
        
        
cpdef object construct_rec_coverage_ext_strand_pad(chrom, pos, bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 'pos': pos,
           'reads_all': 0,
           'reads_fwd': 0,
           'reads_rev': 0,
           'reads_pp': 0,
           'reads_pp_fwd': 0,
           'reads_pp_rev': 0,
           'reads_mate_unmapped': 0,
           'reads_mate_unmapped_fwd': 0,
           'reads_mate_unmapped_rev': 0,
           'reads_mate_other_chr': 0,
           'reads_mate_other_chr_fwd': 0,
           'reads_mate_other_chr_rev': 0,
           'reads_mate_same_strand': 0,
           'reads_mate_same_strand_fwd': 0,
           'reads_mate_same_strand_rev': 0,
           'reads_faceaway': 0,
           'reads_faceaway_fwd': 0,
           'reads_faceaway_rev': 0,
           'reads_softclipped': 0,
           'reads_softclipped_fwd': 0,
           'reads_softclipped_rev': 0,
           'reads_duplicate': 0,
           'reads_duplicate_fwd': 0,
           'reads_duplicate_rev': 0,
           }


def stat_coverage_ext_strand(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup(construct_rec_coverage_ext_strand, construct_rec_coverage_ext_strand_pad, samfile,
                       chrom=chrom, start=start, end=end, one_based=one_based,
                       truncate=truncate, pad=pad, **kwargs)


def write_coverage_ext_strand(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 
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
                  'reads_duplicate',
                  'reads_duplicate_fwd',
                  'reads_duplicate_rev',
                  )
    write_stats(stat_coverage_ext_strand, fieldnames, *args, **kwargs)


def load_coverage_ext_strand(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'),
                     ('reads_fwd', 'i4'),
                     ('reads_rev', 'i4'),
                     ('reads_pp', 'i4'),
                     ('reads_pp_fwd', 'i4'),
                     ('reads_pp_rev', 'i4'),
                     ('reads_mate_unmapped', 'i4'), 
                     ('reads_mate_unmapped_fwd', 'i4'), 
                     ('reads_mate_unmapped_rev', 'i4'), 
                     ('reads_mate_other_chr', 'i4'),
                     ('reads_mate_other_chr_fwd', 'i4'),
                     ('reads_mate_other_chr_rev', 'i4'),
                     ('reads_mate_same_strand', 'i4'),
                     ('reads_mate_same_strand_fwd', 'i4'),
                     ('reads_mate_same_strand_rev', 'i4'),
                     ('reads_faceaway', 'i4'), 
                     ('reads_faceaway_fwd', 'i4'), 
                     ('reads_faceaway_rev', 'i4'), 
                     ('reads_softclipped', 'i4'),
                     ('reads_softclipped_fwd', 'i4'),
                     ('reads_softclipped_rev', 'i4'),
                     ('reads_duplicate', 'i4'),
                     ('reads_duplicate_fwd', 'i4'),
                     ('reads_duplicate_rev', 'i4'),
                    ]
    return load_stats(stat_coverage_ext_strand, default_dtype, *args, **kwargs)
    
    
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
    cdef int reads_pp = 0
    cdef int matches = 0
    cdef int matches_pp = 0
    cdef int mismatches = 0
    cdef int mismatches_pp = 0
    cdef int deletions = 0
    cdef int deletions_pp = 0
    cdef int insertions = 0
    cdef int insertions_pp = 0
    cdef int A = 0
    cdef int A_pp = 0
    cdef int C = 0
    cdef int C_pp = 0
    cdef int T = 0
    cdef int T_pp = 0
    cdef int G = 0
    cdef int G_pp = 0
    cdef int N = 0
    cdef int N_pp = 0
    
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

    return {'chrom': chrom, 'pos': pos, 'ref': refbase,
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


# def stat_variation(Samfile samfile, Fastafile fafile, chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_variation(samfile, fafile, col, one_based)
        
        
def stat_pileup_withref(frec, fpad, Samfile samfile, Fastafile fafile,
                        chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False,
                        **kwargs):
    cdef PileupProxy col
    cdef int curpos
    start, end = normalise_coords(start, end, one_based)
    it = samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate)

    if pad:
        curpos = start
        curchr = chrom
        for col in it:
            curchr = samfile.getrname(col.tid)
            while curpos < col.pos:
                yield fpad(fafile, curchr, curpos, one_based)
                curpos += 1
            yield frec(samfile, fafile, col, one_based)
            curpos = col.pos + 1
        if chrom is not None and end is not None:
            while curpos < end:
                yield fpad(fafile, chrom, curpos, one_based)
                curpos += 1

    else:
        for col in it:
            yield frec(samfile, fafile, col, one_based)


cpdef object construct_rec_variation_pad(Fastafile fafile, chrom, pos, bint one_based=False):
    refbase = fafile.fetch(reference=chrom, start=pos, end=pos+1).upper()
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 'pos': pos, 'ref': refbase,
            'reads_all': 0, 'reads_pp': 0,
            'matches': 0,
            'matches_pp': 0,
            'mismatches': 0,
            'mismatches_pp': 0,
            'deletions': 0,
            'deletions_pp': 0,
            'insertions': 0,
            'insertions_pp': 0,
            'A': 0, 'A_pp': 0,
            'C': 0, 'C_pp': 0,
            'T': 0, 'T_pp': 0,
            'G': 0, 'G_pp': 0,
            'N': 0, 'N_pp': 0}


def stat_variation(Samfile samfile, Fastafile fafile, chrom=None, start=None, end=None,
                   one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup_withref(construct_rec_variation, construct_rec_variation_pad,
                               samfile, fafile,
                               chrom=chrom, start=start, end=end, one_based=one_based,
                               truncate=truncate, pad=pad, **kwargs)


def write_variation(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 'ref', 
                  'reads_all', 'reads_pp',
                  'matches', 'matches_pp',
                  'mismatches', 'mismatches_pp',
                  'deletions', 'deletions_pp',
                  'insertions', 'insertions_pp',
                  'A', 'A_pp', 'C', 'C_pp', 'T', 'T_pp', 'G', 'G_pp', 'N', 'N_pp')
    write_stats(stat_variation, fieldnames, *args, **kwargs)
    
    
def load_variation(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('ref', 'a1'), 
                     ('reads_all', 'i4'), 
                     ('reads_pp', 'i4'),
                     ('matches', 'i4'), 
                     ('matches_pp', 'i4'),
                     ('mismatches', 'i4'), 
                     ('mismatches_pp', 'i4'),
                     ('deletions', 'i4'), 
                     ('deletions_pp', 'i4'),
                     ('insertions', 'i4'), 
                     ('insertions_pp', 'i4'),
                     ('A', 'i4'), 
                     ('A_pp', 'i4'), 
                     ('C', 'i4'), 
                     ('C_pp', 'i4'), 
                     ('T', 'i4'), 
                     ('T_pp', 'i4'), 
                     ('G', 'i4'), 
                     ('G_pp', 'i4'), 
                     ('N', 'i4'), 
                     ('N_pp', 'i4')
                    ]
    return load_stats(stat_variation, default_dtype, *args, **kwargs)

    
#################################
# STRANDED VARIATION STATISTICS #
#################################


cdef struct CountPpStrand:
    int all, pp, fwd, rev, pp_fwd, pp_rev


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

    return {'chrom': chrom, 'pos': pos, 'ref': refbase,
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


# def stat_variation_strand(Samfile samfile, Fastafile fafile, 
#                           chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_variation_strand(samfile, fafile, col, one_based)
        
        
cpdef object construct_rec_variation_strand_pad(Fastafile fafile, chrom, pos, bint one_based=False):
    refbase = fafile.fetch(reference=chrom, start=pos, end=pos+1).upper()
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 'pos': pos, 'ref': refbase,
            'reads_all': 0, 'reads_fwd': 0, 'reads_rev': 0, 
            'reads_pp': 0, 'reads_pp_fwd': 0, 'reads_pp_rev': 0,
            'matches': 0, 'matches_fwd': 0, 'matches_rev': 0, 
            'matches_pp': 0, 'matches_pp_fwd': 0, 'matches_pp_rev': 0,
            'mismatches': 0, 'mismatches_fwd': 0, 'mismatches_rev': 0, 
            'mismatches_pp': 0, 'mismatches_pp_fwd': 0, 'mismatches_pp_rev': 0,
            'deletions': 0, 'deletions_fwd': 0, 'deletions_rev': 0, 
            'deletions_pp': 0, 'deletions_pp_fwd': 0, 'deletions_pp_rev': 0,
            'insertions': 0, 'insertions_fwd': 0, 'insertions_rev': 0, 
            'insertions_pp': 0, 'insertions_pp_fwd': 0, 'insertions_pp_rev': 0,
            'A': 0, 'A_fwd': 0, 'A_rev': 0, 'A_pp': 0, 'A_pp_fwd': 0, 'A_pp_rev': 0,
            'C': 0, 'C_fwd': 0, 'C_rev': 0, 'C_pp': 0, 'C_pp_fwd': 0, 'C_pp_rev': 0,
            'T': 0, 'T_fwd': 0, 'T_rev': 0, 'T_pp': 0, 'T_pp_fwd': 0, 'T_pp_rev': 0,
            'G': 0, 'G_fwd': 0, 'G_rev': 0, 'G_pp': 0, 'G_pp_fwd': 0, 'G_pp_rev': 0,
            'N': 0, 'N_fwd': 0, 'N_rev': 0, 'N_pp': 0, 'N_pp_fwd': 0, 'N_pp_rev': 0,
            }


def stat_variation_strand(Samfile samfile, Fastafile fafile, chrom=None, start=None, end=None,
                          one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup_withref(construct_rec_variation_strand, construct_rec_variation_strand_pad,
                               samfile, fafile,
                               chrom=chrom, start=start, end=end, one_based=one_based,
                               truncate=truncate, pad=pad, **kwargs)


def write_variation_strand(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 'ref', 
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
    
    
def load_variation_strand(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('ref', 'a1'), 
                     ('reads_all', 'i4'), ('reads_fwd', 'i4'), ('reads_rev', 'i4'), 
                     ('reads_pp', 'i4'), ('reads_pp_fwd', 'i4'), ('reads_pp_rev', 'i4'),
                     ('matches', 'i4'), ('matches_fwd', 'i4'), ('matches_rev', 'i4'), 
                     ('matches_pp', 'i4'), ('matches_pp_fwd', 'i4'), ('matches_pp_rev', 'i4'),
                     ('mismatches', 'i4'), ('mismatches_fwd', 'i4'), ('mismatches_rev', 'i4'), 
                     ('mismatches_pp', 'i4'), ('mismatches_pp_fwd', 'i4'), ('mismatches_pp_rev', 'i4'),
                     ('deletions', 'i4'), ('deletions_fwd', 'i4'), ('deletions_rev', 'i4'), 
                     ('deletions_pp', 'i4'), ('deletions_pp_fwd', 'i4'), ('deletions_pp_rev', 'i4'),
                     ('insertions', 'i4'), ('insertions_fwd', 'i4'), ('insertions_rev', 'i4'), 
                     ('insertions_pp', 'i4'), ('insertions_pp_fwd', 'i4'), ('insertions_pp_rev', 'i4'),
                     ('A', 'i4'), ('A_fwd', 'i4'), ('A_rev', 'i4'), ('A_pp', 'i4'), ('A_pp_fwd', 'i4'), ('A_pp_rev', 'i4'),
                     ('C', 'i4'), ('C_fwd', 'i4'), ('C_rev', 'i4'), ('C_pp', 'i4'), ('C_pp_fwd', 'i4'), ('C_pp_rev', 'i4'),
                     ('T', 'i4'), ('T_fwd', 'i4'), ('T_rev', 'i4'), ('T_pp', 'i4'), ('T_pp_fwd', 'i4'), ('T_pp_rev', 'i4'),
                     ('G', 'i4'), ('G_fwd', 'i4'), ('G_rev', 'i4'), ('G_pp', 'i4'), ('G_pp_fwd', 'i4'), ('G_pp_rev', 'i4'),
                     ('N', 'i4'), ('N_fwd', 'i4'), ('N_rev', 'i4'), ('N_pp', 'i4'), ('N_pp_fwd', 'i4'), ('N_pp_rev', 'i4')
                    ]
    return load_stats(stat_variation_strand, default_dtype, *args, **kwargs)
    

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
    cdef int reads_p = 0 # reads "paired", i.e., mate is mapped to same chromosome, so tlen is meaningful
    cdef int reads_pp = 0 # reads "properly paired", as defined by aligner
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
        rms_tlen = rootmean(tlen_p_squared_sum, reads_p)
        variance_tlen = tlen_p_dev_squared_sum * 1. / reads_p
        std_tlen = int(round(sqrt(variance_tlen)))
    else:
        rms_tlen = std_tlen = mean_tlen = 0
    if reads_pp > 0:
        mean_tlen_pp = int(round(tlen_pp_mean))
        rms_tlen_pp = rootmean(tlen_pp_squared_sum, reads_pp)
        variance_tlen_pp = tlen_pp_dev_squared_sum * 1. / reads_pp
        std_tlen_pp = int(round(sqrt(variance_tlen_pp)))
    else:
        rms_tlen_pp = std_tlen_pp = mean_tlen_pp = 0

    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': n, 
            'reads_paired': reads_p,
            'reads_pp': reads_pp,
            'mean_tlen': mean_tlen,
            'mean_tlen_pp': mean_tlen_pp,
            'rms_tlen': rms_tlen,
            'rms_tlen_pp': rms_tlen_pp,
            'std_tlen': std_tlen,
            'std_tlen_pp': std_tlen_pp}


# def stat_tlen(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_tlen(samfile, col, one_based)
        
        
cpdef object construct_rec_tlen_pad(chrom, pos, bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': 0, 
            'reads_paired': 0,
            'reads_pp': 0,
            'mean_tlen': 0,
            'mean_tlen_pp': 0,
            'rms_tlen': 0,
            'rms_tlen_pp': 0,
            'std_tlen': 0,
            'std_tlen_pp': 0,
            }


def stat_tlen(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup(construct_rec_tlen, construct_rec_tlen_pad, samfile,
                       chrom=chrom, start=start, end=end, one_based=one_based,
                       truncate=truncate, pad=pad, **kwargs)


def write_tlen(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 
                  'reads_all', 'reads_paired', 'reads_pp', 
                  'mean_tlen', 'mean_tlen_pp',
                  'rms_tlen', 'rms_tlen_pp',
                  'std_tlen', 'std_tlen_pp')
    write_stats(stat_tlen, fieldnames, *args, **kwargs)
    

def load_tlen(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'), ('reads_paired', 'i4'), ('reads_pp', 'i4'),
                     ('mean_tlen', 'i4'), ('mean_tlen_pp', 'i4'),
                     ('rms_tlen', 'i4'), ('rms_tlen_pp', 'i4'),
                     ('std_tlen', 'i4'), ('std_tlen_pp', 'i4')
                    ]
    return load_stats(stat_tlen, default_dtype, *args, **kwargs)
    
    
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
    cdef int reads_fwd = 0
    cdef int reads_rev = 0
    cdef int reads_p = 0 # reads "paired", i.e., mate is mapped to same chromosome, so tlen is meaningful
    cdef int reads_p_fwd = 0
    cdef int reads_p_rev = 0
    cdef int reads_pp = 0 # reads "properly paired", as defined by aligner
    cdef int reads_pp_fwd = 0
    cdef int reads_pp_rev = 0
    
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
        rms_tlen = rootmean(tlen_p_squared_sum, reads_p)
        variance_tlen = tlen_p_dev_squared_sum * 1. / reads_p
        std_tlen = int(round(sqrt(variance_tlen)))
    else:
        rms_tlen = std_tlen = mean_tlen = 0
    if reads_p_rev > 0:
        mean_tlen_rev = int(round(tlen_p_rev_mean))
        rms_tlen_rev = rootmean(tlen_p_rev_squared_sum, reads_p_rev)
        variance_tlen_rev = tlen_p_rev_dev_squared_sum * 1. / reads_p_rev
        std_tlen_rev = int(round(sqrt(variance_tlen_rev)))
    else:
        rms_tlen_rev = std_tlen_rev = mean_tlen_rev = 0
    if reads_p_fwd > 0:
        mean_tlen_fwd = int(round(tlen_p_fwd_mean))
        rms_tlen_fwd = rootmean(tlen_p_fwd_squared_sum, reads_p_fwd)
        variance_tlen_fwd = tlen_p_fwd_dev_squared_sum * 1. / reads_p_fwd
        std_tlen_fwd = int(round(sqrt(variance_tlen_fwd)))
    else:
        rms_tlen_fwd = std_tlen_fwd = mean_tlen_fwd = 0
    if reads_pp > 0:
        mean_tlen_pp = int(round(tlen_pp_mean))
        rms_tlen_pp = rootmean(tlen_pp_squared_sum, reads_pp)
        variance_tlen_pp = tlen_pp_dev_squared_sum * 1. / reads_pp
        std_tlen_pp = int(round(sqrt(variance_tlen_pp)))
    else:
        rms_tlen_pp = std_tlen_pp = mean_tlen_pp = 0
    if reads_pp_rev > 0:
        mean_tlen_pp_rev = int(round(tlen_pp_rev_mean))
        rms_tlen_pp_rev = rootmean(tlen_pp_rev_squared_sum, reads_pp_rev)
        variance_tlen_pp_rev = tlen_pp_rev_dev_squared_sum * 1. / reads_pp_rev
        std_tlen_pp_rev = int(round(sqrt(variance_tlen_pp_rev)))
    else:
        rms_tlen_pp_rev = std_tlen_pp_rev = mean_tlen_pp_rev = 0
    if reads_pp_fwd > 0:
        mean_tlen_pp_fwd = int(round(tlen_pp_fwd_mean))
        rms_tlen_pp_fwd = rootmean(tlen_pp_fwd_squared_sum, reads_pp_fwd)
        variance_tlen_pp_fwd = tlen_pp_fwd_dev_squared_sum * 1. / reads_pp_fwd
        std_tlen_pp_fwd = int(round(sqrt(variance_tlen_pp_fwd)))
    else:
        rms_tlen_pp_fwd = std_tlen_pp_fwd = mean_tlen_pp_fwd = 0

    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': n, 
            'reads_fwd': reads_fwd, 
            'reads_rev': reads_rev,
            'reads_paired': reads_p,
            'reads_paired_fwd': reads_p_fwd,
            'reads_paired_rev': reads_p_rev,
            'reads_pp': reads_pp,
            'reads_pp_fwd': reads_pp_fwd,
            'reads_pp_rev': reads_pp_rev,
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


# def stat_tlen_strand(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_tlen_strand(samfile, col, one_based)
        
        
cpdef object construct_rec_tlen_strand_pad(chrom, pos, bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': 0, 
            'reads_fwd': 0, 
            'reads_rev': 0,
            'reads_paired': 0,
            'reads_paired_fwd': 0,
            'reads_paired_rev': 0,
            'reads_pp': 0,
            'reads_pp_fwd': 0,
            'reads_pp_rev': 0,
            'mean_tlen': 0,
            'mean_tlen_fwd': 0,
            'mean_tlen_rev': 0,
            'mean_tlen_pp': 0,
            'mean_tlen_pp_fwd': 0,
            'mean_tlen_pp_rev': 0,
            'rms_tlen': 0,
            'rms_tlen_fwd': 0,
            'rms_tlen_rev': 0,
            'rms_tlen_pp': 0,
            'rms_tlen_pp_fwd': 0,
            'rms_tlen_pp_rev': 0,
            'std_tlen': 0,
            'std_tlen_fwd': 0,
            'std_tlen_rev': 0,
            'std_tlen_pp': 0,
            'std_tlen_pp_fwd': 0,
            'std_tlen_pp_rev': 0,
            }


def stat_tlen_strand(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup(construct_rec_tlen_strand, construct_rec_tlen_strand_pad, samfile,
                       chrom=chrom, start=start, end=end, one_based=one_based,
                       truncate=truncate, pad=pad, **kwargs)


def write_tlen_strand(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 
                  'reads_all', 'reads_fwd', 'reads_rev', 
                  'reads_paired', 'reads_paired_fwd', 'reads_paired_rev', 
                  'reads_pp', 'reads_pp_fwd', 'reads_pp_rev', 
                  'mean_tlen', 'mean_tlen_fwd', 'mean_tlen_rev', 
                  'mean_tlen_pp', 'mean_tlen_pp_fwd', 'mean_tlen_pp_rev',
                  'rms_tlen', 'rms_tlen_fwd', 'rms_tlen_rev', 
                  'rms_tlen_pp', 'rms_tlen_pp_fwd', 'rms_tlen_pp_rev',
                  'std_tlen', 'std_tlen_fwd', 'std_tlen_rev', 
                  'std_tlen_pp', 'std_tlen_pp_fwd', 'std_tlen_pp_rev')
    write_stats(stat_tlen_strand, fieldnames, *args, **kwargs)


def load_tlen_strand(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'), ('reads_fwd', 'i4'), ('reads_rev', 'i4'), 
                     ('reads_paired', 'i4'), ('reads_paired_fwd', 'i4'), ('reads_paired_rev', 'i4'), 
                     ('reads_pp', 'i4'), ('reads_pp_fwd', 'i4'), ('reads_pp_rev', 'i4'), 
                     ('mean_tlen', 'i4'), ('mean_tlen_fwd', 'i4'), ('mean_tlen_rev', 'i4'),
                     ('mean_tlen_pp', 'i4'), ('mean_tlen_pp_fwd', 'i4'), ('mean_tlen_pp_rev', 'i4'),
                     ('rms_tlen', 'i4'), ('rms_tlen_fwd', 'i4'), ('rms_tlen_rev', 'i4'),
                     ('rms_tlen_pp', 'i4'), ('rms_tlen_pp_fwd', 'i4'), ('rms_tlen_pp_rev', 'i4'),
                     ('std_tlen', 'i4'), ('std_tlen_fwd', 'i4'), ('std_tlen_rev', 'i4'),
                     ('std_tlen_pp', 'i4'), ('std_tlen_pp_fwd', 'i4'), ('std_tlen_pp_rev', 'i4')
                    ]
    return load_stats(stat_tlen_strand, default_dtype, *args, **kwargs)
        
    
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
    cdef int reads_pp = 0
    cdef int reads_mapq0 = 0
    cdef int reads_mapq0_pp = 0

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
    rms_mapq = rootmean(mapq_squared_sum, n)
    max_mapq = mapq_max
    if reads_pp > 0:
        rms_mapq_pp = rootmean(mapq_pp_squared_sum, reads_pp)
        max_mapq_pp = mapq_pp_max
    else:
        rms_mapq_pp = max_mapq_pp = 0
        
    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': n, 
            'reads_pp': reads_pp,
            'reads_mapq0': reads_mapq0,
            'reads_mapq0_pp': reads_mapq0_pp,
            'rms_mapq': rms_mapq,
            'rms_mapq_pp': rms_mapq_pp,
            'max_mapq': max_mapq,
            'max_mapq_pp': max_mapq_pp}


# def stat_mapq(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_mapq(samfile, col, one_based)
        
        
cpdef object construct_rec_mapq_pad(chrom, pos, bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': 0, 
            'reads_pp': 0,
            'reads_mapq0': 0,
            'reads_mapq0_pp': 0,
            'rms_mapq': 0,
            'rms_mapq_pp': 0,
            'max_mapq': 0,
            'max_mapq_pp': 0,
            }


def stat_mapq(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup(construct_rec_mapq, construct_rec_mapq_pad, samfile,
                       chrom=chrom, start=start, end=end, one_based=one_based,
                       truncate=truncate, pad=pad, **kwargs)


def write_mapq(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 
                  'reads_all', 'reads_pp',
                  'reads_mapq0', 'reads_mapq0_pp',
                  'rms_mapq', 'rms_mapq_pp',
                  'max_mapq', 'max_mapq_pp')
    write_stats(stat_mapq, fieldnames, *args, **kwargs)
    
    
def load_mapq(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'), ('reads_pp', 'i4'),
                     ('reads_mapq0', 'i4'), ('reads_mapq0_pp', 'i4'),
                     ('rms_mapq', 'i4'), ('rms_mapq_pp', 'i4'),
                     ('max_mapq', 'i4'), ('max_mapq_pp', 'i4')
                    ]
    return load_stats(stat_mapq, default_dtype, *args, **kwargs)
        
    
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

    cdef int reads_rev = 0
    cdef int reads_fwd = 0
    cdef int reads_pp = 0
    cdef int reads_pp_rev = 0
    cdef int reads_pp_fwd = 0
    cdef int reads_mapq0 = 0
    cdef int reads_mapq0_fwd = 0
    cdef int reads_mapq0_rev = 0
    cdef int reads_mapq0_pp = 0
    cdef int reads_mapq0_pp_fwd = 0
    cdef int reads_mapq0_pp_rev = 0

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
    rms_mapq = rootmean(mapq_squared_sum, n)
    max_mapq = mapq_max
    if reads_rev > 0:
        rms_mapq_rev = rootmean(mapq_rev_squared_sum, reads_rev)
        max_mapq_rev = mapq_rev_max
    else:
        rms_mapq_rev = max_mapq_rev = 0
    if reads_fwd > 0:
        rms_mapq_fwd = rootmean(mapq_fwd_squared_sum, reads_fwd)
        max_mapq_fwd = mapq_fwd_max
    else:
        rms_mapq_fwd = max_mapq_fwd = 0
    if reads_pp > 0:
        rms_mapq_pp = rootmean(mapq_pp_squared_sum, reads_pp)
        max_mapq_pp = mapq_pp_max
    else:
        rms_mapq_pp = max_mapq_pp = 0
    if reads_pp_fwd > 0:
        rms_mapq_pp_fwd = rootmean(mapq_pp_fwd_squared_sum, reads_pp_fwd)
        max_mapq_pp_fwd = mapq_pp_fwd_max
    else:
        rms_mapq_pp_fwd = max_mapq_pp_fwd = 0
    if reads_pp_rev > 0:
        rms_mapq_pp_rev = rootmean(mapq_pp_rev_squared_sum, reads_pp_rev)
        max_mapq_pp_rev = mapq_pp_rev_max
    else:
        rms_mapq_pp_rev = max_mapq_pp_rev = 0
        
    return {'chrom': chrom, 
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


# def stat_mapq_strand(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_mapq_strand(samfile, col, one_based)
        
        
cpdef object construct_rec_mapq_strand_pad(chrom, pos, bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': 0,
            'reads_fwd': 0, 
            'reads_rev': 0, 
            'reads_pp': 0,
            'reads_pp_fwd': 0,
            'reads_pp_rev': 0,
            'reads_mapq0': 0,
            'reads_mapq0_fwd': 0, 
            'reads_mapq0_rev': 0, 
            'reads_mapq0_pp': 0,
            'reads_mapq0_pp_fwd': 0,
            'reads_mapq0_pp_rev': 0,
            'rms_mapq': 0,
            'rms_mapq_fwd': 0,
            'rms_mapq_rev': 0,
            'rms_mapq_pp': 0,
            'rms_mapq_pp_fwd': 0,
            'rms_mapq_pp_rev': 0,
            'max_mapq': 0,
            'max_mapq_fwd': 0,
            'max_mapq_rev': 0,
            'max_mapq_pp': 0,
            'max_mapq_pp_fwd': 0,
            'max_mapq_pp_rev': 0,
            }


def stat_mapq_strand(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup(construct_rec_mapq_strand, construct_rec_mapq_strand_pad, samfile,
                       chrom=chrom, start=start, end=end, one_based=one_based,
                       truncate=truncate, pad=pad, **kwargs)


def write_mapq_strand(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 
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
    
    
def load_mapq_strand(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'), ('reads_fwd', 'i4'), ('reads_rev', 'i4'), 
                     ('reads_pp', 'i4'), ('reads_pp_fwd', 'i4'), ('reads_pp_rev', 'i4'), 
                     ('reads_mapq0', 'i4'), ('reads_mapq0_fwd', 'i4'), ('reads_mapq0_rev', 'i4'), 
                     ('reads_mapq0_pp', 'i4'), ('reads_mapq0_pp_fwd', 'i4'), ('reads_mapq0_pp_rev', 'i4'), 
                     ('rms_mapq', 'i4'), ('rms_mapq_fwd', 'i4'), ('rms_mapq_rev', 'i4'),
                     ('rms_mapq_pp', 'i4'), ('rms_mapq_pp_fwd', 'i4'), ('rms_mapq_pp_rev', 'i4'),
                     ('max_mapq', 'i4'), ('max_mapq_fwd', 'i4'), ('max_mapq_rev', 'i4'), 
                     ('max_mapq_pp', 'i4'), ('max_mapq_pp_fwd', 'i4'), ('max_mapq_pp_rev', 'i4'), 
                    ]
    return load_stats(stat_mapq_strand, default_dtype, *args, **kwargs)
        
    
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
    cdef int reads_nodel = 0
    cdef int reads_pp = 0
    cdef int reads_pp_nodel = 0
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
        if is_proper_pair:
            reads_pp += 1
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

    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': n, 
            'reads_pp': reads_pp,
            'rms_baseq': rms_baseq,
            'rms_baseq_pp': rms_baseq_pp}


# def stat_baseq(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_baseq(samfile, col, one_based)
        
        
cpdef object construct_rec_baseq_pad(chrom, pos, bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': 0, 
            'reads_pp': 0,
            'rms_baseq': 0,
            'rms_baseq_pp': 0,
            }


def stat_baseq(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup(construct_rec_baseq, construct_rec_baseq_pad, samfile,
                       chrom=chrom, start=start, end=end, one_based=one_based,
                       truncate=truncate, pad=pad, **kwargs)


def write_baseq(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 
                  'reads_all', 
                  'reads_pp', 
                  'rms_baseq', 
                  'rms_baseq_pp',
                  )
    write_stats(stat_baseq, fieldnames, *args, **kwargs)
    
    
def load_baseq(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'),
                     ('reads_pp', 'i4'),
                     ('rms_baseq', 'i4'),
                     ('rms_baseq_pp', 'i4'),
                    ]
    return load_stats(stat_baseq, default_dtype, *args, **kwargs)
        
    
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

    cdef int reads_rev = 0
    cdef int reads_fwd = 0
    cdef int reads_pp = 0
    cdef int reads_pp_rev = 0
    cdef int reads_pp_fwd = 0
    cdef int reads_nodel = 0
    cdef int reads_rev_nodel = 0
    cdef int reads_fwd_nodel = 0
    cdef int reads_pp_nodel = 0
    cdef int reads_pp_rev_nodel = 0
    cdef int reads_pp_fwd_nodel = 0

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
        
    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': n,
            'reads_fwd': reads_fwd, 
            'reads_rev': reads_rev, 
            'reads_pp': reads_pp,
            'reads_pp_fwd': reads_pp_fwd,
            'reads_pp_rev': reads_pp_rev,
            'rms_baseq': rms_baseq,
            'rms_baseq_fwd': rms_baseq_fwd,
            'rms_baseq_rev': rms_baseq_rev,
            'rms_baseq_pp': rms_baseq_pp,
            'rms_baseq_pp_fwd': rms_baseq_pp_fwd,
            'rms_baseq_pp_rev': rms_baseq_pp_rev,
            }


# def stat_baseq_strand(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_baseq_strand(samfile, col, one_based)
        
        
cpdef object construct_rec_baseq_strand_pad(chrom, pos, bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': 0,
            'reads_fwd': 0, 
            'reads_rev': 0, 
            'reads_pp': 0,
            'reads_pp_fwd': 0,
            'reads_pp_rev': 0,
            'rms_baseq': 0,
            'rms_baseq_fwd': 0,
            'rms_baseq_rev': 0,
            'rms_baseq_pp': 0,
            'rms_baseq_pp_fwd': 0,
            'rms_baseq_pp_rev': 0,
            }


def stat_baseq_strand(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup(construct_rec_baseq_strand, construct_rec_baseq_strand_pad, samfile,
                       chrom=chrom, start=start, end=end, one_based=one_based,
                       truncate=truncate, pad=pad, **kwargs)


def write_baseq_strand(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 
                  'reads_all', 'reads_fwd', 'reads_rev', 
                  'reads_pp', 'reads_pp_fwd', 'reads_pp_rev', 
                  'rms_baseq', 'rms_baseq_fwd', 'rms_baseq_rev', 
                  'rms_baseq_pp', 'rms_baseq_pp_fwd', 'rms_baseq_pp_rev', 
                  )
    write_stats(stat_baseq_strand, fieldnames, *args, **kwargs)
    
    
def load_baseq_strand(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'), ('reads_fwd', 'i4'), ('reads_rev', 'i4'), 
                     ('reads_pp', 'i4'), ('reads_pp_fwd', 'i4'), ('reads_pp_rev', 'i4'), 
                     ('rms_baseq', 'i4'), ('rms_baseq_fwd', 'i4'), ('rms_baseq_rev', 'i4'),
                     ('rms_baseq_pp', 'i4'), ('rms_baseq_pp_fwd', 'i4'), ('rms_baseq_pp_rev', 'i4'),
                    ]
    return load_stats(stat_baseq_strand, default_dtype, *args, **kwargs)
        
    
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
    cdef int reads_nodel = 0
    cdef int reads_pp = 0
    cdef int reads_pp_nodel = 0
    cdef int matches = 0
    cdef int matches_pp = 0
    cdef int mismatches = 0
    cdef int mismatches_pp = 0

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
        if is_proper_pair:
            reads_pp += 1
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

    return {'chrom': chrom, 'pos': pos, 'ref': refbase,
            'reads_all': n, 'reads_pp': reads_pp,
            'matches': matches,
            'matches_pp': matches_pp,
            'mismatches': mismatches,
            'mismatches_pp': mismatches_pp,
            'rms_baseq': rms_baseq,
            'rms_baseq_pp': rms_baseq_pp,
            'rms_baseq_matches': rms_baseq_matches,
            'rms_baseq_matches_pp': rms_baseq_matches_pp,
            'rms_baseq_mismatches': rms_baseq_mismatches,
            'rms_baseq_mismatches_pp': rms_baseq_mismatches_pp,
            }


# def stat_baseq_ext(Samfile samfile, Fastafile fafile, 
#                    chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_baseq_ext(samfile, fafile, col, one_based)
        
        
cpdef object construct_rec_baseq_ext_pad(Fastafile fafile, chrom, pos, bint one_based=False):
    refbase = fafile.fetch(reference=chrom, start=pos, end=pos+1).upper()
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 'pos': pos, 'ref': refbase,
            'reads_all': 0, 'reads_pp': 0,
            'matches': 0,
            'matches_pp': 0,
            'mismatches': 0,
            'mismatches_pp': 0,
            'rms_baseq': 0,
            'rms_baseq_pp': 0,
            'rms_baseq_matches': 0,
            'rms_baseq_matches_pp': 0,
            'rms_baseq_mismatches': 0,
            'rms_baseq_mismatches_pp': 0,
            }


def stat_baseq_ext(Samfile samfile, Fastafile fafile, chrom=None, start=None, end=None,
                   one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup_withref(construct_rec_baseq_ext, construct_rec_baseq_ext_pad,
                               samfile, fafile,
                               chrom=chrom, start=start, end=end, one_based=one_based,
                               truncate=truncate, pad=pad, **kwargs)


def write_baseq_ext(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 
                  'ref', 
                  'reads_all', 'reads_pp',
                  'matches', 'matches_pp',
                  'mismatches', 'mismatches_pp',
                  'rms_baseq', 'rms_baseq_pp',
                  'rms_baseq_matches', 'rms_baseq_matches_pp',
                  'rms_baseq_mismatches', 'rms_baseq_mismatches_pp',
                  )
    write_stats(stat_baseq_ext, fieldnames, *args, **kwargs)
    
    
def load_baseq_ext(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('ref', 'a1'), 
                     ('reads_all', 'i4'), ('reads_pp', 'i4'),
                     ('matches', 'i4'), ('matches_pp', 'i4'),
                     ('mismatches', 'i4'), ('mismatches_pp', 'i4'),
                     ('rms_baseq', 'i4'), ('rms_baseq_pp', 'i4'),
                     ('rms_baseq_matches', 'i4'), ('rms_baseq_matches_pp', 'i4'),
                     ('rms_baseq_mismatches', 'i4'), ('rms_baseq_mismatches_pp', 'i4'),
                    ]
    return load_stats(stat_baseq_ext, default_dtype, *args, **kwargs)
        
    
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
    cdef int reads_fwd = 0
    cdef int reads_rev = 0
    cdef int reads_nodel = 0
    cdef int reads_fwd_nodel = 0
    cdef int reads_rev_nodel = 0
    cdef int reads_pp = 0
    cdef int reads_pp_fwd = 0
    cdef int reads_pp_rev = 0
    cdef int reads_pp_nodel = 0
    cdef int reads_pp_fwd_nodel = 0
    cdef int reads_pp_rev_nodel = 0
    cdef int matches = 0
    cdef int matches_fwd = 0
    cdef int matches_rev = 0
    cdef int matches_pp = 0
    cdef int matches_pp_fwd = 0
    cdef int matches_pp_rev = 0
    cdef int mismatches = 0
    cdef int mismatches_fwd = 0
    cdef int mismatches_rev = 0
    cdef int mismatches_pp = 0
    cdef int mismatches_pp_fwd = 0
    cdef int mismatches_pp_rev = 0

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

    return {'chrom': chrom, 'pos': pos, 'ref': refbase,
            'reads_all': n, 'reads_fwd': reads_fwd, 'reads_rev': reads_rev, 
            'reads_pp': reads_pp, 'reads_pp_fwd': reads_pp_fwd, 'reads_pp_rev': reads_pp_rev,
            'matches': matches, 'matches_fwd': matches_fwd, 'matches_rev': matches_rev, 
            'matches_pp': matches_pp, 'matches_pp_fwd': matches_pp_fwd, 'matches_pp_rev': matches_pp_rev,
            'mismatches': mismatches, 'mismatches_fwd': mismatches_fwd, 'mismatches_rev': mismatches_rev, 
            'mismatches_pp': mismatches_pp, 'mismatches_pp_fwd': mismatches_pp_fwd, 'mismatches_pp_rev': mismatches_pp_rev,
            'rms_baseq': rms_baseq, 'rms_baseq_fwd': rms_baseq_fwd, 'rms_baseq_rev': rms_baseq_rev, 
            'rms_baseq_pp': rms_baseq_pp, 'rms_baseq_pp_fwd': rms_baseq_pp_fwd, 'rms_baseq_pp_rev': rms_baseq_pp_rev,
            'rms_baseq_matches': rms_baseq_matches, 'rms_baseq_matches_fwd': rms_baseq_matches_fwd, 'rms_baseq_matches_rev': rms_baseq_matches_rev, 
            'rms_baseq_matches_pp': rms_baseq_matches_pp, 'rms_baseq_matches_pp_fwd': rms_baseq_matches_pp_fwd, 'rms_baseq_matches_pp_rev': rms_baseq_matches_pp_rev,
            'rms_baseq_mismatches': rms_baseq_mismatches, 'rms_baseq_mismatches_fwd': rms_baseq_mismatches_fwd, 'rms_baseq_mismatches_rev': rms_baseq_mismatches_rev, 
            'rms_baseq_mismatches_pp': rms_baseq_mismatches_pp, 'rms_baseq_mismatches_pp_fwd': rms_baseq_mismatches_pp_fwd, 'rms_baseq_mismatches_pp_rev': rms_baseq_mismatches_pp_rev,
            }


# def stat_baseq_ext_strand(Samfile samfile, Fastafile fafile, 
#                           chrom=None, start=None, end=None, one_based=False, truncate=False, **kwargs):
#     start, end = normalise_coords(start, end, one_based)
#     for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
#         yield construct_rec_baseq_ext_strand(samfile, fafile, col, one_based)
        
        
cpdef object construct_rec_baseq_ext_strand_pad(Fastafile fafile, chrom, pos, bint one_based=False):
    refbase = fafile.fetch(reference=chrom, start=pos, end=pos+1).upper()
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 'pos': pos, 'ref': refbase,
            'reads_all': 0, 'reads_fwd': 0, 'reads_rev': 0, 
            'reads_pp': 0, 'reads_pp_fwd': 0, 'reads_pp_rev': 0,
            'matches': 0, 'matches_fwd': 0, 'matches_rev': 0, 
            'matches_pp': 0, 'matches_pp_fwd': 0, 'matches_pp_rev': 0,
            'mismatches': 0, 'mismatches_fwd': 0, 'mismatches_rev': 0, 
            'mismatches_pp': 0, 'mismatches_pp_fwd': 0, 'mismatches_pp_rev': 0,
            'rms_baseq': 0, 'rms_baseq_fwd': 0, 'rms_baseq_rev': 0, 
            'rms_baseq_pp': 0, 'rms_baseq_pp_fwd': 0, 'rms_baseq_pp_rev': 0,
            'rms_baseq_matches': 0, 'rms_baseq_matches_fwd': 0, 'rms_baseq_matches_rev': 0, 
            'rms_baseq_matches_pp': 0, 'rms_baseq_matches_pp_fwd': 0, 'rms_baseq_matches_pp_rev': 0,
            'rms_baseq_mismatches': 0, 'rms_baseq_mismatches_fwd': 0, 'rms_baseq_mismatches_rev': 0, 
            'rms_baseq_mismatches_pp': 0, 'rms_baseq_mismatches_pp_fwd': 0, 'rms_baseq_mismatches_pp_rev': 0,
            }


def stat_baseq_ext_strand(Samfile samfile, Fastafile fafile, chrom=None, start=None, end=None,
                   one_based=False, truncate=False, pad=False, **kwargs):
    return stat_pileup_withref(construct_rec_baseq_ext_strand, construct_rec_baseq_ext_strand_pad,
                               samfile, fafile,
                               chrom=chrom, start=start, end=end, one_based=one_based,
                               truncate=truncate, pad=pad, **kwargs)


def write_baseq_ext_strand(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 'ref', 
                  'reads_all', 'reads_fwd', 'reads_rev', 
                  'reads_pp', 'reads_pp_fwd', 'reads_pp_rev',
                  'matches', 'matches_fwd', 'matches_rev', 
                  'matches_pp', 'matches_pp_fwd', 'matches_pp_rev',
                  'mismatches', 'mismatches_fwd', 'mismatches_rev', 
                  'mismatches_pp', 'mismatches_pp_fwd', 'mismatches_pp_rev', 
                  'rms_baseq', 'rms_baseq_fwd', 'rms_baseq_rev', 
                  'rms_baseq_pp', 'rms_baseq_pp_fwd', 'rms_baseq_pp_rev',
                  'rms_baseq_matches', 'rms_baseq_matches_fwd', 'rms_baseq_matches_rev', 
                  'rms_baseq_matches_pp', 'rms_baseq_matches_pp_fwd', 'rms_baseq_matches_pp_rev',
                  'rms_baseq_mismatches', 'rms_baseq_mismatches_fwd', 'rms_baseq_mismatches_rev', 
                  'rms_baseq_mismatches_pp', 'rms_baseq_mismatches_pp_fwd', 'rms_baseq_mismatches_pp_rev',
                  )
    write_stats(stat_baseq_ext_strand, fieldnames, *args, **kwargs)
    
    
def load_baseq_ext_strand(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('ref', 'a1'), 
                     ('reads_all', 'i4'), ('reads_fwd', 'i4'), ('reads_rev', 'i4'),
                     ('reads_pp', 'i4'), ('reads_pp_fwd', 'i4'), ('reads_pp_rev', 'i4'),
                     ('matches', 'i4'), ('matches_fwd', 'i4'), ('matches_rev', 'i4'),
                     ('matches_pp', 'i4'), ('matches_pp_fwd', 'i4'), ('matches_pp_rev', 'i4'),
                     ('mismatches', 'i4'), ('mismatches_fwd', 'i4'), ('mismatches_rev', 'i4'),
                     ('mismatches_pp', 'i4'), ('mismatches_pp_fwd', 'i4'), ('mismatches_pp_rev', 'i4'),
                     ('rms_baseq', 'i4'), ('rms_baseq_fwd', 'i4'), ('rms_baseq_rev', 'i4'),
                     ('rms_baseq_pp', 'i4'), ('rms_baseq_pp_fwd', 'i4'), ('rms_baseq_pp_rev', 'i4'),
                     ('rms_baseq_matches', 'i4'), ('rms_baseq_matches_fwd', 'i4'), ('rms_baseq_matches_rev', 'i4'),
                     ('rms_baseq_matches_pp', 'i4'), ('rms_baseq_matches_pp_fwd', 'i4'), ('rms_baseq_matches_pp_rev', 'i4'),
                     ('rms_baseq_mismatches', 'i4'), ('rms_baseq_mismatches_fwd', 'i4'), ('rms_baseq_mismatches_rev', 'i4'),
                     ('rms_baseq_mismatches_pp', 'i4'), ('rms_baseq_mismatches_pp_fwd', 'i4'), ('rms_baseq_mismatches_pp_rev', 'i4'),
                    ]
    return load_stats(stat_baseq_ext_strand, default_dtype, *args, **kwargs)
        
    
##############################
# NORMED COVERAGE STATISTICS #
##############################


from bisect import bisect_left


def construct_rec_coverage_normed_pad(chrom, pos, one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 'pos': pos,
            'reads_all': 0,
            'dp_normed_median': 0,
            'dp_normed_mean': 0,
            'dp_percentile': 0}


def stat_coverage_normed(Samfile samfile, chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False,
                         **kwargs):
    start, end = normalise_coords(start, end, one_based)
    
    # first need to load the coverage data into an array, to calculate the median
    it = (col.n for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate))
    a = np.fromiter(it, dtype='u4')
    dp_mean = np.mean(a)
    dp_median = np.median(a)
    dp_percentiles = [np.percentile(a, q) for q in range(101)]

    def construct_rec_coverage_normed(Samfile samfile, PileupProxy col, bint one_based=False):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos
        dp = col.n
        dp_normed_median = dp * 1. / dp_median
        dp_normed_mean = dp * 1. / dp_mean
        dp_percentile = bisect_left(dp_percentiles, dp)
        return {'chrom': chrom, 'pos': pos,
                'reads_all': col.n,
                'dp_normed_median': dp_normed_median,
                'dp_normed_mean': dp_normed_mean,
                'dp_percentile': dp_percentile}

    # then iterate again to generate stats
    return stat_pileup(construct_rec_coverage_normed, construct_rec_coverage_normed_pad,
                       samfile, chrom=chrom, start=start, end=end, one_based=one_based,
                       truncate=truncate, pad=pad, **kwargs)

    # for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
    #     chrom = samfile.getrname(col.tid)
    #     pos = col.pos + 1 if one_based else col.pos
    #     dp = col.n
    #     dp_normed_median = dp * 1. / dp_median
    #     dp_normed_mean = dp * 1. / dp_mean
    #     dp_percentile = bisect_left(dp_percentiles, dp)
    #     yield {'chrom': chrom, 'pos': pos,
    #            'reads_all': col.n,
    #            'dp_normed_median': dp_normed_median,
    #            'dp_normed_mean': dp_normed_mean,
    #            'dp_percentile': dp_percentile}
        
        
def write_coverage_normed(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 
                  'reads_all', 
                  'dp_normed_median', 
                  'dp_normed_mean',
                  'dp_percentile'
                  )
    write_stats(stat_coverage_normed, fieldnames, *args, **kwargs)
    
    
def load_coverage_normed(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'), 
                     ('dp_normed_median', 'f4'), 
                     ('dp_normed_mean', 'f4'),
                     ('dp_percentile', 'u1')
                    ]
    return load_stats(stat_coverage_normed, default_dtype, *args, **kwargs)
        
    
#################################################
# BASIC COVERAGE STATISTICS WITH GC COMPOSITION #
#################################################


from collections import Counter


def stat_coverage_gc(Samfile samfile, Fastafile fafile,
                     chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False,
                     window_size=300, window_offset=None, **kwargs):
    cdef Py_ssize_t i # loop index
    cdef char* seq # sequence window
    cdef int gc_count 
    start, end = normalise_coords(start, end, one_based)
    if window_offset is None:
        window_offset = window_size / 2
        
    def construct_rec_coverage_gc(Samfile samfile, Fastafile fafile, PileupProxy col, bint one_based):
        chrom = samfile.getrname(col.tid)

        ref_window_start = col.pos - window_offset
        ref_window_end = ref_window_start + window_size
        if ref_window_start < 0:
            ref_window_start = 0
        ref_window = fafile.fetch(chrom, ref_window_start, ref_window_end).lower()
        if len(ref_window) == 0:
            gc_percent = -1
        else:
            gc_percent = gc_content(ref_window)

        rec = construct_rec_coverage(samfile, col, one_based)
        rec['gc'] = gc_percent
        return rec

    def construct_rec_coverage_gc_pad(Fastafile fafile, chrom, pos, bint one_based):
        ref_window_start = pos - window_offset
        ref_window_end = ref_window_start + window_size
        if ref_window_start < 0:
            ref_window_start = 0
        ref_window = fafile.fetch(chrom, ref_window_start, ref_window_end).lower()
        if len(ref_window) == 0:
            gc_percent = -1
        else:
            gc_percent = gc_content(ref_window)
        rec = construct_rec_coverage_pad(chrom, pos, one_based)
        rec['gc'] = gc_percent
        return rec

    return stat_pileup_withref(construct_rec_coverage_gc, construct_rec_coverage_gc_pad,
                               samfile, fafile, chrom=chrom, start=start, end=end, one_based=one_based,
                               truncate=truncate, pad=pad, **kwargs)

    # for col in samfile.pileup(reference=chrom, start=start, end=end, truncate=truncate):
    #
    #     chrom = samfile.getrname(col.tid)
    #
    #     if col.pos <= window_offset:
    #         continue # until we get a bit further into the chromosome
    #
    #     ref_window_start = col.pos - window_offset
    #     ref_window_end = ref_window_start + window_size
    #     ref_window = fafile.fetch(chrom, ref_window_start, ref_window_end).lower()
    #
    #     if len(ref_window) == 0:
    #         break # because we've hit the end of the chromosome
    #
    #     gc_percent = gc_content(ref_window)
    #
    #     rec = construct_rec_coverage(samfile, col, one_based)
    #     rec['gc'] = gc_percent
    #     yield rec
        
        
def write_coverage_gc(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 'gc', 'reads_all', 'reads_pp')
    write_stats(stat_coverage_gc, fieldnames, *args, **kwargs)
    

def load_coverage_gc(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('gc', 'u1'),
                     ('reads_all', 'i4'), 
                     ('reads_pp', 'i4'), 
                    ]
    return load_stats(stat_coverage_gc, default_dtype, *args, **kwargs)
        
    
####################################
# COVERAGE STATISTICS NORMED BY GC #
####################################


def stat_coverage_normed_gc(Samfile samfile, Fastafile fafile, 
                            chrom=None, start=None, end=None, one_based=False, truncate=False, pad=False,
                            window_size=300, window_offset=None, **kwargs):
    start, end = normalise_coords(start, end, one_based)
    if window_offset is None:
        window_offset = window_size / 2
    
    # first need to load the coverage data into an array, to calculate the median
    recs = stat_coverage_gc(samfile, fafile, chrom=chrom, start=start, end=end, 
                            one_based=one_based, truncate=truncate, pad=pad,
                            window_size=window_size, window_offset=window_offset)
    it = ((rec['reads_all'], rec['gc']) for rec in recs)
    a = np.fromiter(it, dtype=[('dp', 'u4'), ('gc', 'u1')]).view(np.recarray)
    dp_mean = np.mean(a.dp)
    dp_median = np.median(a.dp)
    dp_percentiles = [np.percentile(a.dp, q) for q in range(101)]    
    dp_mean_bygc = dict()
    dp_median_bygc = dict()
    dp_percentiles_bygc = dict()
    for gc in range(101):
        flt = a.gc == gc
        if np.count_nonzero(flt) > 0:
            b = a[flt].dp
            dp_mean_bygc[gc] = np.mean(b)
            dp_median_bygc[gc] = np.median(b)
            dp_percentiles_bygc[gc] = [np.percentile(b, q) for q in range(101)]
    
    # second pass
    recs = stat_coverage_gc(samfile, fafile, chrom=chrom, start=start, end=end, 
                            one_based=one_based, truncate=truncate, pad=pad,
                            window_size=window_size, window_offset=window_offset)
    for rec in recs:
        dp = rec['reads_all']
        gc = rec['gc']
        dp_normed_median = dp * 1. / dp_median
        dp_normed_mean = dp * 1. / dp_mean
        dp_percentile = bisect_left(dp_percentiles, dp)
        dp_normed_median_bygc = dp * 1. / dp_median_bygc[gc]
        dp_normed_mean_bygc = dp * 1. / dp_mean_bygc[gc]
        dp_percentile_bygc = bisect_left(dp_percentiles_bygc[gc], dp)
        rec['dp_normed_median'] = dp_normed_median
        rec['dp_normed_mean'] = dp_normed_mean
        rec['dp_percentile'] = dp_percentile
        rec['dp_normed_median_gc'] = dp_normed_median_bygc
        rec['dp_normed_mean_gc'] = dp_normed_mean_bygc
        rec['dp_percentile_gc'] = dp_percentile_bygc
        del rec['reads_pp']
        yield rec
        
        
def write_coverage_normed_gc(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 'gc',
                  'reads_all', #'reads_pp',
                  'dp_normed_median', 
                  'dp_normed_mean',
                  'dp_percentile',
                  'dp_normed_median_gc',
                  'dp_normed_mean_gc',
                  'dp_percentile_gc')
    write_stats(stat_coverage_normed_gc, fieldnames, *args, **kwargs)
    

def load_coverage_normed_gc(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('gc', 'u1'),
                      ('reads_all', 'i4'), #('reads_pp', 'i4'),
                      ('dp_normed_median', 'f4'), 
                      ('dp_normed_mean', 'f4'),
                      ('dp_percentile', 'u1'),
                      ('dp_normed_median_gc', 'f4'),
                      ('dp_normed_mean_gc', 'f4'),
                      ('dp_percentile_gc', 'u1')
                    ]
    return load_stats(stat_coverage_normed_gc, default_dtype, *args, **kwargs)
        
    
from itertools import chain


###################
# BINNED COVERAGE #
###################


def stat_coverage_binned(Samfile samfile, Fastafile fastafile, 
                         chrom=None, start=None, end=None, one_based=False,
                         window_size=300, window_offset=None, **kwargs):
    if window_offset is None:
        window_offset = window_size / 2
    if chrom is None:
        it = chain(*[_iter_coverage_binned(samfile, fastafile, chrom, None, None, one_based, window_size, window_offset) 
                     for chrom in sorted(samfile.references)])
    else:
        it = _iter_coverage_binned(samfile, fastafile, chrom, start, end, one_based, window_size, window_offset)
    return it


cdef inline int gc_content(ref_window):
    cdef Py_ssize_t i, n
    cdef char* seq
    cdef int gc_count = 0
    n = len(ref_window)
    seq = ref_window
    for i in range(n):
        if seq[i] == 'g' or seq[i] == 'c':
            gc_count += 1
    gc_percent = int(round(gc_count * 100. / n))
    return gc_percent

        
def _iter_coverage_binned(Samfile samfile, Fastafile fastafile, 
                          chrom, start, end, one_based, 
                          int window_size, int window_offset):
    assert chrom is not None, 'unexpected error: chromosome is None'
    cdef int rtid, rstart, rend, has_coord, bin_start, bin_end
    cdef int reads_all, reads_pp
    cdef bam1_t * b
    cdef uint32_t flag
    cdef bint is_unmapped
    cdef bint is_proper_pair
    cdef IteratorRowRegion it
    start, end = normalise_coords(start, end, one_based)
    has_coord, rtid, rstart, rend = samfile._parseRegion(chrom, start, end, None)
    it = IteratorRowRegion(samfile, rtid, rstart, rend, reopen=False)
    b = it.b
    # setup first bin
    bin_start = rstart
    bin_end = bin_start + window_size
    reads_all = reads_pp = 0

    # iterate over reads
    it.cnext()
    while it.retval > 0:
        while b.core.pos > bin_end: # end of bin, yield record
            # determine %GC
            ref_window = fastafile.fetch(chrom, bin_start, bin_end).lower()
            if len(ref_window) == 0:
                raise StopIteration # because we've hit the end of the chromosome
            gc_percent = gc_content(ref_window)
            # yield record for bin
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'gc': gc_percent, 'reads_all': reads_all, 'reads_pp': reads_pp}
            yield rec
            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size
            reads_all = reads_pp = 0
        # increment counters
        flag = b.core.flag
        is_unmapped = <bint>(flag & BAM_FUNMAP)
        if not is_unmapped:
            reads_all += 1
            is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
            if is_proper_pair:
                reads_pp += 1
        # move iterator on
        it.cnext()

    # deal with last non-empty bin
    ref_window = fastafile.fetch(chrom, bin_start, bin_end).lower()
    if len(ref_window) == 0:
        raise StopIteration # because we've hit the end of the chromosome
    gc_percent = gc_content(ref_window)
    # yield record for bin
    pos = bin_start + window_offset
    if one_based:
        pos += 1
    rec = {'chrom': chrom, 'pos': pos,
           'gc': gc_percent, 'reads_all': reads_all, 'reads_pp': reads_pp}
    yield rec

    # deal with empty bins up to explicit end
    if end is not None:
        while bin_end < end:
            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size
            ref_window = fastafile.fetch(chrom, bin_start, bin_end).lower()
            if len(ref_window) == 0:
                raise StopIteration # because we've hit the end of the chromosome
            gc_percent = gc_content(ref_window)
            # yield record for bin
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'gc': gc_percent, 'reads_all': 0, 'reads_pp': 0}
            yield rec


    
def write_coverage_binned(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 'gc', 'reads_all', 'reads_pp')
    write_stats(stat_coverage_binned, fieldnames, *args, **kwargs)
    

def load_coverage_binned(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('gc', 'u1'),
                     ('reads_all', 'i4'), 
                     ('reads_pp', 'i4'),
                    ]
    return load_stats(stat_coverage_binned, default_dtype, *args, **kwargs)
        
    
############################################
# BINNED COVERAGE WITH EXTENDED PROPERTIES #
############################################


def stat_coverage_ext_binned(Samfile samfile, Fastafile fastafile, 
                         chrom=None, start=None, end=None, one_based=False,
                         window_size=300, window_offset=None, **kwargs):
    if window_offset is None:
        window_offset = window_size / 2
    if chrom is None:
        it = chain(*[_iter_coverage_ext_binned(samfile, fastafile, chrom, None, None, one_based, window_size, window_offset) 
                     for chrom in sorted(samfile.references)])
    else:
        it = _iter_coverage_ext_binned(samfile, fastafile, chrom, start, end, one_based, window_size, window_offset)
    return it
        
        
def _iter_coverage_ext_binned(Samfile samfile, Fastafile fastafile, 
                          chrom, start, end, one_based, 
                          int window_size, int window_offset):
    assert chrom is not None, 'unexpected error: chromosome is None'
    cdef int rtid, rstart, rend, has_coord, bin_start, bin_end
    cdef int reads_all, reads_pp, reads_mate_unmapped, reads_mate_other_chr, reads_mate_same_strand, reads_faceaway, reads_softclipped, reads_duplicate
    cdef bam1_t * b
    cdef uint32_t flag
    cdef bint is_unmapped
    cdef bint is_reverse
    cdef bint is_proper_pair 
    cdef bint is_duplicate
    cdef bint mate_is_unmappped 
    cdef bint mate_is_reverse
    cdef int tlen
    cdef IteratorRowRegion it
    start, end = normalise_coords(start, end, one_based)
    has_coord, rtid, rstart, rend = samfile._parseRegion(chrom, start, end, None)
    it = IteratorRowRegion(samfile, rtid, rstart, rend, reopen=False)
    b = it.b

    # setup first bin
    bin_start = rstart
    bin_end = bin_start + window_size
    reads_all = reads_pp = reads_mate_unmapped = reads_mate_other_chr = reads_mate_same_strand = reads_faceaway = reads_softclipped = reads_duplicate = 0

    # iterate over reads
    it.cnext()
    while it.retval > 0:
        while b.core.pos > bin_end: # end of bin, yield record
            # determine %GC
            ref_window = fastafile.fetch(chrom, bin_start, bin_end).lower()
            if len(ref_window) == 0:
                raise StopIteration # because we've hit the end of the chromosome
            gc_percent = gc_content(ref_window)
            # yield record for bin
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'gc': gc_percent, 'reads_all': reads_all, 'reads_pp': reads_pp,
                   'reads_mate_unmapped': reads_mate_unmapped,
                   'reads_mate_other_chr': reads_mate_other_chr,
                   'reads_mate_same_strand': reads_mate_same_strand,
                   'reads_faceaway': reads_faceaway,
                   'reads_softclipped': reads_softclipped,
                   'reads_duplicate': reads_duplicate}
            yield rec
            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size
            reads_all = reads_pp = reads_mate_unmapped = reads_mate_other_chr = reads_mate_same_strand = reads_faceaway = reads_softclipped = reads_duplicate = 0
        # increment counters
        flag = b.core.flag
        is_unmapped = <bint>(flag & BAM_FUNMAP)
        if not is_unmapped:
            reads_all += 1
            is_reverse = <bint>(flag & BAM_FREVERSE)
            is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
            is_duplicate = <bint>(flag & BAM_FDUP)
            mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
            mate_is_reverse = <bint>(flag & BAM_FMREVERSE)
            tlen = b.core.isize
            if is_duplicate:
                reads_duplicate += 1
            if is_proper_pair:
                reads_pp += 1
            if mate_is_unmapped:
                reads_mate_unmapped += 1
            elif b.core.tid != b.core.mtid:
                reads_mate_other_chr += 1
            elif (is_reverse and mate_is_reverse) or (not is_reverse and not mate_is_reverse):
                reads_mate_same_strand += 1
            elif (is_reverse and tlen > 0) or (not is_reverse and tlen < 0):
                reads_faceaway += 1
            if is_softclipped(b):
                reads_softclipped += 1
        # move iterator on
        it.cnext()

    # deal with last non-empty bin
    ref_window = fastafile.fetch(chrom, bin_start, bin_end).lower()
    if len(ref_window) == 0:
        raise StopIteration # because we've hit the end of the chromosome
    gc_percent = gc_content(ref_window)
    # yield record for bin
    pos = bin_start + window_offset
    if one_based:
        pos += 1
    rec = {'chrom': chrom, 'pos': pos,
           'gc': gc_percent, 'reads_all': reads_all, 'reads_pp': reads_pp,
           'reads_mate_unmapped': reads_mate_unmapped,
           'reads_mate_other_chr': reads_mate_other_chr,
           'reads_mate_same_strand': reads_mate_same_strand,
           'reads_faceaway': reads_faceaway,
           'reads_softclipped': reads_softclipped,
           'reads_duplicate': reads_duplicate}
    yield rec

    # deal with empty bins up to explicit end
    if end is not None:
        while bin_end < end:
            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size
            ref_window = fastafile.fetch(chrom, bin_start, bin_end).lower()
            if len(ref_window) == 0:
                raise StopIteration # because we've hit the end of the chromosome
            gc_percent = gc_content(ref_window)
            # yield record for bin
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'gc': gc_percent, 'reads_all': 0, 'reads_pp': 0,
                   'reads_mate_unmapped': 0,
                   'reads_mate_other_chr': 0,
                   'reads_mate_same_strand': 0,
                   'reads_faceaway': 0,
                   'reads_softclipped': 0,
                   'reads_duplicate': 0}
            yield rec


def write_coverage_ext_binned(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 'gc', 
                  'reads_all', 
                  'reads_pp', 
                  'reads_mate_unmapped', 
                  'reads_mate_other_chr',
                  'reads_mate_same_strand',
                  'reads_faceaway', 
                  'reads_softclipped',
                  'reads_duplicate')
    write_stats(stat_coverage_ext_binned, fieldnames, *args, **kwargs)
    

def load_coverage_ext_binned(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('gc', 'u1'),
                     ('reads_all', 'i4'), 
                     ('reads_pp', 'i4'),
                     ('reads_mate_unmapped', 'i4'), 
                     ('reads_mate_other_chr', 'i4'),
                     ('reads_mate_same_strand', 'i4'),
                     ('reads_faceaway', 'i4'), 
                     ('reads_softclipped', 'i4'),
                     ('reads_duplicate', 'i4')
                    ]
    return load_stats(stat_coverage_ext_binned, default_dtype, *args, **kwargs)
        
    
###############
# BINNED MAPQ #
###############


def stat_mapq_binned(Samfile samfile, 
                     chrom=None, start=None, end=None, one_based=False,
                     window_size=300, window_offset=None, **kwargs):
    if window_offset is None:
        window_offset = window_size / 2
    if chrom is None:
        it = chain(*[_iter_mapq_binned(samfile, chrom, None, None, one_based, window_size, window_offset) 
                     for chrom in sorted(samfile.references)])
    else:
        it = _iter_mapq_binned(samfile, chrom, start, end, one_based, window_size, window_offset)
    return it
        
        
def _iter_mapq_binned(Samfile samfile,  
                          chrom, start, end, one_based, 
                          int window_size, int window_offset):
    assert chrom is not None, 'unexpected error: chromosome is None'
    cdef int rtid, rstart, rend, has_coord, bin_start, bin_end
    cdef int reads_all, reads_mapq0
    cdef bam1_t * b
    cdef uint32_t flag
    cdef bint is_unmapped
    cdef IteratorRowRegion it
    cdef uint64_t mapq
    cdef uint64_t mapq_squared
    cdef uint64_t mapq_squared_sum = 0
    start, end = normalise_coords(start, end, one_based)
    has_coord, rtid, rstart, rend = samfile._parseRegion(chrom, start, end, None)
    it = IteratorRowRegion(samfile, rtid, rstart, rend, reopen=False)
    b = it.b

    # setup first bin
    bin_start = rstart
    bin_end = bin_start + window_size
    reads_all = reads_mapq0 = mapq_squared_sum = 0

    # iterate over reads
    it.cnext()
    while it.retval > 0:
        while b.core.pos > bin_end: # end of bin, yield record
            # yield record for bin
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'reads_all': reads_all,
                   'reads_mapq0': reads_mapq0,
                   'rms_mapq': rootmean(mapq_squared_sum, reads_all)}
            yield rec
            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size
            reads_all = reads_mapq0 = mapq_squared_sum = 0
        # increment counters
        flag = b.core.flag
        is_unmapped = <bint>(flag & BAM_FUNMAP)
        if not is_unmapped:
            reads_all += 1
            mapq = b.core.qual
            mapq_squared = mapq**2
            mapq_squared_sum += mapq_squared
            if mapq == 0:
                reads_mapq0 += 1
        # move iterator on
        it.cnext()

    # deal with last non-empty bin
    pos = bin_start + window_offset
    if one_based:
        pos += 1
    rec = {'chrom': chrom, 'pos': pos,
           'reads_all': reads_all,
           'reads_mapq0': reads_mapq0,
           'rms_mapq': rootmean(mapq_squared_sum, reads_all)}
    yield rec

    # deal with empty bins up to explicit end
    if end is not None:
        while bin_end < end:
            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'reads_all': 0,
                   'reads_mapq0': 0,
                   'rms_mapq': 0}
            yield rec


def write_mapq_binned(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 'reads_all', 'reads_mapq0', 'rms_mapq')
    write_stats(stat_mapq_binned, fieldnames, *args, **kwargs)
    

def load_mapq_binned(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'), 
                     ('reads_mapq0', 'i4'),
                     ('rms_mapq', 'i4'),
                    ]
    return load_stats(stat_mapq_binned, default_dtype, *args, **kwargs)
        
    
################
# BINNED CIGAR #
################


def stat_alignment_binned(Samfile samfile, 
                      chrom=None, start=None, end=None, one_based=False,
                      window_size=300, window_offset=None, **kwargs):
    if window_offset is None:
        window_offset = window_size / 2
    if chrom is None:
        it = chain(*[_iter_alignment_binned(samfile, chrom, None, None, one_based, window_size, window_offset) 
                     for chrom in sorted(samfile.references)])
    else:
        it = _iter_alignment_binned(samfile, chrom, start, end, one_based, window_size, window_offset)
    return it
        
        
def _iter_alignment_binned(Samfile samfile,  
                       chrom, start, end, one_based, 
                       int window_size, int window_offset):
    assert chrom is not None, 'unexpected error: chromosome is None'
    cdef int rtid, rstart, rend, has_coord, bin_start, bin_end
    cdef bam1_t * b
    cdef uint32_t flag
    cdef bint is_unmapped
    cdef bint is_proper_pair
    cdef IteratorRowRegion it
    cdef Py_ssize_t i # loop index
    cdef int reads_all, k, op, l
    cdef int M, I, D, N, S, H, P, EQ, X
    reads_all = M = I = D = N = S = H = P = EQ = X = 0
    start, end = normalise_coords(start, end, one_based)
    has_coord, rtid, rstart, rend = samfile._parseRegion(chrom, start, end, None)
    it = IteratorRowRegion(samfile, rtid, rstart, rend, reopen=False)
    b = it.b

    # setup first bin
    bin_start = rstart
    bin_end = bin_start + window_size
    c = Counter()

    # iterate over reads
    it.cnext()
    while it.retval > 0:
        while b.core.pos > bin_end:  # end of bin, yield record
            # yield record for bin
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos, 'reads_all': reads_all,
                   'M': M, 'I': I, 'D': D, 'N': N, 'S': S, 'H': H, 'P': P, '=': EQ, 'X': X,
                   'bases_all': M + I + S + EQ + X}
            yield rec
            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size
            reads_all = M = I = D = N = S = H = P = EQ = X = 0
        # increment counters
        flag = b.core.flag
        is_unmapped = <bint>(flag & BAM_FUNMAP)
        if not is_unmapped:
            cigar_p = bam1_cigar(b)
            cigar = list()
            for k in range(b.core.n_cigar):
                op = cigar_p[k] & BAM_CIGAR_MASK
                l = cigar_p[k] >> BAM_CIGAR_SHIFT
                cigar.append((op, l))
                if op == BAM_CMATCH:
                    M += l
                elif op == BAM_CINS:
                    I += l
                elif op == BAM_CDEL:
                    D += l
                elif op == BAM_CREF_SKIP:
                    N += l
                elif op == BAM_CSOFT_CLIP:
                    S += l
                elif op == BAM_CHARD_CLIP:
                    H += l
                elif op == BAM_CPAD:
                    P += l
                elif op == BAM_CEQUAL:
                    EQ += l
                elif op == BAM_CDIFF:
                    X += l
            reads_all += 1
        # move iterator on
        it.cnext()

    # deal with last non-empty bin
    pos = bin_start + window_offset
    if one_based:
        pos += 1
    rec = {'chrom': chrom, 'pos': pos, 'reads_all': reads_all,
           'M': M, 'I': I, 'D': D, 'N': N, 'S': S, 'H': H, 'P': P, '=': EQ, 'X': X,
           'bases_all': M + I + S + EQ + X}
    yield rec
    # start new bin
    bin_start = bin_end
    bin_end = bin_start + window_size
    reads_all = M = I = D = N = S = H = P = EQ = X = 0

    # deal with empty bins up to explicit end
    if end is not None:
        while bin_end < end:
            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size
            reads_all = M = I = D = N = S = H = P = EQ = X = 0
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos, 'reads_all': reads_all,
                   'M': M, 'I': I, 'D': D, 'N': N, 'S': S, 'H': H, 'P': P, '=': EQ, 'X': X,
                   'bases_all': M + I + S + EQ + X}
            yield rec

    
def write_alignment_binned(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 'reads_all', 'bases_all', 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X')
    write_stats(stat_alignment_binned, fieldnames, *args, **kwargs)
    

def load_alignment_binned(*args, **kwargs):
    default_dtype = [('chrom', 'a12'), 
                     ('pos', 'i4'),
                     ('reads_all', 'i4'), 
                     ('bases_all', 'i4'), 
                     ('M', 'i4'), 
                     ('I', 'i4'), 
                     ('D', 'i4'), 
                     ('N', 'i4'), 
                     ('S', 'i4'), 
                     ('H', 'i4'), 
                     ('P', 'i4'), 
                     ('=', 'i4'), 
                     ('X', 'i4')
                    ]
    return load_stats(stat_alignment_binned, default_dtype, *args, **kwargs)
        
    
###############
# BINNED TLEN #
###############


def stat_tlen_binned(Samfile samfile,
                     chrom=None, start=None, end=None, one_based=False,
                     window_size=300, window_offset=None, **kwargs):
    if window_offset is None:
        window_offset = window_size / 2
    if chrom is None:
        it = chain(*[_iter_tlen_binned(samfile, chrom, None, None, one_based, window_size, window_offset)
                     for chrom in sorted(samfile.references)])
    else:
        it = _iter_tlen_binned(samfile, chrom, start, end, one_based, window_size, window_offset)
    return it


def _iter_tlen_binned(Samfile samfile,
                          chrom, start, end, one_based,
                          int window_size, int window_offset):
    assert chrom is not None, 'unexpected error: chromosome is None'
    cdef int rtid, rstart, rend, has_coord, bin_start, bin_end
    cdef int reads_all = 0
    cdef int reads_pp = 0
    cdef bam1_t * b
    cdef uint32_t flag
    cdef bint is_unmapped
    cdef bint is_proper_pair
    cdef IteratorRowRegion it
    cdef int64_t tlen
    cdef int64_t tlen_squared
    cdef int64_t tlen_sum = 0
    cdef int64_t tlen_pp_sum = 0
    cdef int64_t tlen_squared_sum = 0
    cdef int64_t tlen_pp_squared_sum = 0
    start, end = normalise_coords(start, end, one_based)
    has_coord, rtid, rstart, rend = samfile._parseRegion(chrom, start, end, None)
    it = IteratorRowRegion(samfile, rtid, rstart, rend, reopen=False)
    b = it.b

    # setup first bin
    bin_start = rstart
    bin_end = bin_start + window_size

    # iterate over reads
    it.cnext()
    while it.retval > 0:
        while b.core.pos > bin_end: # end of bin, yield record
            # yield record for bin
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'reads_all': reads_all, 'reads_pp': reads_pp,
                   'mean_tlen': _mean(tlen_sum, reads_all),
                   'mean_tlen_pp': _mean(tlen_pp_sum, reads_pp),
                   'rms_tlen': rootmean(tlen_squared_sum, reads_all),
                   'rms_tlen_pp': rootmean(tlen_pp_squared_sum, reads_pp)}
            yield rec
            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size
            tlen_sum = tlen_squared_sum = tlen_pp_sum = tlen_pp_squared_sum = reads_all = reads_pp = 0
        # increment counters
        flag = b.core.flag
        is_unmapped = <bint>(flag & BAM_FUNMAP)
        if not is_unmapped:
            reads_all += 1
            tlen = b.core.isize
            tlen_sum += tlen
            tlen_squared = tlen**2
            tlen_squared_sum += tlen_squared
            is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
            if is_proper_pair:
                reads_pp += 1
                tlen_pp_sum += tlen
                tlen_pp_squared_sum += tlen_squared
        # move iterator on
        it.cnext()

    # deal with last non-empty bin
    pos = bin_start + window_offset
    if one_based:
        pos += 1
    rec = {'chrom': chrom, 'pos': pos,
           'reads_all': reads_all, 'reads_pp': reads_pp,
           'mean_tlen': _mean(tlen_sum, reads_all),
           'mean_tlen_pp': _mean(tlen_pp_sum, reads_pp),
           'rms_tlen': rootmean(tlen_squared_sum, reads_all),
           'rms_tlen_pp': rootmean(tlen_pp_squared_sum, reads_pp)}
    yield rec

    # deal with empty bins up to explicit end
    if end is not None:
        while bin_end < end:
            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size
            tlen_sum = tlen_squared_sum = tlen_pp_sum = tlen_pp_squared_sum = reads_all = reads_pp = 0
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec = {'chrom': chrom, 'pos': pos,
                   'reads_all': reads_all, 'reads_pp': reads_pp,
                   'mean_tlen': _mean(tlen_sum, reads_all),
                   'mean_tlen_pp': _mean(tlen_pp_sum, reads_pp),
                   'rms_tlen': rootmean(tlen_squared_sum, reads_all),
                   'rms_tlen_pp': rootmean(tlen_pp_squared_sum, reads_pp)}
            yield rec


def write_tlen_binned(*args, **kwargs):
    fieldnames = ('chrom', 'pos', 'reads_all', 'reads_pp', 'mean_tlen', 'mean_tlen_pp', 'rms_tlen', 'rms_tlen_pp')
    write_stats(stat_tlen_binned, fieldnames, *args, **kwargs)


def load_tlen_binned(*args, **kwargs):
    default_dtype = [('chrom', 'a12'),
                     ('pos', 'i4'),
                     ('reads_all', 'i4'),
                     ('reads_pp', 'i4'),
                     ('mean_tlen', 'i4'),
                     ('mean_tlen_pp', 'i4'),
                     ('rms_tlen', 'i4'),
                     ('rms_tlen_pp', 'i4'),
                    ]
    return load_stats(stat_tlen_binned, default_dtype, *args, **kwargs)


#####################
# UTILITY FUNCTIONS #
#####################  


def normalise_coords(start, end, one_based):
    if one_based:
        start = start - 1 if start is not None else None
        end = end - 1 if end is not None else None
    start = 0 if start is None else start
    return start, end

    
def write_stats(statfun, fieldnames, outfile, samfile, fafile=None,
                dialect=csv.excel_tab, write_header=True, 
                chrom=None, start=None, end=None, 
                one_based=False, progress=None, **kwargs):
    cdef long long counter = 0
    cdef long long modulus
    
    writer = csv.DictWriter(outfile, fieldnames, dialect=dialect)
    
    if write_header:
        writer.writeheader()

    if fafile is None:
        recs = statfun(samfile, chrom=chrom, start=start, end=end, one_based=one_based, **kwargs)
    else:
        recs = statfun(samfile, fafile, chrom=chrom, start=start, end=end, one_based=one_based, **kwargs)

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
    
    
from operator import itemgetter


def flatten(recs, fields):
    getter = itemgetter(*fields)
    rows = (getter(rec) for rec in recs)
    return rows


def flatten_one(recs, field):
    getter = itemgetter(field)
    items = (getter(rec) for rec in recs)
    return items


def load_stats(statfun, default_dtype, *args, **kwargs):
    try:
        fields = kwargs['fields']
        del kwargs['fields']
    except:
        fields = [t[0] for t in default_dtype]
    try:
        dtype_overrides = kwargs['dtype'] # expect dict
        del kwargs['dtype']
        dtype = dict(default_dtype)
        for k in dtype_overrides:
            dtype[k] = dtype_overrides[k]
    except:
        dtype = dict(default_dtype)
    recs = statfun(*args, **kwargs)
    if len(fields) == 1:
        f = fields[0]
        dtype = dtype[f]
        items = flatten_one(recs, f)
    else:
        # trim dtype to selected fields
        dtype = [(f, dtype[f]) for f in fields]
        items = flatten(recs, fields)
    a = np.fromiter(items, dtype=dtype)
    return a.view(np.recarray)
    
                
cdef inline bint is_softclipped(bam1_t * aln):
    cdef int k
    cigar_p = bam1_cigar(aln);
    for k in range(aln.core.n_cigar):
        op = cigar_p[k] & BAM_CIGAR_MASK
        if op == BAM_CSOFT_CLIP:
            return 1
    return 0


#cdef inline object cigar_counter(bam1_t * aln):
#    # TODO optimise here?
#    cdef int k
#    cigar = []
#    cigar_p = bam1_cigar(aln);
#    for k from 0 <= k < aln.core.n_cigar:
#        op = cigar_p[k] & BAM_CIGAR_MASK
#        l = cigar_p[k] >> BAM_CIGAR_SHIFT
#        cigar.append((op, l))
#    return Counter(dict(cigar))    


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


cdef inline int rootmean(uint64_t sqsum, int count):
    if count > 0:
        return int(round(sqrt(sqsum * 1. / count)))
    else:
        return 0
    
    
cdef inline int _mean(int64_t sum, int count):
    if count > 0:
        return int(round(sum * 1. / count))
    else:
        return 0


# SANDBOX


def count_reads(Samfile samfile, chrom=None, start=None, end=None):
    cdef IteratorRowRegion it
    cdef int n = 0
    has_coord, rtid, rstart, rend = samfile._parseRegion(chrom, start, end, None)
    it = IteratorRowRegion(samfile, rtid, rstart, rend, reopen=False)
    while True:
        it.cnext()
        if it.retval > 0:
            n += 1
        else:
            break
    return n
    
    
    
    
    
