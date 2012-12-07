# cython: profile=True


import sys
import numpy as np
#import numpy.ma as ma
cimport numpy as np
import time
import csv
from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t
from libc.math cimport sqrt
from cpython cimport PyBytes_FromStringAndSize
from csamtools cimport Samfile, Fastafile, PileupProxy, bam1_t, bam_pileup1_t, bam1_cigar, bam1_seq


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
    
    
# TODO variation by strand
    
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
    cdef uint64_t tlen_squared
    cdef uint64_t tlen_p_squared_sum = 0
    cdef uint64_t tlen_pp_squared_sum = 0
    
    # initialise variables
    n = col.n
    plp = col.plp
    tlen_p = []
    tlen_pp = []
    
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

        if is_proper_pair:
            reads_p += 1
            tlen_p.append(tlen)
            tlen_p_squared_sum += tlen_squared
            reads_pp += 1
            tlen_pp.append(tlen)
            tlen_pp_squared_sum += tlen_squared
        # N.B. insert size is only meaningful if mate is mapped to same chromosome
        elif not mate_is_unmapped and not mate_other_chr:
            reads_p += 1
            tlen_p.append(tlen)
            tlen_p_squared_sum += tlen_squared

    # calculate output variables
    rms_tlen = sqrt(tlen_p_squared_sum*1. / reads_p)
    std_tlen = np.std(np.array(tlen_p, dtype=np.int))
    rms_tlen_pp = sqrt(tlen_pp_squared_sum*1. / reads_pp)
    std_tlen_pp = np.std(np.array(tlen_pp, dtype=np.int))

    # round values to nearest integer, any finer precision is probably not
    # interesting    
    return {'chr': chrom, 
            'pos': pos, 
            'reads_all': n, 
            'reads_pp': reads_pp,
            'rms_tlen': int(round(rms_tlen)),
            'rms_tlen_pp': int(round(rms_tlen_pp)),
            'std_tlen': int(round(std_tlen)),
            'std_tlen_pp': int(round(std_tlen_pp))}


def stat_tlen(samfile, chrom=None, start=None, end=None, one_based=False):
    start, end = normalise_coords(start, end, one_based)
    for col in samfile.pileup(reference=chrom, start=start, end=end):
        yield construct_rec_tlen(samfile, col, one_based)
        
        
def write_tlen(*args, **kwargs):
    fieldnames = ('chr', 'pos', 'reads_all', 'reads_pp', 
                  'rms_tlen', 'rms_tlen_pp',
                  'std_tlen', 'std_tlen_pp')
    write_stats(stat_tlen, fieldnames, *args, **kwargs)
    
    
# TODO tlen stats by strand


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

