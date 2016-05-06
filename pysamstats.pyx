# cython: profile=False
# cython: embedsignature=True
from __future__ import print_function, division, absolute_import


__version__ = '0.24.3'


import sys as _sys
import itertools as _itertools
import time as _time
import csv as _csv


# PY2/3 compatibility
PY2 = _sys.version_info[0] == 2
if PY2:
    _string_types = basestring,
else:
    _string_types = str,


## These are bits set in the flag.
## have to put these definitions here, in csamtools.pxd they got ignored
## @abstract the read is paired in sequencing, no matter whether it is mapped
#  in a pair */
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


dtype_coverage = [('chrom', 'a12'),
                  ('pos', 'i4'),
                  ('reads_all', 'i4'),
                  ('reads_pp', 'i4')]


fields_coverage = [t[0] for t in dtype_coverage]


cpdef dict _rec_coverage(AlignmentFile alignmentfile, FastaFile fafile,
                         PileupColumn col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i  # loop index
    cdef int reads_all  # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef int reads_pp = 0

    # initialise variables
    n = col.n
    plp = col.plp

    # get chromosome name and position
    chrom = alignmentfile.getrname(col.tid)
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


cpdef dict _rec_coverage_pad(FastaFile fafile, chrom, pos,
                             bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom,
            'pos': pos,
            'reads_all': 0,
            'reads_pp': 0}


def stat_coverage(alignmentfile, **kwargs):
    """Generate coverage statistics per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_coverage, _rec_coverage_pad,
                        alignmentfile, **kwargs)


def load_coverage(*args, **kwargs):
    return _load_stats(stat_coverage, dtype_coverage, *args, **kwargs)
    
    
################################
# STRANDED COVERAGE STATISTICS #
################################


dtype_coverage_strand = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_fwd', 'i4'),
    ('reads_rev', 'i4'),
    ('reads_pp', 'i4'),
    ('reads_pp_fwd', 'i4'),
    ('reads_pp_rev', 'i4'),
]


fields_coverage_strand = [t[0] for t in dtype_coverage_strand]


cpdef dict _rec_coverage_strand(AlignmentFile alignmentfile, FastaFile fafile,
                                PileupColumn col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int reads_all # total number of reads in column
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
    chrom = alignmentfile.getrname(col.tid)
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


cpdef dict _rec_coverage_strand_pad(FastaFile fafile, chrom, pos,
                                    bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom,
            'pos': pos,
            'reads_all': 0,
            'reads_fwd': 0,
            'reads_rev': 0,
            'reads_pp': 0,
            'reads_pp_fwd': 0,
            'reads_pp_rev': 0}


def stat_coverage_strand(alignmentfile, **kwargs):
    """Generate coverage statistics by strand per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_coverage_strand,
                        _rec_coverage_strand_pad,
                        alignmentfile, **kwargs)


def load_coverage_strand(*args, **kwargs):
    return _load_stats(stat_coverage_strand, dtype_coverage_strand,
                      *args, **kwargs)
    
    
################################
# EXTENDED COVERAGE STATISTICS #
################################


dtype_coverage_ext = [
    ('chrom', 'a12'),
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


fields_coverage_ext = [t[0] for t in dtype_coverage_ext]


cpdef dict _rec_coverage_ext(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
                             bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i  # loop index
    cdef int reads_all  # total number of reads in column
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
    chrom = alignmentfile.getrname(tid)
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
        elif (is_reverse and mate_is_reverse) \
                or (not is_reverse and not mate_is_reverse):
            reads_mate_same_strand += 1
        elif (is_reverse and tlen > 0) or (not is_reverse and tlen < 0):
            reads_faceaway += 1
        if _is_softclipped(aln):
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


cpdef dict _rec_coverage_ext_pad(FastaFile fafile, chrom, pos,
                                 bint one_based=False):
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


def stat_coverage_ext(alignmentfile, **kwargs):
    """Generate extended coverage statistics per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_coverage_ext, _rec_coverage_ext_pad, alignmentfile,
                        **kwargs)


def load_coverage_ext(*args, **kwargs):
    return _load_stats(stat_coverage_ext, dtype_coverage_ext, *args, **kwargs)
    
    
##########################################
# EXTENDED COVERAGE STATISTICS BY STRAND #
##########################################


dtype_coverage_ext_strand = [
    ('chrom', 'a12'),
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


fields_coverage_ext_strand = [t[0] for t in dtype_coverage_ext_strand]


cpdef dict _rec_coverage_ext_strand(AlignmentFile alignmentfile, FastaFile fafile,
                                    PileupColumn col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int reads_all # total number of reads in column
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
    chrom = alignmentfile.getrname(tid)
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
        if _is_softclipped(aln):
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


cpdef dict _rec_coverage_ext_strand_pad(FastaFile fafile, chrom, pos,
                                        bint one_based=False):
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


def stat_coverage_ext_strand(alignmentfile, **kwargs):
    """Generate extended coverage statistics by strand per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_coverage_ext_strand,
                        _rec_coverage_ext_strand_pad, alignmentfile, **kwargs)


def load_coverage_ext_strand(*args, **kwargs):
    return _load_stats(stat_coverage_ext_strand, dtype_coverage_ext_strand,
                       *args, **kwargs)
    
    
########################
# VARIATION STATISTICS #
########################


dtype_variation = [
    ('chrom', 'a12'),
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


fields_variation = [t[0] for t in dtype_variation]


cpdef dict _rec_variation(AlignmentFile alignmentfile, FastaFile fafile,
                          PileupColumn col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i  # loop index
    cdef int reads_all  # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef bytes alnbase, refbase_b
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
    cdef int a = 0
    cdef int a_pp = 0
    cdef int c = 0
    cdef int c_pp = 0
    cdef int t = 0
    cdef int t_pp = 0
    cdef int g = 0
    cdef int g_pp = 0
    cdef int n = 0
    cdef int n_pp = 0
    
    # initialise variables
    reads_all = col.n
    plp = col.plp

    # get chromosome name and position
    chrom = alignmentfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # reference base
    refbase = fafile\
        .fetch(reference=chrom, start=col.pos, end=col.pos+1)\
        .upper()
    if not PY2:
        refbase_b = refbase.encode('ascii')
    else:
        refbase_b = refbase
    
    # loop over reads, extract what we need
    for i in range(reads_all):
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
            alnbase = _get_seq_base(aln, read.qpos)
            if alnbase == b'A':
                a += 1
                if is_proper_pair:
                    a_pp += 1
            elif alnbase == b'T':
                t += 1
                if is_proper_pair:
                    t_pp += 1
            elif alnbase == b'C':
                c += 1
                if is_proper_pair:
                    c_pp += 1
            elif alnbase == b'G':
                g += 1
                if is_proper_pair:
                    g_pp += 1
            elif alnbase == b'N':
                n += 1
                if is_proper_pair:
                    n_pp += 1
            if read.indel > 0:
                insertions += 1
                if is_proper_pair:
                    insertions_pp += 1
            if alnbase == refbase_b:
                matches += 1
                if is_proper_pair:
                    matches_pp += 1
            else:
                mismatches += 1
                if is_proper_pair:
                    mismatches_pp += 1

    return {'chrom': chrom, 'pos': pos, 'ref': refbase,
            'reads_all': reads_all, 'reads_pp': reads_pp,
            'matches': matches,
            'matches_pp': matches_pp,
            'mismatches': mismatches,
            'mismatches_pp': mismatches_pp,
            'deletions': deletions,
            'deletions_pp': deletions_pp,
            'insertions': insertions,
            'insertions_pp': insertions_pp,
            'A': a, 'A_pp': a_pp,
            'C': c, 'C_pp': c_pp,
            'T': t, 'T_pp': t_pp,
            'G': g, 'G_pp': g_pp,
            'N': n, 'N_pp': n_pp}


cpdef dict _rec_variation_pad(FastaFile fafile, chrom, pos,
                              bint one_based=False):
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


def stat_variation(alignmentfile, fafile, **kwargs):
    """Generate variation statistics per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    fafile : pysam.FastaFile or string
        FASTA file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_variation, _rec_variation_pad,
                        alignmentfile, fafile=fafile, **kwargs)


def load_variation(*args, **kwargs):
    return _load_stats(stat_variation, dtype_variation, *args, **kwargs)

    
#################################
# STRANDED VARIATION STATISTICS #
#################################


dtype_variation_strand = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('ref', 'a1'),
    ('reads_all', 'i4'),
    ('reads_fwd', 'i4'),
    ('reads_rev', 'i4'),
    ('reads_pp', 'i4'),
    ('reads_pp_fwd', 'i4'),
    ('reads_pp_rev', 'i4'),
    ('matches', 'i4'),
    ('matches_fwd', 'i4'),
    ('matches_rev', 'i4'),
    ('matches_pp', 'i4'),
    ('matches_pp_fwd', 'i4'),
    ('matches_pp_rev', 'i4'),
    ('mismatches', 'i4'),
    ('mismatches_fwd', 'i4'),
    ('mismatches_rev', 'i4'),
    ('mismatches_pp', 'i4'),
    ('mismatches_pp_fwd', 'i4'),
    ('mismatches_pp_rev', 'i4'),
    ('deletions', 'i4'),
    ('deletions_fwd', 'i4'),
    ('deletions_rev', 'i4'),
    ('deletions_pp', 'i4'),
    ('deletions_pp_fwd', 'i4'),
    ('deletions_pp_rev', 'i4'),
    ('insertions', 'i4'),
    ('insertions_fwd', 'i4'),
    ('insertions_rev', 'i4'),
    ('insertions_pp', 'i4'),
    ('insertions_pp_fwd', 'i4'),
    ('insertions_pp_rev', 'i4'),
    ('A', 'i4'), ('A_fwd', 'i4'), ('A_rev', 'i4'),
    ('A_pp', 'i4'), ('A_pp_fwd', 'i4'), ('A_pp_rev', 'i4'),
    ('C', 'i4'), ('C_fwd', 'i4'), ('C_rev', 'i4'),
    ('C_pp', 'i4'), ('C_pp_fwd', 'i4'), ('C_pp_rev', 'i4'),
    ('T', 'i4'), ('T_fwd', 'i4'), ('T_rev', 'i4'),
    ('T_pp', 'i4'), ('T_pp_fwd', 'i4'), ('T_pp_rev', 'i4'),
    ('G', 'i4'), ('G_fwd', 'i4'), ('G_rev', 'i4'),
    ('G_pp', 'i4'), ('G_pp_fwd', 'i4'), ('G_pp_rev', 'i4'),
    ('N', 'i4'), ('N_fwd', 'i4'), ('N_rev', 'i4'),
    ('N_pp', 'i4'), ('N_pp_fwd', 'i4'), ('N_pp_rev', 'i4')
]


fields_variation_strand = [t[0] for t in dtype_variation_strand]


cdef struct _CountPpStrand:
    int all, pp, fwd, rev, pp_fwd, pp_rev


cdef inline _init_pp_strand(_CountPpStrand* c):
    c.all = c.fwd = c.rev = c.pp = c.pp_fwd = c.pp_rev = 0    


cdef inline _incr_pp_strand(_CountPpStrand* c, bint is_reverse,
                           bint is_proper_pair):
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
                

cpdef dict _rec_variation_strand(AlignmentFile alignmentfile, FastaFile fafile,
                                 PileupColumn col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int reads_all # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair, is_reverse
    cdef bytes alnbase, refbase_b
    # counting variables
    cdef _CountPpStrand reads, matches, mismatches, deletions, insertions, \
        A, C, T, G, N
    
    # initialise variables
    n = col.n
    plp = col.plp
    _init_pp_strand(&reads)
    _init_pp_strand(&matches)
    _init_pp_strand(&mismatches)
    _init_pp_strand(&deletions)
    _init_pp_strand(&insertions)
    _init_pp_strand(&A)
    _init_pp_strand(&T)
    _init_pp_strand(&C)
    _init_pp_strand(&G)
    _init_pp_strand(&N)

    # get chromosome name and position
    chrom = alignmentfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # reference base
    refbase = fafile.fetch(reference=chrom, start=col.pos, end=col.pos+1).upper()
    if not PY2:
        refbase_b = refbase.encode('ascii')
    else:
        refbase_b = refbase

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
        _incr_pp_strand(&reads, is_reverse, is_proper_pair)
        if read.is_del:
            _incr_pp_strand(&deletions, is_reverse, is_proper_pair)
        else:
            alnbase = _get_seq_base(aln, read.qpos)
            if alnbase == b'A':
                _incr_pp_strand(&A, is_reverse, is_proper_pair)
            elif alnbase == b'T':
                _incr_pp_strand(&T, is_reverse, is_proper_pair)
            elif alnbase == b'C':
                _incr_pp_strand(&C, is_reverse, is_proper_pair)
            elif alnbase == b'G':
                _incr_pp_strand(&G, is_reverse, is_proper_pair)
            elif alnbase == b'N':
                _incr_pp_strand(&N, is_reverse, is_proper_pair)
            if read.indel > 0:
                _incr_pp_strand(&insertions, is_reverse, is_proper_pair)
            if alnbase == refbase_b:
                _incr_pp_strand(&matches, is_reverse, is_proper_pair)
            else:
                _incr_pp_strand(&mismatches, is_reverse, is_proper_pair)

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


cpdef dict _rec_variation_strand_pad(FastaFile fafile, chrom, pos,
                                     bint one_based=False):
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
            'A': 0, 'A_fwd': 0, 'A_rev': 0,
            'A_pp': 0, 'A_pp_fwd': 0, 'A_pp_rev': 0,
            'C': 0, 'C_fwd': 0, 'C_rev': 0,
            'C_pp': 0, 'C_pp_fwd': 0, 'C_pp_rev': 0,
            'T': 0, 'T_fwd': 0, 'T_rev': 0,
            'T_pp': 0, 'T_pp_fwd': 0, 'T_pp_rev': 0,
            'G': 0, 'G_fwd': 0, 'G_rev': 0,
            'G_pp': 0, 'G_pp_fwd': 0, 'G_pp_rev': 0,
            'N': 0, 'N_fwd': 0, 'N_rev': 0,
            'N_pp': 0, 'N_pp_fwd': 0, 'N_pp_rev': 0,
            }


def stat_variation_strand(alignmentfile, fafile, **kwargs):
    """Generate variation statistics by strand per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    fafile : pysam.FastaFile or string
        FASTA file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_variation_strand,
                        _rec_variation_strand_pad,
                        alignmentfile, fafile=fafile, **kwargs)


def load_variation_strand(*args, **kwargs):
    return _load_stats(stat_variation_strand, dtype_variation_strand,
                       *args, **kwargs)
    

##########################
# INSERT SIZE STATISTICS #
##########################


dtype_tlen = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_paired', 'i4'),
    ('reads_pp', 'i4'),
    ('mean_tlen', 'i4'),
    ('mean_tlen_pp', 'i4'),
    ('rms_tlen', 'i4'),
    ('rms_tlen_pp', 'i4'),
    ('std_tlen', 'i4'),
    ('std_tlen_pp', 'i4')
]


fields_tlen = [t[0] for t in dtype_tlen]


cpdef dict _rec_tlen(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
                     bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i  # loop index
    cdef int reads_all  # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef bint mate_is_unmappped 
    cdef bint mate_other_chr
    # reads "paired", i.e., mate is mapped to same chromosome, so tlen is
    # meaningful
    cdef int reads_p = 0
    # reads "properly paired", as defined by aligner
    cdef int reads_pp = 0
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
    chrom = alignmentfile.getrname(col.tid)
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
        tlen_p_mean = tlen_p_sum / reads_p
    if reads_pp > 0:
        tlen_pp_mean = tlen_pp_sum / reads_pp
        
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
        rms_tlen = _rootmean(tlen_p_squared_sum, reads_p)
        variance_tlen = tlen_p_dev_squared_sum / reads_p
        std_tlen = int(round(sqrt(variance_tlen)))
    else:
        rms_tlen = std_tlen = mean_tlen = median_tlen = 0
    if reads_pp > 0:
        mean_tlen_pp = int(round(tlen_pp_mean))
        rms_tlen_pp = _rootmean(tlen_pp_squared_sum, reads_pp)
        variance_tlen_pp = tlen_pp_dev_squared_sum / reads_pp
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


cpdef dict _rec_tlen_pad(FastaFile fafile, chrom, pos, bint one_based=False):
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


def stat_tlen(alignmentfile, **kwargs):
    """Generate insert size statistics per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_tlen, _rec_tlen_pad, alignmentfile, **kwargs)


def load_tlen(*args, **kwargs):
    return _load_stats(stat_tlen, dtype_tlen, *args, **kwargs)
    
    
####################################
# INSERT SIZE STATISTICS BY STRAND #
####################################


dtype_tlen_strand = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_fwd', 'i4'),
    ('reads_rev', 'i4'),
    ('reads_paired', 'i4'),
    ('reads_paired_fwd', 'i4'),
    ('reads_paired_rev', 'i4'),
    ('reads_pp', 'i4'),
    ('reads_pp_fwd', 'i4'),
    ('reads_pp_rev', 'i4'),
    ('mean_tlen', 'i4'),
    ('mean_tlen_fwd', 'i4'),
    ('mean_tlen_rev', 'i4'),
    ('mean_tlen_pp', 'i4'),
    ('mean_tlen_pp_fwd', 'i4'),
    ('mean_tlen_pp_rev', 'i4'),
    ('rms_tlen', 'i4'),
    ('rms_tlen_fwd', 'i4'),
    ('rms_tlen_rev', 'i4'),
    ('rms_tlen_pp', 'i4'),
    ('rms_tlen_pp_fwd', 'i4'),
    ('rms_tlen_pp_rev', 'i4'),
    ('std_tlen', 'i4'),
    ('std_tlen_fwd', 'i4'),
    ('std_tlen_rev', 'i4'),
    ('std_tlen_pp', 'i4'),
    ('std_tlen_pp_fwd', 'i4'),
    ('std_tlen_pp_rev', 'i4')
]


fields_tlen_strand = [t[0] for t in dtype_tlen_strand]


cpdef dict _rec_tlen_strand(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
                            bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int reads_all # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef bint mate_is_unmappped 
    cdef bint mate_other_chr
    
    # counting variables
    cdef int reads_fwd = 0
    cdef int reads_rev = 0
    # reads "paired", i.e., mate is mapped to same chromosome, so tlen is
    # meaningful
    cdef int reads_p = 0
    cdef int reads_p_fwd = 0
    cdef int reads_p_rev = 0
    # reads "properly paired", as defined by aligner
    cdef int reads_pp = 0
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
    chrom = alignmentfile.getrname(col.tid)
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

        # N.B. insert size is only meaningful if mate is mapped to same
        # chromosome
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
        tlen_p_mean = tlen_p_sum / reads_p
        if reads_p_rev > 0:
            tlen_p_rev_mean = tlen_p_rev_sum / reads_p_rev
        if reads_p_fwd > 0:
            tlen_p_fwd_mean = tlen_p_fwd_sum / reads_p_fwd
    if reads_pp > 0:
        tlen_pp_mean = tlen_pp_sum / reads_pp
        if reads_pp_rev > 0:
            tlen_pp_rev_mean = tlen_pp_rev_sum / reads_pp_rev
        if reads_pp_fwd > 0:
            tlen_pp_fwd_mean = tlen_pp_fwd_sum / reads_pp_fwd
        
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
        # N.B. insert size is only meaningful if mate is mapped to same
        # chromosome
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
        rms_tlen = _rootmean(tlen_p_squared_sum, reads_p)
        variance_tlen = tlen_p_dev_squared_sum / reads_p
        std_tlen = int(round(sqrt(variance_tlen)))
    else:
        rms_tlen = std_tlen = mean_tlen = 0
    if reads_p_rev > 0:
        mean_tlen_rev = int(round(tlen_p_rev_mean))
        rms_tlen_rev = _rootmean(tlen_p_rev_squared_sum, reads_p_rev)
        variance_tlen_rev = tlen_p_rev_dev_squared_sum / reads_p_rev
        std_tlen_rev = int(round(sqrt(variance_tlen_rev)))
    else:
        rms_tlen_rev = std_tlen_rev = mean_tlen_rev = 0
    if reads_p_fwd > 0:
        mean_tlen_fwd = int(round(tlen_p_fwd_mean))
        rms_tlen_fwd = _rootmean(tlen_p_fwd_squared_sum, reads_p_fwd)
        variance_tlen_fwd = tlen_p_fwd_dev_squared_sum / reads_p_fwd
        std_tlen_fwd = int(round(sqrt(variance_tlen_fwd)))
    else:
        rms_tlen_fwd = std_tlen_fwd = mean_tlen_fwd = 0
    if reads_pp > 0:
        mean_tlen_pp = int(round(tlen_pp_mean))
        rms_tlen_pp = _rootmean(tlen_pp_squared_sum, reads_pp)
        variance_tlen_pp = tlen_pp_dev_squared_sum / reads_pp
        std_tlen_pp = int(round(sqrt(variance_tlen_pp)))
    else:
        rms_tlen_pp = std_tlen_pp = mean_tlen_pp = 0
    if reads_pp_rev > 0:
        mean_tlen_pp_rev = int(round(tlen_pp_rev_mean))
        rms_tlen_pp_rev = _rootmean(tlen_pp_rev_squared_sum, reads_pp_rev)
        variance_tlen_pp_rev = tlen_pp_rev_dev_squared_sum / reads_pp_rev
        std_tlen_pp_rev = int(round(sqrt(variance_tlen_pp_rev)))
    else:
        rms_tlen_pp_rev = std_tlen_pp_rev = mean_tlen_pp_rev = 0
    if reads_pp_fwd > 0:
        mean_tlen_pp_fwd = int(round(tlen_pp_fwd_mean))
        rms_tlen_pp_fwd = _rootmean(tlen_pp_fwd_squared_sum, reads_pp_fwd)
        variance_tlen_pp_fwd = tlen_pp_fwd_dev_squared_sum / reads_pp_fwd
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


cpdef dict _rec_tlen_strand_pad(FastaFile fafile, chrom, pos,
                                bint one_based=False):
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


def stat_tlen_strand(alignmentfile, **kwargs):
    """Generate insert size statistics by strand per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_tlen_strand, _rec_tlen_strand_pad, alignmentfile,
                        **kwargs)


def load_tlen_strand(*args, **kwargs):
    return _load_stats(stat_tlen_strand, dtype_tlen_strand, *args, **kwargs)
        
    
##############################
# MAPPING QUALITY STATISTICS #
##############################


dtype_mapq = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_pp', 'i4'),
    ('reads_mapq0', 'i4'),
    ('reads_mapq0_pp', 'i4'),
    ('rms_mapq', 'i4'),
    ('rms_mapq_pp', 'i4'),
    ('max_mapq', 'i4'),
    ('max_mapq_pp', 'i4')
]


fields_mapq = [t[0] for t in dtype_mapq]


cpdef dict _rec_mapq(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
                     bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int reads_all # total number of reads in column
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
    chrom = alignmentfile.getrname(col.tid)
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
    rms_mapq = _rootmean(mapq_squared_sum, n)
    max_mapq = mapq_max
    if reads_pp > 0:
        rms_mapq_pp = _rootmean(mapq_pp_squared_sum, reads_pp)
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


cpdef dict _rec_mapq_pad(FastaFile fafile, chrom, pos, bint one_based=False):
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


def stat_mapq(alignmentfile, **kwargs):
    """Generate mapping quality statistics per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_mapq, _rec_mapq_pad, alignmentfile, **kwargs)


def load_mapq(*args, **kwargs):
    return _load_stats(stat_mapq, dtype_mapq, *args, **kwargs)
        
    
########################################
# MAPPING QUALITY STATISTICS BY STRAND #
########################################


dtype_mapq_strand = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_fwd', 'i4'),
    ('reads_rev', 'i4'),
    ('reads_pp', 'i4'),
    ('reads_pp_fwd', 'i4'),
    ('reads_pp_rev', 'i4'),
    ('reads_mapq0', 'i4'),
    ('reads_mapq0_fwd', 'i4'),
    ('reads_mapq0_rev', 'i4'),
    ('reads_mapq0_pp', 'i4'),
    ('reads_mapq0_pp_fwd', 'i4'),
    ('reads_mapq0_pp_rev', 'i4'),
    ('rms_mapq', 'i4'),
    ('rms_mapq_fwd', 'i4'),
    ('rms_mapq_rev', 'i4'),
    ('rms_mapq_pp', 'i4'),
    ('rms_mapq_pp_fwd', 'i4'),
    ('rms_mapq_pp_rev', 'i4'),
    ('max_mapq', 'i4'),
    ('max_mapq_fwd', 'i4'),
    ('max_mapq_rev', 'i4'),
    ('max_mapq_pp', 'i4'),
    ('max_mapq_pp_fwd', 'i4'),
    ('max_mapq_pp_rev', 'i4'),
]


fields_mapq_strand = [t[0] for t in dtype_mapq_strand]


cpdef dict _rec_mapq_strand(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
                            bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i  # loop index
    cdef int reads_all  # total number of reads in column

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
    chrom = alignmentfile.getrname(col.tid)
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
    rms_mapq = _rootmean(mapq_squared_sum, n)
    max_mapq = mapq_max
    if reads_rev > 0:
        rms_mapq_rev = _rootmean(mapq_rev_squared_sum, reads_rev)
        max_mapq_rev = mapq_rev_max
    else:
        rms_mapq_rev = max_mapq_rev = 0
    if reads_fwd > 0:
        rms_mapq_fwd = _rootmean(mapq_fwd_squared_sum, reads_fwd)
        max_mapq_fwd = mapq_fwd_max
    else:
        rms_mapq_fwd = max_mapq_fwd = 0
    if reads_pp > 0:
        rms_mapq_pp = _rootmean(mapq_pp_squared_sum, reads_pp)
        max_mapq_pp = mapq_pp_max
    else:
        rms_mapq_pp = max_mapq_pp = 0
    if reads_pp_fwd > 0:
        rms_mapq_pp_fwd = _rootmean(mapq_pp_fwd_squared_sum, reads_pp_fwd)
        max_mapq_pp_fwd = mapq_pp_fwd_max
    else:
        rms_mapq_pp_fwd = max_mapq_pp_fwd = 0
    if reads_pp_rev > 0:
        rms_mapq_pp_rev = _rootmean(mapq_pp_rev_squared_sum, reads_pp_rev)
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


cpdef dict _rec_mapq_strand_pad(FastaFile fafile, chrom, pos,
                                bint one_based=False):
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


def stat_mapq_strand(alignmentfile, **kwargs):
    """Generate mapping quality statistics by strand per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_mapq_strand, _rec_mapq_strand_pad, alignmentfile,
                        **kwargs)


def load_mapq_strand(*args, **kwargs):
    return _load_stats(stat_mapq_strand, dtype_mapq_strand, *args, **kwargs)
        
    
###########################
# BASE QUALITY STATISTICS #
###########################


dtype_baseq = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_pp', 'i4'),
    ('rms_baseq', 'i4'),
    ('rms_baseq_pp', 'i4'),
]


fields_baseq = [t[0] for t in dtype_baseq]


cpdef dict _rec_baseq(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
                      bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int reads_all # total number of reads in column
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
    chrom = alignmentfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # loop over reads, extract what we need
    for i in range(n):
        read = &(plp[0][i])
        aln = read.b
        flag = aln.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        if is_proper_pair:
            reads_pp += 1
        # N.B., base quality only makes sense if the aligned read is not a
        # deletion
        if not read.is_del:
            reads_nodel += 1
            baseq = pysam_bam_get_qual(aln)[read.qpos]
            baseq_squared = baseq**2
            baseq_squared_sum += baseq_squared
            if is_proper_pair:
                reads_pp_nodel += 1
                baseq_pp_squared_sum += baseq_squared

    # output variables
    rms_baseq = _rootmean(baseq_squared_sum, reads_nodel)
    rms_baseq_pp = _rootmean(baseq_pp_squared_sum, reads_pp_nodel)

    return {'chrom': chrom,
            'pos': pos, 
            'reads_all': n, 
            'reads_pp': reads_pp,
            'rms_baseq': rms_baseq,
            'rms_baseq_pp': rms_baseq_pp}


cpdef dict _rec_baseq_pad(FastaFile fafile, chrom, pos, bint one_based=False):
    pos = pos + 1 if one_based else pos
    return {'chrom': chrom, 
            'pos': pos, 
            'reads_all': 0, 
            'reads_pp': 0,
            'rms_baseq': 0,
            'rms_baseq_pp': 0,
            }


def stat_baseq(alignmentfile, **kwargs):
    """Generate base quality statistics per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_baseq, _rec_baseq_pad, alignmentfile, **kwargs)


def load_baseq(*args, **kwargs):
    return _load_stats(stat_baseq, dtype_baseq, *args, **kwargs)
        
    
#####################################
# BASE QUALITY STATISTICS BY STRAND #
#####################################


dtype_baseq_strand = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_fwd', 'i4'),
    ('reads_rev', 'i4'),
    ('reads_pp', 'i4'),
    ('reads_pp_fwd', 'i4'),
    ('reads_pp_rev', 'i4'),
    ('rms_baseq', 'i4'),
    ('rms_baseq_fwd', 'i4'),
    ('rms_baseq_rev', 'i4'),
    ('rms_baseq_pp', 'i4'),
    ('rms_baseq_pp_fwd', 'i4'),
    ('rms_baseq_pp_rev', 'i4'),
]


fields_baseq_strand = [t[0] for t in dtype_baseq_strand]


cpdef dict _rec_baseq_strand(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
                             bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int reads_all # total number of reads in column

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
    chrom = alignmentfile.getrname(col.tid)
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
            baseq = pysam_bam_get_qual(aln)[read.qpos]
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
    rms_baseq = _rootmean(baseq_squared_sum, reads_nodel)
    rms_baseq_rev = _rootmean(baseq_rev_squared_sum, reads_rev_nodel)
    rms_baseq_fwd = _rootmean(baseq_fwd_squared_sum, reads_fwd_nodel)
    rms_baseq_pp = _rootmean(baseq_pp_squared_sum, reads_pp_nodel)
    rms_baseq_pp_fwd = _rootmean(baseq_pp_fwd_squared_sum, reads_pp_fwd_nodel)
    rms_baseq_pp_rev = _rootmean(baseq_pp_rev_squared_sum, reads_pp_rev_nodel)
        
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


cpdef dict _rec_baseq_strand_pad(FastaFile fafile, chrom, pos,
                                 bint one_based=False):
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


def stat_baseq_strand(alignmentfile, **kwargs):
    """Generate base quality statistics by strand per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_baseq_strand, _rec_baseq_strand_pad,
                        alignmentfile, **kwargs)


def load_baseq_strand(*args, **kwargs):
    return _load_stats(stat_baseq_strand, dtype_baseq_strand,
                       *args, **kwargs)
        
    
####################################
# EXTENDED BASE QUALITY STATISTICS #
####################################


dtype_baseq_ext = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('ref', 'a1'),
    ('reads_all', 'i4'),
    ('reads_pp', 'i4'),
    ('matches', 'i4'),
    ('matches_pp', 'i4'),
    ('mismatches', 'i4'),
    ('mismatches_pp', 'i4'),
    ('rms_baseq', 'i4'),
    ('rms_baseq_pp', 'i4'),
    ('rms_baseq_matches', 'i4'),
    ('rms_baseq_matches_pp', 'i4'),
    ('rms_baseq_mismatches', 'i4'),
    ('rms_baseq_mismatches_pp', 'i4'),
]


fields_baseq_ext = [t[0] for t in dtype_baseq_ext]


cpdef dict _rec_baseq_ext(AlignmentFile alignmentfile, FastaFile fafile,
                          PileupColumn col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int reads_all # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef bytes refbase_b, alnbase
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
    chrom = alignmentfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # reference base
    refbase = fafile.fetch(reference=chrom, start=col.pos,
                           end=col.pos + 1).upper()
    if not PY2:
        refbase_b = refbase.encode('ascii')
    else:
        refbase_b = refbase

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
            baseq = pysam_bam_get_qual(aln)[read.qpos]
            baseq_squared = baseq**2
            baseq_squared_sum += baseq_squared
            if is_proper_pair:
                reads_pp_nodel += 1
                baseq_pp_squared_sum += baseq_squared
            alnbase = _get_seq_base(aln, read.qpos)
            if alnbase == refbase_b:
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
    rms_baseq = _rootmean(baseq_squared_sum, reads_nodel)
    rms_baseq_pp = _rootmean(baseq_pp_squared_sum, reads_pp_nodel)
    rms_baseq_matches = _rootmean(baseq_matches_squared_sum, matches)
    rms_baseq_matches_pp = _rootmean(baseq_matches_pp_squared_sum, matches_pp)
    rms_baseq_mismatches = _rootmean(baseq_mismatches_squared_sum, mismatches)
    rms_baseq_mismatches_pp = _rootmean(baseq_mismatches_pp_squared_sum,
                                        mismatches_pp)

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


cpdef dict _rec_baseq_ext_pad(FastaFile fafile, chrom, pos,
                              bint one_based=False):
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


def stat_baseq_ext(alignmentfile, fafile, **kwargs):
    """Generate extended base quality statistics per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    fafile : pysam.FastaFile or string
        FASTA file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_baseq_ext, _rec_baseq_ext_pad, alignmentfile,
                        fafile=fafile, **kwargs)


def load_baseq_ext(*args, **kwargs):
    return _load_stats(stat_baseq_ext, dtype_baseq_ext, *args, **kwargs)
        
    
##############################################
# EXTENDED BASE QUALITY STATISTICS BY STRAND #
##############################################


dtype_baseq_ext_strand = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('ref', 'a1'),
    ('reads_all', 'i4'),
    ('reads_fwd', 'i4'),
    ('reads_rev', 'i4'),
    ('reads_pp', 'i4'),
    ('reads_pp_fwd', 'i4'),
    ('reads_pp_rev', 'i4'),
    ('matches', 'i4'),
    ('matches_fwd', 'i4'),
    ('matches_rev', 'i4'),
    ('matches_pp', 'i4'),
    ('matches_pp_fwd', 'i4'),
    ('matches_pp_rev', 'i4'),
    ('mismatches', 'i4'),
    ('mismatches_fwd', 'i4'),
    ('mismatches_rev', 'i4'),
    ('mismatches_pp', 'i4'),
    ('mismatches_pp_fwd', 'i4'),
    ('mismatches_pp_rev', 'i4'),
    ('rms_baseq', 'i4'),
    ('rms_baseq_fwd', 'i4'),
    ('rms_baseq_rev', 'i4'),
    ('rms_baseq_pp', 'i4'),
    ('rms_baseq_pp_fwd', 'i4'),
    ('rms_baseq_pp_rev', 'i4'),
    ('rms_baseq_matches', 'i4'),
    ('rms_baseq_matches_fwd', 'i4'),
    ('rms_baseq_matches_rev', 'i4'),
    ('rms_baseq_matches_pp', 'i4'),
    ('rms_baseq_matches_pp_fwd', 'i4'),
    ('rms_baseq_matches_pp_rev', 'i4'),
    ('rms_baseq_mismatches', 'i4'),
    ('rms_baseq_mismatches_fwd', 'i4'),
    ('rms_baseq_mismatches_rev', 'i4'),
    ('rms_baseq_mismatches_pp', 'i4'),
    ('rms_baseq_mismatches_pp_fwd', 'i4'),
    ('rms_baseq_mismatches_pp_rev', 'i4')
]


fields_baseq_ext_strand = [t[0] for t in dtype_baseq_ext_strand]


cpdef dict _rec_baseq_ext_strand(AlignmentFile alignmentfile, FastaFile fafile,
                                 PileupColumn col, bint one_based=False):

    # statically typed variables
    cdef bam_pileup1_t ** plp
    cdef bam_pileup1_t * read
    cdef bam1_t * aln
    cdef int i # loop index
    cdef int reads_all # total number of reads in column
    cdef uint32_t flag
    cdef bint is_proper_pair
    cdef bint is_reverse
    cdef bytes alnbase, refbase_b
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
    chrom = alignmentfile.getrname(col.tid)
    pos = col.pos + 1 if one_based else col.pos
    
    # reference base
    refbase = fafile.fetch(reference=chrom, start=col.pos, end=col.pos+1).upper()
    if not PY2:
        refbase_b = refbase.encode('ascii')
    else:
        refbase_b = refbase

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
            baseq = pysam_bam_get_qual(aln)[read.qpos]
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
            alnbase = _get_seq_base(aln, read.qpos)
            if alnbase == refbase_b:
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
    rms_baseq = _rootmean(baseq_squared_sum, reads_nodel)
    rms_baseq_fwd = _rootmean(baseq_fwd_squared_sum, reads_fwd_nodel)
    rms_baseq_rev = _rootmean(baseq_rev_squared_sum, reads_rev_nodel)
    rms_baseq_pp = _rootmean(baseq_pp_squared_sum, reads_pp_nodel)
    rms_baseq_pp_fwd = _rootmean(baseq_pp_fwd_squared_sum, reads_pp_fwd_nodel)
    rms_baseq_pp_rev = _rootmean(baseq_pp_rev_squared_sum, reads_pp_rev_nodel)
    rms_baseq_matches = _rootmean(baseq_matches_squared_sum, matches)
    rms_baseq_matches_fwd = _rootmean(baseq_matches_fwd_squared_sum,
                                      matches_fwd)
    rms_baseq_matches_rev = _rootmean(baseq_matches_rev_squared_sum,
                                      matches_rev)
    rms_baseq_matches_pp = _rootmean(baseq_matches_pp_squared_sum, matches_pp)
    rms_baseq_matches_pp_fwd = _rootmean(baseq_matches_pp_fwd_squared_sum,
                                         matches_pp_fwd)
    rms_baseq_matches_pp_rev = _rootmean(baseq_matches_pp_rev_squared_sum,
                                         matches_pp_rev)
    rms_baseq_mismatches = _rootmean(baseq_mismatches_squared_sum, mismatches)
    rms_baseq_mismatches_fwd = _rootmean(baseq_mismatches_fwd_squared_sum,
                                         mismatches_fwd)
    rms_baseq_mismatches_rev = _rootmean(baseq_mismatches_rev_squared_sum,
                                         mismatches_rev)
    rms_baseq_mismatches_pp = _rootmean(baseq_mismatches_pp_squared_sum,
                                        mismatches_pp)
    rms_baseq_mismatches_pp_fwd = _rootmean(baseq_mismatches_pp_fwd_squared_sum,
                                            mismatches_pp_fwd)
    rms_baseq_mismatches_pp_rev = _rootmean(baseq_mismatches_pp_rev_squared_sum,
                                            mismatches_pp_rev)

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


cpdef dict _rec_baseq_ext_strand_pad(FastaFile fafile, chrom, pos,
                                     bint one_based=False):
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


def stat_baseq_ext_strand(alignmentfile, fafile, **kwargs):
    """Generate extended base quality statistics by strand per genome position.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    fafile : pysam.FastaFile or string
        FASTA file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column

    Returns
    -------

    recs : iterator
        record generator

    """

    return _iter_pileup(_rec_baseq_ext_strand, _rec_baseq_ext_strand_pad,
                        alignmentfile, fafile=fafile, **kwargs)


def load_baseq_ext_strand(*args, **kwargs):
    return _load_stats(stat_baseq_ext_strand, dtype_baseq_ext_strand,
                       *args, **kwargs)

    
#################################################
# BASIC COVERAGE STATISTICS WITH GC COMPOSITION #
#################################################


dtype_coverage_gc = [('chrom', 'a12'),
                     ('pos', 'i4'),
                     ('gc', 'u1'),
                     ('reads_all', 'i4'),
                     ('reads_pp', 'i4')]


fields_coverage_gc = [t[0] for t in dtype_coverage_gc]


def stat_coverage_gc(alignmentfile, fafile, chrom=None, start=None, end=None,
                     one_based=False, truncate=False, pad=False, max_depth=8000,
                     window_size=300, window_offset=None, **kwargs):
    """Generate coverage statistics per genome position with reference genome
    %GC composition in surrounding window.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    fafile : pysam.FastaFile or string
        FASTA file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column
    window_size : int
        size of surrounding window in base pairs
    window_offset :
        distance from window start to record position

    Returns
    -------

    recs : iterator
        record generator

    """


    if window_offset is None:
        window_offset = window_size / 2
        
    def _rec_coverage_gc(AlignmentFile alignmentfile, FastaFile fafile,
                         PileupColumn col, bint one_based):
        chrom = alignmentfile.getrname(col.tid)

        ref_window_start = col.pos - window_offset
        ref_window_end = ref_window_start + window_size
        if ref_window_start < 0:
            ref_window_start = 0
        ref_window = fafile\
            .fetch(chrom, ref_window_start, ref_window_end)\
            .lower()
        if len(ref_window) == 0:
            gc_percent = -1
        else:
            gc_percent = _gc_content(ref_window)

        rec = _rec_coverage(alignmentfile, fafile, col, one_based)
        rec['gc'] = gc_percent
        return rec

    def _rec_coverage_gc_pad(FastaFile fafile, chrom, pos,
                             bint one_based):
        ref_window_start = pos - window_offset
        ref_window_end = ref_window_start + window_size
        if ref_window_start < 0:
            ref_window_start = 0
        ref_window = fafile\
            .fetch(chrom, ref_window_start, ref_window_end)\
            .lower()
        if len(ref_window) == 0:
            gc_percent = -1
        else:
            gc_percent = _gc_content(ref_window)
        rec = _rec_coverage_pad(fafile, chrom, pos, one_based)
        rec['gc'] = gc_percent
        return rec

    return _iter_pileup(_rec_coverage_gc,
                        _rec_coverage_gc_pad,
                        alignmentfile, fafile=fafile, chrom=chrom, start=start,
                        end=end, one_based=one_based,
                        truncate=truncate, pad=pad,
                        max_depth=max_depth,
                        **kwargs)

        
def load_coverage_gc(*args, **kwargs):
    return _load_stats(stat_coverage_gc, dtype_coverage_gc, *args, **kwargs)


###################
# BINNED COVERAGE #
###################


dtype_coverage_binned = [('chrom', 'a12'),
                         ('pos', 'i4'),
                         ('gc', 'u1'),
                         ('reads_all', 'i4'),
                         ('reads_pp', 'i4')]


fields_coverage_binned = [t[0] for t in dtype_coverage_binned]


def stat_coverage_binned(alignmentfile, fafile, **kwargs):
    """Generate binned coverage statistics.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    fafile : pysam.FastaFile or string
        FASTA file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column
    window_size : int
        size of surrounding window in base pairs
    window_offset :
        distance from window start to record position

    Returns
    -------

    recs : iterator
        record generator

    """

    stat = _CoverageBinned()
    return _iter_binned(stat, alignmentfile=alignmentfile, fafile=fafile, **kwargs)


cdef int _gc_content(ref_window) except -1:
    cdef Py_ssize_t i, n
    cdef char* seq
    cdef int gc_count = 0
    n = len(ref_window)
    if not PY2:
        # TODO optimise?
        ref_window = ref_window.encode('ascii')
    seq = ref_window
    for i in range(n):
        if seq[i] == b'g' or seq[i] == b'c':
            gc_count += 1
    gc_percent = int(round(gc_count * 100. / n))
    return gc_percent


cdef class _StatBinned(object):

    cdef dict rec(self, chrom, bin_start, bin_end, FastaFile fafile):
        return dict()

    cdef recv(self, bam1_t * b):
        pass


cdef class _CoverageBinned(_StatBinned):

    cdef int reads_all, reads_pp

    def __cinit__(self):
        self.reads_all = self.reads_pp = 0

    cdef dict rec(self, chrom, bin_start, bin_end, FastaFile fafile):

        # determine %GC
        ref_window = fafile.fetch(chrom, bin_start, bin_end).lower()
        gc_percent = _gc_content(ref_window)

        # make record for bin
        rec = {'gc': gc_percent,
               'reads_all': self.reads_all,
               'reads_pp': self.reads_pp}

        # reset counters
        self.reads_all = self.reads_pp = 0

        return rec

    cdef recv(self, bam1_t * b):
        cdef uint32_t flag
        cdef bint is_unmapped
        cdef bint is_proper_pair
        flag = b.core.flag
        is_unmapped = <bint>(flag & BAM_FUNMAP)
        if not is_unmapped:
            self.reads_all += 1
            is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
            if is_proper_pair:
                self.reads_pp += 1


def load_coverage_binned(*args, **kwargs):
    return _load_stats(stat_coverage_binned, dtype_coverage_binned,
                       *args, **kwargs)
        
    
############################################
# BINNED COVERAGE WITH EXTENDED PROPERTIES #
############################################


dtype_coverage_ext_binned = [('chrom', 'a12'),
                             ('pos', 'i4'),
                             ('gc', 'u1'),
                             ('reads_all', 'i4'),
                             ('reads_pp', 'i4'),
                             ('reads_mate_unmapped', 'i4'),
                             ('reads_mate_other_chr', 'i4'),
                             ('reads_mate_same_strand', 'i4'),
                             ('reads_faceaway', 'i4'),
                             ('reads_softclipped', 'i4'),
                             ('reads_duplicate', 'i4')]


fields_coverage_ext_binned = [t[0] for t in dtype_coverage_ext_binned]


def stat_coverage_ext_binned(alignmentfile, fafile, **kwargs):
    """Generate binned extended coverage statistics.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    fafile : pysam.FastaFile or string
        FASTA file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column
    window_size : int
        size of surrounding window in base pairs
    window_offset :
        distance from window start to record position

    Returns
    -------

    recs : iterator
        record generator

    """

    stat = _CoverageExtBinned()
    return _iter_binned(stat, alignmentfile=alignmentfile, fafile=fafile,
                        **kwargs)


cdef class _CoverageExtBinned(_StatBinned):

    cdef int reads_all, reads_pp, reads_mate_unmapped, reads_mate_other_chr, \
        reads_mate_same_strand, reads_faceaway, reads_softclipped, \
        reads_duplicate

    def __cinit__(self):
        self.reads_all = self.reads_pp = self.reads_mate_unmapped \
            = self.reads_mate_other_chr = self.reads_mate_same_strand\
            = self.reads_faceaway = self.reads_softclipped \
            = self.reads_duplicate = 0

    cdef dict rec(self, chrom, bin_start, bin_end, FastaFile fafile):

        # determine %GC
        ref_window = fafile.fetch(chrom, bin_start, bin_end).lower()
        gc_percent = _gc_content(ref_window)

        # make record for bin
        rec = {'gc': gc_percent,
               'reads_all': self.reads_all,
               'reads_pp': self.reads_pp,
               'reads_mate_unmapped': self.reads_mate_unmapped,
               'reads_mate_other_chr': self.reads_mate_other_chr,
               'reads_mate_same_strand': self.reads_mate_same_strand,
               'reads_faceaway': self.reads_faceaway,
               'reads_softclipped': self.reads_softclipped,
               'reads_duplicate': self.reads_duplicate}

        # reset counters
        self.reads_all = self.reads_pp = self.reads_mate_unmapped \
            = self.reads_mate_other_chr = self.reads_mate_same_strand\
            = self.reads_faceaway = self.reads_softclipped \
            = self.reads_duplicate = 0

        return rec

    cdef recv(self, bam1_t * b):
        cdef uint32_t flag
        cdef bint is_unmapped
        cdef bint is_reverse
        cdef bint is_proper_pair
        cdef bint is_duplicate
        cdef bint mate_is_unmappped
        cdef bint mate_is_reverse
        cdef int tlen
        flag = b.core.flag
        is_unmapped = <bint>(flag & BAM_FUNMAP)
        if not is_unmapped:
            self.reads_all += 1
            is_reverse = <bint>(flag & BAM_FREVERSE)
            is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
            is_duplicate = <bint>(flag & BAM_FDUP)
            mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
            mate_is_reverse = <bint>(flag & BAM_FMREVERSE)
            tlen = b.core.isize
            if is_duplicate:
                self.reads_duplicate += 1
            if is_proper_pair:
                self.reads_pp += 1
            if mate_is_unmapped:
                self.reads_mate_unmapped += 1
            elif b.core.tid != b.core.mtid:
                self.reads_mate_other_chr += 1
            elif (is_reverse and mate_is_reverse) \
                    or (not is_reverse and not mate_is_reverse):
                self.reads_mate_same_strand += 1
            elif (is_reverse and tlen > 0) or (not is_reverse and tlen < 0):
                self.reads_faceaway += 1
            if _is_softclipped(b):
                self.reads_softclipped += 1


def load_coverage_ext_binned(*args, **kwargs):
    return _load_stats(stat_coverage_ext_binned, dtype_coverage_ext_binned,
                       *args, **kwargs)
        
    
###############
# BINNED MAPQ #
###############


dtype_mapq_binned = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_mapq0', 'i4'),
    ('rms_mapq', 'i4'),
]


fields_mapq_binned = [t[0] for t in dtype_mapq_binned]


def stat_mapq_binned(alignmentfile, **kwargs):
    """Generate binned mapping quality statistics.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column
    window_size : int
        size of surrounding window in base pairs
    window_offset :
        distance from window start to record position

    Returns
    -------

    recs : iterator
        record generator

    """

    stat = _MapqBinned()
    return _iter_binned(stat, alignmentfile=alignmentfile, **kwargs)


cdef class _MapqBinned(_StatBinned):

    cdef int reads_all, reads_mapq0
    cdef uint64_t mapq, mapq_squared_sum

    def __cinit__(self):
        self.reads_all = self.reads_mapq0 = self.mapq_squared_sum = 0

    cdef dict rec(self, chrom, bin_start, bin_end, FastaFile fafile):

        # make record for bin
        rec = {'reads_all': self.reads_all,
               'reads_mapq0': self.reads_mapq0,
               'rms_mapq': _rootmean(self.mapq_squared_sum, self.reads_all)}

        # reset counters
        self.reads_all = self.reads_mapq0 = self.mapq_squared_sum = 0

        return rec

    cdef recv(self, bam1_t * b):
        cdef uint32_t flag
        cdef bint is_unmapped
        cdef uint64_t mapq
        cdef uint64_t mapq_squared
        flag = b.core.flag
        is_unmapped = <bint>(flag & BAM_FUNMAP)
        if not is_unmapped:
            self.reads_all += 1
            mapq = b.core.qual
            mapq_squared = mapq**2
            self.mapq_squared_sum += mapq_squared
            if mapq == 0:
                self.reads_mapq0 += 1


def load_mapq_binned(*args, **kwargs):
    return _load_stats(stat_mapq_binned, dtype_mapq_binned, *args, **kwargs)
        
    
################
# BINNED CIGAR #
################


dtype_alignment_binned = [
    ('chrom', 'a12'),
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


fields_alignment_binned = [t[0] for t in dtype_alignment_binned]


def stat_alignment_binned(alignmentfile, **kwargs):
    """Generate binned alignment statistics.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    fafile : pysam.FastaFile or string
        FASTA file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column
    window_size : int
        size of surrounding window in base pairs
    window_offset :
        distance from window start to record position

    Returns
    -------

    recs : iterator
        record generator

    """

    stat = _AlignmentBinned()
    return _iter_binned(stat, alignmentfile=alignmentfile, **kwargs)


cdef class _AlignmentBinned(_StatBinned):

    cdef int reads_all, M, I, D, N, S, H, P, EQ, X

    def __cinit__(self):
        self.reads_all = self.M = self.I = self.D = self.N = self.S = self.H\
            = self.P = self.EQ = self.X = 0

    cdef dict rec(self, chrom, bin_start, bin_end, FastaFile fafile):

        # make record for bin
        rec = {'reads_all': self.reads_all,
               'M': self.M, 'I': self.I, 'D': self.D, 'N': self.N, 'S': self.S,
               'H': self.H, 'P': self.P, '=': self.EQ, 'X': self.X,
               'bases_all': self.M + self.I + self.S + self.EQ + self.X}

        # reset counters
        self.reads_all = self.M = self.I = self.D = self.N = self.S = self.H\
            = self.P = self.EQ = self.X = 0

        return rec

    cdef recv(self, bam1_t * b):
        cdef uint32_t flag
        cdef bint is_unmapped
        cdef int k, op, l
        flag = b.core.flag
        is_unmapped = <bint>(flag & BAM_FUNMAP)
        if not is_unmapped:
            cigar_p = pysam_bam_get_cigar(b)
            cigar = list()
            for k in range(b.core.n_cigar):
                op = cigar_p[k] & BAM_CIGAR_MASK
                l = cigar_p[k] >> BAM_CIGAR_SHIFT
                cigar.append((op, l))
                if op == BAM_CMATCH:
                    self.M += l
                elif op == BAM_CINS:
                    self.I += l
                elif op == BAM_CDEL:
                    self.D += l
                elif op == BAM_CREF_SKIP:
                    self.N += l
                elif op == BAM_CSOFT_CLIP:
                    self.S += l
                elif op == BAM_CHARD_CLIP:
                    self.H += l
                elif op == BAM_CPAD:
                    self.P += l
                elif op == BAM_CEQUAL:
                    self.EQ += l
                elif op == BAM_CDIFF:
                    self.X += l
            self.reads_all += 1

    
def load_alignment_binned(*args, **kwargs):
    return _load_stats(stat_alignment_binned, dtype_alignment_binned,
                       *args, **kwargs)
        
    
###############
# BINNED TLEN #
###############


dtype_tlen_binned = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_pp', 'i4'),
    ('mean_tlen', 'i4'),
    ('mean_tlen_pp', 'i4'),
    ('rms_tlen', 'i4'),
    ('rms_tlen_pp', 'i4'),
]


fields_tlen_binned = [t[0] for t in dtype_tlen_binned]


def stat_tlen_binned(alignmentfile, **kwargs):
    """Generate binned insert size statistics.

    Parameters
    ----------

    alignmentfile : pysam.AlignmentFile or string
        SAM or BAM file or file path
    chrom : string
        chromosome/contig
    start : int
        start position
    end : int
        end position
    one_based : bool
        coordinate system
    truncate : bool
        if True, truncate output to selected region
    pad : bool
        if True, emit records for every position, even if no reads are aligned
    max_depth : int
        maximum depth to allow in pileup column
    window_size : int
        size of surrounding window in base pairs
    window_offset :
        distance from window start to record position

    Returns
    -------

    recs : iterator
        record generator

    """

    stat = _TlenBinned()
    return _iter_binned(stat, alignmentfile=alignmentfile, **kwargs)


cdef class _TlenBinned(_StatBinned):

    cdef int reads_all
    cdef int reads_pp
    cdef int64_t tlen_sum
    cdef int64_t tlen_pp_sum
    cdef int64_t tlen_squared_sum
    cdef int64_t tlen_pp_squared_sum

    def __cinit__(self):
        self.tlen_sum = self.tlen_squared_sum = self.tlen_pp_sum = \
            self.tlen_pp_squared_sum = self.reads_all = self.reads_pp = 0

    cdef dict rec(self, chrom, bin_start, bin_end, FastaFile fafile):

        # make record for bin
        rec = {'reads_all': self.reads_all, 'reads_pp': self.reads_pp,
               'mean_tlen': _mean(self.tlen_sum, self.reads_all),
               'mean_tlen_pp': _mean(self.tlen_pp_sum, self.reads_pp),
               'rms_tlen': _rootmean(self.tlen_squared_sum, self.reads_all),
               'rms_tlen_pp': _rootmean(self.tlen_pp_squared_sum,
                                        self.reads_pp)}

        # reset counters
        self.tlen_sum = self.tlen_squared_sum = self.tlen_pp_sum = \
            self.tlen_pp_squared_sum = self.reads_all = self.reads_pp = 0

        return rec

    cdef recv(self, bam1_t * b):
        cdef uint32_t flag
        cdef bint is_unmapped
        cdef bint is_proper_pair
        cdef int64_t tlen
        cdef int64_t tlen_squared
        flag = b.core.flag
        is_unmapped = <bint>(flag & BAM_FUNMAP)
        if not is_unmapped:
            self.reads_all += 1
            tlen = b.core.isize
            self.tlen_sum += tlen
            tlen_squared = tlen**2
            self.tlen_squared_sum += tlen_squared
            is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
            if is_proper_pair:
                self.reads_pp += 1
                self.tlen_pp_sum += tlen
                self.tlen_pp_squared_sum += tlen_squared


def load_tlen_binned(*args, **kwargs):
    return _load_stats(stat_tlen_binned, dtype_tlen_binned, *args, **kwargs)


#####################
# UTILITY FUNCTIONS #
#####################


def _iter_pileup(frec, fpad, alignmentfile, fafile=None, chrom=None, start=None,
                 end=None, one_based=False, truncate=False, pad=False,
                 max_depth=8000, **kwargs):
    """General purpose function to generate statistics, where one record is
    generated for each genome position in the selected range, based on a
    pileup column.

    """

    if isinstance(alignmentfile, _string_types):
        alignmentfile = AlignmentFile(alignmentfile)

    if isinstance(fafile, _string_types):
        fafile = FastaFile(fafile)

    if pad:
        return _iter_pileup_padded(frec, fpad, alignmentfile,
                                   fafile=fafile, chrom=chrom,
                                   start=start, end=end, one_based=one_based,
                                   truncate=truncate, max_depth=max_depth,
                                   **kwargs)
    else:
        return _iter_pileup_default(frec, alignmentfile, fafile=fafile, chrom=chrom,
                                    start=start, end=end, one_based=one_based,
                                    truncate=truncate, max_depth=max_depth,
                                    **kwargs)


def _iter_pileup_default(frec, alignmentfile, fafile=None,
                         chrom=None, start=None, end=None,
                         one_based=False, truncate=False,
                         max_depth=8000, **kwargs):
    start, end = _normalise_coords(alignmentfile, chrom, start, end, one_based)
    it = alignmentfile.pileup(reference=chrom, start=start, end=end,
                        truncate=truncate, max_depth=max_depth)
    for col in it:
        yield frec(alignmentfile, fafile, col, one_based)


def _iter_pileup_padded(frec, fpad, alignmentfile, fafile=None, chrom=None,
                        start=None, end=None, one_based=False,
                        truncate=False, max_depth=8000, **kwargs):

    if chrom is not None:
        it = _iter_pileup_padded_chrom(frec, fpad, alignmentfile, fafile=fafile,
                                       chrom=chrom, start=start, end=end,
                                       one_based=one_based,
                                       truncate=truncate,
                                       max_depth=max_depth, **kwargs)
    else:
        its = list()
        for chrom in alignmentfile.references:
            itc = _iter_pileup_padded_chrom(frec, fpad, alignmentfile, fafile=fafile,
                                            chrom=chrom, start=None,
                                            end=None,
                                            one_based=one_based,
                                            truncate=truncate,
                                            max_depth=max_depth,
                                            **kwargs)
            its.append(itc)
        it = _itertools.chain(*its)
    return it


def _iter_pileup_padded_chrom(frec, fpad, alignmentfile, fafile, chrom, start=None,
                              end=None, one_based=False, truncate=False,
                              max_depth=8000, **kwargs):
    cdef PileupColumn col
    cdef int curpos
    assert chrom is not None, 'chromosome is None'
    start, end = _normalise_coords(alignmentfile, chrom, start, end, one_based)
    it = alignmentfile.pileup(reference=chrom, start=start, end=end,
                        truncate=truncate, max_depth=max_depth)
    curpos = start
    for col in it:
        while curpos < col.pos:
            yield fpad(fafile, chrom, curpos, one_based)
            curpos += 1
        yield frec(alignmentfile, fafile, col, one_based)
        curpos = col.pos + 1
    while curpos < end:
        yield fpad(fafile, chrom, curpos, one_based)
        curpos += 1


def _iter_binned(stat, alignmentfile, fafile=None, chrom=None, start=None, end=None,
                 one_based=False, window_size=300, window_offset=None,
                 **kwargs):

    if isinstance(alignmentfile, _string_types):
        alignmentfile = AlignmentFile(alignmentfile)

    if isinstance(fafile, _string_types):
        fafile = FastaFile(fafile)

    if window_offset is None:
        window_offset = window_size / 2

    if chrom is None:
        its = list()
        for chrom in alignmentfile.references:
            itc = _iter_binned_chrom(stat, alignmentfile=alignmentfile,
                                     fafile=fafile, chrom=chrom, start=None,
                                     end=None, one_based=one_based,
                                     window_size=window_size,
                                     window_offset=window_offset)
            its.append(itc)
        it = _itertools.chain(*its)

    else:
        it = _iter_binned_chrom(stat, alignmentfile=alignmentfile, fafile=fafile,
                                chrom=chrom, start=start, end=end,
                                one_based=one_based, window_size=window_size,
                                window_offset=window_offset)

    return it


def _iter_binned_chrom(_StatBinned stat, AlignmentFile alignmentfile, FastaFile fafile,
                       chrom, start=None, end=None, one_based=False,
                       int window_size=300, int window_offset=150, **kwargs):

    # setup

    assert chrom is not None, 'chromosome is None'
    start, end = _normalise_coords(alignmentfile, chrom, start, end, one_based)

    cdef int rtid, rstart, rend, has_coord, bin_start, bin_end
    has_coord, rtid, rstart, rend = \
        alignmentfile.parse_region(chrom, start, end, None)

    cdef IteratorRowRegion it
    cdef bam1_t * b
    it = IteratorRowRegion(alignmentfile, rtid, rstart, rend,
                           multiple_iterators=False)
    b = it.b

    # setup first bin
    bin_start = rstart
    bin_end = bin_start + window_size

    # iterate over reads
    it.cnext()
    while it.retval > 0:
        while b.core.pos > bin_end:  # end of bin

            # yield record for bin
            rec = stat.rec(chrom, bin_start, bin_end, fafile)
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec['chrom'] = chrom
            rec['pos'] = pos
            yield rec

            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size

        # process current read
        stat.recv(b)

        # move iterator on
        it.cnext()

    # deal with last non-empty bin
    rec = stat.rec(chrom, bin_start, bin_end, fafile)
    pos = bin_start + window_offset
    if one_based:
        pos += 1
    rec['chrom'] = chrom
    rec['pos'] = pos
    yield rec

    # deal with empty bins up to explicit end
    if end is not None:
        while bin_end < end:
            # start new bin
            bin_start = bin_end
            bin_end = bin_start + window_size

            # yield record
            rec = stat.rec(chrom, bin_start, bin_end, fafile)
            pos = bin_start + window_offset
            if one_based:
                pos += 1
            rec['chrom'] = chrom
            rec['pos'] = pos
            yield rec


def _normalise_coords(AlignmentFile alignmentfile, chrom, start, end, one_based):
    """Convert to zero-based coordinates and deal with unspecified start
    and/or end.

    """

    if chrom is None:
        return None, None

    else:
        assert chrom in alignmentfile.references, \
            'chromosome not in SAM references: %s' % chrom

        if one_based:
            start = start - 1 if start is not None else None
            end = end - 1 if end is not None else None

        chrlen = alignmentfile.lengths[alignmentfile.references.index(chrom)]
        if start is None:
            start = 0
        if end is None:
            end = chrlen
        if end > chrlen:
            end = chrlen

        return start, end

    
def write_csv(stats_type, outfile, alignmentfile, fields=None, dialect='excel-tab',
              write_header=True, progress=None, **kwargs):
    """Write statistics output to a CSV file.

    Parameters
    ----------

    stats_type : string
        statistics type, one of 'coverage', 'coverage_ext', ...
    outfile : file-like
        output file to write to
    alignmentfile : pysam.AlignmentFile or string
        input BAM or SAM file or file path
    fields : list of strings
        list of field names to output (all by default)
    dialect : string
        CSV dialect
    write_header : bool
        if True write a header row
    progress : int
        log progress to stderr every N rows
    **kwargs
        passed through to the statistics function

    """

    cdef long counter, modulus

    # lookup stats function
    stats_function = globals()['stat_' + stats_type]

    # determine field names
    if not fields:
        fields = globals()['fields_' + stats_type]

    # setup record generator
    recs = stats_function(alignmentfile, **kwargs)

    # flatten records to rows
    rows = flatten(recs, *fields)

    # initialise writer
    writer = _csv.writer(outfile, dialect=dialect)

    # write header row
    if write_header:
        writer.writerow(fields)

    if progress is None:
        # N.B., don't use writer.writerows(recs)!
        for row in rows:
            writer.writerow(row)

    else:
        counter = 0
        modulus = progress
        before = _time.time()
        before_all = before
        for row in rows:
            counter += 1
            writer.writerow(row)
            if counter % modulus == 0:
                after = _time.time()
                elapsed = after - before_all
                batch_elapsed = after - before
                msg = '[pysamstats] %s rows in %.2fs (%d rows/s); batch in ' \
                      '%.2fs (%d rows/s)' \
                      % (counter, elapsed, counter / elapsed, batch_elapsed,
                         progress / batch_elapsed)
                print(msg, file=_sys.stderr)
                before = after
        after_all = _time.time()
        elapsed_all = after_all - before_all
        msg = '[pysamstats] %s rows in %.2fs (%d rows/s)' \
              % (counter, elapsed_all, counter / elapsed_all)
        print(msg, file=_sys.stderr)
    

def write_hdf5(stats_type, outfile, alignmentfile, fields=None, progress=None,
               hdf5_group='/', hdf5_dataset='data', hdf5_complevel=5,
               hdf5_complib='zlib', hdf5_shuffle=True,
               hdf5_fletcher32=False, hdf5_chunksize=2**20, **kwargs):
    """Write statistics output to an HDF5 file. Requires PyTables.

    Parameters
    ----------

    stats_type : string
        statistics type, one of 'coverage', 'coverage_ext', ...
    outfile : string
        output file path
    alignmentfile : pysam.AlignmentFile or string
        input BAM or SAM file or file path
    fields : list of strings
        list of field names to output (all by default)
    progress : int
        log progress to stderr approximately every N rows
    hdf5_group : string
        group to write new dataset to
    hdf5_dataset : string
        name of dataset to create
    hdf5_chunksize : int
        size of chunks in number of bytes
    dtype : dict
        override dtype
    **kwargs
        passed through to the statistics function

    Notes
    -----

    The length of the chunks in number of items is calculated by dividing the
    chunk size in number of bytes by the size of each row in number of bytes as
    determined from the dtype.

    """

    import tables
    import numpy as np
    h5file = None

    # lookup stats function
    stats_function = globals()['stat_' + stats_type]

    # determine field names
    if not fields:
        fields = globals()['fields_' + stats_type]

    # determine dtype
    default_dtype = globals()['dtype_' + stats_type]
    dtype = dict(default_dtype)
    dtype_overrides = kwargs.pop('dtype', None)
    if dtype_overrides:
        # expect dict
        dtype.update(dtype_overrides)
    if len(fields) == 1:
        dtype = dtype[fields[0]]
    else:
        dtype = [(f, dtype[f]) for f in fields]
    dtype = np.dtype(dtype)

    # setup record generator
    recs = stats_function(alignmentfile, **kwargs)

    # flatten records to rows
    rows = flatten(recs, *fields)

    try:

        # open output file
        h5file = tables.open_file(outfile, mode='a')

        # determine chunk shape
        hdf5_chunklen = int(hdf5_chunksize/dtype.itemsize)
        hdf5_chunkshape = (hdf5_chunklen,)

        # replace any existing node at that location
        try:
            h5file.remove_node(hdf5_group, hdf5_dataset)
        except tables.NoSuchNodeError:
            pass

        # create dataset
        h5table = h5file.create_table(
            hdf5_group, hdf5_dataset, dtype,
            title=stats_type,
            filters=tables.Filters(complevel=hdf5_complevel,
                                   complib=hdf5_complib,
                                   shuffle=hdf5_shuffle,
                                   fletcher32=hdf5_fletcher32),
            createparents=True,
            chunkshape=hdf5_chunkshape)

        # record initial time
        counter = 0
        counter_before = 0
        before = _time.time()
        before_all = before

        # load data in batches of size `hdf5_chunklen`
        chunk = list(_itertools.islice(rows, hdf5_chunklen))

        # load chunk at a time
        while chunk:

            # write chunk
            h5table.append(chunk)
            h5table.flush()

            # keep track of number of records loaded
            n = len(chunk)  # may be shorter than chunklen if final batch
            counter += n

            # log progress
            if progress and (counter % progress) < hdf5_chunklen:
                after = _time.time()
                elapsed = after - before_all
                batch_elapsed = after - before
                batch_size = counter - counter_before
                msg = '[pysamstats] %s rows in %.2fs (%d rows/s); last %s ' \
                      'rows in %.2fs (%d rows/s)' \
                      % (counter, elapsed, counter / elapsed,
                         batch_size, batch_elapsed, batch_size / batch_elapsed)
                print(msg, file=_sys.stderr)
                before = after
                counter_before = counter

            # load next batch
            chunk = list(_itertools.islice(rows, hdf5_chunklen))

        if progress:
            after_all = _time.time()
            elapsed_all = after_all - before_all
            msg = '[pysamstats] %s rows in %.2fs (%d rows/s)' \
                  % (counter, elapsed_all, counter / elapsed_all)
            print(msg, file=_sys.stderr)

    finally:
        if h5file is not None:
            h5file.close()


from operator import itemgetter


def flatten(recs, *fields):
    """Convert a record (dict) iterator to a row (tuple) iterator.

    Parameters
    ----------

    recs : iterator of dicts
        records generator
    fields : list of strings
        names of fields to select

    Returns
    -------

    rows : iterator of tuples
        rows generator

    """

    getter = itemgetter(*fields)
    it = (getter(rec) for rec in recs)
    return it


def tabulate(stat, *args, **kwargs):
    """Tabulate statistics.

    Parameters
    ----------

    stat : string
        statistics type
    *args
        passed through to statistics function
    fields : list of strings
        names of fields to select
    **args
        passed through to statistics function

    Returns
    -------

    table : row container

    """

    return _StatsTable(stat, *args, **kwargs)


class _StatsTable(object):

    def __init__(self, stats_type, *args, **kwargs):
        try:
            self.stats_function = globals()['stat_' + stats_type]
            self.fields = kwargs.pop('fields', None)
            if self.fields is None:
                self.fields = globals()['fields_' + stats_type]
            self.args = args
            self.kwargs = kwargs
        except KeyError:
            raise Exception('statistics type not found: %r' % stats_type)

    def __iter__(self):
        recs = self.stats_function(*self.args, **self.kwargs)
        fields = tuple(self.fields)
        rows = flatten(recs, *fields)
        yield fields
        for row in rows:
            yield row


def _load_stats(statfun, default_dtype, *args, **kwargs):

    import numpy as np

    # determine fields to load
    fields = kwargs.pop('fields', None)
    if not fields:
        fields = [t[0] for t in default_dtype]

    # determine dtype
    dtype = dict(default_dtype)
    dtype_overrides = kwargs.pop('dtype', None)
    if dtype_overrides:
        # expect dict
        dtype.update(dtype_overrides)
    if len(fields) == 1:
        dtype = dtype[fields[0]]
    else:
        dtype = [(f, dtype[f]) for f in fields]

    # setup record generator
    recs = statfun(*args, **kwargs)

    # flatten records
    it = flatten(recs, *fields)

    # load into a Numpy array
    a = np.fromiter(it, dtype=dtype)

    # view as recarray for convenience
    if len(fields) > 1:
        a = a.view(np.recarray)

    return a
    
                
cdef inline bint _is_softclipped(bam1_t * aln):
    cdef int k
    cigar_p = pysam_bam_get_cigar(aln)
    for k in range(aln.core.n_cigar):
        op = cigar_p[k] & BAM_CIGAR_MASK
        if op == BAM_CSOFT_CLIP:
            return 1
    return 0


cdef inline object _get_seq_base(bam1_t *src, uint32_t k):
    cdef uint8_t * p
    cdef char * s

    if not src.core.l_qseq:
        return None

    seq = PyBytes_FromStringAndSize(NULL, 1)
    s   = <char*>seq
    p   = pysam_bam_get_seq(src)

    s[0] = bam_nt16_rev_table[p[k//2] >> 4 * (1 - k%2) & 0xf]

    return seq


cdef inline int _rootmean(uint64_t sqsum, int count):
    if count > 0:
        return int(round(sqrt(sqsum / count)))
    else:
        return 0
    
    
cdef inline int _mean(int64_t total, int count):
    if count > 0:
        return int(round(total / count))
    else:
        return 0


# SANDBOX


def count_reads(AlignmentFile alignmentfile, chrom=None, start=None, end=None):
    cdef IteratorRowRegion it
    cdef int n = 0
    has_coord, rtid, rstart, rend = alignmentfile.parse_region(chrom, start,
                                                               end, None)
    it = IteratorRowRegion(alignmentfile, rtid, rstart, rend,
                           multiple_iterators=False)
    while True:
        it.cnext()
        if it.retval > 0:
            n += 1
        else:
            break
    return n
