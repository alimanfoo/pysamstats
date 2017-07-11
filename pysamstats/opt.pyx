# cython: profile=False
# cython: embedsignature=True
from __future__ import print_function, division, absolute_import


from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t
from libc.math cimport sqrt
from cpython cimport PyBytes_FromStringAndSize
from pysam.libchtslib cimport bam1_t, bam_pileup1_t
from pysam.libcfaidx cimport FastaFile
from pysam.libcalignmentfile cimport AlignmentFile, IteratorRowRegion
from pysam.libcalignedsegment cimport pysam_bam_get_cigar, pysam_bam_get_seq, \
    pysam_bam_get_qual, PileupColumn


import sys as _sys
import itertools


# PY2/3 compatibility
PY2 = _sys.version_info[0] == 2
if PY2:
    # noinspection PyUnresolvedReferences
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

cdef char* bam_nt16_rev_table = "=ACMGRSVTWYHKDBN"


cdef class PileupStat(object):

    cdef dict rec(self, chrom, pos, FastaFile fafile):
        return dict()

    cdef void recv(self, PileupColumn col, bam1_t* b):
        pass


#############################
# BASIC COVERAGE STATISTICS #
#############################


# noinspection PyAttributeOutsideInit
cdef class Coverage(PileupStat):

    cdef:
        int reads_all
        int reads_pp

    def __init__(self):
        self.reset()

    def reset(self):
        self.reads_all = 0
        self.reads_pp = 0

    cdef void recv(self, PileupColumn col, bam1_t* b):
        cdef:
            bint is_proper_pair

        self.reads_all += 1
        is_proper_pair = <bint>(b.core.flag & BAM_FPROPER_PAIR)
        if is_proper_pair:
            self.reads_pp += 1

    cdef dict rec(self, chrom, pos, FastaFile fafile):

        # make record
        rec = {'reads_all': self.reads_all,
               'reads_pp': self.reads_pp}

        # reset counters
        self.reset()

        return rec


################################
# STRANDED COVERAGE STATISTICS #
################################


# noinspection PyAttributeOutsideInit
cdef class CoverageStrand(PileupStat):

    cdef:
        int reads_all
        int reads_fwd
        int reads_rev
        int reads_pp
        int reads_pp_fwd
        int reads_pp_rev

    def __init__(self):
        self.reset()

    def reset(self):
        self.reads_all = 0
        self.reads_pp = 0
        self.reads_fwd = 0
        self.reads_rev = 0
        self.reads_pp = 0
        self.reads_pp_fwd = 0
        self.reads_pp_rev = 0

    cdef void recv(self, PileupColumn col, bam1_t* b):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bint is_reverse

        self.reads_all += 1
        flag = b.core.flag
        is_reverse = <bint>(flag & BAM_FREVERSE)
        if is_reverse:
            self.reads_rev += 1
        else:
            self.reads_fwd += 1
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        if is_proper_pair:
            self.reads_pp += 1
            if is_reverse:
                self.reads_pp_rev += 1
            else:
                self.reads_pp_fwd += 1

    cdef dict rec(self, chrom, pos, FastaFile fafile):

        # make record
        rec = {'reads_all': self.reads_all,
               'reads_fwd': self.reads_fwd,
               'reads_rev': self.reads_rev,
               'reads_pp': self.reads_pp,
               'reads_pp_fwd': self.reads_pp_fwd,
               'reads_pp_rev': self.reads_pp_rev}

        # reset counters
        self.reset()

        return rec


################################
# EXTENDED COVERAGE STATISTICS #
################################


# noinspection PyAttributeOutsideInit
cdef class CoverageExt(PileupStat):

    cdef:
        int reads_all
        int reads_pp
        int reads_mate_unmapped
        int reads_mate_other_chr
        int reads_mate_same_strand
        int reads_faceaway
        int reads_softclipped
        int reads_duplicate

    def __init__(self):
        self.reset()

    def reset(self):
        self.reads_all = 0
        self.reads_pp = 0
        self.reads_mate_unmapped = 0
        self.reads_mate_other_chr = 0
        self.reads_mate_same_strand = 0
        self.reads_faceaway = 0
        self.reads_softclipped = 0
        self.reads_duplicate = 0

    cdef void recv(self, PileupColumn col, bam1_t* b):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bint is_reverse
            bint is_duplicate
            bint mate_is_unmappped
            bint mate_is_reverse
            int tlen

        self.reads_all += 1
        flag = b.core.flag
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
        elif col.tid != b.core.mtid:
            self.reads_mate_other_chr += 1
        elif (is_reverse and mate_is_reverse) or (not is_reverse and not mate_is_reverse):
            self.reads_mate_same_strand += 1
        elif (is_reverse and tlen > 0) or (not is_reverse and tlen < 0):
            self.reads_faceaway += 1
        if is_softclipped(b):
            self.reads_softclipped += 1

    cdef dict rec(self, chrom, pos, FastaFile fafile):

        # make record
        rec = {'reads_all': self.reads_all,
               'reads_pp': self.reads_pp,
               'reads_mate_unmapped': self.reads_mate_unmapped,
               'reads_mate_other_chr': self.reads_mate_other_chr,
               'reads_mate_same_strand': self.reads_mate_same_strand,
               'reads_faceaway': self.reads_faceaway,
               'reads_softclipped': self.reads_softclipped,
               'reads_duplicate': self.reads_duplicate}

        # reset counters
        self.reset()

        return rec


# cpdef dict rec_coverage_ext(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
#                             bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i, n  # loop index
#     cdef int reads_all  # total number of reads in column
#     cdef bint is_reverse
#     cdef bint is_proper_pair
#     cdef bint is_duplicate
#     cdef bint mate_is_unmappped
#     cdef bint mate_is_reverse
#     cdef int tlen
#     # counting variables
#     cdef int reads_pp = 0
#     cdef int reads_mate_unmapped = 0
#     cdef int reads_mate_other_chr = 0
#     cdef int reads_mate_same_strand = 0
#     cdef int reads_faceaway = 0
#     cdef int reads_softclipped = 0
#     cdef int reads_duplicate = 0
#
#     # initialise variables
#     n = col.n
#     plp = col.plp
#
#     # get chromosome name and position
#     tid = col.tid
#     chrom = alignmentfile.getrname(tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # loop over reads, extract what we need
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#         is_reverse = <bint>(flag & BAM_FREVERSE)
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         is_duplicate = <bint>(flag & BAM_FDUP)
#         mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
#         mate_is_reverse = <bint>(flag & BAM_FMREVERSE)
#         tlen = aln.core.isize
#         if is_duplicate:
#             reads_duplicate += 1
#         if is_proper_pair:
#             reads_pp += 1
#         if mate_is_unmapped:
#             reads_mate_unmapped += 1
#         elif tid != aln.core.mtid:
#             reads_mate_other_chr += 1
#         elif (is_reverse and mate_is_reverse) \
#                 or (not is_reverse and not mate_is_reverse):
#             reads_mate_same_strand += 1
#         elif (is_reverse and tlen > 0) or (not is_reverse and tlen < 0):
#             reads_faceaway += 1
#         if is_softclipped(aln):
#             reads_softclipped += 1
#
#     return {'chrom': chrom, 'pos': pos,
#             'reads_all': n,
#             'reads_pp': reads_pp,
#             'reads_mate_unmapped': reads_mate_unmapped,
#             'reads_mate_other_chr': reads_mate_other_chr,
#             'reads_mate_same_strand': reads_mate_same_strand,
#             'reads_faceaway': reads_faceaway,
#             'reads_softclipped': reads_softclipped,
#             'reads_duplicate': reads_duplicate}
#
#
# cpdef dict rec_coverage_ext_pad(FastaFile fafile, chrom, pos,
#                                 bint one_based=False):
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom, 'pos': pos,
#             'reads_all': 0,
#             'reads_pp': 0,
#             'reads_mate_unmapped': 0,
#             'reads_mate_other_chr': 0,
#             'reads_mate_same_strand': 0,
#             'reads_faceaway': 0,
#             'reads_softclipped': 0,
#             'reads_duplicate': 0}
#
#
# ##########################################
# # EXTENDED COVERAGE STATISTICS BY STRAND #
# ##########################################
#
#
# cpdef dict rec_coverage_ext_strand(AlignmentFile alignmentfile, FastaFile fafile,
#                                    PileupColumn col, bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i, n # loop index
#     cdef int reads_all # total number of reads in column
#     cdef bint is_reverse
#     cdef bint is_proper_pair
#     cdef bint is_duplicate
#     cdef bint mate_is_unmappped
#     cdef bint mate_is_reverse
#     cdef int tlen
#     # counting variables
#     cdef int reads_rev = 0
#     cdef int reads_fwd = 0
#     cdef int reads_pp = 0
#     cdef int reads_pp_rev = 0
#     cdef int reads_pp_fwd = 0
#     cdef int reads_mate_unmapped = 0
#     cdef int reads_mate_unmapped_rev = 0
#     cdef int reads_mate_unmapped_fwd = 0
#     cdef int reads_mate_other_chr = 0
#     cdef int reads_mate_other_chr_rev = 0
#     cdef int reads_mate_other_chr_fwd = 0
#     cdef int reads_mate_same_strand = 0
#     cdef int reads_mate_same_strand_rev = 0
#     cdef int reads_mate_same_strand_fwd = 0
#     cdef int reads_faceaway = 0
#     cdef int reads_faceaway_rev = 0
#     cdef int reads_faceaway_fwd = 0
#     cdef int reads_softclipped = 0
#     cdef int reads_softclipped_rev = 0
#     cdef int reads_softclipped_fwd = 0
#     cdef int reads_duplicate = 0
#     cdef int reads_duplicate_rev = 0
#     cdef int reads_duplicate_fwd = 0
#
#     # initialise variables
#     n = col.n
#     plp = col.plp
#
#     # get chromosome name and position
#     tid = col.tid
#     chrom = alignmentfile.getrname(tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # loop over reads, extract what we need
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#         is_reverse = <bint>(flag & BAM_FREVERSE)
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         is_duplicate = <bint>(flag & BAM_FDUP)
#         mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
#         mate_is_reverse = <bint>(flag & BAM_FMREVERSE)
#         tlen = aln.core.isize
#         if is_reverse:
#             reads_rev += 1
#         else:
#             reads_fwd += 1
#         if is_proper_pair:
#             reads_pp += 1
#             if is_reverse:
#                 reads_pp_rev += 1
#             else:
#                 reads_pp_fwd += 1
#         if mate_is_unmapped:
#             reads_mate_unmapped += 1
#             if is_reverse:
#                 reads_mate_unmapped_rev += 1
#             else:
#                 reads_mate_unmapped_fwd += 1
#         elif tid != aln.core.mtid:
#             reads_mate_other_chr += 1
#             if is_reverse:
#                 reads_mate_other_chr_rev += 1
#             else:
#                 reads_mate_other_chr_fwd += 1
#         elif is_reverse and mate_is_reverse:
#             reads_mate_same_strand += 1
#             reads_mate_same_strand_rev += 1
#         elif not is_reverse and not mate_is_reverse:
#             reads_mate_same_strand += 1
#             reads_mate_same_strand_fwd += 1
#         elif (is_reverse and tlen > 0) or (not is_reverse and tlen < 0):
#             reads_faceaway += 1
#             if is_reverse:
#                 reads_faceaway_rev += 1
#             else:
#                 reads_faceaway_fwd += 1
#         if is_softclipped(aln):
#             reads_softclipped += 1
#             if is_reverse:
#                 reads_softclipped_rev += 1
#             else:
#                 reads_softclipped_fwd += 1
#         if is_duplicate:
#             reads_duplicate += 1
#             if is_reverse:
#                 reads_duplicate_rev += 1
#             else:
#                 reads_duplicate_fwd += 1
#
#     return {'chrom': chrom, 'pos': pos,
#            'reads_all': n,
#            'reads_fwd': reads_fwd,
#            'reads_rev': reads_rev,
#            'reads_pp': reads_pp,
#            'reads_pp_fwd': reads_pp_fwd,
#            'reads_pp_rev': reads_pp_rev,
#            'reads_mate_unmapped': reads_mate_unmapped,
#            'reads_mate_unmapped_fwd': reads_mate_unmapped_fwd,
#            'reads_mate_unmapped_rev': reads_mate_unmapped_rev,
#            'reads_mate_other_chr': reads_mate_other_chr,
#            'reads_mate_other_chr_fwd': reads_mate_other_chr_fwd,
#            'reads_mate_other_chr_rev': reads_mate_other_chr_rev,
#            'reads_mate_same_strand': reads_mate_same_strand,
#            'reads_mate_same_strand_fwd': reads_mate_same_strand_fwd,
#            'reads_mate_same_strand_rev': reads_mate_same_strand_rev,
#            'reads_faceaway': reads_faceaway,
#            'reads_faceaway_fwd': reads_faceaway_fwd,
#            'reads_faceaway_rev': reads_faceaway_rev,
#            'reads_softclipped': reads_softclipped,
#            'reads_softclipped_fwd': reads_softclipped_fwd,
#            'reads_softclipped_rev': reads_softclipped_rev,
#            'reads_duplicate': reads_duplicate,
#            'reads_duplicate_fwd': reads_duplicate_fwd,
#            'reads_duplicate_rev': reads_duplicate_rev,
#            }
#
#
# cpdef dict rec_coverage_ext_strand_pad(FastaFile fafile, chrom, pos,
#                                        bint one_based=False):
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom, 'pos': pos,
#            'reads_all': 0,
#            'reads_fwd': 0,
#            'reads_rev': 0,
#            'reads_pp': 0,
#            'reads_pp_fwd': 0,
#            'reads_pp_rev': 0,
#            'reads_mate_unmapped': 0,
#            'reads_mate_unmapped_fwd': 0,
#            'reads_mate_unmapped_rev': 0,
#            'reads_mate_other_chr': 0,
#            'reads_mate_other_chr_fwd': 0,
#            'reads_mate_other_chr_rev': 0,
#            'reads_mate_same_strand': 0,
#            'reads_mate_same_strand_fwd': 0,
#            'reads_mate_same_strand_rev': 0,
#            'reads_faceaway': 0,
#            'reads_faceaway_fwd': 0,
#            'reads_faceaway_rev': 0,
#            'reads_softclipped': 0,
#            'reads_softclipped_fwd': 0,
#            'reads_softclipped_rev': 0,
#            'reads_duplicate': 0,
#            'reads_duplicate_fwd': 0,
#            'reads_duplicate_rev': 0,
#            }
#
#
# ########################
# # VARIATION STATISTICS #
# ########################
#
#
# cpdef dict rec_variation(AlignmentFile alignmentfile, FastaFile fafile,
#                          PileupColumn col, bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i  # loop index
#     cdef int reads_all  # total number of reads in column
#     cdef uint32_t flag
#     cdef bint is_proper_pair
#     cdef bytes alnbase, refbase_b
#     # counting variables
#     cdef int reads_pp = 0
#     cdef int matches = 0
#     cdef int matches_pp = 0
#     cdef int mismatches = 0
#     cdef int mismatches_pp = 0
#     cdef int deletions = 0
#     cdef int deletions_pp = 0
#     cdef int insertions = 0
#     cdef int insertions_pp = 0
#     cdef int a = 0
#     cdef int a_pp = 0
#     cdef int c = 0
#     cdef int c_pp = 0
#     cdef int t = 0
#     cdef int t_pp = 0
#     cdef int g = 0
#     cdef int g_pp = 0
#     cdef int n = 0
#     cdef int n_pp = 0
#
#     # initialise variables
#     reads_all = col.n
#     plp = col.plp
#
#     # get chromosome name and position
#     chrom = alignmentfile.getrname(col.tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # reference base
#     refbase = fafile\
#         .fetch(reference=chrom, start=col.pos, end=col.pos+1)\
#         .upper()
#     if not PY2:
#         refbase_b = refbase.encode('ascii')
#     else:
#         refbase_b = refbase
#
#     # loop over reads, extract what we need
#     for i in range(reads_all):
#         read = &(plp[0][i])
#         # read.qpos
#         # read.is_del
#         # read.indel
#         aln = read.b
#         flag = aln.core.flag
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         if is_proper_pair:
#             reads_pp += 1
#         if read.is_del:
#             deletions += 1
#             if is_proper_pair:
#                 deletions_pp += 1
#         else:
# #            alnbase = get_seq_range(aln, 0, aln.core.l_qseq)[read.qpos]
#             alnbase = get_seq_base(aln, read.qpos)
#             if alnbase == b'A':
#                 a += 1
#                 if is_proper_pair:
#                     a_pp += 1
#             elif alnbase == b'T':
#                 t += 1
#                 if is_proper_pair:
#                     t_pp += 1
#             elif alnbase == b'C':
#                 c += 1
#                 if is_proper_pair:
#                     c_pp += 1
#             elif alnbase == b'G':
#                 g += 1
#                 if is_proper_pair:
#                     g_pp += 1
#             elif alnbase == b'N':
#                 n += 1
#                 if is_proper_pair:
#                     n_pp += 1
#             if read.indel > 0:
#                 insertions += 1
#                 if is_proper_pair:
#                     insertions_pp += 1
#             if alnbase == refbase_b:
#                 matches += 1
#                 if is_proper_pair:
#                     matches_pp += 1
#             else:
#                 mismatches += 1
#                 if is_proper_pair:
#                     mismatches_pp += 1
#
#     return {'chrom': chrom, 'pos': pos, 'ref': refbase,
#             'reads_all': reads_all, 'reads_pp': reads_pp,
#             'matches': matches,
#             'matches_pp': matches_pp,
#             'mismatches': mismatches,
#             'mismatches_pp': mismatches_pp,
#             'deletions': deletions,
#             'deletions_pp': deletions_pp,
#             'insertions': insertions,
#             'insertions_pp': insertions_pp,
#             'A': a, 'A_pp': a_pp,
#             'C': c, 'C_pp': c_pp,
#             'T': t, 'T_pp': t_pp,
#             'G': g, 'G_pp': g_pp,
#             'N': n, 'N_pp': n_pp}
#
#
# cpdef dict rec_variation_pad(FastaFile fafile, chrom, pos,
#                              bint one_based=False):
#     refbase = fafile.fetch(reference=chrom, start=pos, end=pos+1).upper()
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom, 'pos': pos, 'ref': refbase,
#             'reads_all': 0, 'reads_pp': 0,
#             'matches': 0,
#             'matches_pp': 0,
#             'mismatches': 0,
#             'mismatches_pp': 0,
#             'deletions': 0,
#             'deletions_pp': 0,
#             'insertions': 0,
#             'insertions_pp': 0,
#             'A': 0, 'A_pp': 0,
#             'C': 0, 'C_pp': 0,
#             'T': 0, 'T_pp': 0,
#             'G': 0, 'G_pp': 0,
#             'N': 0, 'N_pp': 0}
#
#
# #################################
# # STRANDED VARIATION STATISTICS #
# #################################
#
#
# cdef struct CountPpStrand:
#     int all, pp, fwd, rev, pp_fwd, pp_rev
#
#
# cdef inline init_pp_strand(CountPpStrand* c):
#     c.all = c.fwd = c.rev = c.pp = c.pp_fwd = c.pp_rev = 0
#
#
# cdef inline incr_pp_strand(CountPpStrand* c, bint is_reverse,
#                            bint is_proper_pair):
#     c.all += 1
#     if is_reverse:
#         c.rev += 1
#         if is_proper_pair:
#             c.pp += 1
#             c.pp_rev += 1
#     else:
#         c.fwd += 1
#         if is_proper_pair:
#             c.pp += 1
#             c.pp_fwd += 1
#
#
# cpdef dict rec_variation_strand(AlignmentFile alignmentfile, FastaFile fafile,
#                                 PileupColumn col, bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i, n # loop index
#     cdef int reads_all # total number of reads in column
#     cdef uint32_t flag
#     cdef bint is_proper_pair, is_reverse
#     cdef bytes alnbase, refbase_b
#     # counting variables
#     cdef CountPpStrand reads, matches, mismatches, deletions, insertions, \
#         A, C, T, G, N
#
#     # initialise variables
#     n = col.n
#     plp = col.plp
#     init_pp_strand(&reads)
#     init_pp_strand(&matches)
#     init_pp_strand(&mismatches)
#     init_pp_strand(&deletions)
#     init_pp_strand(&insertions)
#     init_pp_strand(&A)
#     init_pp_strand(&T)
#     init_pp_strand(&C)
#     init_pp_strand(&G)
#     init_pp_strand(&N)
#
#     # get chromosome name and position
#     chrom = alignmentfile.getrname(col.tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # reference base
#     refbase = fafile.fetch(reference=chrom, start=col.pos, end=col.pos+1).upper()
#     if not PY2:
#         refbase_b = refbase.encode('ascii')
#     else:
#         refbase_b = refbase
#
#     # loop over reads, extract what we need
#     for i in range(n):
#         read = &(plp[0][i])
#         # read.qpos
#         # read.is_del
#         # read.indel
#         aln = read.b
#         flag = aln.core.flag
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         is_reverse = <bint>(flag & BAM_FREVERSE)
#         incr_pp_strand(&reads, is_reverse, is_proper_pair)
#         if read.is_del:
#             incr_pp_strand(&deletions, is_reverse, is_proper_pair)
#         else:
#             alnbase = get_seq_base(aln, read.qpos)
#             if alnbase == b'A':
#                 incr_pp_strand(&A, is_reverse, is_proper_pair)
#             elif alnbase == b'T':
#                 incr_pp_strand(&T, is_reverse, is_proper_pair)
#             elif alnbase == b'C':
#                 incr_pp_strand(&C, is_reverse, is_proper_pair)
#             elif alnbase == b'G':
#                 incr_pp_strand(&G, is_reverse, is_proper_pair)
#             elif alnbase == b'N':
#                 incr_pp_strand(&N, is_reverse, is_proper_pair)
#             if read.indel > 0:
#                 incr_pp_strand(&insertions, is_reverse, is_proper_pair)
#             if alnbase == refbase_b:
#                 incr_pp_strand(&matches, is_reverse, is_proper_pair)
#             else:
#                 incr_pp_strand(&mismatches, is_reverse, is_proper_pair)
#
#     return {'chrom': chrom, 'pos': pos, 'ref': refbase,
#             'reads_all': n, 'reads_fwd': reads.fwd, 'reads_rev': reads.rev,
#             'reads_pp': reads.pp, 'reads_pp_fwd': reads.pp_fwd, 'reads_pp_rev': reads.pp_rev,
#             'matches': matches.all, 'matches_fwd': matches.fwd, 'matches_rev': matches.rev,
#             'matches_pp': matches.pp, 'matches_pp_fwd': matches.pp_fwd, 'matches_pp_rev': matches.pp_rev,
#             'mismatches': mismatches.all, 'mismatches_fwd': mismatches.fwd, 'mismatches_rev': mismatches.rev,
#             'mismatches_pp': mismatches.pp, 'mismatches_pp_fwd': mismatches.pp_fwd, 'mismatches_pp_rev': mismatches.pp_rev,
#             'deletions': deletions.all, 'deletions_fwd': deletions.fwd, 'deletions_rev': deletions.rev,
#             'deletions_pp': deletions.pp, 'deletions_pp_fwd': deletions.pp_fwd, 'deletions_pp_rev': deletions.pp_rev,
#             'insertions': insertions.all, 'insertions_fwd': insertions.fwd, 'insertions_rev': insertions.rev,
#             'insertions_pp': insertions.pp, 'insertions_pp_fwd': insertions.pp_fwd, 'insertions_pp_rev': insertions.pp_rev,
#             'A': A.all, 'A_fwd': A.fwd, 'A_rev': A.rev, 'A_pp': A.pp, 'A_pp_fwd': A.pp_fwd, 'A_pp_rev': A.pp_rev,
#             'C': C.all, 'C_fwd': C.fwd, 'C_rev': C.rev, 'C_pp': C.pp, 'C_pp_fwd': C.pp_fwd, 'C_pp_rev': C.pp_rev,
#             'T': T.all, 'T_fwd': T.fwd, 'T_rev': T.rev, 'T_pp': T.pp, 'T_pp_fwd': T.pp_fwd, 'T_pp_rev': T.pp_rev,
#             'G': G.all, 'G_fwd': G.fwd, 'G_rev': G.rev, 'G_pp': G.pp, 'G_pp_fwd': G.pp_fwd, 'G_pp_rev': G.pp_rev,
#             'N': N.all, 'N_fwd': N.fwd, 'N_rev': N.rev, 'N_pp': N.pp, 'N_pp_fwd': N.pp_fwd, 'N_pp_rev': N.pp_rev,
#             }
#
#
# cpdef dict rec_variation_strand_pad(FastaFile fafile, chrom, pos,
#                                     bint one_based=False):
#     refbase = fafile.fetch(reference=chrom, start=pos, end=pos+1).upper()
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom, 'pos': pos, 'ref': refbase,
#             'reads_all': 0, 'reads_fwd': 0, 'reads_rev': 0,
#             'reads_pp': 0, 'reads_pp_fwd': 0, 'reads_pp_rev': 0,
#             'matches': 0, 'matches_fwd': 0, 'matches_rev': 0,
#             'matches_pp': 0, 'matches_pp_fwd': 0, 'matches_pp_rev': 0,
#             'mismatches': 0, 'mismatches_fwd': 0, 'mismatches_rev': 0,
#             'mismatches_pp': 0, 'mismatches_pp_fwd': 0, 'mismatches_pp_rev': 0,
#             'deletions': 0, 'deletions_fwd': 0, 'deletions_rev': 0,
#             'deletions_pp': 0, 'deletions_pp_fwd': 0, 'deletions_pp_rev': 0,
#             'insertions': 0, 'insertions_fwd': 0, 'insertions_rev': 0,
#             'insertions_pp': 0, 'insertions_pp_fwd': 0, 'insertions_pp_rev': 0,
#             'A': 0, 'A_fwd': 0, 'A_rev': 0,
#             'A_pp': 0, 'A_pp_fwd': 0, 'A_pp_rev': 0,
#             'C': 0, 'C_fwd': 0, 'C_rev': 0,
#             'C_pp': 0, 'C_pp_fwd': 0, 'C_pp_rev': 0,
#             'T': 0, 'T_fwd': 0, 'T_rev': 0,
#             'T_pp': 0, 'T_pp_fwd': 0, 'T_pp_rev': 0,
#             'G': 0, 'G_fwd': 0, 'G_rev': 0,
#             'G_pp': 0, 'G_pp_fwd': 0, 'G_pp_rev': 0,
#             'N': 0, 'N_fwd': 0, 'N_rev': 0,
#             'N_pp': 0, 'N_pp_fwd': 0, 'N_pp_rev': 0,
#             }
#
#
# ##########################
# # INSERT SIZE STATISTICS #
# ##########################
#
#
# cpdef dict rec_tlen(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
#                     bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i, n  # loop index
#     cdef int reads_all  # total number of reads in column
#     cdef uint32_t flag
#     cdef bint is_proper_pair
#     cdef bint mate_is_unmappped
#     cdef bint mate_other_chr
#     # reads "paired", i.e., mate is mapped to same chromosome, so tlen is meaningful
#     cdef int reads_p = 0
#     # reads "properly paired", as defined by aligner
#     cdef int reads_pp = 0
#     cdef int64_t tlen
#     cdef int64_t tlen_squared
#     cdef int64_t tlen_p_sum = 0
#     cdef double tlen_p_mean = 0
#     cdef double tlen_p_dev_squared
#     cdef double tlen_p_dev_squared_sum = 0
#     cdef int64_t tlen_p_squared_sum = 0
#     cdef int64_t tlen_pp_sum = 0
#     cdef double tlen_pp_mean = 0
#     cdef double tlen_pp_dev_squared
#     cdef double tlen_pp_dev_squared_sum = 0
#     cdef int64_t tlen_pp_squared_sum = 0
#
#     # initialise variables
#     n = col.n
#     plp = col.plp
#
#     # get chromosome name and position
#     chrom = alignmentfile.getrname(col.tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # loop over reads
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
#         mate_other_chr = <bint>(aln.core.tid != aln.core.mtid)
#
#         # N.B., pysam exposes this property as 'tlen' rather than 'isize' so we
#         # follow their naming convention
#         tlen = aln.core.isize
#         tlen_squared = tlen**2
#
#         # N.B. insert size is only meaningful if mate is mapped to same chromosome
#         if not mate_is_unmapped and not mate_other_chr:
#             reads_p += 1
#             tlen_p_sum += tlen
#             tlen_p_squared_sum += tlen_squared
#             if is_proper_pair:
#                 reads_pp += 1
#                 tlen_pp_sum += tlen
#                 tlen_pp_squared_sum += tlen_squared
#
#     # calculate intermediate variables
#     if reads_p > 0:
#         tlen_p_mean = tlen_p_sum / reads_p
#     if reads_pp > 0:
#         tlen_pp_mean = tlen_pp_sum / reads_pp
#
#     # loop over reads again to calculate variance (and hence std)
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
#         mate_other_chr = <bint>(aln.core.tid != aln.core.mtid)
#         tlen = aln.core.isize
#         # N.B. insert size is only meaningful if mate is mapped to same chromosome
#         if not mate_is_unmapped and not mate_other_chr:
#             tlen_p_dev_squared = (tlen - tlen_p_mean)**2
#             tlen_p_dev_squared_sum += tlen_p_dev_squared
#             if is_proper_pair:
#                 tlen_pp_dev_squared = (tlen - tlen_pp_mean)**2
#                 tlen_pp_dev_squared_sum += tlen_pp_dev_squared
#
#     # calculate output variables
#     # N.B. round values to nearest integer, any finer precision is probably not
#     # interesting
#     if reads_p > 0:
#         mean_tlen = int(round(tlen_p_mean))
#         rms_tlen = rootmean(tlen_p_squared_sum, reads_p)
#         variance_tlen = tlen_p_dev_squared_sum / reads_p
#         std_tlen = int(round(sqrt(variance_tlen)))
#     else:
#         rms_tlen = std_tlen = mean_tlen = median_tlen = 0
#     if reads_pp > 0:
#         mean_tlen_pp = int(round(tlen_pp_mean))
#         rms_tlen_pp = rootmean(tlen_pp_squared_sum, reads_pp)
#         variance_tlen_pp = tlen_pp_dev_squared_sum / reads_pp
#         std_tlen_pp = int(round(sqrt(variance_tlen_pp)))
#     else:
#         rms_tlen_pp = std_tlen_pp = mean_tlen_pp = 0
#
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': n,
#             'reads_paired': reads_p,
#             'reads_pp': reads_pp,
#             'mean_tlen': mean_tlen,
#             'mean_tlen_pp': mean_tlen_pp,
#             'rms_tlen': rms_tlen,
#             'rms_tlen_pp': rms_tlen_pp,
#             'std_tlen': std_tlen,
#             'std_tlen_pp': std_tlen_pp}
#
#
# cpdef dict rec_tlen_pad(FastaFile fafile, chrom, pos, bint one_based=False):
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': 0,
#             'reads_paired': 0,
#             'reads_pp': 0,
#             'mean_tlen': 0,
#             'mean_tlen_pp': 0,
#             'rms_tlen': 0,
#             'rms_tlen_pp': 0,
#             'std_tlen': 0,
#             'std_tlen_pp': 0,
#             }
#
#
# ####################################
# # INSERT SIZE STATISTICS BY STRAND #
# ####################################
#
#
# cpdef dict rec_tlen_strand(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
#                            bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i, n # loop index
#     cdef int reads_all # total number of reads in column
#     cdef uint32_t flag
#     cdef bint is_proper_pair
#     cdef bint mate_is_unmappped
#     cdef bint mate_other_chr
#
#     # counting variables
#     cdef int reads_fwd = 0
#     cdef int reads_rev = 0
#     # reads "paired", i.e., mate is mapped to same chromosome, so tlen is
#     # meaningful
#     cdef int reads_p = 0
#     cdef int reads_p_fwd = 0
#     cdef int reads_p_rev = 0
#     # reads "properly paired", as defined by aligner
#     cdef int reads_pp = 0
#     cdef int reads_pp_fwd = 0
#     cdef int reads_pp_rev = 0
#
#     cdef int64_t tlen
#     cdef int64_t tlen_squared
#
#     cdef int64_t tlen_p_sum = 0
#     cdef double tlen_p_mean = 0
#     cdef double tlen_p_dev_squared
#     cdef double tlen_p_dev_squared_sum = 0
#     cdef int64_t tlen_p_squared_sum = 0
#     cdef int64_t tlen_p_fwd_sum = 0
#     cdef double tlen_p_fwd_mean = 0
#     cdef double tlen_p_fwd_dev_squared
#     cdef double tlen_p_fwd_dev_squared_sum = 0
#     cdef int64_t tlen_p_fwd_squared_sum = 0
#     cdef int64_t tlen_p_rev_sum = 0
#     cdef double tlen_p_rev_mean = 0
#     cdef double tlen_p_rev_dev_squared
#     cdef double tlen_p_rev_dev_squared_sum = 0
#     cdef int64_t tlen_p_rev_squared_sum = 0
#
#     cdef int64_t tlen_pp_sum = 0
#     cdef double tlen_pp_mean = 0
#     cdef double tlen_pp_dev_squared
#     cdef double tlen_pp_dev_squared_sum = 0
#     cdef int64_t tlen_pp_squared_sum = 0
#     cdef int64_t tlen_pp_fwd_sum = 0
#     cdef double tlen_pp_fwd_mean = 0
#     cdef double tlen_pp_fwd_dev_squared
#     cdef double tlen_pp_fwd_dev_squared_sum = 0
#     cdef int64_t tlen_pp_fwd_squared_sum = 0
#     cdef int64_t tlen_pp_rev_sum = 0
#     cdef double tlen_pp_rev_mean = 0
#     cdef double tlen_pp_rev_dev_squared
#     cdef double tlen_pp_rev_dev_squared_sum = 0
#     cdef int64_t tlen_pp_rev_squared_sum = 0
#
#     # initialise variables
#     n = col.n
#     plp = col.plp
#
#     # get chromosome name and position
#     chrom = alignmentfile.getrname(col.tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # loop over reads
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         is_reverse = <bint>(flag & BAM_FREVERSE)
#         mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
#         mate_other_chr = <bint>(aln.core.tid != aln.core.mtid)
#
#         # not sure these are really needed
#         if is_reverse:
#             reads_rev += 1
#         else:
#             reads_fwd += 1
#
#         # N.B., pysam exposes this property as 'tlen' rather than 'isize' so we
#         # follow their naming convention
#         tlen = aln.core.isize
#         tlen_squared = tlen**2
#
#         # N.B. insert size is only meaningful if mate is mapped to same
#         # chromosome
#         if not mate_is_unmapped and not mate_other_chr:
#             reads_p += 1
#             tlen_p_sum += tlen
#             tlen_p_squared_sum += tlen_squared
#             if is_reverse:
#                 reads_p_rev += 1
#                 tlen_p_rev_sum += tlen
#                 tlen_p_rev_squared_sum += tlen_squared
#             else:
#                 reads_p_fwd += 1
#                 tlen_p_fwd_sum += tlen
#                 tlen_p_fwd_squared_sum += tlen_squared
#
#             if is_proper_pair:
#                 reads_pp += 1
#                 tlen_pp_sum += tlen
#                 tlen_pp_squared_sum += tlen_squared
#                 if is_reverse:
#                     reads_pp_rev += 1
#                     tlen_pp_rev_sum += tlen
#                     tlen_pp_rev_squared_sum += tlen_squared
#                 else:
#                     reads_pp_fwd += 1
#                     tlen_pp_fwd_sum += tlen
#                     tlen_pp_fwd_squared_sum += tlen_squared
#
#     # calculate intermediate variables
#     if reads_p > 0:
#         tlen_p_mean = tlen_p_sum / reads_p
#         if reads_p_rev > 0:
#             tlen_p_rev_mean = tlen_p_rev_sum / reads_p_rev
#         if reads_p_fwd > 0:
#             tlen_p_fwd_mean = tlen_p_fwd_sum / reads_p_fwd
#     if reads_pp > 0:
#         tlen_pp_mean = tlen_pp_sum / reads_pp
#         if reads_pp_rev > 0:
#             tlen_pp_rev_mean = tlen_pp_rev_sum / reads_pp_rev
#         if reads_pp_fwd > 0:
#             tlen_pp_fwd_mean = tlen_pp_fwd_sum / reads_pp_fwd
#
#     # loop over reads again to calculate variance (and hence std)
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         is_reverse = <bint>(flag & BAM_FREVERSE)
#         mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
#         mate_other_chr = <bint>(aln.core.tid != aln.core.mtid)
#         tlen = aln.core.isize
#         # N.B. insert size is only meaningful if mate is mapped to same
#         # chromosome
#         if not mate_is_unmapped and not mate_other_chr:
#             tlen_p_dev_squared = (tlen - tlen_p_mean)**2
#             tlen_p_dev_squared_sum += tlen_p_dev_squared
#             if is_reverse:
#                 tlen_p_rev_dev_squared = (tlen - tlen_p_rev_mean)**2
#                 tlen_p_rev_dev_squared_sum += tlen_p_rev_dev_squared
#             else:
#                 tlen_p_fwd_dev_squared = (tlen - tlen_p_fwd_mean)**2
#                 tlen_p_fwd_dev_squared_sum += tlen_p_fwd_dev_squared
#             if is_proper_pair:
#                 tlen_pp_dev_squared = (tlen - tlen_pp_mean)**2
#                 tlen_pp_dev_squared_sum += tlen_pp_dev_squared
#                 if is_reverse:
#                     tlen_pp_rev_dev_squared = (tlen - tlen_pp_rev_mean)**2
#                     tlen_pp_rev_dev_squared_sum += tlen_pp_rev_dev_squared
#                 else:
#                     tlen_pp_fwd_dev_squared = (tlen - tlen_pp_fwd_mean)**2
#                     tlen_pp_fwd_dev_squared_sum += tlen_pp_fwd_dev_squared
#
#     # calculate output variables
#     # N.B. round values to nearest integer, any finer precision is probably not
#     # interesting
#     if reads_p > 0:
#         mean_tlen = int(round(tlen_p_mean))
#         rms_tlen = rootmean(tlen_p_squared_sum, reads_p)
#         variance_tlen = tlen_p_dev_squared_sum / reads_p
#         std_tlen = int(round(sqrt(variance_tlen)))
#     else:
#         rms_tlen = std_tlen = mean_tlen = 0
#     if reads_p_rev > 0:
#         mean_tlen_rev = int(round(tlen_p_rev_mean))
#         rms_tlen_rev = rootmean(tlen_p_rev_squared_sum, reads_p_rev)
#         variance_tlen_rev = tlen_p_rev_dev_squared_sum / reads_p_rev
#         std_tlen_rev = int(round(sqrt(variance_tlen_rev)))
#     else:
#         rms_tlen_rev = std_tlen_rev = mean_tlen_rev = 0
#     if reads_p_fwd > 0:
#         mean_tlen_fwd = int(round(tlen_p_fwd_mean))
#         rms_tlen_fwd = rootmean(tlen_p_fwd_squared_sum, reads_p_fwd)
#         variance_tlen_fwd = tlen_p_fwd_dev_squared_sum / reads_p_fwd
#         std_tlen_fwd = int(round(sqrt(variance_tlen_fwd)))
#     else:
#         rms_tlen_fwd = std_tlen_fwd = mean_tlen_fwd = 0
#     if reads_pp > 0:
#         mean_tlen_pp = int(round(tlen_pp_mean))
#         rms_tlen_pp = rootmean(tlen_pp_squared_sum, reads_pp)
#         variance_tlen_pp = tlen_pp_dev_squared_sum / reads_pp
#         std_tlen_pp = int(round(sqrt(variance_tlen_pp)))
#     else:
#         rms_tlen_pp = std_tlen_pp = mean_tlen_pp = 0
#     if reads_pp_rev > 0:
#         mean_tlen_pp_rev = int(round(tlen_pp_rev_mean))
#         rms_tlen_pp_rev = rootmean(tlen_pp_rev_squared_sum, reads_pp_rev)
#         variance_tlen_pp_rev = tlen_pp_rev_dev_squared_sum / reads_pp_rev
#         std_tlen_pp_rev = int(round(sqrt(variance_tlen_pp_rev)))
#     else:
#         rms_tlen_pp_rev = std_tlen_pp_rev = mean_tlen_pp_rev = 0
#     if reads_pp_fwd > 0:
#         mean_tlen_pp_fwd = int(round(tlen_pp_fwd_mean))
#         rms_tlen_pp_fwd = rootmean(tlen_pp_fwd_squared_sum, reads_pp_fwd)
#         variance_tlen_pp_fwd = tlen_pp_fwd_dev_squared_sum / reads_pp_fwd
#         std_tlen_pp_fwd = int(round(sqrt(variance_tlen_pp_fwd)))
#     else:
#         rms_tlen_pp_fwd = std_tlen_pp_fwd = mean_tlen_pp_fwd = 0
#
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': n,
#             'reads_fwd': reads_fwd,
#             'reads_rev': reads_rev,
#             'reads_paired': reads_p,
#             'reads_paired_fwd': reads_p_fwd,
#             'reads_paired_rev': reads_p_rev,
#             'reads_pp': reads_pp,
#             'reads_pp_fwd': reads_pp_fwd,
#             'reads_pp_rev': reads_pp_rev,
#             'mean_tlen': mean_tlen,
#             'mean_tlen_fwd': mean_tlen_fwd,
#             'mean_tlen_rev': mean_tlen_rev,
#             'mean_tlen_pp': mean_tlen_pp,
#             'mean_tlen_pp_fwd': mean_tlen_pp_fwd,
#             'mean_tlen_pp_rev': mean_tlen_pp_rev,
#             'rms_tlen': rms_tlen,
#             'rms_tlen_fwd': rms_tlen_fwd,
#             'rms_tlen_rev': rms_tlen_rev,
#             'rms_tlen_pp': rms_tlen_pp,
#             'rms_tlen_pp_fwd': rms_tlen_pp_fwd,
#             'rms_tlen_pp_rev': rms_tlen_pp_rev,
#             'std_tlen': std_tlen,
#             'std_tlen_fwd': std_tlen_fwd,
#             'std_tlen_rev': std_tlen_rev,
#             'std_tlen_pp': std_tlen_pp,
#             'std_tlen_pp_fwd': std_tlen_pp_fwd,
#             'std_tlen_pp_rev': std_tlen_pp_rev}
#
#
# cpdef dict rec_tlen_strand_pad(FastaFile fafile, chrom, pos, bint one_based=False):
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': 0,
#             'reads_fwd': 0,
#             'reads_rev': 0,
#             'reads_paired': 0,
#             'reads_paired_fwd': 0,
#             'reads_paired_rev': 0,
#             'reads_pp': 0,
#             'reads_pp_fwd': 0,
#             'reads_pp_rev': 0,
#             'mean_tlen': 0,
#             'mean_tlen_fwd': 0,
#             'mean_tlen_rev': 0,
#             'mean_tlen_pp': 0,
#             'mean_tlen_pp_fwd': 0,
#             'mean_tlen_pp_rev': 0,
#             'rms_tlen': 0,
#             'rms_tlen_fwd': 0,
#             'rms_tlen_rev': 0,
#             'rms_tlen_pp': 0,
#             'rms_tlen_pp_fwd': 0,
#             'rms_tlen_pp_rev': 0,
#             'std_tlen': 0,
#             'std_tlen_fwd': 0,
#             'std_tlen_rev': 0,
#             'std_tlen_pp': 0,
#             'std_tlen_pp_fwd': 0,
#             'std_tlen_pp_rev': 0,
#             }
#
#
# ##############################
# # MAPPING QUALITY STATISTICS #
# ##############################
#
#
# cpdef dict rec_mapq(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
#                     bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i, n # loop index
#     cdef int reads_all # total number of reads in column
#     cdef uint32_t flag
#     cdef uint64_t mapq
#     cdef uint64_t mapq_max = 0
#     cdef uint64_t mapq_pp_max = 0
#     cdef uint64_t mapq_squared
#     cdef uint64_t mapq_squared_sum = 0
#     cdef uint64_t mapq_pp_squared_sum = 0
#     cdef bint is_proper_pair
#     cdef int reads_pp = 0
#     cdef int reads_mapq0 = 0
#     cdef int reads_mapq0_pp = 0
#
#     # initialise variables
#     n = col.n
#     plp = col.plp
#
#     # get chromosome name and position
#     chrom = alignmentfile.getrname(col.tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # loop over reads, extract what we need
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         mapq = aln.core.qual
#         mapq_squared = mapq**2
#         mapq_squared_sum += mapq_squared
#         if mapq == 0:
#             reads_mapq0 += 1
#         if mapq > mapq_max:
#             mapq_max = mapq
#         if is_proper_pair:
#             reads_pp += 1
#             mapq_pp_squared_sum += mapq_squared
#             if mapq > mapq_pp_max:
#                 mapq_pp_max = mapq
#             if mapq == 0:
#                 reads_mapq0_pp += 1
#
#     # construct output variables
#     rms_mapq = rootmean(mapq_squared_sum, n)
#     max_mapq = mapq_max
#     if reads_pp > 0:
#         rms_mapq_pp = rootmean(mapq_pp_squared_sum, reads_pp)
#         max_mapq_pp = mapq_pp_max
#     else:
#         rms_mapq_pp = max_mapq_pp = 0
#
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': n,
#             'reads_pp': reads_pp,
#             'reads_mapq0': reads_mapq0,
#             'reads_mapq0_pp': reads_mapq0_pp,
#             'rms_mapq': rms_mapq,
#             'rms_mapq_pp': rms_mapq_pp,
#             'max_mapq': max_mapq,
#             'max_mapq_pp': max_mapq_pp}
#
#
# cpdef dict rec_mapq_pad(FastaFile fafile, chrom, pos, bint one_based=False):
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': 0,
#             'reads_pp': 0,
#             'reads_mapq0': 0,
#             'reads_mapq0_pp': 0,
#             'rms_mapq': 0,
#             'rms_mapq_pp': 0,
#             'max_mapq': 0,
#             'max_mapq_pp': 0,
#             }
#
#
# ########################################
# # MAPPING QUALITY STATISTICS BY STRAND #
# ########################################
#
#
# cpdef dict rec_mapq_strand(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
#                            bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i, n  # loop index
#     cdef int reads_all  # total number of reads in column
#
#     cdef uint32_t flag
#     cdef uint64_t mapq
#     cdef uint64_t mapq_squared
#     cdef bint is_proper_pair
#     cdef bint is_reverse
#
#     cdef uint64_t mapq_max = 0
#     cdef uint64_t mapq_rev_max = 0
#     cdef uint64_t mapq_fwd_max = 0
#     cdef uint64_t mapq_pp_max = 0
#     cdef uint64_t mapq_pp_rev_max = 0
#     cdef uint64_t mapq_pp_fwd_max = 0
#     cdef uint64_t mapq_squared_sum = 0
#     cdef uint64_t mapq_fwd_squared_sum = 0
#     cdef uint64_t mapq_rev_squared_sum = 0
#     cdef uint64_t mapq_pp_squared_sum = 0
#     cdef uint64_t mapq_pp_fwd_squared_sum = 0
#     cdef uint64_t mapq_pp_rev_squared_sum = 0
#
#     cdef int reads_rev = 0
#     cdef int reads_fwd = 0
#     cdef int reads_pp = 0
#     cdef int reads_pp_rev = 0
#     cdef int reads_pp_fwd = 0
#     cdef int reads_mapq0 = 0
#     cdef int reads_mapq0_fwd = 0
#     cdef int reads_mapq0_rev = 0
#     cdef int reads_mapq0_pp = 0
#     cdef int reads_mapq0_pp_fwd = 0
#     cdef int reads_mapq0_pp_rev = 0
#
#     # initialise variables
#     n = col.n
#     plp = col.plp
#
#     # get chromosome name and position
#     chrom = alignmentfile.getrname(col.tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # loop over reads, extract what we need
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         is_reverse = <bint>(flag & BAM_FREVERSE)
#         mapq = aln.core.qual
#         mapq_squared = mapq**2
#
#         mapq_squared_sum += mapq_squared
#         if mapq > mapq_max:
#             mapq_max = mapq
#         if is_reverse:
#             reads_rev += 1
#             mapq_rev_squared_sum += mapq_squared
#             if mapq > mapq_rev_max:
#                 mapq_rev_max = mapq
#         else:
#             reads_fwd += 1
#             mapq_fwd_squared_sum += mapq_squared
#             if mapq > mapq_fwd_max:
#                 mapq_fwd_max = mapq
#
#         if mapq == 0:
#             reads_mapq0 += 1
#             if is_reverse:
#                 reads_mapq0_rev += 1
#             else:
#                 reads_mapq0_fwd += 1
#
#         if is_proper_pair:
#             reads_pp += 1
#             mapq_pp_squared_sum += mapq_squared
#             if mapq > mapq_pp_max:
#                 mapq_pp_max = mapq
#             if is_reverse:
#                 reads_pp_rev += 1
#                 mapq_pp_rev_squared_sum += mapq_squared
#                 if mapq > mapq_pp_rev_max:
#                     mapq_pp_rev_max = mapq
#             else:
#                 reads_pp_fwd += 1
#                 mapq_pp_fwd_squared_sum += mapq_squared
#                 if mapq > mapq_pp_fwd_max:
#                     mapq_pp_fwd_max = mapq
#             if mapq == 0:
#                 reads_mapq0_pp += 1
#                 if is_reverse:
#                     reads_mapq0_pp_rev += 1
#                 else:
#                     reads_mapq0_pp_fwd += 1
#
#     # construct output variables
#     rms_mapq = rootmean(mapq_squared_sum, n)
#     max_mapq = mapq_max
#     if reads_rev > 0:
#         rms_mapq_rev = rootmean(mapq_rev_squared_sum, reads_rev)
#         max_mapq_rev = mapq_rev_max
#     else:
#         rms_mapq_rev = max_mapq_rev = 0
#     if reads_fwd > 0:
#         rms_mapq_fwd = rootmean(mapq_fwd_squared_sum, reads_fwd)
#         max_mapq_fwd = mapq_fwd_max
#     else:
#         rms_mapq_fwd = max_mapq_fwd = 0
#     if reads_pp > 0:
#         rms_mapq_pp = rootmean(mapq_pp_squared_sum, reads_pp)
#         max_mapq_pp = mapq_pp_max
#     else:
#         rms_mapq_pp = max_mapq_pp = 0
#     if reads_pp_fwd > 0:
#         rms_mapq_pp_fwd = rootmean(mapq_pp_fwd_squared_sum, reads_pp_fwd)
#         max_mapq_pp_fwd = mapq_pp_fwd_max
#     else:
#         rms_mapq_pp_fwd = max_mapq_pp_fwd = 0
#     if reads_pp_rev > 0:
#         rms_mapq_pp_rev = rootmean(mapq_pp_rev_squared_sum, reads_pp_rev)
#         max_mapq_pp_rev = mapq_pp_rev_max
#     else:
#         rms_mapq_pp_rev = max_mapq_pp_rev = 0
#
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': n,
#             'reads_fwd': reads_fwd,
#             'reads_rev': reads_rev,
#             'reads_pp': reads_pp,
#             'reads_pp_fwd': reads_pp_fwd,
#             'reads_pp_rev': reads_pp_rev,
#             'reads_mapq0': reads_mapq0,
#             'reads_mapq0_fwd': reads_mapq0_fwd,
#             'reads_mapq0_rev': reads_mapq0_rev,
#             'reads_mapq0_pp': reads_mapq0_pp,
#             'reads_mapq0_pp_fwd': reads_mapq0_pp_fwd,
#             'reads_mapq0_pp_rev': reads_mapq0_pp_rev,
#             'rms_mapq': rms_mapq,
#             'rms_mapq_fwd': rms_mapq_fwd,
#             'rms_mapq_rev': rms_mapq_rev,
#             'rms_mapq_pp': rms_mapq_pp,
#             'rms_mapq_pp_fwd': rms_mapq_pp_fwd,
#             'rms_mapq_pp_rev': rms_mapq_pp_rev,
#             'max_mapq': max_mapq,
#             'max_mapq_fwd': max_mapq_fwd,
#             'max_mapq_rev': max_mapq_rev,
#             'max_mapq_pp': max_mapq_pp,
#             'max_mapq_pp_fwd': max_mapq_pp_fwd,
#             'max_mapq_pp_rev': max_mapq_pp_rev,
#             }
#
#
# cpdef dict rec_mapq_strand_pad(FastaFile fafile, chrom, pos, bint one_based=False):
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': 0,
#             'reads_fwd': 0,
#             'reads_rev': 0,
#             'reads_pp': 0,
#             'reads_pp_fwd': 0,
#             'reads_pp_rev': 0,
#             'reads_mapq0': 0,
#             'reads_mapq0_fwd': 0,
#             'reads_mapq0_rev': 0,
#             'reads_mapq0_pp': 0,
#             'reads_mapq0_pp_fwd': 0,
#             'reads_mapq0_pp_rev': 0,
#             'rms_mapq': 0,
#             'rms_mapq_fwd': 0,
#             'rms_mapq_rev': 0,
#             'rms_mapq_pp': 0,
#             'rms_mapq_pp_fwd': 0,
#             'rms_mapq_pp_rev': 0,
#             'max_mapq': 0,
#             'max_mapq_fwd': 0,
#             'max_mapq_rev': 0,
#             'max_mapq_pp': 0,
#             'max_mapq_pp_fwd': 0,
#             'max_mapq_pp_rev': 0,
#             }
#
#
# ###########################
# # BASE QUALITY STATISTICS #
# ###########################
#
#
# cpdef dict rec_baseq(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
#                      bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i, n # loop index
#     cdef int reads_all # total number of reads in column
#     cdef uint32_t flag
#     cdef bint is_proper_pair
#     cdef int reads_nodel = 0
#     cdef int reads_pp = 0
#     cdef int reads_pp_nodel = 0
#     cdef uint64_t baseq, baseq_squared
#     cdef uint64_t baseq_squared_sum = 0
#     cdef uint64_t baseq_pp_squared_sum = 0
#
#     # initialise variables
#     n = col.n
#     plp = col.plp
#
#     # get chromosome name and position
#     chrom = alignmentfile.getrname(col.tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # loop over reads, extract what we need
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         if is_proper_pair:
#             reads_pp += 1
#         # N.B., base quality only makes sense if the aligned read is not a
#         # deletion
#         if not read.is_del:
#             reads_nodel += 1
#             baseq = pysam_bam_get_qual(aln)[read.qpos]
#             baseq_squared = baseq**2
#             baseq_squared_sum += baseq_squared
#             if is_proper_pair:
#                 reads_pp_nodel += 1
#                 baseq_pp_squared_sum += baseq_squared
#
#     # output variables
#     rms_baseq = rootmean(baseq_squared_sum, reads_nodel)
#     rms_baseq_pp = rootmean(baseq_pp_squared_sum, reads_pp_nodel)
#
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': n,
#             'reads_pp': reads_pp,
#             'rms_baseq': rms_baseq,
#             'rms_baseq_pp': rms_baseq_pp}
#
#
# cpdef dict rec_baseq_pad(FastaFile fafile, chrom, pos, bint one_based=False):
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': 0,
#             'reads_pp': 0,
#             'rms_baseq': 0,
#             'rms_baseq_pp': 0,
#             }
#
#
# #####################################
# # BASE QUALITY STATISTICS BY STRAND #
# #####################################
#
#
# cpdef dict rec_baseq_strand(AlignmentFile alignmentfile, FastaFile fafile, PileupColumn col,
#                             bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i, n # loop index
#     cdef int reads_all # total number of reads in column
#
#     cdef uint32_t flag
#     cdef bint is_proper_pair
#     cdef bint is_reverse
#     cdef uint64_t baseq
#     cdef uint64_t baseq_squared
#
#     cdef uint64_t baseq_squared_sum = 0
#     cdef uint64_t baseq_fwd_squared_sum = 0
#     cdef uint64_t baseq_rev_squared_sum = 0
#     cdef uint64_t baseq_pp_squared_sum = 0
#     cdef uint64_t baseq_pp_fwd_squared_sum = 0
#     cdef uint64_t baseq_pp_rev_squared_sum = 0
#
#     cdef int reads_rev = 0
#     cdef int reads_fwd = 0
#     cdef int reads_pp = 0
#     cdef int reads_pp_rev = 0
#     cdef int reads_pp_fwd = 0
#     cdef int reads_nodel = 0
#     cdef int reads_rev_nodel = 0
#     cdef int reads_fwd_nodel = 0
#     cdef int reads_pp_nodel = 0
#     cdef int reads_pp_rev_nodel = 0
#     cdef int reads_pp_fwd_nodel = 0
#
#     # initialise variables
#     n = col.n
#     plp = col.plp
#
#     # get chromosome name and position
#     chrom = alignmentfile.getrname(col.tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # loop over reads, extract what we need
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         is_reverse = <bint>(flag & BAM_FREVERSE)
#
#         if is_reverse:
#             reads_rev += 1
#         else:
#             reads_fwd += 1
#         if is_proper_pair:
#             reads_pp += 1
#             if is_reverse:
#                 reads_pp_rev += 1
#             else:
#                 reads_pp_fwd += 1
#
#         # N.B., baseq only makes sense if not a deletion
#         if not read.is_del:
#             reads_nodel += 1
#             baseq = pysam_bam_get_qual(aln)[read.qpos]
#             baseq_squared = baseq**2
#             baseq_squared_sum += baseq_squared
#             if is_reverse:
#                 reads_rev_nodel += 1
#                 baseq_rev_squared_sum += baseq_squared
#             else:
#                 reads_fwd_nodel += 1
#                 baseq_fwd_squared_sum += baseq_squared
#             if is_proper_pair:
#                 reads_pp_nodel += 1
#                 baseq_pp_squared_sum += baseq_squared
#                 if is_reverse:
#                     reads_pp_rev_nodel += 1
#                     baseq_pp_rev_squared_sum += baseq_squared
#                 else:
#                     reads_pp_fwd_nodel += 1
#                     baseq_pp_fwd_squared_sum += baseq_squared
#
#     # construct output variables
#     rms_baseq = rootmean(baseq_squared_sum, reads_nodel)
#     rms_baseq_rev = rootmean(baseq_rev_squared_sum, reads_rev_nodel)
#     rms_baseq_fwd = rootmean(baseq_fwd_squared_sum, reads_fwd_nodel)
#     rms_baseq_pp = rootmean(baseq_pp_squared_sum, reads_pp_nodel)
#     rms_baseq_pp_fwd = rootmean(baseq_pp_fwd_squared_sum, reads_pp_fwd_nodel)
#     rms_baseq_pp_rev = rootmean(baseq_pp_rev_squared_sum, reads_pp_rev_nodel)
#
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': n,
#             'reads_fwd': reads_fwd,
#             'reads_rev': reads_rev,
#             'reads_pp': reads_pp,
#             'reads_pp_fwd': reads_pp_fwd,
#             'reads_pp_rev': reads_pp_rev,
#             'rms_baseq': rms_baseq,
#             'rms_baseq_fwd': rms_baseq_fwd,
#             'rms_baseq_rev': rms_baseq_rev,
#             'rms_baseq_pp': rms_baseq_pp,
#             'rms_baseq_pp_fwd': rms_baseq_pp_fwd,
#             'rms_baseq_pp_rev': rms_baseq_pp_rev,
#             }
#
#
# cpdef dict rec_baseq_strand_pad(FastaFile fafile, chrom, pos, bint one_based=False):
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom,
#             'pos': pos,
#             'reads_all': 0,
#             'reads_fwd': 0,
#             'reads_rev': 0,
#             'reads_pp': 0,
#             'reads_pp_fwd': 0,
#             'reads_pp_rev': 0,
#             'rms_baseq': 0,
#             'rms_baseq_fwd': 0,
#             'rms_baseq_rev': 0,
#             'rms_baseq_pp': 0,
#             'rms_baseq_pp_fwd': 0,
#             'rms_baseq_pp_rev': 0,
#             }
#
#
# ####################################
# # EXTENDED BASE QUALITY STATISTICS #
# ####################################
#
#
# cpdef dict rec_baseq_ext(AlignmentFile alignmentfile, FastaFile fafile,
#                          PileupColumn col, bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i, n # loop index
#     cdef int reads_all # total number of reads in column
#     cdef uint32_t flag
#     cdef bint is_proper_pair
#     cdef bytes refbase_b, alnbase
#     # counting variables
#     cdef int reads_nodel = 0
#     cdef int reads_pp = 0
#     cdef int reads_pp_nodel = 0
#     cdef int matches = 0
#     cdef int matches_pp = 0
#     cdef int mismatches = 0
#     cdef int mismatches_pp = 0
#
#     cdef uint64_t baseq
#     cdef uint64_t baseq_squared
#
#     cdef uint64_t baseq_squared_sum = 0
#     cdef uint64_t baseq_pp_squared_sum = 0
#     cdef uint64_t baseq_matches_squared_sum = 0
#     cdef uint64_t baseq_matches_pp_squared_sum = 0
#     cdef uint64_t baseq_mismatches_squared_sum = 0
#     cdef uint64_t baseq_mismatches_pp_squared_sum = 0
#
#     # initialise variables
#     n = col.n
#     plp = col.plp
#
#     # get chromosome name and position
#     chrom = alignmentfile.getrname(col.tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # reference base
#     refbase = fafile.fetch(reference=chrom, start=col.pos,
#                            end=col.pos + 1).upper()
#     if not PY2:
#         refbase_b = refbase.encode('ascii')
#     else:
#         refbase_b = refbase
#
#     # loop over reads, extract what we need
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         if is_proper_pair:
#             reads_pp += 1
#         if not read.is_del:
#             reads_nodel += 1
#             baseq = pysam_bam_get_qual(aln)[read.qpos]
#             baseq_squared = baseq**2
#             baseq_squared_sum += baseq_squared
#             if is_proper_pair:
#                 reads_pp_nodel += 1
#                 baseq_pp_squared_sum += baseq_squared
#             alnbase = get_seq_base(aln, read.qpos)
#             if alnbase == refbase_b:
#                 matches += 1
#                 baseq_matches_squared_sum += baseq_squared
#                 if is_proper_pair:
#                     matches_pp += 1
#                     baseq_matches_pp_squared_sum += baseq_squared
#             else:
#                 mismatches += 1
#                 baseq_mismatches_squared_sum += baseq_squared
#                 if is_proper_pair:
#                     mismatches_pp += 1
#                     baseq_mismatches_pp_squared_sum += baseq_squared
#
#     # construct output variables
#     rms_baseq = rootmean(baseq_squared_sum, reads_nodel)
#     rms_baseq_pp = rootmean(baseq_pp_squared_sum, reads_pp_nodel)
#     rms_baseq_matches = rootmean(baseq_matches_squared_sum, matches)
#     rms_baseq_matches_pp = rootmean(baseq_matches_pp_squared_sum, matches_pp)
#     rms_baseq_mismatches = rootmean(baseq_mismatches_squared_sum, mismatches)
#     rms_baseq_mismatches_pp = rootmean(baseq_mismatches_pp_squared_sum,
#                                        mismatches_pp)
#
#     return {'chrom': chrom, 'pos': pos, 'ref': refbase,
#             'reads_all': n, 'reads_pp': reads_pp,
#             'matches': matches,
#             'matches_pp': matches_pp,
#             'mismatches': mismatches,
#             'mismatches_pp': mismatches_pp,
#             'rms_baseq': rms_baseq,
#             'rms_baseq_pp': rms_baseq_pp,
#             'rms_baseq_matches': rms_baseq_matches,
#             'rms_baseq_matches_pp': rms_baseq_matches_pp,
#             'rms_baseq_mismatches': rms_baseq_mismatches,
#             'rms_baseq_mismatches_pp': rms_baseq_mismatches_pp,
#             }
#
#
# cpdef dict rec_baseq_ext_pad(FastaFile fafile, chrom, pos, bint one_based=False):
#     refbase = fafile.fetch(reference=chrom, start=pos, end=pos+1).upper()
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom, 'pos': pos, 'ref': refbase,
#             'reads_all': 0, 'reads_pp': 0,
#             'matches': 0,
#             'matches_pp': 0,
#             'mismatches': 0,
#             'mismatches_pp': 0,
#             'rms_baseq': 0,
#             'rms_baseq_pp': 0,
#             'rms_baseq_matches': 0,
#             'rms_baseq_matches_pp': 0,
#             'rms_baseq_mismatches': 0,
#             'rms_baseq_mismatches_pp': 0,
#             }
#
#
# ##############################################
# # EXTENDED BASE QUALITY STATISTICS BY STRAND #
# ##############################################
#
#
# cpdef dict rec_baseq_ext_strand(AlignmentFile alignmentfile, FastaFile fafile,
#                                 PileupColumn col, bint one_based=False):
#
#     # statically typed variables
#     cdef bam_pileup1_t** plp
#     cdef bam_pileup1_t* read
#     cdef bam1_t* aln
#     cdef int i, n # loop index
#     cdef int reads_all # total number of reads in column
#     cdef uint32_t flag
#     cdef bint is_proper_pair
#     cdef bint is_reverse
#     cdef bytes alnbase, refbase_b
#     # counting variables
#     cdef int reads_fwd = 0
#     cdef int reads_rev = 0
#     cdef int reads_nodel = 0
#     cdef int reads_fwd_nodel = 0
#     cdef int reads_rev_nodel = 0
#     cdef int reads_pp = 0
#     cdef int reads_pp_fwd = 0
#     cdef int reads_pp_rev = 0
#     cdef int reads_pp_nodel = 0
#     cdef int reads_pp_fwd_nodel = 0
#     cdef int reads_pp_rev_nodel = 0
#     cdef int matches = 0
#     cdef int matches_fwd = 0
#     cdef int matches_rev = 0
#     cdef int matches_pp = 0
#     cdef int matches_pp_fwd = 0
#     cdef int matches_pp_rev = 0
#     cdef int mismatches = 0
#     cdef int mismatches_fwd = 0
#     cdef int mismatches_rev = 0
#     cdef int mismatches_pp = 0
#     cdef int mismatches_pp_fwd = 0
#     cdef int mismatches_pp_rev = 0
#
#     cdef uint64_t baseq
#     cdef uint64_t baseq_squared
#
#     cdef uint64_t baseq_squared_sum = 0
#     cdef uint64_t baseq_fwd_squared_sum = 0
#     cdef uint64_t baseq_rev_squared_sum = 0
#     cdef uint64_t baseq_pp_squared_sum = 0
#     cdef uint64_t baseq_pp_fwd_squared_sum = 0
#     cdef uint64_t baseq_pp_rev_squared_sum = 0
#     cdef uint64_t baseq_matches_squared_sum = 0
#     cdef uint64_t baseq_matches_fwd_squared_sum = 0
#     cdef uint64_t baseq_matches_rev_squared_sum = 0
#     cdef uint64_t baseq_matches_pp_squared_sum = 0
#     cdef uint64_t baseq_matches_pp_fwd_squared_sum = 0
#     cdef uint64_t baseq_matches_pp_rev_squared_sum = 0
#     cdef uint64_t baseq_mismatches_squared_sum = 0
#     cdef uint64_t baseq_mismatches_fwd_squared_sum = 0
#     cdef uint64_t baseq_mismatches_rev_squared_sum = 0
#     cdef uint64_t baseq_mismatches_pp_squared_sum = 0
#     cdef uint64_t baseq_mismatches_pp_fwd_squared_sum = 0
#     cdef uint64_t baseq_mismatches_pp_rev_squared_sum = 0
#
#     # initialise variables
#     n = col.n
#     plp = col.plp
#
#     # get chromosome name and position
#     chrom = alignmentfile.getrname(col.tid)
#     pos = col.pos + 1 if one_based else col.pos
#
#     # reference base
#     refbase = fafile.fetch(reference=chrom, start=col.pos, end=col.pos+1).upper()
#     if not PY2:
#         refbase_b = refbase.encode('ascii')
#     else:
#         refbase_b = refbase
#
#     # loop over reads, extract what we need
#     for i in range(n):
#         read = &(plp[0][i])
#         aln = read.b
#         flag = aln.core.flag
#         is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
#         is_reverse = <bint>(flag & BAM_FREVERSE)
#
#         if is_reverse:
#             reads_rev += 1
#         else:
#             reads_fwd += 1
#         if is_proper_pair:
#             reads_pp += 1
#             if is_reverse:
#                 reads_pp_rev += 1
#             else:
#                 reads_pp_fwd += 1
#
#         if not read.is_del:
#             reads_nodel += 1
#             baseq = pysam_bam_get_qual(aln)[read.qpos]
#             baseq_squared = baseq**2
#             baseq_squared_sum += baseq_squared
#             if is_reverse:
#                 reads_rev_nodel += 1
#                 baseq_rev_squared_sum += baseq_squared
#             else:
#                 reads_fwd_nodel += 1
#                 baseq_fwd_squared_sum += baseq_squared
#             if is_proper_pair:
#                 reads_pp_nodel += 1
#                 baseq_pp_squared_sum += baseq_squared
#                 if is_reverse:
#                     reads_pp_rev_nodel += 1
#                     baseq_pp_rev_squared_sum += baseq_squared
#                 else:
#                     reads_pp_fwd_nodel += 1
#                     baseq_pp_fwd_squared_sum += baseq_squared
#             alnbase = get_seq_base(aln, read.qpos)
#             if alnbase == refbase_b:
#                 matches += 1
#                 baseq_matches_squared_sum += baseq_squared
#                 if is_reverse:
#                     matches_rev += 1
#                     baseq_matches_rev_squared_sum += baseq_squared
#                 else:
#                     matches_fwd += 1
#                     baseq_matches_fwd_squared_sum += baseq_squared
#
#                 if is_proper_pair:
#                     matches_pp += 1
#                     baseq_matches_pp_squared_sum += baseq_squared
#                     if is_reverse:
#                         matches_pp_rev += 1
#                         baseq_matches_pp_rev_squared_sum += baseq_squared
#                     else:
#                         matches_pp_fwd += 1
#                         baseq_matches_pp_fwd_squared_sum += baseq_squared
#             else:
#                 mismatches += 1
#                 baseq_mismatches_squared_sum += baseq_squared
#                 if is_reverse:
#                     mismatches_rev += 1
#                     baseq_mismatches_rev_squared_sum += baseq_squared
#                 else:
#                     mismatches_fwd += 1
#                     baseq_mismatches_fwd_squared_sum += baseq_squared
#
#                 if is_proper_pair:
#                     mismatches_pp += 1
#                     baseq_mismatches_pp_squared_sum += baseq_squared
#                     if is_reverse:
#                         mismatches_pp_rev += 1
#                         baseq_mismatches_pp_rev_squared_sum += baseq_squared
#                     else:
#                         mismatches_pp_fwd += 1
#                         baseq_mismatches_pp_fwd_squared_sum += baseq_squared
#
#     # construct output variables
#     rms_baseq = rootmean(baseq_squared_sum, reads_nodel)
#     rms_baseq_fwd = rootmean(baseq_fwd_squared_sum, reads_fwd_nodel)
#     rms_baseq_rev = rootmean(baseq_rev_squared_sum, reads_rev_nodel)
#     rms_baseq_pp = rootmean(baseq_pp_squared_sum, reads_pp_nodel)
#     rms_baseq_pp_fwd = rootmean(baseq_pp_fwd_squared_sum, reads_pp_fwd_nodel)
#     rms_baseq_pp_rev = rootmean(baseq_pp_rev_squared_sum, reads_pp_rev_nodel)
#     rms_baseq_matches = rootmean(baseq_matches_squared_sum, matches)
#     rms_baseq_matches_fwd = rootmean(baseq_matches_fwd_squared_sum,
#                                      matches_fwd)
#     rms_baseq_matches_rev = rootmean(baseq_matches_rev_squared_sum,
#                                      matches_rev)
#     rms_baseq_matches_pp = rootmean(baseq_matches_pp_squared_sum, matches_pp)
#     rms_baseq_matches_pp_fwd = rootmean(baseq_matches_pp_fwd_squared_sum,
#                                         matches_pp_fwd)
#     rms_baseq_matches_pp_rev = rootmean(baseq_matches_pp_rev_squared_sum,
#                                         matches_pp_rev)
#     rms_baseq_mismatches = rootmean(baseq_mismatches_squared_sum, mismatches)
#     rms_baseq_mismatches_fwd = rootmean(baseq_mismatches_fwd_squared_sum,
#                                         mismatches_fwd)
#     rms_baseq_mismatches_rev = rootmean(baseq_mismatches_rev_squared_sum,
#                                         mismatches_rev)
#     rms_baseq_mismatches_pp = rootmean(baseq_mismatches_pp_squared_sum,
#                                        mismatches_pp)
#     rms_baseq_mismatches_pp_fwd = rootmean(baseq_mismatches_pp_fwd_squared_sum,
#                                            mismatches_pp_fwd)
#     rms_baseq_mismatches_pp_rev = rootmean(baseq_mismatches_pp_rev_squared_sum,
#                                            mismatches_pp_rev)
#
#     return {'chrom': chrom, 'pos': pos, 'ref': refbase,
#             'reads_all': n, 'reads_fwd': reads_fwd, 'reads_rev': reads_rev,
#             'reads_pp': reads_pp, 'reads_pp_fwd': reads_pp_fwd, 'reads_pp_rev': reads_pp_rev,
#             'matches': matches, 'matches_fwd': matches_fwd, 'matches_rev': matches_rev,
#             'matches_pp': matches_pp, 'matches_pp_fwd': matches_pp_fwd, 'matches_pp_rev': matches_pp_rev,
#             'mismatches': mismatches, 'mismatches_fwd': mismatches_fwd, 'mismatches_rev': mismatches_rev,
#             'mismatches_pp': mismatches_pp, 'mismatches_pp_fwd': mismatches_pp_fwd, 'mismatches_pp_rev': mismatches_pp_rev,
#             'rms_baseq': rms_baseq, 'rms_baseq_fwd': rms_baseq_fwd, 'rms_baseq_rev': rms_baseq_rev,
#             'rms_baseq_pp': rms_baseq_pp, 'rms_baseq_pp_fwd': rms_baseq_pp_fwd, 'rms_baseq_pp_rev': rms_baseq_pp_rev,
#             'rms_baseq_matches': rms_baseq_matches, 'rms_baseq_matches_fwd': rms_baseq_matches_fwd, 'rms_baseq_matches_rev': rms_baseq_matches_rev,
#             'rms_baseq_matches_pp': rms_baseq_matches_pp, 'rms_baseq_matches_pp_fwd': rms_baseq_matches_pp_fwd, 'rms_baseq_matches_pp_rev': rms_baseq_matches_pp_rev,
#             'rms_baseq_mismatches': rms_baseq_mismatches, 'rms_baseq_mismatches_fwd': rms_baseq_mismatches_fwd, 'rms_baseq_mismatches_rev': rms_baseq_mismatches_rev,
#             'rms_baseq_mismatches_pp': rms_baseq_mismatches_pp, 'rms_baseq_mismatches_pp_fwd': rms_baseq_mismatches_pp_fwd, 'rms_baseq_mismatches_pp_rev': rms_baseq_mismatches_pp_rev,
#             }
#
#
# cpdef dict rec_baseq_ext_strand_pad(FastaFile fafile, chrom, pos, bint one_based=False):
#     refbase = fafile.fetch(reference=chrom, start=pos, end=pos+1).upper()
#     pos = pos + 1 if one_based else pos
#     return {'chrom': chrom, 'pos': pos, 'ref': refbase,
#             'reads_all': 0, 'reads_fwd': 0, 'reads_rev': 0,
#             'reads_pp': 0, 'reads_pp_fwd': 0, 'reads_pp_rev': 0,
#             'matches': 0, 'matches_fwd': 0, 'matches_rev': 0,
#             'matches_pp': 0, 'matches_pp_fwd': 0, 'matches_pp_rev': 0,
#             'mismatches': 0, 'mismatches_fwd': 0, 'mismatches_rev': 0,
#             'mismatches_pp': 0, 'mismatches_pp_fwd': 0, 'mismatches_pp_rev': 0,
#             'rms_baseq': 0, 'rms_baseq_fwd': 0, 'rms_baseq_rev': 0,
#             'rms_baseq_pp': 0, 'rms_baseq_pp_fwd': 0, 'rms_baseq_pp_rev': 0,
#             'rms_baseq_matches': 0, 'rms_baseq_matches_fwd': 0, 'rms_baseq_matches_rev': 0,
#             'rms_baseq_matches_pp': 0, 'rms_baseq_matches_pp_fwd': 0, 'rms_baseq_matches_pp_rev': 0,
#             'rms_baseq_mismatches': 0, 'rms_baseq_mismatches_fwd': 0, 'rms_baseq_mismatches_rev': 0,
#             'rms_baseq_mismatches_pp': 0, 'rms_baseq_mismatches_pp_fwd': 0, 'rms_baseq_mismatches_pp_rev': 0,
#             }
#
#
# #################################################
# # BASIC COVERAGE STATISTICS WITH GC COMPOSITION #
# #################################################
#
#
# def frecs_coverage_gc(window_size=300, window_offset=None):
#
#     if window_offset is None:
#         window_offset = window_size / 2
#
#     def rec_coverage_gc(AlignmentFile alignmentfile, FastaFile fafile,
#                          PileupColumn col, bint one_based):
#         chrom = alignmentfile.getrname(col.tid)
#
#         ref_window_start = col.pos - window_offset
#         ref_window_end = ref_window_start + window_size
#         if ref_window_start < 0:
#             ref_window_start = 0
#         ref_window = fafile\
#             .fetch(chrom, ref_window_start, ref_window_end)\
#             .lower()
#         if len(ref_window) == 0:
#             gc = -1
#         else:
#             gc = gc_content(ref_window)
#
#         rec = rec_coverage(alignmentfile, fafile, col, one_based)
#         rec['gc'] = gc
#         return rec
#
#     def rec_coverage_gc_pad(FastaFile fafile, chrom, pos,
#                              bint one_based):
#         ref_window_start = pos - window_offset
#         ref_window_end = ref_window_start + window_size
#         if ref_window_start < 0:
#             ref_window_start = 0
#         ref_window = fafile\
#             .fetch(chrom, ref_window_start, ref_window_end)\
#             .lower()
#         if len(ref_window) == 0:
#             gc = -1
#         else:
#             gc = gc_content(ref_window)
#         rec = rec_coverage_pad(fafile, chrom, pos, one_based)
#         rec['gc'] = gc
#         return rec
#
#     return rec_coverage_gc, rec_coverage_gc_pad
#

###################
# BINNED COVERAGE #
###################


cdef int gc_content(ref_window) except -1:
    cdef Py_ssize_t i, n
    cdef char* seq
    cdef int gc_count = 0
    n = len(ref_window)
    if not PY2:
        ref_window = ref_window.encode('ascii')
    seq = ref_window
    for i in range(n):
        if seq[i] == b'g' or seq[i] == b'c':
            gc_count += 1
    gc = int(round(gc_count * 100. / n))
    return gc


cdef class BinnedStat(object):

    cdef dict rec(self, chrom, bin_start, bin_end, FastaFile fafile):
        return dict()

    cdef void recv(self, bam1_t* b):
        pass


cdef class CoverageBinned(BinnedStat):

    cdef int reads_all, reads_pp

    def __init__(self):
        self.reads_all = self.reads_pp = 0

    cdef dict rec(self, chrom, bin_start, bin_end, FastaFile fafile):

        # determine %GC
        ref_window = fafile.fetch(chrom, bin_start, bin_end).lower()
        gc = gc_content(ref_window)

        # make record for bin
        rec = {'gc': gc,
               'reads_all': self.reads_all,
               'reads_pp': self.reads_pp}

        # reset counters
        self.reads_all = self.reads_pp = 0

        return rec

    cdef void recv(self, bam1_t* b):
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


############################################
# BINNED COVERAGE WITH EXTENDED PROPERTIES #
############################################


cdef class CoverageExtBinned(BinnedStat):

    cdef int reads_all, reads_pp, reads_mate_unmapped, reads_mate_other_chr, \
        reads_mate_same_strand, reads_faceaway, reads_softclipped, \
        reads_duplicate

    def __init__(self):
        self.reads_all = self.reads_pp = self.reads_mate_unmapped \
            = self.reads_mate_other_chr = self.reads_mate_same_strand\
            = self.reads_faceaway = self.reads_softclipped \
            = self.reads_duplicate = 0

    cdef dict rec(self, chrom, bin_start, bin_end, FastaFile fafile):

        # determine %GC
        ref_window = fafile.fetch(chrom, bin_start, bin_end).lower()
        gc = gc_content(ref_window)

        # make record for bin
        rec = {'gc': gc,
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

    cdef void recv(self, bam1_t* b):
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
            if is_softclipped(b):
                self.reads_softclipped += 1


###############
# BINNED MAPQ #
###############


cdef class MapqBinned(BinnedStat):

    cdef int reads_all, reads_mapq0
    cdef uint64_t mapq, mapq_squared_sum

    def __init__(self):
        self.reads_all = self.reads_mapq0 = self.mapq_squared_sum = 0

    cdef dict rec(self, chrom, bin_start, bin_end, FastaFile fafile):

        # make record for bin
        rec = {'reads_all': self.reads_all,
               'reads_mapq0': self.reads_mapq0,
               'rms_mapq': rootmean(self.mapq_squared_sum, self.reads_all)}

        # reset counters
        self.reads_all = self.reads_mapq0 = self.mapq_squared_sum = 0

        return rec

    cdef void recv(self, bam1_t* b):
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


################
# BINNED CIGAR #
################


cdef class AlignmentBinned(BinnedStat):

    cdef int reads_all, M, I, D, N, S, H, P, EQ, X

    def __init__(self):
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

    cdef void recv(self, bam1_t* b):
        cdef uint32_t flag
        cdef bint is_unmapped
        cdef int k, op, l
        flag = b.core.flag
        is_unmapped = <bint>(flag & BAM_FUNMAP)
        if not is_unmapped:
            cigar_p = pysam_bam_get_cigar(b)
            # cigar = list()
            for k in range(b.core.n_cigar):
                op = cigar_p[k] & BAM_CIGAR_MASK
                l = cigar_p[k] >> BAM_CIGAR_SHIFT
                # cigar.append((op, l))
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


###############
# BINNED TLEN #
###############


cdef class TlenBinned(BinnedStat):

    cdef int reads_all
    cdef int reads_pp
    cdef int64_t tlen_sum
    cdef int64_t tlen_pp_sum
    cdef int64_t tlen_squared_sum
    cdef int64_t tlen_pp_squared_sum

    def __init__(self):
        self.tlen_sum = self.tlen_squared_sum = self.tlen_pp_sum = \
            self.tlen_pp_squared_sum = self.reads_all = self.reads_pp = 0

    cdef dict rec(self, chrom, bin_start, bin_end, FastaFile fafile):

        # make record for bin
        rec = {'reads_all': self.reads_all, 'reads_pp': self.reads_pp,
               'mean_tlen': mean(self.tlen_sum, self.reads_all),
               'mean_tlen_pp': mean(self.tlen_pp_sum, self.reads_pp),
               'rms_tlen': rootmean(self.tlen_squared_sum, self.reads_all),
               'rms_tlen_pp': rootmean(self.tlen_pp_squared_sum,
                                       self.reads_pp)}

        # reset counters
        self.tlen_sum = self.tlen_squared_sum = self.tlen_pp_sum = \
            self.tlen_pp_squared_sum = self.reads_all = self.reads_pp = 0

        return rec

    cdef void recv(self, bam1_t* b):
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


#####################
# UTILITY FUNCTIONS #
#####################


def iter_pileup(stat, alignmentfile, fafile, pad, **kwargs):
    """General purpose function to generate statistics, where one record is
    generated for each genome position in the selected range, based on a
    pileup column.

    """

    if isinstance(alignmentfile, _string_types):
        alignmentfile = AlignmentFile(alignmentfile)

    if isinstance(fafile, _string_types):
        fafile = FastaFile(fafile)

    if pad:
        return iter_pileup_padded(stat, alignmentfile=alignmentfile, fafile=fafile, **kwargs)
    else:
        return iter_pileup_default(stat, alignmentfile=alignmentfile, fafile=fafile, **kwargs)


def iter_pileup_default(stat, alignmentfile, fafile, chrom, start, end, one_based, truncate,
                        max_depth, int min_mapq, int min_baseq, bint no_del, bint no_dup):
    cdef:
        PileupColumn col

    # obtain pileup iterator
    start, end = normalise_coords(alignmentfile, chrom, start, end, one_based)
    it = alignmentfile.pileup(reference=chrom, start=start, end=end, truncate=truncate,
                              max_depth=max_depth)

    # iterate over pileup columns
    for col in it:

        rec = stat_pileup(stat, col, fafile=fafile, chrom=chrom, one_based=one_based,
                          min_mapq=min_mapq, min_baseq=min_baseq, no_del=no_del,
                          no_dup=no_dup)

        yield rec


def stat_pileup(PileupStat stat, PileupColumn col, fafile, chrom, bint one_based, int min_mapq,
                int min_baseq, bint no_del, bint no_dup):
    cdef:
        bam_pileup1_t** plp
        bam_pileup1_t* read
        bam1_t* b
        int i, n

    n = col.n
    plp = col.plp

    # iterate over reads in the column
    for i in range(n):
        read = &(plp[0][i])
        b = read.b

        # read filters
        if min_mapq > 0:
            if b.core.qual < min_mapq:
                continue
        if min_baseq > 0:
            if not read.is_del:
                if pysam_bam_get_qual(read.b)[read.qpos] < min_baseq:
                    continue
        if no_del:
            if read.is_del:
                continue
        if no_dup:
            if b.core.flag & BAM_FDUP:
                continue

        # accumulate statistics
        stat.recv(col, b)

    # construct record
    # TODO check this doesn't foul-up GC window centre
    pos = col.pos + 1 if one_based else col.pos
    rec = stat.rec(chrom, pos, fafile)
    rec['chrom'] = chrom
    rec['pos'] = pos
    return rec


def iter_pileup_padded(stat, alignmentfile, fafile, chrom, **kwargs):
    if chrom is not None:
        it = iter_pileup_padded_chrom(stat, alignmentfile=alignmentfile, fafile=fafile,
                                      chrom=chrom, **kwargs)
    else:
        its = list()
        for chrom in alignmentfile.references:
            itc = iter_pileup_padded_chrom(stat, alignmentfile=alignmentfile, fafile=fafile,
                                           chrom=chrom, **kwargs)
            its.append(itc)
        it = itertools.chain(*its)
    return it


def iter_pileup_padded_chrom(stat, alignmentfile, fafile, chrom, start, end,
                             one_based, truncate, max_depth, min_mapq, min_baseq, no_del, no_dup):
    cdef:
        PileupColumn col
        int curpos

    # obtain pileup iterator
    assert chrom is not None, 'chromosome is None'
    start, end = normalise_coords(alignmentfile, chrom, start, end, one_based)
    it = alignmentfile.pileup(reference=chrom, start=start, end=end, truncate=truncate,
                              max_depth=max_depth)

    # keep track of current position
    curpos = start

    # iterate over pileup columns
    for col in it:

        # pad up to next pileup column
        while curpos < col.pos:
            # construct record
            rec = stat.rec(chrom, curpos, fafile)
            rec['chrom'] = chrom
            rec['pos'] = curpos
            yield rec
            curpos += 1

        # construct and yield statistics record for current pileup column
        rec = stat_pileup(stat, col, fafile=fafile, chrom=chrom, one_based=one_based,
                          min_mapq=min_mapq, min_baseq=min_baseq, no_del=no_del,
                          no_dup=no_dup)
        yield rec

        # advance current position
        curpos = col.pos + 1

    while curpos < end:

        # pad out to end position
        rec = stat.rec(chrom, curpos, fafile)
        rec['chrom'] = chrom
        rec['pos'] = curpos
        yield rec
        curpos += 1


def iter_binned(stat, alignmentfile, fafile, chrom, window_size=300, window_offset=None,
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
            itc = iter_binned_chrom(stat, alignmentfile=alignmentfile, fafile=fafile,
                                    chrom=chrom, window_size=window_size,
                                    window_offset=window_offset, **kwargs)
            its.append(itc)
        it = itertools.chain(*its)

    else:
        it = iter_binned_chrom(stat, alignmentfile=alignmentfile, fafile=fafile,
                               chrom=chrom, window_size=window_size, window_offset=window_offset,
                               **kwargs)

    return it


def iter_binned_chrom(BinnedStat stat, AlignmentFile alignmentfile, FastaFile fafile,
                      chrom, start, end, one_based, int window_size, int window_offset,
                      int min_mapq, int no_dup):

    # setup

    assert chrom is not None, 'chromosome is None'
    start, end = normalise_coords(alignmentfile, chrom, start, end, one_based)

    cdef int rtid, rstart, rend, has_coord, bin_start, bin_end
    has_coord, rtid, rstart, rend = alignmentfile.parse_region(chrom, start, end, None)

    cdef IteratorRowRegion it
    cdef bam1_t* b
    it = IteratorRowRegion(alignmentfile, rtid, rstart, rend, multiple_iterators=False)
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
        if min_mapq > 0:
            if b.core.qual < min_mapq:
                it.cnext()
                continue
        if no_dup:
            if <bint>(b.core.flag & BAM_FDUP):
                it.cnext()
                continue
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


def normalise_coords(AlignmentFile alignmentfile, chrom, start, end, one_based):
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


cdef inline bint is_softclipped(bam1_t* aln):
    cdef int k
    cigar_p = pysam_bam_get_cigar(aln)
    for k in range(aln.core.n_cigar):
        op = cigar_p[k] & BAM_CIGAR_MASK
        if op == BAM_CSOFT_CLIP:
            return 1
    return 0


cdef inline object get_seq_base(bam1_t *src, uint32_t k):
    cdef uint8_t* p
    cdef char* s

    if not src.core.l_qseq:
        return None

    seq = PyBytes_FromStringAndSize(NULL, 1)
    s   = <char*>seq
    p   = pysam_bam_get_seq(src)

    s[0] = bam_nt16_rev_table[p[k//2] >> 4 * (1 - k%2) & 0xf]

    return seq


cdef inline int rootmean(uint64_t sqsum, int count):
    if count > 0:
        return int(round(sqrt(sqsum / count)))
    else:
        return 0


cdef inline int mean(int64_t total, int count):
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
