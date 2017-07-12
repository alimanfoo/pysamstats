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

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):
        return dict()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
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

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            bint is_proper_pair

        # convenience variables
        is_proper_pair = <bint>(read.b.core.flag & BAM_FPROPER_PAIR)

        # do the counting
        self.reads_all += 1
        if is_proper_pair:
            self.reads_pp += 1

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

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

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bint is_reverse

        # convenience variables
        flag = read.b.core.flag
        is_reverse = <bint>(flag & BAM_FREVERSE)
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)

        # do the counting
        self.reads_all += 1
        if is_reverse:
            self.reads_rev += 1
        else:
            self.reads_fwd += 1
        if is_proper_pair:
            self.reads_pp += 1
            if is_reverse:
                self.reads_pp_rev += 1
            else:
                self.reads_pp_fwd += 1

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

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

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bint is_reverse
            bint is_duplicate
            bint mate_is_unmappped
            bint mate_is_reverse
            int tlen

        # convenience variables
        flag = read.b.core.flag
        is_reverse = <bint>(flag & BAM_FREVERSE)
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_duplicate = <bint>(flag & BAM_FDUP)
        mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
        mate_is_reverse = <bint>(flag & BAM_FMREVERSE)
        tlen = read.b.core.isize

        # do the counting
        self.reads_all += 1
        if is_duplicate:
            self.reads_duplicate += 1
        if is_proper_pair:
            self.reads_pp += 1
        if mate_is_unmapped:
            self.reads_mate_unmapped += 1
        elif col.tid != read.b.core.mtid:
            self.reads_mate_other_chr += 1
        elif (is_reverse and mate_is_reverse) or (not is_reverse and not mate_is_reverse):
            self.reads_mate_same_strand += 1
        elif (is_reverse and tlen > 0) or (not is_reverse and tlen < 0):
            self.reads_faceaway += 1
        if is_softclipped(read.b):
            self.reads_softclipped += 1

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

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


##########################################
# EXTENDED COVERAGE STATISTICS BY STRAND #
##########################################


# noinspection PyAttributeOutsideInit
cdef class CoverageExtStrand(PileupStat):

    cdef:
        int reads_all
        int reads_rev
        int reads_fwd
        int reads_pp
        int reads_pp_rev
        int reads_pp_fwd
        int reads_mate_unmapped
        int reads_mate_unmapped_rev
        int reads_mate_unmapped_fwd
        int reads_mate_other_chr
        int reads_mate_other_chr_rev
        int reads_mate_other_chr_fwd
        int reads_mate_same_strand
        int reads_mate_same_strand_rev
        int reads_mate_same_strand_fwd
        int reads_faceaway
        int reads_faceaway_rev
        int reads_faceaway_fwd
        int reads_softclipped
        int reads_softclipped_rev
        int reads_softclipped_fwd
        int reads_duplicate
        int reads_duplicate_rev
        int reads_duplicate_fwd

    def __init__(self):
        self.reset()

    def reset(self):
        self.reads_all = 0
        self.reads_rev = 0
        self.reads_fwd = 0
        self.reads_pp = 0
        self.reads_pp_rev = 0
        self.reads_pp_fwd = 0
        self.reads_mate_unmapped = 0
        self.reads_mate_unmapped_rev = 0
        self.reads_mate_unmapped_fwd = 0
        self.reads_mate_other_chr = 0
        self.reads_mate_other_chr_rev = 0
        self.reads_mate_other_chr_fwd = 0
        self.reads_mate_same_strand = 0
        self.reads_mate_same_strand_rev = 0
        self.reads_mate_same_strand_fwd = 0
        self.reads_faceaway = 0
        self.reads_faceaway_rev = 0
        self.reads_faceaway_fwd = 0
        self.reads_softclipped = 0
        self.reads_softclipped_rev = 0
        self.reads_softclipped_fwd = 0
        self.reads_duplicate = 0
        self.reads_duplicate_rev = 0
        self.reads_duplicate_fwd = 0

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bint is_reverse
            bint is_duplicate
            bint mate_is_unmappped
            bint mate_is_reverse
            int tlen

        # convenience variables
        flag = read.b.core.flag
        is_reverse = <bint>(flag & BAM_FREVERSE)
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_duplicate = <bint>(flag & BAM_FDUP)
        mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
        mate_is_reverse = <bint>(flag & BAM_FMREVERSE)
        tlen = read.b.core.isize

        # do the counting
        self.reads_all += 1
        if is_reverse:
            self.reads_rev += 1
        else:
            self.reads_fwd += 1
        if is_proper_pair:
            self.reads_pp += 1
            if is_reverse:
                self.reads_pp_rev += 1
            else:
                self.reads_pp_fwd += 1
        if mate_is_unmapped:
            self.reads_mate_unmapped += 1
            if is_reverse:
                self.reads_mate_unmapped_rev += 1
            else:
                self.reads_mate_unmapped_fwd += 1
        elif col.tid != read.b.core.mtid:
            self.reads_mate_other_chr += 1
            if is_reverse:
                self.reads_mate_other_chr_rev += 1
            else:
                self.reads_mate_other_chr_fwd += 1
        elif is_reverse and mate_is_reverse:
            self.reads_mate_same_strand += 1
            self.reads_mate_same_strand_rev += 1
        elif not is_reverse and not mate_is_reverse:
            self.reads_mate_same_strand += 1
            self.reads_mate_same_strand_fwd += 1
        elif (is_reverse and tlen > 0) or (not is_reverse and tlen < 0):
            self.reads_faceaway += 1
            if is_reverse:
                self.reads_faceaway_rev += 1
            else:
                self.reads_faceaway_fwd += 1
        if is_softclipped(read.b):
            self.reads_softclipped += 1
            if is_reverse:
                self.reads_softclipped_rev += 1
            else:
                self.reads_softclipped_fwd += 1
        if is_duplicate:
            self.reads_duplicate += 1
            if is_reverse:
                self.reads_duplicate_rev += 1
            else:
                self.reads_duplicate_fwd += 1

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        rec = {'reads_all': self.reads_all,
               'reads_fwd': self.reads_fwd,
               'reads_rev': self.reads_rev,
               'reads_pp': self.reads_pp,
               'reads_pp_fwd': self.reads_pp_fwd,
               'reads_pp_rev': self.reads_pp_rev,
               'reads_mate_unmapped': self.reads_mate_unmapped,
               'reads_mate_unmapped_fwd': self.reads_mate_unmapped_fwd,
               'reads_mate_unmapped_rev': self.reads_mate_unmapped_rev,
               'reads_mate_other_chr': self.reads_mate_other_chr,
               'reads_mate_other_chr_fwd': self.reads_mate_other_chr_fwd,
               'reads_mate_other_chr_rev': self.reads_mate_other_chr_rev,
               'reads_mate_same_strand': self.reads_mate_same_strand,
               'reads_mate_same_strand_fwd': self.reads_mate_same_strand_fwd,
               'reads_mate_same_strand_rev': self.reads_mate_same_strand_rev,
               'reads_faceaway': self.reads_faceaway,
               'reads_faceaway_fwd': self.reads_faceaway_fwd,
               'reads_faceaway_rev': self.reads_faceaway_rev,
               'reads_softclipped': self.reads_softclipped,
               'reads_softclipped_fwd': self.reads_softclipped_fwd,
               'reads_softclipped_rev': self.reads_softclipped_rev,
               'reads_duplicate': self.reads_duplicate,
               'reads_duplicate_fwd': self.reads_duplicate_fwd,
               'reads_duplicate_rev': self.reads_duplicate_rev}

        # reset counters
        self.reset()

        return rec


########################
# VARIATION STATISTICS #
########################


# noinspection PyAttributeOutsideInit
cdef class Variation(PileupStat):

    cdef:
        int reads_all
        int reads_pp
        int matches
        int matches_pp
        int mismatches
        int mismatches_pp
        int deletions
        int deletions_pp
        int insertions
        int insertions_pp
        int A
        int A_pp
        int C
        int C_pp
        int T
        int T_pp
        int G
        int G_pp
        int N
        int N_pp

    def __init__(self):
        self.reset()

    def reset(self):
        self.reads_all = 0
        self.reads_pp = 0
        self.matches = 0
        self.matches_pp = 0
        self.mismatches = 0
        self.mismatches_pp = 0
        self.deletions = 0
        self.deletions_pp = 0
        self.insertions = 0
        self.insertions_pp = 0
        self.A = 0
        self.A_pp = 0
        self.C = 0
        self.C_pp = 0
        self.T = 0
        self.T_pp = 0
        self.G = 0
        self.G_pp = 0
        self.N = 0
        self.N_pp = 0

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bytes alnbase

        # convenience variables
        flag = read.b.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)

        # do the counting
        self.reads_all += 1
        if is_proper_pair:
            self.reads_pp += 1
        if read.is_del:
            self.deletions += 1
            if is_proper_pair:
                self.deletions_pp += 1
        else:
            alnbase = get_seq_base(read.b, read.qpos)
            if alnbase == b'A':
                self.A += 1
                if is_proper_pair:
                    self.A_pp += 1
            elif alnbase == b'T':
                self.T += 1
                if is_proper_pair:
                    self.T_pp += 1
            elif alnbase == b'C':
                self.C += 1
                if is_proper_pair:
                    self.C_pp += 1
            elif alnbase == b'G':
                self.G += 1
                if is_proper_pair:
                    self.G_pp += 1
            elif alnbase == b'N':
                self.N += 1
                if is_proper_pair:
                    self.N_pp += 1
            if read.indel > 0:
                self.insertions += 1
                if is_proper_pair:
                    self.insertions_pp += 1
            if alnbase == refbase:
                self.matches += 1
                if is_proper_pair:
                    self.matches_pp += 1
            else:
                self.mismatches += 1
                if is_proper_pair:
                    self.mismatches_pp += 1

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        if not PY2:
            ref = str(refbase, 'ascii')
        else:
            ref = refbase
        rec = {'ref': ref,
               'reads_all': self.reads_all,
               'reads_pp': self.reads_pp,
               'matches': self.matches,
               'matches_pp': self.matches_pp,
               'mismatches': self.mismatches,
               'mismatches_pp': self.mismatches_pp,
               'deletions': self.deletions,
               'deletions_pp': self.deletions_pp,
               'insertions': self.insertions,
               'insertions_pp': self.insertions_pp,
               'A': self.A, 'A_pp': self.A_pp,
               'C': self.C, 'C_pp': self.C_pp,
               'T': self.T, 'T_pp': self.T_pp,
               'G': self.G, 'G_pp': self.G_pp,
               'N': self.N, 'N_pp': self.N_pp}

        # reset counters
        self.reset()

        return rec


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


# noinspection PyAttributeOutsideInit
cdef class VariationStrand(PileupStat):

    cdef:
        CountPpStrand reads
        CountPpStrand matches
        CountPpStrand mismatches
        CountPpStrand deletions
        CountPpStrand insertions
        CountPpStrand A
        CountPpStrand C
        CountPpStrand T
        CountPpStrand G
        CountPpStrand N

    def __init__(self):
        self.reset()

    def reset(self):
        init_pp_strand(&self.reads)
        init_pp_strand(&self.matches)
        init_pp_strand(&self.mismatches)
        init_pp_strand(&self.deletions)
        init_pp_strand(&self.insertions)
        init_pp_strand(&self.A)
        init_pp_strand(&self.T)
        init_pp_strand(&self.C)
        init_pp_strand(&self.G)
        init_pp_strand(&self.N)

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bytes alnbase

        # convenience variables
        flag = read.b.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)

        # do the counting
        incr_pp_strand(&self.reads, is_reverse, is_proper_pair)
        if read.is_del:
            incr_pp_strand(&self.deletions, is_reverse, is_proper_pair)
        else:
            alnbase = get_seq_base(read.b, read.qpos)
            if alnbase == b'A':
                incr_pp_strand(&self.A, is_reverse, is_proper_pair)
            elif alnbase == b'T':
                incr_pp_strand(&self.T, is_reverse, is_proper_pair)
            elif alnbase == b'C':
                incr_pp_strand(&self.C, is_reverse, is_proper_pair)
            elif alnbase == b'G':
                incr_pp_strand(&self.G, is_reverse, is_proper_pair)
            elif alnbase == b'N':
                incr_pp_strand(&self.N, is_reverse, is_proper_pair)
            if read.indel > 0:
                incr_pp_strand(&self.insertions, is_reverse, is_proper_pair)
            if alnbase == refbase:
                incr_pp_strand(&self.matches, is_reverse, is_proper_pair)
            else:
                incr_pp_strand(&self.mismatches, is_reverse, is_proper_pair)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        if not PY2:
            ref = str(refbase, 'ascii')
        else:
            ref = refbase
        rec = {'ref': ref,
               'reads_all': self.reads.all,
               'reads_fwd': self.reads.fwd,
               'reads_rev': self.reads.rev,
               'reads_pp': self.reads.pp,
               'reads_pp_fwd': self.reads.pp_fwd,
               'reads_pp_rev': self.reads.pp_rev,
               'matches': self.matches.all,
               'matches_fwd': self.matches.fwd,
               'matches_rev': self.matches.rev,
               'matches_pp': self.matches.pp,
               'matches_pp_fwd': self.matches.pp_fwd,
               'matches_pp_rev': self.matches.pp_rev,
               'mismatches': self.mismatches.all,
               'mismatches_fwd': self.mismatches.fwd,
               'mismatches_rev': self.mismatches.rev,
               'mismatches_pp': self.mismatches.pp,
               'mismatches_pp_fwd': self.mismatches.pp_fwd,
               'mismatches_pp_rev': self.mismatches.pp_rev,
               'deletions': self.deletions.all,
               'deletions_fwd': self.deletions.fwd,
               'deletions_rev': self.deletions.rev,
               'deletions_pp': self.deletions.pp,
               'deletions_pp_fwd': self.deletions.pp_fwd,
               'deletions_pp_rev': self.deletions.pp_rev,
               'insertions': self.insertions.all,
               'insertions_fwd': self.insertions.fwd,
               'insertions_rev': self.insertions.rev,
               'insertions_pp': self.insertions.pp,
               'insertions_pp_fwd': self.insertions.pp_fwd,
               'insertions_pp_rev': self.insertions.pp_rev,
               'A': self.A.all, 'A_fwd': self.A.fwd, 'A_rev': self.A.rev,
               'A_pp': self.A.pp, 'A_pp_fwd': self.A.pp_fwd, 'A_pp_rev': self.A.pp_rev,
               'C': self.C.all, 'C_fwd': self.C.fwd, 'C_rev': self.C.rev,
               'C_pp': self.C.pp, 'C_pp_fwd': self.C.pp_fwd, 'C_pp_rev': self.C.pp_rev,
               'T': self.T.all, 'T_fwd': self.T.fwd, 'T_rev': self.T.rev,
               'T_pp': self.T.pp, 'T_pp_fwd': self.T.pp_fwd, 'T_pp_rev': self.T.pp_rev,
               'G': self.G.all, 'G_fwd': self.G.fwd, 'G_rev': self.G.rev,
               'G_pp': self.G.pp, 'G_pp_fwd': self.G.pp_fwd, 'G_pp_rev': self.G.pp_rev,
               'N': self.N.all, 'N_fwd': self.N.fwd, 'N_rev': self.N.rev,
               'N_pp': self.N.pp, 'N_pp_fwd': self.N.pp_fwd, 'N_pp_rev': self.N.pp_rev}

        # reset counters
        self.reset()

        return rec


##########################
# INSERT SIZE STATISTICS #
##########################


# noinspection PyAttributeOutsideInit
cdef class OnlineStats:

    cdef:
        long n
        double m
        double m2
        double d
        double d2
        int64_t s  # sum
        int64_t sq  # sum of squares

    def __init__(self):
        self.reset()

    def reset(self):
        self.n = 0
        self.m = 0
        self.m2 = 0
        self.d = 0
        self.d2 = 0
        self.s = 0
        self.sq = 0

    cdef void update(self, int64_t x):
        self.n += 1
        self.s += x
        self.sq += x ** 2
        self.d = x - self.m
        self.m += self.d / self.n
        self.d2 = x - self.m
        self.m2 += self.d * self.d2

    def variance(self):
        if self.n < 2:
            return 0
        else:
            return int(round(self.m2 / (self.n - 1)))

    def std(self):
        if self.n < 2:
            return 0
        else:
            return int(round(sqrt(self.m2 / (self.n - 1))))

    def mean(self):
        return mean(self.s, self.n)

    def rms(self):
        return rootmean(self.sq, self.n)


# noinspection PyAttributeOutsideInit
cdef class Tlen(PileupStat):

    cdef:
        int reads_all
        OnlineStats tlen
        OnlineStats tlen_pp

    def __init__(self):
        self.tlen = OnlineStats()
        self.tlen_pp = OnlineStats()
        self.reset()

    def reset(self):
        self.reads_all = 0
        self.tlen.reset()
        self.tlen_pp.reset()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bint mate_is_unmappped
            bint mate_other_chr
            int64_t tlen

        # convenience variables
        flag = read.b.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
        mate_other_chr = <bint>(read.b.core.tid != read.b.core.mtid)
        tlen = read.b.core.isize

        # do the counting
        self.reads_all += 1
        # N.B. insert size is only meaningful if mate is mapped to same chromosome
        if not mate_is_unmapped and not mate_other_chr:
            self.tlen.update(tlen)
        if is_proper_pair:
            self.tlen_pp.update(tlen)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        rec = {'reads_all': self.reads_all,
               'reads_paired': self.tlen.n,
               'reads_pp': self.tlen_pp.n,
               'mean_tlen': self.tlen.mean(),
               'mean_tlen_pp': self.tlen_pp.mean(),
               'rms_tlen': self.tlen.rms(),
               'rms_tlen_pp': self.tlen_pp.rms(),
               'std_tlen': self.tlen.std(),
               'std_tlen_pp': self.tlen_pp.std()}

        # reset counters
        self.reset()

        return rec


####################################
# INSERT SIZE STATISTICS BY STRAND #
####################################


# noinspection PyAttributeOutsideInit
cdef class TlenStrand(PileupStat):

    cdef:
        int reads_all
        int reads_fwd
        int reads_rev
        OnlineStats tlen
        OnlineStats tlen_fwd
        OnlineStats tlen_rev
        OnlineStats tlen_pp
        OnlineStats tlen_pp_fwd
        OnlineStats tlen_pp_rev

    def __init__(self):
        self.tlen = OnlineStats()
        self.tlen_fwd = OnlineStats()
        self.tlen_rev = OnlineStats()
        self.tlen_pp = OnlineStats()
        self.tlen_pp_fwd = OnlineStats()
        self.tlen_pp_rev = OnlineStats()
        self.reset()

    def reset(self):
        self.reads_all = 0
        self.reads_fwd = 0
        self.reads_rev = 0
        self.tlen.reset()
        self.tlen_fwd.reset()
        self.tlen_rev.reset()
        self.tlen_pp.reset()
        self.tlen_pp_fwd.reset()
        self.tlen_pp_rev.reset()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bint is_reverse
            bint mate_is_unmappped
            bint mate_other_chr
            int64_t tlen

        # convenience variables
        flag = read.b.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)
        mate_is_unmapped = <bint>(flag & BAM_FMUNMAP)
        mate_other_chr = <bint>(read.b.core.tid != read.b.core.mtid)
        tlen = read.b.core.isize

        # do the counting
        self.reads_all += 1
        if is_reverse:
            self.reads_rev += 1
        else:
            self.reads_fwd += 1

        # N.B. insert size is only meaningful if mate is mapped to same chromosome
        if not mate_is_unmapped and not mate_other_chr:
            self.tlen.update(tlen)
            if is_reverse:
                self.tlen_rev.update(tlen)
            else:
                self.tlen_fwd.update(tlen)
        if is_proper_pair:
            self.tlen_pp.update(tlen)
            if is_reverse:
                self.tlen_pp_rev.update(tlen)
            else:
                self.tlen_pp_fwd.update(tlen)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        rec = {
            'reads_all': self.reads_all,
            'reads_fwd': self.reads_fwd,
            'reads_rev': self.reads_rev,
            'reads_paired': self.tlen.n,
            'reads_paired_fwd': self.tlen_fwd.n,
            'reads_paired_rev': self.tlen_rev.n,
            'reads_pp': self.tlen_pp.n,
            'reads_pp_fwd': self.tlen_pp_fwd.n,
            'reads_pp_rev': self.tlen_pp_rev.n,
            'mean_tlen': self.tlen.mean(),
            'mean_tlen_fwd': self.tlen_fwd.mean(),
            'mean_tlen_rev': self.tlen_rev.mean(),
            'mean_tlen_pp': self.tlen_pp.mean(),
            'mean_tlen_pp_fwd': self.tlen_pp_fwd.mean(),
            'mean_tlen_pp_rev': self.tlen_pp_rev.mean(),
            'rms_tlen': self.tlen.rms(),
            'rms_tlen_fwd': self.tlen_fwd.rms(),
            'rms_tlen_rev': self.tlen_rev.rms(),
            'rms_tlen_pp': self.tlen_pp.rms(),
            'rms_tlen_pp_fwd': self.tlen_pp_fwd.rms(),
            'rms_tlen_pp_rev': self.tlen_pp_rev.rms(),
            'std_tlen': self.tlen.std(),
            'std_tlen_fwd': self.tlen_fwd.std(),
            'std_tlen_rev': self.tlen_rev.std(),
            'std_tlen_pp': self.tlen_pp.std(),
            'std_tlen_pp_fwd': self.tlen_pp_fwd.std(),
            'std_tlen_pp_rev': self.tlen_pp_rev.std(),
        }

        # reset counters
        self.reset()

        return rec


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

    # setup intermediate variables
    n = col.n
    plp = col.plp

    # setup reference base
    if fafile is not None:
        refbase = fafile.fetch(reference=chrom, start=col.pos, end=col.pos+1).upper()
        if not PY2:
            refbase = refbase.encode('ascii')
    else:
        refbase = None

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
        stat.recv(read, col, refbase)

    # construct record
    # TODO check this doesn't foul-up GC window centre
    pos = col.pos + 1 if one_based else col.pos
    rec = stat.rec(chrom, pos, fafile, refbase)
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


def rootmean(uint64_t sqsum, int count):
    if count > 0:
        return int(round(sqrt(sqsum / count)))
    else:
        return 0


def mean(int64_t total, int count):
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
