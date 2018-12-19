# cython: profile=False
# cython: embedsignature=True
from __future__ import print_function, division, absolute_import


from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t
from libc.math cimport sqrt
from numpy import around as round
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
cdef class CountPp:
    cdef:
        int all, pp

    def __init__(self):
        self.reset()

    def reset(self):
        self.all = self.pp = 0

    cdef void incr(self, bint is_proper_pair):
        self.all += 1
        if is_proper_pair:
            self.pp += 1


# noinspection PyAttributeOutsideInit
cdef class Coverage(PileupStat):

    cdef:
        CountPp reads

    def __init__(self):
        self.reads = CountPp()
        self.reset()

    def reset(self):
        self.reads.reset()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            bint is_proper_pair

        # convenience variables
        is_proper_pair = <bint>(read.b.core.flag & BAM_FPROPER_PAIR)

        # do the counting
        self.reads.incr(is_proper_pair)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        rec = {'reads_all': self.reads.all,
               'reads_pp': self.reads.pp}

        # reset counters
        self.reset()

        return rec


################################
# STRANDED COVERAGE STATISTICS #
################################


# noinspection PyAttributeOutsideInit
cdef class CountPpStrand:
    cdef:
        int all, pp, fwd, rev, pp_fwd, pp_rev

    def __init__(self):
        self.reset()

    def reset(self):
        self.all = self.fwd = self.rev = self.pp = self.pp_fwd = self.pp_rev = 0

    cdef void incr(self, bint is_reverse, bint is_proper_pair):
        self.all += 1
        if is_reverse:
            self.rev += 1
            if is_proper_pair:
                self.pp += 1
                self.pp_rev += 1
        else:
            self.fwd += 1
            if is_proper_pair:
                self.pp += 1
                self.pp_fwd += 1


# noinspection PyAttributeOutsideInit
cdef class CoverageStrand(PileupStat):

    cdef:
        CountPpStrand reads

    def __init__(self):
        self.reads = CountPpStrand()
        self.reset()

    def reset(self):
        self.reads.reset()

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
        self.reads.incr(is_reverse, is_proper_pair)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        rec = {'reads_all': self.reads.all,
               'reads_fwd': self.reads.fwd,
               'reads_rev': self.reads.rev,
               'reads_pp': self.reads.pp,
               'reads_pp_fwd': self.reads.pp_fwd,
               'reads_pp_rev': self.reads.pp_rev}

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
cdef class CountStrand:
    cdef:
        int all, fwd, rev

    def __init__(self):
        self.reset()

    def reset(self):
        self.all = self.fwd = self.rev = 0

    cdef void incr(self, bint is_reverse):
        self.all += 1
        if is_reverse:
            self.rev += 1
        else:
            self.fwd += 1


# noinspection PyAttributeOutsideInit
cdef class CoverageExtStrand(PileupStat):

    cdef:
        CountStrand all
        CountStrand pp
        CountStrand mate_unmapped
        CountStrand mate_other_chr
        CountStrand same_strand
        CountStrand faceaway
        CountStrand softclipped
        CountStrand duplicate

    def __init__(self):
        self.all = CountStrand()
        self.pp = CountStrand()
        self.mate_unmapped = CountStrand()
        self.mate_other_chr = CountStrand()
        self.same_strand = CountStrand()
        self.faceaway = CountStrand()
        self.softclipped = CountStrand()
        self.duplicate = CountStrand()
        self.reset()

    def reset(self):
        self.all.reset()
        self.pp.reset()
        self.mate_unmapped.reset()
        self.mate_other_chr.reset()
        self.same_strand.reset()
        self.faceaway.reset()
        self.softclipped.reset()
        self.duplicate.reset()

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
        self.all.incr(is_reverse)
        if is_proper_pair:
            self.pp.incr(is_reverse)
        if mate_is_unmapped:
            self.mate_unmapped.incr(is_reverse)
        elif col.tid != read.b.core.mtid:
            self.mate_other_chr.incr(is_reverse)
        elif is_reverse and mate_is_reverse:
            self.same_strand.incr(is_reverse)
        elif not is_reverse and not mate_is_reverse:
            self.same_strand.incr(is_reverse)
        elif (is_reverse and tlen > 0) or (not is_reverse and tlen < 0):
            self.faceaway.incr(is_reverse)
        if is_softclipped(read.b):
            self.softclipped.incr(is_reverse)
        if is_duplicate:
            self.duplicate.incr(is_reverse)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        rec = {'reads_all': self.all.all,
               'reads_fwd': self.all.fwd,
               'reads_rev': self.all.rev,
               'reads_pp': self.pp.all,
               'reads_pp_fwd': self.pp.fwd,
               'reads_pp_rev': self.pp.rev,
               'reads_mate_unmapped': self.mate_unmapped.all,
               'reads_mate_unmapped_fwd': self.mate_unmapped.fwd,
               'reads_mate_unmapped_rev': self.mate_unmapped.rev,
               'reads_mate_other_chr': self.mate_other_chr.all,
               'reads_mate_other_chr_fwd': self.mate_other_chr.fwd,
               'reads_mate_other_chr_rev': self.mate_other_chr.rev,
               'reads_mate_same_strand': self.same_strand.all,
               'reads_mate_same_strand_fwd': self.same_strand.fwd,
               'reads_mate_same_strand_rev': self.same_strand.rev,
               'reads_faceaway': self.faceaway.all,
               'reads_faceaway_fwd': self.faceaway.fwd,
               'reads_faceaway_rev': self.faceaway.rev,
               'reads_softclipped': self.softclipped.all,
               'reads_softclipped_fwd': self.softclipped.fwd,
               'reads_softclipped_rev': self.softclipped.rev,
               'reads_duplicate': self.duplicate.all,
               'reads_duplicate_fwd': self.duplicate.fwd,
               'reads_duplicate_rev': self.duplicate.rev}

        # reset counters
        self.reset()

        return rec


########################
# VARIATION STATISTICS #
########################


# noinspection PyAttributeOutsideInit
cdef class Variation(PileupStat):

    cdef:
        CountPp reads
        CountPp matches
        CountPp mismatches
        CountPp deletions
        CountPp insertions
        CountPp A
        CountPp C
        CountPp T
        CountPp G
        CountPp N

    def __init__(self):
        self.reads = CountPp()
        self.matches = CountPp()
        self.mismatches = CountPp()
        self.deletions = CountPp()
        self.insertions = CountPp()
        self.A = CountPp()
        self.C = CountPp()
        self.T = CountPp()
        self.G = CountPp()
        self.N = CountPp()
        self.reset()

    def reset(self):
        self.reads.reset()
        self.matches.reset()
        self.mismatches.reset()
        self.deletions.reset()
        self.insertions.reset()
        self.A.reset()
        self.C.reset()
        self.T.reset()
        self.G.reset()
        self.N.reset()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bytes alnbase

        # convenience variables
        flag = read.b.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)

        # do the counting
        self.reads.incr(is_proper_pair)
        if read.is_del:
            if not read.is_refskip:
                self.deletions.incr(is_proper_pair)
        else:
            alnbase = get_seq_base(read.b, read.qpos)
            if alnbase == b'A':
                self.A.incr(is_proper_pair)
            elif alnbase == b'T':
                self.T.incr(is_proper_pair)
            elif alnbase == b'C':
                self.C.incr(is_proper_pair)
            elif alnbase == b'G':
                self.G.incr(is_proper_pair)
            elif alnbase == b'N':
                self.N.incr(is_proper_pair)
            if read.indel > 0:
                self.insertions.incr(is_proper_pair)
            if alnbase == refbase:
                self.matches.incr(is_proper_pair)
            else:
                self.mismatches.incr(is_proper_pair)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        if PY2:
            ref = refbase
        else:
            ref = str(refbase, 'ascii')
        rec = {'ref': ref,
               'reads_all': self.reads.all,
               'reads_pp': self.reads.pp,
               'matches': self.matches.all,
               'matches_pp': self.matches.pp,
               'mismatches': self.mismatches.all,
               'mismatches_pp': self.mismatches.pp,
               'deletions': self.deletions.all,
               'deletions_pp': self.deletions.pp,
               'insertions': self.insertions.all,
               'insertions_pp': self.insertions.pp,
               'A': self.A.all, 'A_pp': self.A.pp,
               'C': self.C.all, 'C_pp': self.C.pp,
               'T': self.T.all, 'T_pp': self.T.pp,
               'G': self.G.all, 'G_pp': self.G.pp,
               'N': self.N.all, 'N_pp': self.N.pp}

        # reset counters
        self.reset()

        return rec


#################################
# STRANDED VARIATION STATISTICS #
#################################


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
        self.reads = CountPpStrand()
        self.matches = CountPpStrand()
        self.mismatches = CountPpStrand()
        self.deletions = CountPpStrand()
        self.insertions = CountPpStrand()
        self.A = CountPpStrand()
        self.C = CountPpStrand()
        self.T = CountPpStrand()
        self.G = CountPpStrand()
        self.N = CountPpStrand()
        self.reset()

    def reset(self):
        self.reads.reset()
        self.matches.reset()
        self.mismatches.reset()
        self.deletions.reset()
        self.insertions.reset()
        self.A.reset()
        self.T.reset()
        self.C.reset()
        self.G.reset()
        self.N.reset()

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
        self.reads.incr(is_reverse, is_proper_pair)
        if read.is_del:
            if not read.is_refskip:
                self.deletions.incr(is_reverse, is_proper_pair)
        else:
            alnbase = get_seq_base(read.b, read.qpos)
            if alnbase == b'A':
                self.A.incr(is_reverse, is_proper_pair)
            elif alnbase == b'T':
                self.T.incr(is_reverse, is_proper_pair)
            elif alnbase == b'C':
                self.C.incr(is_reverse, is_proper_pair)
            elif alnbase == b'G':
                self.G.incr(is_reverse, is_proper_pair)
            elif alnbase == b'N':
                self.N.incr(is_reverse, is_proper_pair)
            if read.indel > 0:
                self.insertions.incr(is_reverse, is_proper_pair)
            if alnbase == refbase:
                self.matches.incr(is_reverse, is_proper_pair)
            else:
                self.mismatches.incr(is_reverse, is_proper_pair)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        if PY2:
            ref = refbase
        else:
            ref = str(refbase, 'ascii')
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
cdef class TlenHelper:

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
        TlenHelper tlen
        TlenHelper tlen_pp

    def __init__(self):
        self.tlen = TlenHelper()
        self.tlen_pp = TlenHelper()
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
        TlenHelper tlen
        TlenHelper tlen_fwd
        TlenHelper tlen_rev
        TlenHelper tlen_pp
        TlenHelper tlen_pp_fwd
        TlenHelper tlen_pp_rev

    def __init__(self):
        self.tlen = TlenHelper()
        self.tlen_fwd = TlenHelper()
        self.tlen_rev = TlenHelper()
        self.tlen_pp = TlenHelper()
        self.tlen_pp_fwd = TlenHelper()
        self.tlen_pp_rev = TlenHelper()
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


##############################
# MAPPING QUALITY STATISTICS #
##############################


# noinspection PyAttributeOutsideInit
cdef class MapqHelper:
    cdef:
        int n
        int nz
        uint64_t max
        uint64_t sqsum

    def __init__(self):
        self.reset()

    def reset(self):
        self.n = 0
        self.nz = 0
        self.max = 0
        self.sqsum = 0

    cdef void update(self, uint64_t mapq):
        self.n += 1
        if mapq == 0:
            self.nz += 1
        if mapq > self.max:
            self.max = mapq
        self.sqsum += mapq**2

    def rms(self):
        return rootmean(self.sqsum, self.n)


# noinspection PyAttributeOutsideInit
cdef class Mapq(PileupStat):

    cdef:
        MapqHelper all
        MapqHelper pp

    def __init__(self):
        self.all = MapqHelper()
        self.pp = MapqHelper()
        self.reset()

    def reset(self):
        self.all.reset()
        self.pp.reset()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            uint64_t mapq
            uint64_t mapq_squared

        # convenience variables
        flag = read.b.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)

        # do the counting
        mapq = read.b.core.qual
        self.all.update(mapq)
        if is_proper_pair:
            self.pp.update(mapq)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        rec = {'reads_all': self.all.n,
               'reads_pp': self.pp.n,
               'reads_mapq0': self.all.nz,
               'reads_mapq0_pp': self.pp.nz,
               'rms_mapq': self.all.rms(),
               'rms_mapq_pp': self.pp.rms(),
               'max_mapq': self.all.max,
               'max_mapq_pp': self.pp.max}

        # reset counters
        self.reset()

        return rec


########################################
# MAPPING QUALITY STATISTICS BY STRAND #
########################################


# noinspection PyAttributeOutsideInit
cdef class MapqStrand(PileupStat):

    cdef:
        MapqHelper all
        MapqHelper fwd
        MapqHelper rev
        MapqHelper pp
        MapqHelper pp_fwd
        MapqHelper pp_rev

    def __init__(self):
        self.all = MapqHelper()
        self.fwd = MapqHelper()
        self.rev = MapqHelper()
        self.pp = MapqHelper()
        self.pp_fwd = MapqHelper()
        self.pp_rev = MapqHelper()
        self.reset()

    def reset(self):
        self.all.reset()
        self.fwd.reset()
        self.rev.reset()
        self.pp.reset()
        self.pp_fwd.reset()
        self.pp_rev.reset()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bint is_reverse
            uint64_t mapq
            uint64_t mapq_squared

        # convenience variables
        flag = read.b.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)

        # do the counting
        mapq = read.b.core.qual
        self.all.update(mapq)
        if is_reverse:
            self.rev.update(mapq)
        else:
            self.fwd.update(mapq)
        if is_proper_pair:
            self.pp.update(mapq)
            if is_reverse:
                self.pp_rev.update(mapq)
            else:
                self.pp_fwd.update(mapq)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        rec = {
            'reads_all': self.all.n,
            'reads_fwd': self.fwd.n,
            'reads_rev': self.rev.n,
            'reads_pp': self.pp.n,
            'reads_pp_fwd': self.pp_fwd.n,
            'reads_pp_rev': self.pp_rev.n,
            'reads_mapq0': self.all.nz,
            'reads_mapq0_fwd': self.fwd.nz,
            'reads_mapq0_rev': self.rev.nz,
            'reads_mapq0_pp': self.pp.nz,
            'reads_mapq0_pp_fwd': self.pp_fwd.nz,
            'reads_mapq0_pp_rev': self.pp_rev.nz,
            'rms_mapq': self.all.rms(),
            'rms_mapq_fwd': self.fwd.rms(),
            'rms_mapq_rev': self.rev.rms(),
            'rms_mapq_pp': self.pp.rms(),
            'rms_mapq_pp_fwd': self.pp_fwd.rms(),
            'rms_mapq_pp_rev': self.pp_rev.rms(),
            'max_mapq': self.all.max,
            'max_mapq_fwd': self.fwd.max,
            'max_mapq_rev': self.rev.max,
            'max_mapq_pp': self.pp.max,
            'max_mapq_pp_fwd': self.pp_fwd.max,
            'max_mapq_pp_rev': self.pp_rev.max,
        }

        # reset counters
        self.reset()

        return rec


###########################
# BASE QUALITY STATISTICS #
###########################


# noinspection PyAttributeOutsideInit
cdef class BaseqHelper:

    cdef:
        int n
        int n_nodel
        uint64_t sqsum

    def __init__(self):
        self.reset()

    def reset(self):
        self.n = 0
        self.n_nodel = 0
        self.sqsum = 0

    cdef void update(self, int64_t baseq_squared):
        self.n += 1
        if baseq_squared >= 0:
            self.n_nodel += 1
            self.sqsum += baseq_squared

    def rms(self):
        return rootmean(self.sqsum, self.n_nodel)


# noinspection PyAttributeOutsideInit
cdef class BaseqPpHelper:

    cdef:
        BaseqHelper all
        BaseqHelper pp

    def __init__(self):
        self.all = BaseqHelper()
        self.pp = BaseqHelper()
        self.reset()

    def reset(self):
        self.all.reset()
        self.pp.reset()

    cdef void update(self, int64_t baseq_squared, bint is_proper_pair):
        self.all.update(baseq_squared)
        if is_proper_pair:
            self.pp.update(baseq_squared)


# noinspection PyAttributeOutsideInit
cdef class Baseq(PileupStat):

    cdef:
        BaseqPpHelper helper

    def __init__(self):
        self.helper = BaseqPpHelper()
        self.reset()

    def reset(self):
        self.helper.reset()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            int64_t baseq = -1
            int64_t baseq_squared = -1

        # convenience variables
        flag = read.b.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)

        # do the counting
        if not read.is_del:
            baseq = pysam_bam_get_qual(read.b)[read.qpos]
            baseq_squared = baseq**2

        self.helper.update(baseq_squared, is_proper_pair)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        rec = {
            'reads_all': self.helper.all.n,
            'reads_pp': self.helper.pp.n,
            'rms_baseq': self.helper.all.rms(),
            'rms_baseq_pp': self.helper.pp.rms(),
        }

        # reset counters
        self.reset()

        return rec


# noinspection PyAttributeOutsideInit
cdef class BaseqStrandPpHelper:

    cdef:
        BaseqHelper all
        BaseqHelper fwd
        BaseqHelper rev
        BaseqHelper pp
        BaseqHelper pp_fwd
        BaseqHelper pp_rev

    def __init__(self):
        self.all = BaseqHelper()
        self.fwd = BaseqHelper()
        self.rev = BaseqHelper()
        self.pp = BaseqHelper()
        self.pp_fwd = BaseqHelper()
        self.pp_rev = BaseqHelper()
        self.reset()

    def reset(self):
        self.all.reset()
        self.fwd.reset()
        self.rev.reset()
        self.pp.reset()
        self.pp_fwd.reset()
        self.pp_rev.reset()

    cdef void update(self, int64_t baseq_squared, bint is_proper_pair, bint is_reverse):
        self.all.update(baseq_squared)
        if is_reverse:
            self.rev.update(baseq_squared)
        else:
            self.fwd.update(baseq_squared)
        if is_proper_pair:
            self.pp.update(baseq_squared)
            if is_reverse:
                self.pp_rev.update(baseq_squared)
            else:
                self.pp_fwd.update(baseq_squared)


# noinspection PyAttributeOutsideInit
cdef class BaseqStrand(PileupStat):

    cdef:
        BaseqStrandPpHelper helper

    def __init__(self):
        self.helper = BaseqStrandPpHelper()
        self.reset()

    def reset(self):
        self.helper.reset()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bint is_reverse
            int64_t baseq = -1
            int64_t baseq_squared = -1

        # convenience variables
        flag = read.b.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)

        # do the counting
        if not read.is_del:
            baseq = pysam_bam_get_qual(read.b)[read.qpos]
            baseq_squared = baseq**2

        self.helper.update(baseq_squared, is_proper_pair, is_reverse)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        rec = {
            'reads_all': self.helper.all.n,
            'reads_fwd': self.helper.fwd.n,
            'reads_rev': self.helper.rev.n,
            'reads_pp': self.helper.pp.n,
            'reads_pp_fwd': self.helper.pp_fwd.n,
            'reads_pp_rev': self.helper.pp_rev.n,
            'rms_baseq': self.helper.all.rms(),
            'rms_baseq_fwd': self.helper.fwd.rms(),
            'rms_baseq_rev': self.helper.rev.rms(),
            'rms_baseq_pp': self.helper.pp.rms(),
            'rms_baseq_pp_fwd': self.helper.pp_fwd.rms(),
            'rms_baseq_pp_rev': self.helper.pp_rev.rms(),
        }

        # reset counters
        self.reset()

        return rec


####################################
# EXTENDED BASE QUALITY STATISTICS #
####################################


# noinspection PyAttributeOutsideInit
cdef class BaseqExt(PileupStat):

    cdef:
        BaseqPpHelper all
        BaseqPpHelper matches
        BaseqPpHelper mismatches

    def __init__(self):
        self.all = BaseqPpHelper()
        self.matches = BaseqPpHelper()
        self.mismatches = BaseqPpHelper()
        self.reset()

    def reset(self):
        self.all.reset()
        self.matches.reset()
        self.mismatches.reset()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bint is_reverse
            int64_t baseq = -1
            int64_t baseq_squared = -1
            bytes alnbase

        # convenience variables
        flag = read.b.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)

        # do the counting
        if not read.is_del:
            baseq = pysam_bam_get_qual(read.b)[read.qpos]
            baseq_squared = baseq**2

        self.all.update(baseq_squared, is_proper_pair)
        if not read.is_del:
            alnbase = get_seq_base(read.b, read.qpos)
            if alnbase == refbase:
                self.matches.update(baseq_squared, is_proper_pair)
            else:
                self.mismatches.update(baseq_squared, is_proper_pair)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        if PY2:
            ref = refbase
        else:
            ref = str(refbase, 'ascii')
        rec = {
            'ref': ref,
            'reads_all': self.all.all.n,
            'reads_pp': self.all.pp.n,
            'matches': self.matches.all.n,
            'matches_pp': self.matches.pp.n,
            'mismatches': self.mismatches.all.n,
            'mismatches_pp': self.mismatches.pp.n,
            'rms_baseq': self.all.all.rms(),
            'rms_baseq_pp': self.all.pp.rms(),
            'rms_baseq_matches': self.matches.all.rms(),
            'rms_baseq_matches_pp': self.matches.pp.rms(),
            'rms_baseq_mismatches': self.mismatches.all.rms(),
            'rms_baseq_mismatches_pp': self.mismatches.pp.rms(),
        }

        # reset counters
        self.reset()

        return rec


# noinspection PyAttributeOutsideInit
cdef class BaseqExtStrand(PileupStat):

    cdef:
        BaseqStrandPpHelper all
        BaseqStrandPpHelper matches
        BaseqStrandPpHelper mismatches

    def __init__(self):
        self.all = BaseqStrandPpHelper()
        self.matches = BaseqStrandPpHelper()
        self.mismatches = BaseqStrandPpHelper()
        self.reset()

    def reset(self):
        self.all.reset()
        self.matches.reset()
        self.mismatches.reset()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            uint32_t flag
            bint is_proper_pair
            bint is_reverse
            int64_t baseq = -1
            int64_t baseq_squared = -1
            bytes alnbase

        # convenience variables
        flag = read.b.core.flag
        is_proper_pair = <bint>(flag & BAM_FPROPER_PAIR)
        is_reverse = <bint>(flag & BAM_FREVERSE)

        # do the counting
        if not read.is_del:
            baseq = pysam_bam_get_qual(read.b)[read.qpos]
            baseq_squared = baseq**2

        self.all.update(baseq_squared, is_proper_pair, is_reverse)
        if not read.is_del:
            alnbase = get_seq_base(read.b, read.qpos)
            if alnbase == refbase:
                self.matches.update(baseq_squared, is_proper_pair, is_reverse)
            else:
                self.mismatches.update(baseq_squared, is_proper_pair, is_reverse)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # make record
        if PY2:
            ref = refbase
        else:
            ref = str(refbase, 'ascii')
        rec = {
            'ref': ref,
            'reads_all': self.all.all.n,
            'reads_fwd': self.all.fwd.n,
            'reads_rev': self.all.rev.n,
            'reads_pp': self.all.pp.n,
            'reads_pp_fwd': self.all.pp_fwd.n,
            'reads_pp_rev': self.all.pp_rev.n,
            'matches': self.matches.all.n,
            'matches_fwd': self.matches.fwd.n,
            'matches_rev': self.matches.rev.n,
            'matches_pp': self.matches.pp.n,
            'matches_pp_fwd': self.matches.pp_fwd.n,
            'matches_pp_rev': self.matches.pp_rev.n,
            'mismatches': self.mismatches.all.n,
            'mismatches_fwd': self.mismatches.fwd.n,
            'mismatches_rev': self.mismatches.rev.n,
            'mismatches_pp': self.mismatches.pp.n,
            'mismatches_pp_fwd': self.mismatches.pp_fwd.n,
            'mismatches_pp_rev': self.mismatches.pp_rev.n,
            'rms_baseq': self.all.all.rms(),
            'rms_baseq_fwd': self.all.fwd.rms(),
            'rms_baseq_rev': self.all.rev.rms(),
            'rms_baseq_pp': self.all.pp.rms(),
            'rms_baseq_pp_fwd': self.all.pp_fwd.rms(),
            'rms_baseq_pp_rev': self.all.pp_rev.rms(),
            'rms_baseq_matches': self.matches.all.rms(),
            'rms_baseq_matches_fwd': self.matches.fwd.rms(),
            'rms_baseq_matches_rev': self.matches.rev.rms(),
            'rms_baseq_matches_pp': self.matches.pp.rms(),
            'rms_baseq_matches_pp_fwd': self.matches.pp_fwd.rms(),
            'rms_baseq_matches_pp_rev': self.matches.pp_rev.rms(),
            'rms_baseq_mismatches': self.mismatches.all.rms(),
            'rms_baseq_mismatches_fwd': self.mismatches.fwd.rms(),
            'rms_baseq_mismatches_rev': self.mismatches.rev.rms(),
            'rms_baseq_mismatches_pp': self.mismatches.pp.rms(),
            'rms_baseq_mismatches_pp_fwd': self.mismatches.pp_fwd.rms(),
            'rms_baseq_mismatches_pp_rev': self.mismatches.pp_rev.rms(),
        }

        # reset counters
        self.reset()

        return rec


#################################################
# BASIC COVERAGE STATISTICS WITH GC COMPOSITION #
#################################################


# noinspection PyAttributeOutsideInit
cdef class CoverageGC(PileupStat):

    cdef:
        CountPp reads
        int window_size, window_offset

    def __init__(self, window_size, window_offset):
        self.window_size = window_size
        if window_offset is None:
            window_offset = window_size // 2
        self.window_offset = window_offset
        self.reads = CountPp()
        self.reset()

    def reset(self):
        self.reads.reset()

    cdef void recv(self, bam_pileup1_t* read, PileupColumn col, bytes refbase):
        cdef:
            bint is_proper_pair

        # convenience variables
        is_proper_pair = <bint>(read.b.core.flag & BAM_FPROPER_PAIR)

        # do the counting
        self.reads.incr(is_proper_pair)

    cdef dict rec(self, chrom, pos, FastaFile fafile, bytes refbase):

        # compute GC
        ref_window_start = pos - self.window_offset
        ref_window_end = ref_window_start + self.window_size
        if ref_window_start < 0:
            ref_window_start = 0
        ref_window = fafile.fetch(chrom, ref_window_start, ref_window_end).lower()
        if len(ref_window) == 0:
            gc = -1
        else:
            gc = gc_content(ref_window)

        # make record
        rec = {'reads_all': self.reads.all,
               'reads_pp': self.reads.pp,
               'gc': gc}

        # reset counters
        self.reset()

        return rec


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


def iter_pileup_default(stat, alignmentfile, fafile, chrom, start, end, one_based, truncate, stepper,
                        max_depth, int min_mapq, int min_baseq, bint no_del, bint no_dup):
    cdef:
        PileupColumn col

    # obtain pileup iterator
    start, end = normalise_coords(alignmentfile, chrom, start, end, one_based)
    it = alignmentfile.pileup(reference=chrom, start=start, end=end, truncate=truncate, stepper=stepper,
                              max_depth=max_depth)

    # iterate over pileup columns
    for col in it:

        rec = stat_pileup(stat, col, alignmentfile=alignmentfile, fafile=fafile,
                          one_based=one_based, min_mapq=min_mapq, min_baseq=min_baseq,
                          no_del=no_del, no_dup=no_dup)

        yield rec


def stat_pileup(PileupStat stat,
                PileupColumn col,
                AlignmentFile alignmentfile,
                FastaFile fafile,
                bint one_based,
                int min_mapq,
                int min_baseq,
                bint no_del,
                bint no_dup):
    cdef:
        bam_pileup1_t** plp
        bam_pileup1_t* read
        bam1_t* b
        int i, n

    # setup intermediate variables
    n = col.n
    plp = col.plp

    # get current chromosome
    chrom = alignmentfile.getrname(col.tid)

    # setup reference base
    refbase = get_refbase(fafile, chrom, col.pos)

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
    rec = stat.rec(chrom, col.pos, fafile, refbase)
    pos = col.pos + 1 if one_based else col.pos
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


def iter_pileup_padded_chrom(PileupStat stat, alignmentfile, fafile, chrom, start, end,
                             one_based, truncate, stepper, max_depth, min_mapq, min_baseq, no_del, no_dup):
    cdef:
        PileupColumn col
        int curpos

    # obtain pileup iterator
    assert chrom is not None, 'chromosome is None'
    start, end = normalise_coords(alignmentfile, chrom, start, end, one_based)
    it = alignmentfile.pileup(reference=chrom, start=start, end=end, truncate=truncate, stepper=stepper,
                              max_depth=max_depth)

    # keep track of current position
    curpos = start

    # iterate over pileup columns
    for col in it:

        # pad up to next pileup column
        while curpos < col.pos:
            # setup reference base
            refbase = get_refbase(fafile, chrom, curpos)
            # construct record
            rec = stat.rec(chrom, curpos, fafile, refbase)
            rec['chrom'] = chrom
            pos = curpos + 1 if one_based else curpos
            rec['pos'] = pos
            yield rec
            curpos += 1

        # construct and yield statistics record for current pileup column
        rec = stat_pileup(stat, col, alignmentfile=alignmentfile, fafile=fafile,
                          one_based=one_based, min_mapq=min_mapq, min_baseq=min_baseq,
                          no_del=no_del, no_dup=no_dup)
        yield rec

        # advance current position
        curpos = col.pos + 1

    while curpos < end:

        # pad out to end position
        refbase = get_refbase(fafile, chrom, curpos)
        rec = stat.rec(chrom, curpos, fafile, refbase)
        rec['chrom'] = chrom
        pos = curpos + 1 if one_based else curpos
        rec['pos'] = pos
        yield rec
        curpos += 1


def get_refbase(fafile, chrom, pos):
    if fafile is not None:
        refbase = fafile.fetch(reference=chrom, start=pos, end=pos+1).upper()
        if not PY2:
            refbase = refbase.encode('ascii')
    else:
        refbase = None
    return refbase


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
