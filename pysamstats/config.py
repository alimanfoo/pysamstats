# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division


stats_types_noref = ('coverage',
                     'coverage_strand',
                     'coverage_ext',
                     'coverage_ext_strand',
                     'tlen',
                     'tlen_strand',
                     'mapq',
                     'mapq_strand',
                     'baseq',
                     'baseq_strand',
                     'mapq_binned',
                     'alignment_binned',
                     'tlen_binned')

stats_types_withref = ('variation',
                       'variation_strand',
                       'baseq_ext',
                       'baseq_ext_strand',
                       'coverage_gc',
                       'coverage_binned',
                       'coverage_ext_binned')

stats_types = sorted(stats_types_noref + stats_types_withref)

stepper_types = ('nofilter',
                 'samtools',
                 'all')

dtype_coverage = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_pp', 'i4')
]

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

dtype_baseq = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_pp', 'i4'),
    ('rms_baseq', 'i4'),
    ('rms_baseq_pp', 'i4'),
]

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

dtype_coverage_gc = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('gc', 'u1'),
    ('reads_all', 'i4'),
    ('reads_pp', 'i4')
]

dtype_coverage_binned = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('gc', 'u1'),
    ('reads_all', 'i4'),
    ('reads_pp', 'i4')
]

dtype_coverage_ext_binned = [
    ('chrom', 'a12'),
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

dtype_mapq_binned = [
    ('chrom', 'a12'),
    ('pos', 'i4'),
    ('reads_all', 'i4'),
    ('reads_mapq0', 'i4'),
    ('rms_mapq', 'i4'),
]

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
