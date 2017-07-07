from __future__ import print_function, division, absolute_import


"""
Tests for the pysamstats module.

The strategy here is to compare the outputs of the functions under test with
unoptimised, pure-python reference implementations of the same functions, over
an example dataset.

"""


# TODO break up tests into modules
# TODO setup CI


import sys
from pysam import Samfile, Fastafile
import numpy as np
import logging
import pysamstats


# # pileup_functions = [
# #     (pysamstats.load_coverage, 0),
# #     (pysamstats.load_coverage_strand, 0),
# #     (pysamstats.load_coverage_ext, 0),
# #     (pysamstats.load_coverage_ext_strand, 0),
# #     (pysamstats.load_variation, 1),
# #     (pysamstats.load_variation_strand, 1),
# #     (pysamstats.load_tlen, 0),
# #     (pysamstats.load_tlen_strand, 0),
# #     (pysamstats.load_mapq, 0),
# #     (pysamstats.load_mapq_strand, 0),
# #     (pysamstats.load_baseq, 0),
# #     (pysamstats.load_baseq_strand, 0),
# #     (pysamstats.load_baseq_ext, 1),
# #     (pysamstats.load_baseq_ext_strand, 1),
# #     (pysamstats.load_coverage_gc, 1),
# # ]
#
#
# def test_pileup_truncate():
#     kwargs_notrunc = {'chrom': 'Pf3D7_01_v3',
#                       'start': 2000,
#                       'end': 2100,
#                       'one_based': False,
#                       'truncate': False}
#     kwargs_trunc = {'chrom': 'Pf3D7_01_v3',
#                     'start': 2000,
#                     'end': 2100,
#                     'one_based': False,
#                     'truncate': True}
#     for f, needs_ref in pileup_functions:
#         debug(f.__name__)
#         # test no truncate
#         if needs_ref:
#             a = f(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'),
#                   **kwargs_notrunc)
#         else:
#             a = f(Samfile('fixture/test.bam'), **kwargs_notrunc)
#         debug(a[:5])
#         eq_(1952, a['pos'][0])
#         eq_(2154, a['pos'][-1])
#         # test truncate
#         if needs_ref:
#             a = f(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'),
#                   **kwargs_trunc)
#         else:
#             a = f(Samfile('fixture/test.bam'), **kwargs_trunc)
#         eq_(2000, a['pos'][0])
#         eq_(2099, a['pos'][-1])
#
#
# def test_pileup_pad():
#     kwargs_nopad = {'chrom': 'Pf3D7_01_v3',
#                     'start': 0,
#                     'end': 20000,
#                     'one_based': False,
#                     'pad': False}
#     kwargs_pad = {'chrom': 'Pf3D7_01_v3',
#                   'start': 0,
#                   'end': 20000,
#                   'one_based': False,
#                   'pad': True}
#     for f, needs_ref in pileup_functions:
#         debug(f.__name__)
#         # test no pad
#         if needs_ref:
#             a = f(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'),
#                   **kwargs_nopad)
#         else:
#             a = f(Samfile('fixture/test.bam'), **kwargs_nopad)
#         eq_(924, a['pos'][0])
#         eq_(9935, a['pos'][-1])
#         # test pad
#         if needs_ref:
#             a = f(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'),
#                   **kwargs_pad)
#         else:
#             a = f(Samfile('fixture/test.bam'), **kwargs_pad)
#         eq_(0, a['pos'][0])
#         eq_(19999, a['pos'][-1])
#
#
# def test_pileup_pad_wg():
#     # whole genome
#     expected = stat_coverage_refimpl(Samfile('fixture/test.bam'))
#     actual = pysamstats.stat_coverage(Samfile('fixture/test.bam'))
#     _compare_iterators(expected, actual)
#     kwargs_nopad = {'pad': False}
#     kwargs_pad = {'pad': True}
#     for f, needs_ref in pileup_functions:
#         debug(f.__name__)
#         # test no pad
#         if needs_ref:
#             a = f(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'),
#                   **kwargs_nopad)
#         else:
#             a = f(Samfile('fixture/test.bam'), **kwargs_nopad)
#         eq_(sorted(set(a['chrom'])), [b'Pf3D7_01_v3', b'Pf3D7_02_v3'])
#         eq_(924, a[a['chrom'] == b'Pf3D7_01_v3']['pos'][0])
#         eq_(9935, a[a['chrom'] == b'Pf3D7_01_v3']['pos'][-1])
#         eq_(926, a[a['chrom'] == b'Pf3D7_02_v3']['pos'][0])
#         eq_(10074, a[a['chrom'] == b'Pf3D7_02_v3']['pos'][-1])
#         # test pad
#         if needs_ref:
#             a = f(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'),
#                   **kwargs_pad)
#         else:
#             a = f(Samfile('fixture/test.bam'), **kwargs_pad)
#         eq_(sorted(set(a['chrom'])),
#             [b'Pf3D7_01_v3', b'Pf3D7_02_v3', b'Pf3D7_03_v3'])
#         eq_(0, a[a['chrom'] == b'Pf3D7_01_v3']['pos'][0])
#         eq_(50000, a[a['chrom'] == b'Pf3D7_01_v3']['pos'][-1])
#         eq_(0, a[a['chrom'] == b'Pf3D7_02_v3']['pos'][0])
#         eq_(60000, a[a['chrom'] == b'Pf3D7_02_v3']['pos'][-1])
#         eq_(0, a[a['chrom'] == b'Pf3D7_03_v3']['pos'][0])
#         eq_(70000, a[a['chrom'] == b'Pf3D7_03_v3']['pos'][-1])
#
#
# # binned_functions = [
# #     (pysamstats.load_coverage_binned, 1),
# #     (pysamstats.load_coverage_ext_binned, 1),
# #     (pysamstats.load_mapq_binned, 0),
# #     (pysamstats.load_alignment_binned, 0),
# #     (pysamstats.load_tlen_binned, 0),
# # ]
#
#
# def test_binned_pad_region():
#     kwargs = {'chrom': 'Pf3D7_01_v3',
#               'start': 1000,
#               'end': 20000,
#               'one_based': False,
#               'window_size': 200,
#               'window_offset': 100}
#     for f, needs_ref in binned_functions:
#         debug(f.__name__)
#         if needs_ref:
#             a = f(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'),
#                   **kwargs)
#         else:
#             a = f(Samfile('fixture/test.bam'), **kwargs)
#         assert set(a['chrom']) == {b'Pf3D7_01_v3'}
#         eq_(1100, a['pos'][0])
#         eq_(19900, a['pos'][-1])
#
#
# def test_binned_pad_wg():
#     expected = stat_coverage_binned_refimpl(
#         Samfile('fixture/test.bam'),
#         Fastafile('fixture/ref.fa'))
#
#     actual = pysamstats.stat_coverage_binned(Samfile('fixture/test.bam'),
#                                              Fastafile('fixture/ref.fa'))
#     _compare_iterators(expected, actual)
#     kwargs = {'window_size': 200,
#               'window_offset': 100}
#     for f, needs_ref in binned_functions:
#         debug(f.__name__)
#         if needs_ref:
#             a = f(Samfile('fixture/test.bam'), Fastafile('fixture/ref.fa'),
#                   **kwargs)
#         else:
#             a = f(Samfile('fixture/test.bam'), **kwargs)
#         assert sorted(set(a['chrom'])) == [b'Pf3D7_01_v3', b'Pf3D7_02_v3',
#                                            b'Pf3D7_03_v3']
#         eq_(100, a[a['chrom'] == b'Pf3D7_01_v3']['pos'][0])
#         eq_(50100, a[a['chrom'] == b'Pf3D7_01_v3']['pos'][-1])
#         eq_(100, a[a['chrom'] == b'Pf3D7_02_v3']['pos'][0])
#         eq_(60100, a[a['chrom'] == b'Pf3D7_02_v3']['pos'][-1])
#         eq_(100, a[a['chrom'] == b'Pf3D7_03_v3']['pos'][0])
#         eq_(70100, a[a['chrom'] == b'Pf3D7_03_v3']['pos'][-1])
#
#
# def test_pileup_limit():
#
#     for f, needs_ref in pileup_functions:
#         debug(f.__name__)
#
#         # test with effectively no limit
#         kwargs = dict(fields=['reads_all'], max_depth=1000000)
#         if needs_ref:
#             a = f(Samfile('fixture/deep.bam'), Fastafile('fixture/ref.fa'),
#                   **kwargs)
#         else:
#             a = f(Samfile('fixture/deep.bam'), **kwargs)
#         eq_(26169, a[70])
#
#         # test with specific limit
#         kwargs = dict(fields=['reads_all'], max_depth=12000)
#         if needs_ref:
#             a = f(Samfile('fixture/deep.bam'), Fastafile('fixture/ref.fa'),
#                   **kwargs)
#         else:
#             a = f(Samfile('fixture/deep.bam'), **kwargs)
#         eq_(12046, a[70])  # no idea why limit is not exact
#
#         # test with default limit
#         kwargs = dict(fields=['reads_all'])
#         if needs_ref:
#             a = f(Samfile('fixture/deep.bam'), Fastafile('fixture/ref.fa'),
#                   **kwargs)
#         else:
#             a = f(Samfile('fixture/deep.bam'), **kwargs)
#         eq_(8052, a[70])  # no idea why limit is not exact
