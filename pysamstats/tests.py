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
