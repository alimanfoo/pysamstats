# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import logging

from pysam import Samfile
from pysamstats import load_coverage, stat_coverage

from pysamstats.io import write_hdf5
import tables

logger = logging.getLogger(__name__)
debug = logger.debug


def write_hdf5_dtype(arg):

    # testing auto dtype determination.
    dtype, alignment, result, label = arg
    import tempfile

    # use auto
    tmp = tempfile.mktemp(suffix=".h5")
    write_hdf5("coverage", tmp, alignment, chrom=label, dtype=dtype)
    with tables.open_file(tmp, mode="r") as h5file:
        assert result == h5file.root.data.dtype["chrom"].itemsize


def test_write_hdf5():

    contig_label = "AS2_scf7180000696055"
    bampath = "fixture/longcontignames.bam"

    dtypes = [None, {"chrom": "a20"}, {"chrom": "a20"}]
    alignments = [Samfile(bampath), Samfile(bampath), bampath]
    results = [len(contig_label), 20, 20]
    labels = [contig_label, contig_label, contig_label]

    for arg in zip(dtypes, alignments, results, labels):
        yield (write_hdf5_dtype, arg)
