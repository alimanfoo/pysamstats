# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import logging

from pysam import Samfile

from pysamstats.io import write_hdf5
import tables

logger = logging.getLogger(__name__)
debug = logger.debug


def check_write_hdf5_chrom_dtype(arg):

    # testing auto dtype determination.
    dtype, alignment, result, label = arg
    import tempfile

    # use auto
    with tempfile.NamedTemporaryFile(suffix=".h5") as tmp:

        write_hdf5("coverage", tmp.name, alignment, chrom=label, dtype=dtype)

        with tables.open_file(tmp.name, mode="r") as h5file:
            return result == h5file.root.data.dtype["chrom"].itemsize


def test_write_hdf5_chrom_dtype():

    contig_label = "AS2_scf7180000696055"
    bampath = "fixture/longcontignames.bam"

    dtypes = [None, {"chrom": "a20"}, {"chrom": "a20"}]
    alignments = [Samfile(bampath), Samfile(bampath), bampath]
    results = [len(contig_label), 20, 20]
    labels = [contig_label, contig_label, contig_label]

    for arg in zip(dtypes, alignments, results, labels):
        assert check_write_hdf5_chrom_dtype(arg)
