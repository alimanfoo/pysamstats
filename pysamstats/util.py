# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, division
import sys
import csv
import time
import itertools
from operator import itemgetter


def write_csv(stats_type, outfile, alignmentfile, fields=None, dialect='excel-tab',
              write_header=True, progress=None, **kwargs):
    """Write statistics output to a CSV file.

    Parameters
    ----------

    stats_type : string
        Statistics type, one of 'coverage', 'coverage_ext', etc.
    outfile : file-like
        Output file to write to.
    alignmentfile : pysam.AlignmentFile or string
        Input BAM or SAM file or file path.
    fields : list of strings
        List of field names to output (all by default).
    dialect : string
        CSV dialect.
    write_header : bool
        If True write a header row.
    progress : int
        Log progress to stderr every N rows.
    **kwargs
        Passed through to the statistics function.

    """

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
    writer = csv.writer(outfile, dialect=dialect)

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
        before = time.time()
        before_all = before
        for row in rows:
            counter += 1
            writer.writerow(row)
            if counter % modulus == 0:
                after = time.time()
                elapsed = after - before_all
                batch_elapsed = after - before
                msg = '[pysamstats] %s rows in %.2fs (%d rows/s); batch in ' \
                      '%.2fs (%d rows/s)' \
                      % (counter, elapsed, counter / elapsed, batch_elapsed,
                         progress / batch_elapsed)
                print(msg, file=sys.stderr)
                before = after
        after_all = time.time()
        elapsed_all = after_all - before_all
        msg = '[pysamstats] %s rows in %.2fs (%d rows/s)' \
              % (counter, elapsed_all, counter / elapsed_all)
        print(msg, file=sys.stderr)


def write_hdf5(stats_type, outfile, alignmentfile, fields=None, progress=None, hdf5_group='/',
               hdf5_dataset='data', hdf5_complevel=1, hdf5_complib='zlib', hdf5_shuffle=True,
               hdf5_fletcher32=False, hdf5_chunksize=2**20, dtype=None, **kwargs):
    """Write statistics output to an HDF5 file. Requires PyTables.

    Parameters
    ----------
    stats_type : string
        Statistics type, one of 'coverage', 'coverage_ext', etc.
    outfile : string
        Output file path.
    alignmentfile : pysam.AlignmentFile or string
        Input BAM or SAM file or file path.
    fields : list of strings
        List of field names to output (all by default).
    progress : int
        Log progress to stderr approximately every N rows.
    hdf5_group : string
        Group to write new dataset to.
    hdf5_dataset : string
        Name of dataset to create.
    hdf5_complib : string
        Name of compression library (defaults to 'zlib').
    hdf5_complevel : int
        Compression level.
    hdf5_chunksize : int
        Size of chunks in number of bytes.
    hdf5_shuffle : bool
        If True, use byte shuffle filter.
    hdf5_fletcher32 : bool
        If True, use fletcher 32 filter.
    dtype : dict
        Override dtype.
    **kwargs
        Passed through to the statistics function.

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
    default_dtype = dict(globals()['dtype_' + stats_type])
    if dtype is None:
        dtype = default_dtype
    else:
        default_dtype.update(dtype)
        dtype = default_dtype
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
        before = time.time()
        before_all = before

        # load data in batches of size `hdf5_chunklen`
        chunk = list(itertools.islice(rows, hdf5_chunklen))

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
                after = time.time()
                elapsed = after - before_all
                batch_elapsed = after - before
                batch_size = counter - counter_before
                msg = '[pysamstats] %s rows in %.2fs (%d rows/s); last %s ' \
                      'rows in %.2fs (%d rows/s)' \
                      % (counter, elapsed, counter / elapsed,
                         batch_size, batch_elapsed, batch_size / batch_elapsed)
                print(msg, file=sys.stderr)
                before = after
                counter_before = counter

            # load next batch
            chunk = list(itertools.islice(rows, hdf5_chunklen))

        if progress:
            after_all = time.time()
            elapsed_all = after_all - before_all
            msg = '[pysamstats] %s rows in %.2fs (%d rows/s)' \
                  % (counter, elapsed_all, counter / elapsed_all)
            print(msg, file=sys.stderr)

    finally:
        if h5file is not None:
            h5file.close()


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


def load_stats(statfun, default_dtype, *args, **kwargs):

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
