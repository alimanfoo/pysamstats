pysamstats
==========

A small Python utility for calculating statistics per genome position
based on pileups from a SAM or BAM file.

* Source: https://gihub.com/alimanfoo/pysamstats
* Download: http://pypi.python.org/pypi/pysamstats

Installation
------------

```
$ pip install --upgrade pysamstats
```

N.B., pysamstats depends on [pysam](http://code.google.com/p/pysam/)
and [numpy](http://www.numpy.org/). These *should* be installed
automatically if you run the command above, but if you have any
problems, you might try installing pysam and numpy separately first.

Alternatively, clone the git repo and build in-place (requires
cython):

```
$ git clone git://github.com/alimanfoo/pysamstats.git
$ cd pysamstats
$ python setup_dev.py build_ext --inplace 
```

Usage
-----

From the command line:

```
$ pysamstats --help
Usage: pysamstats [options] FILE

Calculate statistics per genome position based on pileups from a SAM or BAM
file and print them to stdout.

Options:
  -h, --help            show this help message and exit
  -t TYPE, --type=TYPE  type of statistics to print: coverage,
                        coverage_strand, coverage_ext, coverage_ext_strand,
                        coverage_normed, coverage_gc, coverage_normed_gc,
                        variation, variation_strand, tlen, tlen_strand, mapq,
                        mapq_strand, baseq, baseq_strand, baseq_ext,
                        baseq_ext_strand
  -c CHROMOSOME, --chromosome=CHROMOSOME
                        chromosome name
  -s START, --start=START
                        start position (1-based)
  -e END, --end=END     end position (1-based)
  -z, --zero-based      use zero-based coordinates (default is false, i.e.,
                        use one-based coords)
  -f FASTA, --fasta=FASTA
                        reference sequence file, only required for some
                        statistics
  --gc-window-length=N  size of window to use for %GC calculations [300]
  --gc-window-offset=N  window offset to use for deciding which genome
                        position to report %GC calculations against [150]
  -o, --omit-header     omit header row from output
  -p N, --progress=N    report progress every N rows

Supported statistics types:

    * coverage            - number of reads aligned to each genome position 
                            (total and properly paired)
    * coverage_strand     - as coverage but with forward/reverse strand counts
    * coverage_ext        - various additional coverage metrics, including 
                            coverage for reads not properly paired (mate 
                            unmapped, mate on other chromosome, ...)
    * coverage_ext_strand - as coverage_ext but with forward/reverse strand counts 
    * coverage_normed     - depth of coverage normalised by median or mean
    * coverage_gc         - as coverage but also includes a column for %GC
    * coverage_normed_gc  - as coverage_normed but also includes columns for normalisation
                            by %GC      
    * variation           - numbers of matches, mismatches, deletions, 
                            insertions, etc.
    * variation_strand    - as variation but with forward/reverse strand counts
    * tlen                - insert size statistics
    * tlen_strand         - as tlen but with statistics by forward/reverse strand
    * mapq                - mapping quality statistics
    * mapq_strand         - as mapq but with statistics by forward/reverse strand
    * baseq               - baseq quality statistics
    * baseq_strand        - as baseq but with statistics by forward/reverse strand
    * baseq_ext           - extended base quality statistics, including qualities
                            of bases matching and mismatching reference
    * baseq_ext_strand    - as baseq_ext but with statistics by forward/reverse strand
    
Examples:

    pysamstats --type coverage example.bam > example.coverage.txt
    pysamstats --type coverage --chromosome Pf3D7_v3_01 --start 100000 --end 200000 example.bam > example.coverage.txt
```

From Python:

```python
import pysam
import pysamstats

mybam = pysam.Samfile('/path/to/your/bamfile.bam')
for rec in pysamstats.stat_coverage(mybam, chrom='Pf3D7_01_v3', start=10000, end=20000):
    print rec['chrom'], rec['pos'], rec['reads_all'], rec['reads_pp']
    ...

```

Field Definitions
-----------------

The suffix **_fwd** means the field is restricted to reads mapped to
the forward strand, and **_rev** means the field is restricted to
reads mapped to the reverse strand. E.g., **reads_fwd** means the
number of reads mapped to the forward strand.

The suffix **_pp** means the field is restricted to reads flagged as
properly paired. 

* **chrom** - Chromosome name.  

* **pos** - Position within chromosome. One-based by default when
    using the command line, zero-based by default when using the
    python API.

* **reads_all** - Number of reads aligned at the position. N.b., this
    is really the total, i.e., includes reads where the mate is
    unmapped or otherwise not properly paired.

* **reads_pp** - Number of reads flagged as properly paired by the
    aligner.

* **reads_mate_unmapped** - Number of reads where the mate is
    unmapped.

* **reads_mate_other_chr** - Number of reads where the mate is mapped
    to another chromosome.

* **reads_mate_same_strand** - Number of reads where the mate is
    mapped to the same strand.

* **reads_faceaway** - Number of reads where the read and its mate are
    mapped facing away from each other.

* **reads_softclipped** - Number of reads where there is some
    softclipping at some point in the read's alignment (not
    necessarily at this position).

* **reads_duplicate** - Number of reads that are flagged as duplicate.

* **dp_normed_median** - Number of reads divided by the median number
    of reads over all positions in the specified region, or whole
    genome if no region specified.

* **dp_normed_mean** - Number of reads divided by the mean number of
    reads over all positions in the specified region, or whole genome
    if no region specified.

* **dp_percentile** - Percentile within which the number of reads falls
    considering all positions in the specified region, or whole genome
    if no region specified.

* **gc** - Percentage GC content in the reference at this position
    (depends on window length and offset specified).

* **dp_normed_median_gc** - As *dp_normed_median* but normalised by
    positions with the same percent GC composition.

* **dp_normed_mean_gc** - As *dp_normed_mean* but normalised by
    positions with the same percent GC composition.

* **dp_percentile_gc** - As *dp_percentile* but only considering
    positions with the same percent GC composition.

* **matches** - Number of reads where the aligned base matches the
    reference.

* **mismatches** - Number of reads where the aligned base does not
    match the reference (but is not a deletion).

* **deletions** - Number of reads where there is a deletion in the
    alignment at this position.

* **insertions** - Number of reads where there is an insertion in the
    alignment at this position.

* **A/C/T/G/N** - Number of reads where the aligned base is an A/C/T/G/N.

* **mean_tlen** - Mean value of outer distance between reads and their
    mates for paired reads aligned at this position. N.B., leftmost
    reads in a pair have a positive tlen, rightmost reads have a
    negative tlen, so if there is no strand bias, this value should be
    0.

* **rms_tlen** - Root-mean-square value of outer distance between
    reads and their mates for paired reads aligned at this position.

* **std_tlen** - Standard deviation of outer distance between reads
    and their mates for paired reads aligned at this position.

* **reads_mapq0** - Number of reads where mapping quality is zero.

* **rms_mapq** - Root-mean-square mapping quality for reads aligned at
    this position.

* **max_mapq** - Maximum value of mapping quality for reads aligned at
    this position.

* **rms_baseq** - Root-mean-square value of base qualities for bases
    aligned at this position.

* **rms_baseq_matches** - Root-mean-square value of base qualities for
    bases aligned at this position where the base matches the
    reference.

* **rms_baseq_mismatches** - Root-mean-square value of base qualities
    for bases aligned at this position where the base does not match
    the reference.

