pysamstats
==========

A small Python utility for calculating statistics per genome position
based on pileups from a SAM or BAM file.

* Source: https://gihub.com/alimanfoo/pysamstats
* Download: http://pypi.python.org/pypi/pysamstats (TODO)

Installation
------------

```
$ pip install --upgrade pysam pysamstats
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
