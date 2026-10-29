import subprocess
import sys

from .util import fixture_path


def test_library_validation_survives_optimized_python():
    script = """
import sys

from pysam import AlignmentFile
import pysamstats
from pysamstats import opt

alignment = AlignmentFile(sys.argv[1])
checks = [
    lambda: next(pysamstats.stat_coverage(alignment, chrom="missing")),
    lambda: next(opt.iter_pileup_padded_chrom(
        opt.Coverage(), alignment, None, None, None, None, False, False,
        "all", 8000, 0, 0, False, False,
    )),
    lambda: next(opt.iter_binned_chrom(
        opt.CoverageBinned(), alignment, None, None, None, None, False,
        300, 150, 0, False,
    )),
]

for check in checks:
    try:
        check()
    except ValueError as error:
        print(error)
    else:
        raise RuntimeError("validation did not raise ValueError")
"""

    result = subprocess.run(
        [sys.executable, "-O", "-c", script, fixture_path("test.bam")],
        capture_output=True,
        check=False,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == [
        "chromosome not in SAM references: missing",
        "chromosome is None",
        "chromosome is None",
    ]
