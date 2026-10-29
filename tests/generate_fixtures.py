from pathlib import Path

import numpy as np
from pysam import AlignmentFile

import pysamstats
from pysamstats.config import stats_types, stats_types_withref


FIXTURE_DIR = Path(__file__).with_name("fixture")


def generate_fixtures():
    bampath = str(FIXTURE_DIR / "test.bam")
    fastapath = str(FIXTURE_DIR / "ref.fa")
    archive = FIXTURE_DIR / "regression.npz"
    assert not archive.is_file()

    data = {}
    for stats_type in stats_types:
        load = getattr(pysamstats, "load_" + stats_type)
        if stats_type in stats_types_withref:
            data[stats_type] = load(AlignmentFile(bampath), fafile=fastapath)
        else:
            data[stats_type] = load(AlignmentFile(bampath))

    np.savez_compressed(archive, **data)


if __name__ == "__main__":
    generate_fixtures()
