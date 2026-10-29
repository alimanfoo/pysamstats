import csv
import io
from pathlib import Path

import pytest

from pysamstats.cli import main


FIXTURE_DIR = Path(__file__).parent / "fixture"
ALIGNMENT_PATH = FIXTURE_DIR / "test.bam"
REFERENCE_PATH = FIXTURE_DIR / "ref.fa"

COVERAGE_FIELDS = ["chrom", "pos", "reads_all", "reads_pp"]
BINNED_COVERAGE_FIELDS = ["chrom", "pos", "gc", "reads_all", "reads_pp"]

COVERAGE_CASE = pytest.param(
    ["--type", "coverage"],
    COVERAGE_FIELDS,
    17956,
    1553308,
    id="coverage",
)
BINNED_COVERAGE_CASE = pytest.param(
    ["--type", "coverage_binned", f"--fasta={REFERENCE_PATH}"],
    BINNED_COVERAGE_FIELDS,
    602,
    20698,
    id="binned-coverage",
)


def assert_text_rows(rows, expected_fields, expected_row_count, expected_reads_all):
    assert rows[0] == expected_fields
    assert len(rows) - 1 == expected_row_count
    reads_all_index = expected_fields.index("reads_all")
    assert sum(int(row[reads_all_index]) for row in rows[1:]) == expected_reads_all


@pytest.mark.parametrize(
    ("arguments", "expected_fields", "expected_row_count", "expected_reads_all"),
    [
        COVERAGE_CASE,
        pytest.param(
            ["--type", "coverage", "--fields=pos,reads_all"],
            ["pos", "reads_all"],
            17956,
            1553308,
            id="selected-fields",
        ),
        BINNED_COVERAGE_CASE,
        pytest.param(
            [
                "--type=coverage",
                "--min-mapq=27",
                "--min-baseq=17",
                "--no-dup",
                "--no-del",
            ],
            COVERAGE_FIELDS,
            17956,
            205577,
            id="filtered-coverage",
        ),
        pytest.param(
            [
                "--type=coverage_binned",
                f"--fasta={REFERENCE_PATH}",
                "--min-mapq=27",
                "--no-dup",
            ],
            BINNED_COVERAGE_FIELDS,
            602,
            2981,
            id="filtered-binned-coverage",
        ),
    ],
)
def test_main_writes_text_to_stdout(
    arguments, expected_fields, expected_row_count, expected_reads_all, capsys
):
    assert main([*arguments, str(ALIGNMENT_PATH)]) is None

    rows = list(csv.reader(io.StringIO(capsys.readouterr().out), dialect="excel-tab"))
    assert_text_rows(rows, expected_fields, expected_row_count, expected_reads_all)


@pytest.mark.parametrize(
    ("arguments", "expected_fields", "expected_row_count", "expected_reads_all"),
    [COVERAGE_CASE, BINNED_COVERAGE_CASE],
)
def test_main_writes_text_file(
    arguments,
    expected_fields,
    expected_row_count,
    expected_reads_all,
    tmp_path,
    capsys,
):
    output_path = tmp_path / "output.txt"

    assert main(
        [*arguments, f"--output={output_path}", str(ALIGNMENT_PATH)]
    ) is None

    assert capsys.readouterr().out == ""
    with output_path.open(newline="") as output:
        rows = list(csv.reader(output, dialect="excel-tab"))
    assert_text_rows(rows, expected_fields, expected_row_count, expected_reads_all)


@pytest.mark.parametrize(
    ("arguments", "expected_fields", "expected_row_count", "expected_reads_all"),
    [COVERAGE_CASE, BINNED_COVERAGE_CASE],
)
def test_main_writes_hdf5_file(
    arguments, expected_fields, expected_row_count, expected_reads_all, tmp_path
):
    tables = pytest.importorskip("tables")
    output_path = tmp_path / "output.h5"

    assert main(
        [
            *arguments,
            f"--output={output_path}",
            "--format=hdf5",
            str(ALIGNMENT_PATH),
        ]
    ) is None

    with tables.open_file(output_path) as h5file:
        assert h5file.root.data.colnames == expected_fields
        assert h5file.root.data.nrows == expected_row_count
        assert h5file.root.data.col("reads_all").sum() == expected_reads_all
