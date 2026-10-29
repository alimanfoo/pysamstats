# Reviving pysamstats: roadmap

The design is `specs/2026-10-29-revival/design.md`. Each stage below is one pull request, except stages marked **manual**, which a maintainer does by hand outside a pull request.

Every session that implements a stage should correct the specifications as details change, and say on its pull request what it corrected. It may change how it delivers the stage when earlier work or implementation discoveries justify that. It should raise any change to the roadmap's stages, boundaries, order or dependencies with the user instead of restructuring the roadmap itself.

## 1. Switch off the noise (manual)

**Delivers:**
- The pyup GitHub app is uninstalled from the repository. Deleting `.pyup.yml` alone isn't enough, because pyup falls back to its defaults.
- Bot PRs #194 (pyup) and #138 (dependabot) are closed, each with a one-line comment saying dependencies are moving to `pyproject.toml`.

**Proof:** the pyup app no longer appears in the repository's installed GitHub apps, and #194 and #138 are closed.

**Leaves for later:** deleting `.pyup.yml` itself, which is stage 3.

## 2. Merge #190 and #192 (manual)

**Delivers:** @DrYak's two PRs, merged into the default branch as written, so the fixes carry their author's credit:
- #190 changes the `'a'` dtype codes to `'S'` in `config.py`, `io.py`, `util.py` and the tests.
- #192 makes the pileup array pointer in `opt.pyx` `const`.

**Proof:** Bioconda's pysamstats build 16 has applied exactly these two patches as recipe patches since August 2026, and builds and passes its tests with them. There's no working CI to prove them in this repo yet. Stage 3 adds it.

**Leaves for later:** the rest of the `const` change, which is stage 9.

## 3. Modern build and green CI

**Delivers:** the package builds and its tests run on the current stack (pysam 0.24.1, Cython ≥3.1, numpy ≥2, CPython 3.11+), and CI proves it. Specifically:
- `pyproject.toml` is rewritten to the design's metadata, build requirements, runtime pins, `hdf5` extra and `test` dependency group.
- `setup.py` declares a plain `Extension` on `pysamstats/opt.pyx` with pysam's include directories and defines. There's no `cythonize()`, no pysam version check and no `distutils`. For now it keeps `scripts=['scripts/pysamstats']` and package discovery for the flat layout.
- `opt.pyx` gains a `# cython: language_level=3` header and nothing else.
- `pysamstats/opt.c` is deleted.
- A top-level MIT `LICENSE` (copyright Alistair Miles) is added.
- Tests move from nose to pytest helpers. The `Samfile` and `Fastafile` aliases stay for now.
- `test_io.py` imports PyTables with `pytest.importorskip`.
- The pileup reference comparisons cover the whole of `fixture/test.bam`, not just `Pf3D7_01_v3:0-2000`.
- `fixture/regression.npz` is regenerated on pysam 0.24.1 / htslib 1.24 with the existing `generate_fixtures()`. The commit message records the pysam and htslib versions and cites htslib commit `222971491`.
- `.github/workflows/python-package.yml` is replaced by a PR/push workflow on Linux x86_64 and macOS arm64, on one Python version. It installs with uv, builds the extension in place as the old CI did, and runs pytest from the repo root.
- `.travis.yml`, `.pyup.yml`, `requirements.txt` and `requirements_dev.txt` are deleted.

**Proof:** CI is green on both platforms, and that includes `test_against_fixtures`, which fails on the default branch today.

**Leaves for later:**
- the src layout and `MANIFEST.in` (stage 4);
- the console script (stage 7);
- the full Python matrix and the floors (stage 6);
- every other change inside `opt.pyx` (stages 8, 9 and 11);
- moving `generate_fixtures()` out of the test module (stage 4).

## 4. Src layout

**Delivers:**
- The package moves to `src/pysamstats/`.
- The test modules and fixtures move to a top-level `tests/` folder, and tests find their fixtures relative to their own file.
- `generate_fixtures()` moves into a script under `tests/` that pytest doesn't collect.
- `fixture/gc.txt` and `MANIFEST.in` are deleted.
- `setup.py`'s extension source path and package discovery follow the move.
- CI installs the package (not editable, not in place) and runs pytest against the installed copy.

**Proof:** CI is green with a non-editable install. A CI step, not a test (`uv sync` installs editable for developers), confirms that `pysamstats.__file__` resolves under site-packages, which shows the tests ran against the installed package.

**Leaves for later:** the console script (stage 7). `scripts=` moves with the layout and still points at `scripts/pysamstats`.

## 5. Release pipeline

**Delivers:** the release workflow described in the design, working as a dry run:
- It runs as a dry run on every pull request and every push to main, can also be dispatched by hand, and publishes only on `v*` tags.
- It checks the tag against `__version__` on tag runs.
- It builds the sdist, then uses cibuildwheel to build wheels from that sdist on Linux x86_64, Linux aarch64 and macOS arm64, for CPython 3.11–3.15.
- It pins pysam to the floor (0.24.1) as a build-only constraint, so the test install resolves pysam normally.
- It installs the `test` group in each wheel's test environment and runs the full suite against every wheel.
- A dry run stops after building and testing, and keeps the sdist and wheels as workflow artifacts.
- The publish job (PyPI Trusted Publishing through a `pypi` environment, then a GitHub release) exists, but only runs on tags.

**Proof:**
- The workflow's automatic dry run on this PR is green on all three runner types and all five Pythons, with 15 wheels and one sdist as artifacts. Manual dispatch can't prove anything before merge, because GitHub only dispatches workflows that exist on the default branch.
- The build logs show pysam 0.24.1 in the build environment and the newest 0.24.x in the test environment.

**Leaves for later:** the publish job's first real run, at stage 14. It can't run before then without publishing.

**Open, settled here:**
- The exact mechanism for the build-only pysam constraint in cibuildwheel. The requirement is that it applies to the build and not to the test install.
- Whether the `test` group installs on CPython 3.15 inside the manylinux containers. PyTables' dependencies numexpr and ndindex have no cp315 wheels yet and build from source. They built fine on macOS, but nobody has tried them in manylinux.

## 6. Floors matrix

**Delivers:**
- The PR/push workflow grows to CPython 3.11–3.15 on Linux x86_64 and macOS arm64.
- One more row runs on Python 3.11. It pins the build dependencies to their floors (setuptools 77, Cython 3.1) with a build constraint, and installs the runtime dependencies with uv's `--resolution lowest-direct` (pysam 0.24.1, numpy 2.0, tables 3.10.1).

**Proof:** every row is green, and the floors row's log shows those exact versions installed.

**Leaves for later:** nothing. This is the final shape of the matrix.

**Open, settled here:** whether the declared floors actually install and pass together. If one doesn't, this stage raises it in `pyproject.toml` and in the design, and says so on the PR.

## 7. Console script

**Delivers:**
- `scripts/pysamstats` becomes `src/pysamstats/cli.py` with `main(argv=None)`, registered as the console script `pysamstats = "pysamstats.cli:main"`.
- `scripts=` and the `scripts/` folder are deleted.
- The `OptionParser.format_epilog` monkeypatch becomes a small optparse subclass.
- Options, help text, output and exit behaviour are unchanged, quirks included.
- The command lines from the old CI become pytest tests that call `main(argv)` in-process. Each one asserts a zero exit, the expected header fields and the expected row count on `test.bam`. The HDF5 ones skip when PyTables is missing.

**Proof:** the new CLI tests pass, and `pysamstats --help` runs from the console script in CI.

**Leaves for later:** the HDF5 usage error for a missing `--output` (stage 11).

## 8. Internal tidy

**Delivers:**
- Remove the Python 2 code and the `__future__` imports from the package and the tests (`PY2`, `basestring`, the `izip_longest` shim).
- Replace the `DEF` constants in `opt.pyx` with `cdef enum` or the matching htslib declarations.
- `getrname` becomes `get_reference_name`.
- In the tests, `Samfile` and `Fastafile` become `AlignmentFile` and `FastaFile`.

There's no change in behaviour.

**Proof:** CI is green, `regression.npz` is unchanged, and the Cython build shows no deprecation warnings.

**Leaves for later:** the `const` signatures (stage 9) and the `assert` changes (stage 11).

## 9. Const fix finished

**Delivers:** in `opt.pyx`, `read` in `stat_pileup` and every `recv(...)` signature in the pileup stat classes take `const bam_pileup1_t*`, and so does any helper they pass it to.

**Proof:** on this PR, CI compiles the extension with discarded-qualifier warnings turned into errors (`-Werror=discarded-qualifiers` for gcc, `-Werror=incompatible-pointer-types-discards-qualifiers` for clang), and stays green with `regression.npz` unchanged. The design asks for a warning-free build, not a standing guard, so the flag needn't stay in the workflow after this PR.

**Leaves for later:** nothing.

## 10. #92 line endings

**Delivers:** `io.write_csv` sets the line terminator to `\n` explicitly, so TSV and CSV rows no longer end in `\r\n`.

**Proof:** a test asserts that `write_csv` output contains `\n` line endings and no `\r`.

**Leaves for later:** the changelog line (stage 13).

## 11. #127 asserts become errors

**Delivers:**
- In `opt.pyx`, the input-validation `assert`s (in `normalise_coords` and the `chrom is None` guards) raise `ValueError` with the same messages.
- In `cli.py`, HDF5 output without `--output` becomes a `parser.error`, which exits with status 2.
- `io.write_hdf5` raises an `ImportError` naming the `pysamstats[hdf5]` extra when PyTables is missing.

**Proof:** tests cover each of these. One test runs the library checks under `python -O`, in a subprocess, to show they no longer depend on `assert`.

**Leaves for later:** the changelog line (stage 13).

## 12. Housekeeping

**Delivers:**
- The README install section is rewritten: pip or uv, Bioconda, the BioContainers image, and the sentence about installing pysam and pysamstats from the same place.
- The `Dockerfile`, `release.txt`, `sandbox.ipynb` and `performance.py` are deleted.

**Proof:** CI is green, and the README's install commands match `pyproject.toml` and the release workflow.

**Leaves for later:** the release notes (stage 13).

## 13. Changelog and version

**Delivers:**
- A 1.2.0 entry in the README's release notes, covering everything in the design's list for it.
- `__version__` becomes `1.2.0`.

**Proof:** CI is green, including the release workflow's dry run on this PR.

**Leaves for later:** the tag (stage 14).

## 14. Publish 1.2.0 (manual)

**Delivers:**
- Configure the PyPI Trusted Publisher for `alimanfoo/pysamstats`, naming the release workflow file and the `pypi` environment.
- Create the `pypi` environment on GitHub.
- Push tag `v1.2.0`.

If the publish job fails on configuration, fix the setting on PyPI and re-run that job. The built wheels are reused.

**Proof:** PyPI shows 1.2.0 with 15 wheels and an sdist, and `pip install pysamstats==1.2.0` works on a clean machine, both Linux and macOS arm64.

**Leaves for later:** Bioconda (stage 15).

## 15. Bioconda, co-maintainer, issues (manual)

**Delivers:**
- @DrYak is invited as co-maintainer.
- The Bioconda autobump PR for 1.2.0 is edited into the design's target recipe, ideally together with @DrYak.
- The issues the design lists as resolved by 1.2.0 are closed, each with a pointer to the release.

**Proof:** Bioconda publishes pysamstats 1.2.0 with no patches, and `conda install -c conda-forge -c bioconda pysamstats=1.2.0` imports and runs `pysamstats --help`.

**Leaves for later:** nothing.
