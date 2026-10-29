# Reviving pysamstats: design

## What we're building

pysamstats 1.2.0. The library and command-line tool stay as they are, with the same functions, the same statistics types and the same output columns. What changes is everything around them. The release installs with `pip install pysamstats` (or `uv add pysamstats`) on CPython 3.11–3.15, from prebuilt wheels, with no compiler needed. It also installs from Bioconda without the recipe patches Bioconda carries today. CI goes green and stays green, and it tests the same things that get shipped. Cutting the next release becomes "push a tag".

Where things stand (all checked on 2026-10-29): `pip install pysamstats` fails on any modern Python. The 1.1.2 sdist pins `pysam<0.16`, which forces a source build of pysam 0.15.4, and that dies with `No module named 'pkg_resources'`. Installing from the default branch fails the same way. Bioconda works, but only because the recipe (build 16) rebuilds the 2018 tarball with PR #190 and PR #192 applied as patches. CI has failed on every run since at least mid-2026. With #190 and #192 applied, the Cython core builds against pysam 0.24.1 with Cython 3.3.0 on Python 3.12 and 3.14 (tried on macOS arm64). It then matches the pure-Python reference implementations in the test suite exactly, for all 15 pileup statistics types, with and without read filters, at every position of `fixture/test.bam`.

## How it works

### One case, end to end

The maintainer bumps `__version__` to `1.2.0` in `src/pysamstats/__init__.py`, merges, and pushes tag `v1.2.0`. The release workflow starts by checking that the tag matches `__version__`, and stops if it doesn't. Next it builds the sdist. Then cibuildwheel builds wheels from that sdist on three runner types (Linux x86_64, Linux aarch64, macOS arm64) for CPython 3.11–3.15. Each build installs pysam at the floor of the supported range (0.24.1) as a build dependency, runs Cython over `opt.pyx`, and compiles the extension against the pysam headers. cibuildwheel installs each finished wheel into a clean environment, along with the wheel's own runtime dependencies, and runs the test suite from `tests/` against that installed copy. If every wheel passes, the workflow publishes them all to PyPI through Trusted Publishing and creates a GitHub release.

A user on a Linux cluster then runs `pip install pysamstats` under Python 3.13. pip picks our `cp313` manylinux wheel, and the newest pysam wheel inside `>=0.24.1,<0.25`, and numpy ≥2.0. When `import pysamstats` loads the compiled `opt` module, Cython checks pysam's extension types at import time. That passes because pysam is inside the tested range. `pysamstats --type coverage x.bam` runs through the `pysamstats` console script.

A day or so later, Bioconda's autobump bot opens a recipe PR for 1.2.0. That PR fails to build until someone edits it into the target recipe described under "Bioconda" below.

### The pysam coupling (the central constraint)

`opt.pyx` cimports pysam's Cython declarations and reads C struct fields directly: `PileupColumn.plp`, and `IteratorRowRegion.b`, `.cnext()` and `.retval`. When the module is imported, the generated code looks up 24 pysam extension types, even though only about 5 are used, and checks their sizes. If any of those types is missing, the import fails with `AttributeError`, as in #157. If one is smaller than at build time, it fails with "wrong size", as in #94 and #45. If one is larger, the import warns, as in #112.

On the other hand, the built extension links no pysam or htslib symbols at all (`nm -u` lists only libc). Everything it takes from htslib comes through static-inline helpers and struct layouts. So a wheel holds nothing but our own module. There are no shared libraries to bundle, no rpath to manage, and nothing for auditwheel or delocate to vendor.

Rules that follow from this:

- **Build against the floor, run anywhere in the range.** Cython tolerates types that grew and rejects types that shrank or vanished, so the safe direction is built-old, run-new. Tested: built on 0.23.3 and run on 0.24.1 passes, built on 0.24.1 and run on 0.23.3 passes, and built on 0.24.1 and run on 0.22.1 fails on import. The floor is 0.24.1 rather than 0.24.0 because 0.24.0 has no CPython 3.15 wheels. Building against it on 3.15 would mean compiling pysam from source.
- **The runtime requirement has an upper bound,** `pysam>=0.24.1,<0.25`. This is the one case where an upper bound on a library is justified. Without it, pip will happily install a combination that fails on import. The build requirement uses the same range. Wheel builds pin to the floor through a build constraint.
- **Widening the range is a deliberate act.** When pysam ships a new minor version (GitHub's release notification for pysam is the trigger), the maintainer opens a PR that widens the upper bound, for example to `<0.26`. The release workflow's dry run runs on it automatically. That builds against the floor and tests every wheel against the newest pysam the widened range allows, on every platform and Python. If it passes, the branch becomes a patch release. If it fails, the next release moves the floor up and rebuilds against it.
- **Mixed installs are the user's responsibility.** pysam from conda with pysamstats from pip, `--no-deps`, and similar setups get no special handling beyond one sentence in the README install section: install pysam and pysamstats from the same place.

### Source layout

It moves to the src layout:

- `src/pysamstats/` holds `__init__.py`, `opt.pyx`, `pileup.py`, `binned.py`, `config.py`, `io.py`, `util.py` and a new `cli.py`.
- `tests/` sits at the top level and holds the test modules and their fixtures (moved from `fixture/`).
- Tests find fixtures relative to their own file, not the current directory.

The point is that the only importable `pysamstats` is the installed one. Today, running pytest from the repo root imports the source folder, which has no compiled extension in it. That would defeat cibuildwheel's test-the-installed-wheel step. As a side effect, the test modules stop shipping inside wheels (`find_packages()` includes `pysamstats/test` today). The fixtures were never shipped, and that doesn't change. `fixture/gc.txt` isn't referenced anywhere, so it's deleted rather than moved.

### Build

`pyproject.toml` holds the standard project metadata (PEP 621):

- name, description, `readme = "README.md"`, authors, project URLs and classifiers (without `readme`, PyPI shows no description);
- `license = "MIT"` as a PEP 639 expression, with a top-level `LICENSE` file, copyright Alistair Miles;
- `requires-python = ">=3.11"`;
- runtime dependencies `pysam>=0.24.1,<0.25` and `numpy>=2.0`;
- an optional extra `hdf5 = ["tables>=3.10.1"]` (3.10.1 is the first PyTables release that accepts numpy 2);
- a console script `pysamstats = "pysamstats.cli:main"`;
- `version` as a dynamic field, which setuptools reads statically from `pysamstats.__version__`.

The build system is setuptools, with build requirements `setuptools>=77` (the first version that understands a PEP 639 license expression), `cython>=3.1` and `pysam>=0.24.1,<0.25`.

`setup.py` shrinks to one job: declaring the `pysamstats.opt` extension as a plain `Extension` whose source is `src/pysamstats/opt.pyx`, with pysam's include directories and defines. It does **not** wrap the extension in `cythonize()`. setuptools then Cythonizes the `.pyx` at build time, and, just as importantly, puts the `.pyx` in the sdist. When `cythonize()` is used, the sdist carries a generated `opt.c` and no `.pyx`, and building from it fails with `'src/pysamstats/opt.pyx' doesn't match any files` (reproduced). Cython directives such as `language_level=3` go in a `# cython:` header comment at the top of `opt.pyx`, not in `setup.py`. The import-time pysam version check, the `distutils.LooseVersion` call and the duplicated metadata all go. The generated `opt.c` is no longer committed.

The version stays a plain string in `__init__.py`, not something derived from git. Bioconda builds from GitHub's archive tarball, which has no `.git` folder for a tag-based tool to read.

### Runtime and CLI

The package API is unchanged. The `scripts/pysamstats` file becomes `src/pysamstats/cli.py`:
- Its body is wrapped in `main(argv=None)`, so tests can call it in-process.
- It keeps optparse and the same options, help text and exit behaviour, quirks included. For example, `--format` is case-sensitive, and an unknown format writes nothing and exits 0. The only exceptions are the two fixes below.
- The `OptionParser.format_epilog` monkeypatch, which today modifies optparse's class for the whole process, moves into a small subclass, so importing `cli` doesn't change optparse for anyone else.

Two long-standing issues are fixed in 1.2.0, both with changelog lines:
- **#92.** TSV and CSV output end rows with `\n`, not `\r\n`. `io.write_csv` sets the line terminator explicitly instead of inheriting the `excel` dialects' `\r\n`. This changes output bytes. Normal CSV readers won't notice.
- **#127.** Input validation no longer uses `assert`, so it survives `python -O`.
  - In the library, the checks in `opt.pyx`'s `normalise_coords` (chromosome not in the BAM's references) and the `chrom is None` guards raise `ValueError` with the same messages.
  - In the CLI, HDF5 output without `--output` becomes a `parser.error`, which exits 2 with a message like the CLI's other usage errors.

Asking for HDF5 output without PyTables fails inside `io.write_hdf5`, not in the CLI, so library users get the same message. It's an `ImportError` that names the `pysamstats[hdf5]` extra, and from the CLI it ends with a non-zero exit like any other uncaught error.

Two changes from open PRs make the code work with current dependencies:
- #190 changes the `'a'` dtype codes to `'S'`, which numpy ≥2.5 requires. It touches `config.py`, `io.py`, `util.py` and the tests, not `opt.pyx`.
- #192 marks the pileup array pointer in `opt.pyx` `const`, which recent htslib with gcc requires. It leaves a discarded-qualifier warning on `read = &(plp[0][i])`. This work finishes the job: `read` and every `recv(...)` signature in the pileup stat classes take a `const bam_pileup1_t*`, and so does any helper they pass it to. All of those only read from it, so the build is warning-free. That matters because gcc 14 already turned the neighbouring incompatible-pointer warning into an error, which is what broke Bioconda in the first place.

The `opt.pyx` internals are also brought up to date, with no change in behaviour. The aim is to leave no deprecated or dead code that a future pysam, Cython or compiler release would force us to revisit:
- Remove the Python 2 branches (`PY2`, `basestring`).
- Replace the `DEF` constants, which Cython 3 deprecates, with `cdef enum` or the matching htslib declarations.
- Replace `getrname` with `get_reference_name`.

The Python modules and tests also lose their `__future__` imports and Python 2 shims.

### Tests

The tests move from nose to pytest:
- `eq_` becomes `assert ==`;
- `assert_raises` becomes `pytest.raises`;
- `assert_almost_equal` becomes `pytest.approx`;
- the old pysam aliases `Samfile` and `Fastafile` become `AlignmentFile` and `FastaFile`.

There are three layers:
1. **Reference comparisons.** Each pileup statistics type is compared against a pure-Python reference implementation built on pysam's own pileup API. These are widened from `Pf3D7_01_v3:0-2000` to the whole of `test.bam`, which costs a few seconds. The binned-statistics tests stay as they are.
2. **Regression fixture.** This already exists: `test_regression.py::test_against_fixtures` checks every statistics type, binned ones included, over the whole BAM against `regression.npz`. It has been failing since htslib 1.13. The existing `generate_fixtures()` helper moves out of the test module into a small script under `tests/` that pytest doesn't collect. It locates its files relative to itself, and records the pysam and htslib versions in the commit that regenerates the fixture. The fixture is regenerated once, on pysam 0.24.1 / htslib 1.24.
3. **CLI smoke tests.** The command-line invocations from the old CI move into pytest, calling `main(argv)` in-process. Each one asserts a zero exit, a header row with the expected fields, and the expected number of rows on `test.bam`.

PyTables is optional throughout. The HDF5 tests, including `test_io.py`, which today imports `tables` at module level and so breaks collection when it's missing, use `pytest.importorskip("tables")`. Every CI and cibuildwheel test environment installs the `test` dependency group (pytest plus tables), so HDF5 always runs in CI.

The regression fixture has a job the reference comparisons can't do. They run on the same htslib as the code under test, so when htslib changes behaviour, both drift together and stay green. Only the fixture notices when the numbers users get have changed. That's exactly how htslib's April 2021 change to overlapping-pair quality handling showed up (commit `222971491`, "Make mpileup's overlap removal choose a random sequence"). Against today's fixture, the change shows up only in `baseq_strand`, `baseq_ext` and `baseq_ext_strand`:
- the strand-split `rms_baseq*` fields differ at up to about 3,600 of 17,956 positions;
- the strand-split `*_mismatches_*` fields differ at about 70 positions;
- the unsplit `rms_baseq_matches*` and `rms_baseq_mismatches*` in `baseq_ext` differ at 1–5 positions each.

Counts, coverage, mapq, tlen and variation are all unchanged.

The rule for the fixture: when it fails after a pysam bump, find the upstream cause, regenerate the fixture, and add a changelog line. Never loosen the comparison.

### CI

There are two workflows:

1. **Pull requests and pushes to main.** Build from source and run the tests on Linux x86_64 for CPython 3.11–3.15, and on macOS arm64 for 3.11 and 3.15 only. GitHub Free allows five macOS jobs at once across the whole account, and the release dry run already builds and tests a macOS wheel for every Python on every PR, so the middle macOS rows only duplicated it and slowed the queue. One extra row, on Python 3.11, uses the floor of every declared direct dependency. Runtime and test dependencies come from uv's `--resolution lowest-direct` (pysam 0.24.1, numpy 2.0, pytest 7.2.1, tables 3.10.1). Build dependencies are pinned to their floors with a build constraint (setuptools 77.0.1, Cython 3.1.0, pysam 0.24.1). That tests every declared direct dependency floor.
2. **Release.** Runs as a dry run, which builds and tests everything but doesn't publish, on every pull request, every push to main and manual dispatch. It publishes only on `v*` tags. That's the shape cibuildwheel's own example uses. The build is a single small extension, so running it on every PR costs a few minutes, and it avoids a path filter that would need maintaining. In order, it does the following:
   - on tags, checks the tag against `__version__`;
   - builds the sdist;
   - builds the wheels **from that sdist**, not from the git checkout, so a broken sdist fails the release instead of shipping. cibuildwheel can do this directly;
   - publishes to PyPI through Trusted Publishing (a `pypi` environment, no stored token);
   - creates the GitHub release.

   The pysam floor pin applies to the build step only, as a build constraint. cibuildwheel uses pip as its build frontend so pip receives that constraint, while the test install receives no constraint and resolves pysam normally to the newest version inside the range.

Both workflows share a concurrency group per workflow and ref for pull requests. A new push to a pull request cancels that PR's older run, so superseded runs don't hold the queue. Every run on `main` or a tag has a unique group, so neither an in-progress nor a queued run can be cancelled by a newer run.

There's no separate pysam-compatibility workflow, because the release dry run on a widening PR does that job (see "The pysam coupling"). The cost: trying an arbitrary pysam version, such as an older one or a pre-release, means editing `pyproject.toml` on a branch.

There's no scheduled run either. pysam releases two or three times a year, and GitHub disables scheduled workflows after 60 days with no repository activity, which on this repo would turn the job off exactly when it was needed.

### Platforms and Python versions

Wheels are built for manylinux_2_28 x86_64 and aarch64, and for macOS arm64.

There are no macOS x86_64 wheels. GitHub's last Intel runner (`macos-15-intel`) goes away around August 2027, and a cross-compiled wheel would be one nobody has run. Intel Mac users can install from Bioconda, which builds osx-64, or build from the sdist, which needs Xcode's command-line tools. Adding the target back while the runner exists would cost one more runner type in the release workflow.

Python versions are CPython 3.11–3.15. A version is dropped at CPython's own end of life (3.11 in October 2027), not on SPEC 0's earlier date, so that conda pipeline environments on older Pythons keep receiving releases. There is no musllinux, no Windows (pysam has no Windows wheels), and no free-threaded build (pysam has no free-threaded wheels). The sdist covers everything else.

### Bioconda

This is the target recipe for 1.2.0. It's applied by hand to the autobump PR after the release, as a maintainer step that's outside the code work:
- Delete both patches, since the code they fix is now in the source.
- Delete `build.sh` and replace it with `script: $PYTHON -m pip install . --no-deps --no-build-isolation -vv` in `meta.yaml`. Its `C_INCLUDE_PATH` export was only there for the zlib headers.
- Host requirements: `python`, `pip`, `setuptools >=77`, `cython >=3.1`, `pysam >=0.24.1,<0.25`.
- Drop `zlib` from host and run. The extension doesn't include or link zlib, and pysam brings its own.
- Runtime pysam pin: `pin_compatible('pysam', max_pin='x.x')`, replacing `x.x.x`, which is what freezes conda users on pysam 0.24.0 today.
- Runtime dependencies keep `numpy` and `pytables`.

The import test, the `pysamstats --help` test and the additional platforms stay as they are. BioContainers builds a Docker image from the recipe automatically, so the repository's own Dockerfile goes.

### Repository and project hygiene

Some of these steps need maintainer access and can't be done from a pull request: setting up the PyPI Trusted Publisher and the `pypi` environment, uninstalling pyup, editing the Bioconda recipe, closing issues and PRs, and inviting a co-maintainer.

- Add `LICENSE` (MIT, Alistair Miles). This closes #132.
- Rewrite the README install section: pip/uv, Bioconda, the BioContainers image, and the same-channel sentence. Release notes stay in the README. The 1.2.0 entry covers the new Python and pysam support, the htslib overlap effect on base-quality statistics (`baseq_strand`, `baseq_ext`, `baseq_ext_strand`), the `\n` line endings (#92), `ValueError` instead of `AssertionError` for bad input (#127), and the console script replacing `scripts/pysamstats`.
- Turn off pyup. Deleting `.pyup.yml` isn't enough, because pyup falls back to its defaults when the file is missing. The pyup GitHub app has to be uninstalled from the repo.
- Close bot PRs #194 and #138. Merge (or supersede, with credit) #190 and #192.
- After the release, close the issues that 1.2.0 resolves: #45, #56, #58, #83, #92, #94, #95, #112, #122, #125, #127, #157. Each closing comment points to 1.2.0.
- Invite @DrYak as co-maintainer once 1.2.0 is out. The Bioconda recipe update is a natural first piece of work to do together.

## What changes, and what goes away

- **`pysamstats/opt.c`** (65k lines generated by Cython 0.29.23). Deleting it ends the stale-generated-C class of bug (#122) and removes a file every pysam bump would otherwise need regenerated. Cost: Cython becomes a hard build requirement, which every route (pip, cibuildwheel, Bioconda) already provides.
- **The pysam version check, LooseVersion and duplicated metadata in `setup.py`; the inert `[tool.pysamstats]` table in `pyproject.toml`.** This buys one source of truth for metadata, and a build that works on Python 3.12+, where distutils is gone.
- **`scripts/pysamstats`.** Replaced by a console entry point that works the same on every platform and under every installer.
- **`fixture/` at the root and `pysamstats/test/` inside the package.** Moved under `tests/`, for the import-isolation reason given under Source layout. `fixture/gc.txt` (298 KB) is deleted. Nothing references it.
- **`.travis.yml`.** Travis has been dead for this repo since 2021.
- **`.pyup.yml`, `requirements.txt`, `requirements_dev.txt`.** Dependencies are declared once, in `pyproject.toml`. The pinned dev requirements are what pyup kept bumping, producing 74 bot PRs and a red CI every month. Development environments come from `uv sync` against a `test` dependency group (pytest, tables).
- **`MANIFEST.in`.** With the src layout and a plain `Extension` whose source is the `.pyx`, setuptools puts the `.pyx` in the sdist by itself, and the fixtures stay out. The release workflow's build-wheels-from-the-sdist step checks this on every release.
- **`Dockerfile`** (miniconda 4.5.4 from 2018). BioContainers covers this. Closes #58.
- **`release.txt`** (`setup.py register sdist upload`). Replaced by the release workflow.
- **`sandbox.ipynb` and `performance.py`.** Scratch files. They remain in git history if anyone wants them, and rewriting the benchmark would be about 70 lines.
- **nose**, which doesn't run on Python 3.12+.
- **Python 2 compatibility code throughout.** Dead since `requires-python >= 3.11`.

## Why it's this and not something else

- **Bioconda only (the status quo).** Rejected. It shuts out pip and uv users, and it depends on third parties patching the recipe. The PyPI wheels that were impractical in 2018 are practical now, because pysam itself ships wheels for every platform we target.
- **An exact pysam pin per release (Bioconda's model, copied to PyPI).** Rejected. It would force a pysamstats release for every pysam patch release, and the cross-version tests show the range works.
- **Removing the pysam coupling** (link htslib ourselves, or rewrite in Rust on noodles or rust-htslib). Rejected for now. It would end the ABI dance for good, but it's a rewrite that has to reproduce pysam's pileup steppers and htslib's overlap handling exactly, and it means shipping our own htslib build. Worth reconsidering only if minor pysam releases keep breaking the coupling.
- **Pure Python on pysam's public API.** Rejected. Measured 12× slower for `coverage` and 31–32× slower for `variation` and `baseq_ext_strand`, mostly because per-read flags force a Python object per read.
- **A scheduled canary, a separate pysam-compatibility workflow, and an import-time guard for mixed installs.** All cut. The canary disables itself on a quiet repo. The compatibility workflow duplicated what a release dry run on a widening PR already tests. The guard was machinery for a mistake the user can avoid.
- **setuptools-scm.** Rejected because Bioconda's tarball builds have no git metadata.
- **meson-python or scikit-build-core.** No advantage for a single Cython extension with nothing to link.

## What's still open

- How the work is split into pull requests. That's being planned separately.
- Whether the `lowest-direct` combination (numpy 2.0, tables 3.10.1, pysam 0.24.1 on Python 3.11) actually installs and passes. Its first CI run will show. If it doesn't, the floors move up to whatever does.
- Linux has not been built or tested yet. Everything above was verified on macOS arm64. The first CI run on Linux is the first real evidence there.
