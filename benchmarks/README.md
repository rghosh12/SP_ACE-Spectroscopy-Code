# Reference-star benchmarks

This directory is reserved for **curated regression cases** (standard stars with
published parameters and spectra you are allowed to use).

## Intended layout (future)

- ``README.md`` — table of star name, literature Teff/logg/[M/H], spectrum URL,
  citation, and tolerance used in tests.
- ``cases/<name>/`` — ``space.par`` (or generator script), optional small
  spectrum excerpt, and **expected** ``*_TGM_ABD.dat`` line (or scalar checks).

## Why empty for now

Full GCOG libraries and proprietary spectra cannot be invented here. Add your
first case once you have:

1. A working Fortran ``space`` binary (see ``Fortran-95_code/Makefile``).
2. Legal rights to store the spectrum (or a script that downloads it).
3. Documented tolerances (pipeline vs literature systematics).

CI can later run ``pytest -m benchmark`` only when ``BENCHMARK_DATA=1`` is set.
