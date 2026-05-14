# Fortran SP_Ace (`PROGRAM space`)

This directory contains the **reference** Fortran 95 implementation (Corrado Boeche, GPLv3).

## Prerequisites

- **gfortran** (GCC Fortran). On Ubuntu/Debian:

  ```bash
  sudo apt-get install gfortran liblapack-dev libblas-dev
  ```

- **LAPACK + BLAS** — `uncertains2.f95` calls `dgelss` / `ilaenv`. The root `Makefile` links `-llapack -lblas` after the object sources.

## Build (single command)

From the **repository root**:

```bash
make -C Fortran-95_code
```

Or from this directory:

```bash
make
```

Produces the executable:

| Platform   | Output   |
|-----------|----------|
| Linux/macOS | `space` |
| Windows (MinGW) | `space.exe` |

Compiler flags default to `-O2 -ffree-line-length-none` (override with `make FFLAGS=...`).

## Run

The program reads the parameter file from the first command-line argument (default `space.par` if omitted in the Fortran `getarg` handling — see `space.f95`).

```bash
./space /path/to/run.par
```

Set **`SP_ACE_EXE`** to this binary when using the Python launcher in `scripts/run_fortran_sp_ace.py` (see `scripts/README.md`).

## CI

GitHub Actions builds this target on every push (`.github/workflows/ci.yml`).

## Upstream

Original SVN layout and history:  
http://svn.ari.uni-heidelberg.de/svn/gavo/hdinputs/sp_ace/src/
