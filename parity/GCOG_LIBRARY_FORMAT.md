# GCOG library file layout (parity notes — [issue #3](https://github.com/rghosh12/SP_ACE-Spectroscopy-Code/issues/3))

This documents behaviour encoded in **`Fortran-95_code/read_GCOG.f95`** and **`ConvertedPython_code/read_GCOG.py`**. Any Python parity test must match these conventions **exactly**.

## Quick-load file: `space_6degpoly.dat`

- Path: `{GCOGlib}/space_6degpoly.dat` (see `gcog_quick_poly_path` in `fortran_tools/gcog_filenames.py`).
- One text line per **selected** line in the line list (same mask as `read_sp_ll` / `llist_mask`).
- Each data line must contain **at least 84** whitespace-separated floats (6th-degree polynomial coefficients for the quick model).

## 4th-degree grid file per (Teff, log g, [M/H])

- Filename pattern (must match Fortran string logic — covered by `tests/test_gcog_filenames.py`):

  `{met_label}{temp_label}{logg_label}`

  Example: `ap00-t5750-g45-4degpoly.dat`

- Each selected line list entry has **one row** with **at least 70** floats (4th-degree coefficients).

## In-memory targets

- `load_GCOG_4deg` fills `share.coeff_4deg[:, i]` for selected lines `i`.
- `load_GCOG_4deg_quick` fills `share.coeff_4deg_quick[:, i]`.

## Regression / golden tests (next step)

When you can ship a **tiny snippet** of a real GCOG tree (copyright permitting), add:

1. A Fortran run or one-off program that dumps one row of `coeff_4deg` for a fixed grid point.
2. A pytest that runs `read_GCOG.load_GCOG_4deg` under the Python port with the same inputs and compares arrays at **machine epsilon** or a documented tolerance.

Until then, structural tests (`test_gcog_filenames.py`, `test_gcog_line_tokens.py`) only guard **naming and row width**, not physics.
