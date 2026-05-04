# SP_Ace — Python Port

Fortran-to-Python translation of the SP_Ace stellar spectral analysis tool.
Original code by Corrado Boeche (GPLv3).  
Source: http://svn.ari.uni-heidelberg.de/svn/gavo/hdinputs/sp_ace/src/

## File mapping

| Python file        | Fortran original   | Description |
|--------------------|--------------------|-------------|
| `num_type.py`      | `num_type.f95`     | NumPy dtype constants (DP, I4B, …) |
| `data_lib.py`      | `data_lib.f95`     | Stellar atmosphere grids & coefficients |
| `share.py`         | `share.f95`        | Global shared variables |
| `stats.py`         | `stats.f95`        | Mean / variance / std-dev routines |
| `func_poly.py`     | `func_poly.f95`    | Polynomial basis & Voigt profile |
| `space_pars.py`    | `space_pars.f95`   | Parameter-file (`space.par`) reader |
| `error.py`         | `error.f95`        | Error messages and clean exit |
| `read_sp_ll.py`    | `read_sp_ll.f95`   | Spectrum & line-list I/O |
| `read_GCOG.py`     | `read_GCOG.f95`    | GCOG library loader |
| `fit_cont.py`      | `fit_cont.f95`     | Pseudo-continuum fitting |
| `make_model.py`    | `make_model.f95`   | Spectral model construction |
| `utils.py`         | `utils.f95`        | Utility routines (S/N, grids, etc.) |
| `minimize.py`      | `minimize.f95`     | χ² functions + Levenberg-Marquardt |
| `uncertains2.py`   | `uncertains2.f95`  | Parameter uncertainty estimation |
| `write_res.py`     | `write_res.f95`    | Output writer |
| `interfaces.py`    | `interfaces.f95`   | Stub (Fortran interface declarations) |
| `space.py`         | `space.f95`        | Main program entry point |

## Dependencies

```
numpy
scipy        # for scipy.linalg.lstsq in uncertains2.py
```

Install with (from repository root):
```bash
pip install -r requirements.txt
```

Or: `pip install numpy scipy`

## Usage

```bash
python space.py space.par
```

The parameter file `space.par` follows the same keyword format as the
original Fortran version.  Required keywords:

```
obs_sp_file   /path/to/spectrum.dat
GCOGlib       /path/to/GCOG/library/
fwhm          2.0
wave_lims     4500 6800
```

## Translation notes

- Fortran `MODULE` variables → Python module-level globals in `share.py`
- `REAL(DP)` / `INTEGER(I4B)` → `np.float64` / `np.int32`
- `ALLOCATABLE` arrays → `None` until allocated by `read_files()`
- `WHERE(mask) arr=val` → `arr[mask] = val` (NumPy boolean indexing)
- `DOT_PRODUCT` → `np.dot`
- `MINLOC`/`MAXLOC` → `np.argmin`/`np.argmax`
- `PACK(arr, mask)` → `arr[mask]`
- `UNPACK(packed, mask, fill)` → reconstructed with NumPy indexing
- LAPACK `DGELSS` → `scipy.linalg.lstsq`
- The Levenberg-Marquardt solver (`minimize.py`) is a direct port of the
  MINPACK routines included in the Fortran source.
