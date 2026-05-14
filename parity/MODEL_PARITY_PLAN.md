# Model spectrum parity plan ([issue #4](https://github.com/rghosh12/SP_ACE-Spectroscopy-Code/issues/4))

Goal: prove **`ConvertedPython_code/make_model.py`**, **`func_poly.py`**, and continuum helpers reproduce **`make_model.f95`**, **`func_poly.f95`**, **`fit_cont.f95`** for identical inputs (same `TGM`, `ABD`, masks, wavelength grid, FWHM / σ).

## Recommended sequence

1. **Freeze inputs** — save `w_sp`, `f_sp`, masks, `coeff_*`, and `TGM` / `ABD` vectors from one Fortran run (small wavelength window).
2. **Dump Fortran model** — write `f_model` (and optionally `cont`) to a NumPy `.npz` or ASCII from Fortran.
3. **Python-only call** — after the port loads the same shared state, call the Python `make_model_*` entry points and compare arrays.
4. **Tolerance** — document `rtol`/`atol` (e.g. `1e-10` relative on flux) and any known ordering differences.

## Voigt / broadening / continuum

Document in the golden-test PR: which `func_poly` entry points are exercised, how `sigma` / FWHM enter `make_model`, and how `fit_cont` interacts with the normalisation loop in `space.f95` / `space.py`.

This repository does **not** yet include bundled golden vectors (GCOG + spectrum IP). Track downloads under `benchmarks/` when rights allow.
