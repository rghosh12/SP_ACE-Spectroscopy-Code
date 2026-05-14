# Reference star benchmark template ([issue #7](https://github.com/rghosh12/SP_ACE-Spectroscopy-Code/issues/7))

Use one subdirectory per star when you add regression data (spectra + `space.par` + expected outputs).

## `cases/<star_id>/README.md` (suggested fields)

| Field | Example |
|-------|---------|
| Star ID / name | HD 12345 |
| Spectrum source | URL or publication DOI |
| Citation | BibTeX key or ADS bibcode |
| Literature Teff / log g / [M/H] | 5777 K, 4.44 dex, 0.0 |
| SP_Ace version / commit | Fortran `space` built from commit `…` |
| Tolerances on Teff / [M/H] | ±50 K, ±0.05 dex (justify vs pipeline) |

## Optional `expected_TGM_ABD.txt`

Paste the second line of `*_TGM_ABD.dat` (values) for scripted diff, or store scalar checks only.

## Automation

- Set environment variable **`BENCHMARK_DATA=1`** when `cases/` contains at least one complete example; then enable `pytest -m benchmark` (see `tests/test_benchmark_regression.py`).
- Do **not** commit large proprietary GCOG trees; use download scripts + `.gitignore`.
