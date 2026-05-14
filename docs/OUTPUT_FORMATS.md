# SP_Ace ASCII output formats

Column layouts match **`ConvertedPython_code/write_res.py`** and **`Fortran-95_code/write_res.f95`**.

## `{stem}_model.dat`

One row per wavelength pixel, **7** whitespace-separated columns:

| Col | Name (conceptual) | Description |
|-----|-------------------|-------------|
| 1 | `w_sp` | Wavelength (Å) |
| 2 | `f_sp` | Raw observed flux |
| 3 | `f_sp_norm` | Normalised flux |
| 4 | `f_model` | Model flux |
| 5 | `cont` | Continuum fit |
| 6 | `weights` | Pixel weight (0–1) |
| 7 | `sn_var` | Stored as integer (`NINT` in Fortran) |

Load in Python: `fortran_tools.parse_outputs.load_model_dat`.

## `{stem}_TGM_ABD.dat`

**Two lines**: first line = header tokens, second line = value tokens (same token count when parsing succeeds). Keys include `conv`, `RV`, `FWHM`, `S/N`, `chisq`, `Teff`, metallicity and per-element columns depending on `ele2write` / `flag_alpha`.

Load: `fortran_tools.parse_outputs.read_tgm_abd_table` or `tgm_abd_as_dict`.

## `{stem}_ew_meas.dat`

Three columns per selected line: wavelength (Å), element id, EW (mÅ).

Load: `fortran_tools.parse_outputs.load_ew_meas` (`numpy` array).

## Stem rule

If the parameter file is `myrun.par`, outputs use stem `myrun` in the **same directory** as the parameter path stored in `share.space_params_file` / Fortran’s `space_params_file`.
