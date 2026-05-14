# Scripts (Stages 1–2)

All scripts expect the **repository root** on `PYTHONPATH` (they insert it automatically).

## `run_fortran_sp_ace.py` — Fortran driver ([issue #1](https://github.com/rghosh12/SP_ACE-Spectroscopy-Code/issues/1))

1. Build Fortran: `make -C Fortran-95_code` (see `Fortran-95_code/README.md`).
2. Point **`SP_ACE_EXE`** at the `space` / `space.exe` binary.
3. Run with an existing parameter file, or generate a minimal `space.par` and run.

```powershell
$env:SP_ACE_EXE = "C:\path\to\Fortran-95_code\space.exe"
python scripts/run_fortran_sp_ace.py C:\work\myrun.par
```

Generate `space.par` and run (paths must match your machine; `GCOGlib` must contain `linelist.dat`, `space_6degpoly.dat`, and grid files):

```powershell
python scripts/run_fortran_sp_ace.py --write-par C:\work\myrun.par `
  --obs-sp C:\path\to\SP_ACE-Spectroscopy-Code\fast_spectrum.txt `
  --gcoglib C:\path\to\GCOG `
  --fwhm 2.0 --wave-lims 4500 6800
```

**`--dry-run`** — only write the parameter file (or validate paths); do not invoke the binary (useful when the executable is not built yet).

Stdout/stderr and the process exit code are those of the Fortran process (`fortran_tools.launcher.run_sp_ace`).

## `plot_sp_ace_diagnostics.py` / `plot_sp_ace_results.py` — plots ([issue #2](https://github.com/rghosh12/SP_ACE-Spectroscopy-Code/issues/2))

After a successful SP_Ace run, plot observed vs model and residuals from `*_model.dat`:

```bash
pip install -r requirements.txt -r requirements-dev.txt
python scripts/plot_sp_ace_diagnostics.py path\to\myrun_model.dat -o diagnostic.png
```

`plot_sp_ace_results.py` is the same workflow under the name used in the roadmap issue.

## `print_sp_ace_summary.py` — text summary ([issue #2](https://github.com/rghosh12/SP_ACE-Spectroscopy-Code/issues/2))

Print parsed header/value tokens from `*_TGM_ABD.dat` and basic stats from `*_model.dat`:

```bash
python scripts/print_sp_ace_summary.py path\to\myrun.par
```

(Expects `myrun_TGM_ABD.dat` and `myrun_model.dat` next to the same stem as `myrun.par`.)

## Output column layouts

See `docs/OUTPUT_FORMATS.md` (derived from `write_res.f95` / `write_res.py`).
