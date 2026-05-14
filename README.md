# SP_ACE-Spectroscopy-Code

SP_Ace (Stellar Parameters And Chemical abundances Estimator) derives stellar
parameters — effective temperature (Teff), surface gravity (log g), and
metallicity ([M/H]) — as well as element abundances from optical stellar
spectra, assuming Local Thermodynamic Equilibrium (LTE) and 1D stellar
atmosphere models.

**Original author:** Corrado Boeche (2016–2020)  
**License:** GNU General Public License v3 — see `LICENSE` for the full text.  
**Original source:** http://svn.ari.uni-heidelberg.de/svn/gavo/hdinputs/sp_ace/src/

---

## Roadmap (Fortran as reference; Python parity first)

The authors kept the production code in **Fortran** for good reasons. A Python rewrite is only trustworthy if **GCOG I/O, polynomial basis, coefficient layout, masks, damping/profile assumptions, and abundances** match the Fortran exactly; otherwise results can look plausible but be wrong (e.g. metallicity).

Suggested order of work:

| Stage | Goal |
|--------|------|
| **1** | Thin **Python wrapper** (or small UI) around the **Fortran executable**: build `Fortran-95_code/`, fixed `space.par`, capture stdout and output files. |
| **2** | **Parse SP_Ace outputs** and diagnostic plots (residuals, selected lines, convergence). |
| **3** | **Python GCOG reader** (`read_GCOG.f95` / `read_GCOG.py`) tested **line-by-line** against Fortran-loaded in-memory arrays or logged dumps. |
| **4** | **Python model spectrum** (`make_model.f95`, `func_poly.f95`, continuum pieces) tested against Fortran models on identical inputs. |
| **5** | **Optimizer in Python** (e.g. `scipy.optimize.least_squares`) only **after** stages 3–4 match Fortran χ² and spectra. |

The current `ConvertedPython_code/` tree is a **line-by-line port** intended for parity work—not a drop-in replacement until validated against Fortran on **known standard stars**.

---

## Repository structure

```
SP_ACE-Spectroscopy-Code/
├── LICENSE
├── requirements.txt          # numpy, scipy (ConvertedPython_code)
├── requirements-dev.txt      # pytest, matplotlib (CI + tooling)
├── fortran_tools/            # Python wrapper, parsers, diagnostics (no physics)
├── scripts/                  # CLI entry points
├── tests/                    # pytest (GCOG paths, parsers, launcher)
├── parity/                   # Parity workflow notes vs Fortran
├── benchmarks/               # Placeholder for reference-star regression data
├── .github/workflows/ci.yml
├── Fortran-95_code/          # Fortran 95 + Makefile → ``space`` binary
│   ├── Makefile
│   ├── space.f95             # PROGRAM space (main)
│   └── …
│
└── ConvertedPython_code/     # Python 3 port (must match Fortran numerically)
    └── …
```

## Fortran executable (reference)

Build the official program (``PROGRAM space``). The Makefile forces **gfortran**
(GNU make’s default ``FC`` is legacy ``f77``, which breaks this project) and
compiles sources in **module dependency order**.

```bash
make -C Fortran-95_code
```

On Linux, install LAPACK/BLAS before linking (``uncertains2`` calls ``dgelss``):

```bash
sudo apt-get install liblapack-dev libblas-dev
```

The Makefile appends ``-llapack -lblas`` to the link line. On Windows/MSYS2, use a
build that provides those libraries (e.g. OpenBLAS) and adjust ``LIBS`` if needed.

On Windows, use [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries) / MinGW and ``mingw32-make`` or ``make`` from MSYS2, then set:

```bash
set SP_ACE_EXE=C:\path\to\Fortran-95_code\space.exe
```

## Python launcher (Stage 1 — I/O only)

Generate a minimal ``space.par`` and run the Fortran binary (requires ``SP_ACE_EXE`` and a real GCOG tree):

```bash
pip install -r requirements.txt
python scripts/run_fortran_sp_ace.py --write-par work\run.par ^
  --obs-sp fast_spectrum.txt --gcoglib C:\path\to\GCOG ^
  --fwhm 2.0 --wave-lims 4500 6800
```

Or run with an existing parameter file:

```bash
python scripts/run_fortran_sp_ace.py path\to\space.par
```

## Output parsing and plots (Stage 2)

After a successful run, SP_Ace writes ``<par_stem>_model.dat``, ``*_TGM_ABD.dat``, ``*_ew_meas.dat``. Load/plot from Python:

```bash
pip install -r requirements.txt -r requirements-dev.txt
python scripts/plot_sp_ace_diagnostics.py path\to\run_model.dat -o diagnostic.png
```

Use ``fortran_tools.parse_outputs`` in your own scripts for structured access.

## Parity and benchmarks

- **``parity/README.md``** — how to compare GCOG rows and model outputs with Fortran.
- **``benchmarks/README.md``** — how to add reference-star cases when data and libraries are available.

## Tests (CI)

```bash
pip install -r requirements.txt -r requirements-dev.txt
python -m pytest -q
```

A **SciPy-based replacement for the LM loop** is intentionally **not** wired in; see ``fortran_tools/optimizer_placeholder.py`` and roadmap Stage 5.

---

## Running the Python port (ConvertedPython_code)

```bash
pip install -r requirements.txt
cd ConvertedPython_code
python space.py space.par
```

See `ConvertedPython_code/README.md` for full usage details.
