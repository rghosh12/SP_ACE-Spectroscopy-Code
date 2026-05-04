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
├── LICENSE                   # GNU GPL version 3 (full text)
├── requirements.txt          # numpy, scipy for ConvertedPython_code
├── Fortran-95_code/          # Original Fortran 95 source files (.f95)
│   ├── data_lib.f95
│   ├── error.f95
│   ├── fit_cont.f95
│   ├── func_poly.f95
│   ├── interfaces.f95
│   ├── make_model.f95
│   ├── minimize.f95
│   ├── num_type.f95
│   ├── read_GCOG.f95
│   ├── read_sp_ll.f95
│   ├── share.f95
│   ├── space.f95
│   ├── space_pars.f95
│   ├── stats.f95
│   ├── uncertains2.f95
│   ├── utils.f95
│   └── write_res.f95
│
└── ConvertedPython_code/     # Python 3 translations of each Fortran module
    ├── num_type.py
    ├── data_lib.py
    ├── share.py
    ├── stats.py
    ├── func_poly.py
    ├── space_pars.py
    ├── error.py
    ├── read_sp_ll.py
    ├── read_GCOG.py
    ├── fit_cont.py
    ├── make_model.py
    ├── utils.py
    ├── minimize.py
    ├── uncertains2.py
    ├── write_res.py
    ├── interfaces.py
    ├── space.py              # Main entry point
    └── README.md             # Python port documentation
```

## Running the Python version

```bash
pip install -r requirements.txt
cd ConvertedPython_code
python space.py space.par
```

See `ConvertedPython_code/README.md` for full usage details.
