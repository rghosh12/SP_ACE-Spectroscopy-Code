# SP_ACE-Spectroscopy-Code

SP_Ace (Stellar Parameters And Chemical abundances Estimator) derives stellar
parameters — effective temperature (Teff), surface gravity (log g), and
metallicity ([M/H]) — as well as element abundances from optical stellar
spectra, assuming Local Thermodynamic Equilibrium (LTE) and 1D stellar
atmosphere models.

**Original author:** Corrado Boeche (2016–2020)  
**License:** GNU General Public License v3  
**Original source:** http://svn.ari.uni-heidelberg.de/svn/gavo/hdinputs/sp_ace/src/

---

## Repository structure

```
SP_ACE-Spectroscopy-Code/
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
cd ConvertedPython_code
pip install numpy scipy
python space.py space.par
```

See `ConvertedPython_code/README.md` for full usage details.
