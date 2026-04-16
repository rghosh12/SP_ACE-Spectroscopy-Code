"""
share.py — SP_Ace shared global variables
Translated from share.f95 (Corrado Boeche, 2016, GPLv3)

All Fortran module-level ALLOCATABLE arrays are initialised to None;
they are allocated (assigned NumPy arrays) by read_sp_ll.read_files()
and other routines before use.
"""

import numpy as np
from num_type import DP, I4B, I2B

# ── File paths ────────────────────────────────────────────────────────────────
obs_sp_file: str = ""
space_params_file: str = ""

# ── Spectrum & model arrays (allocated dynamically) ───────────────────────────
w_sp = None          # wavelength grid of observed spectrum
f_sp = None          # observed flux
f_sp_norm = None     # normalised observed flux
f_model = None       # model spectrum
cont0 = None         # initial continuum
cont = None          # continuum
weights = None       # pixel weights
f_discrep = None     # flux discrepancy

# ── Line-list arrays ──────────────────────────────────────────────────────────
wave_ll = None       # line-list wavelengths
ele_ll = None        # atomic number of each line
ew = None            # equivalent widths
disp = None          # spectral dispersion per pixel
Ex_inf = None        # excitation potential

# ── Matrix for GCOG polynomial coefficients ───────────────────────────────────
X_abd = None         # shape (70, dim_ele) — abundance transform matrix

# ── Integer arrays ────────────────────────────────────────────────────────────
wave_center_ll = None   # central pixel index for each line
flag_lines = None       # line type flags (0=normal,1=H,2=Na)

# ── Rejection arrays ──────────────────────────────────────────────────────────
wave_rej = None
rad_rej = None
sn_var = None
sig_noise = None

# ── Polynomial coefficient arrays ─────────────────────────────────────────────
coeff_4deg = None        # shape (70, dim_ll)
coeff_4deg_quick = None  # shape (84, dim_ll)

# ── Scalar dimensions ─────────────────────────────────────────────────────────
dimsp: int = 0
n_weig: int = 0
dim_ll: int = 0
dim_rej: int = 0
n_ele_symb: int = 0
dim_ele: int = 0
dim_ele_dy: int = 0

# ── Element index arrays ──────────────────────────────────────────────────────
ele2meas = None   # atomic numbers of elements to measure
ele2write = None  # atomic numbers of elements to write

# ── Pixel radius array for continuum normalisation ────────────────────────────
rad_pix = None

# ── Stellar parameter vector [Teff, logg, [M/H], sigma, RV, norm] ─────────────
TGM = np.zeros(6, dtype=DP)

# ── Abundance arrays ──────────────────────────────────────────────────────────
ABD = None
up_ABD = None
lo_ABD = None
ABD_old = None
residABD = None

# ── Scalar observables ────────────────────────────────────────────────────────
sn: float = 0.0       # signal-to-noise ratio
sigma: float = 0.4    # line width (Gaussian sigma, Å)
fwhm: float = 0.0     # FWHM (Å)
norm_rad: float = 30. # normalisation radius (Å)

# ── Parameter space bounds (user-supplied) ────────────────────────────────────
temp_infS: float = 0.0
temp_supS: float = 0.0
logg_infS: float = 0.0
logg_supS: float = 0.0
met_infS: float = 0.0
met_supS: float = 0.0

# ── Closest grid point to current TGM estimate ───────────────────────────────
TGM_prox = np.zeros(3, dtype=DP)

# ── Boolean masks ─────────────────────────────────────────────────────────────
select_ll_mask = None   # which lines are selected for measurement
cosmic_mask = None      # pixels flagged as cosmic rays
ABD_mask = None         # which abundances are being fitted
write_ABD_mask = None   # which abundances to write to output
alpha_mask = None       # alpha-element mask
TGM_mask = np.array([True, True, True, True, True, True])  # 6-element mask
