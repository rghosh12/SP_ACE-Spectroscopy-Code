"""
utils.py — SP_Ace utility routines
Translated from utils.f95 (Corrado Boeche, 2016, GPLv3)
"""

import numpy as np
import share
import data_lib
from num_type import DP
from stats import avg_std_

# ── Module-level masks (mirrors Fortran USE share,ONLY:ABD_mask,write_ABD_mask)
# These are references; the actual arrays live in share.
ABD_mask = None
write_ABD_mask = None


def _sync():
    """Keep local references in sync with share module."""
    global ABD_mask, write_ABD_mask
    ABD_mask = share.ABD_mask
    write_ABD_mask = share.write_ABD_mask


# ── Cosmic-ray rejection ──────────────────────────────────────────────────────

def cosmic_rej():
    """
    Flag pixels as cosmic rays when the observed flux exceeds the model
    by more than 4σ (residual), or when flux ≤ 0.01.
    Updates *share.cosmic_mask* in place.
    """
    f_sp = share.f_sp
    f_model = share.f_model
    weights = share.weights
    cosmic_mask = share.cosmic_mask

    mask = weights > 1e-3
    resid = f_sp - f_model
    _, _, sig_r = avg_std_(resid, mask)

    cosmic_mask[f_sp > (1.0 + 4.0 * sig_r)] = True
    cosmic_mask[f_sp <= 0.01] = True
    share.cosmic_mask = cosmic_mask


# ── Signal-to-noise computation ───────────────────────────────────────────────

def new_sn(sn):
    """
    Estimate the global S/N from (f_sp_norm − f_model) residuals.

    Parameters
    ----------
    sn : float — current S/N estimate (may be overwritten)

    Returns
    -------
    sn : float — updated S/N (−1 if residuals vanish)
    """
    import read_sp_ll as rsl
    add_pix = rsl.add_pix_sp
    f_sp_norm = share.f_sp_norm
    f_model = share.f_model
    weights = share.weights
    dimsp = share.dimsp

    resid = f_sp_norm - f_model
    mask1 = np.array(weights > 0.001, dtype=bool)
    mask1[:add_pix] = False
    mask1[dimsp - add_pix:] = False

    mean, _, sdev1 = avg_std_(resid, mask1)
    mask2 = np.abs(resid - mean) < (3.0 * sdev1)
    mask3 = mask1 & mask2

    if not np.any(mask3):
        return -1.0

    _, _, sdev2 = avg_std_(resid, mask3)
    if sdev2 < 1e-4:
        sdev2 = 10.0
    return 1.0 / sdev2


def find_sn_var(f_sp_norm, f_model, sn_var_out):
    """
    Compute local S/N at each pixel within a ±25-pixel window.

    Parameters
    ----------
    f_sp_norm  : ndarray, shape (N,)
    f_model    : ndarray, shape (N,)
    sn_var_out : ndarray, shape (N,) — updated in place

    Returns
    -------
    sn_var_out : ndarray
    """
    weights = share.weights
    dimsp = share.dimsp
    resid = f_sp_norm - f_model
    int_pix = 25

    for i in range(dimsp):
        if weights[i] > 1e-3:
            iinf = max(0, i - int_pix)
            isup = min(dimsp, i + int_pix + 1)
            mask = weights[iinf:isup] > 1e-3

            avg_r, _, sig_r = avg_std_(resid[iinf:isup], mask)
            clip = np.abs(resid[iinf:isup] - avg_r) <= (3 * sig_r)
            mask2 = mask & clip

            if np.count_nonzero(mask2) > 2:
                _, _, sig_r = avg_std_(resid[iinf:isup], mask2)
                if sig_r < 1e-4:
                    sig_r = 10.0
                sn_var_out[i] = 1.0 / sig_r
            else:
                sn_var_out[i] = 0.1
        else:
            sn_var_out[i] = 0.1

    return sn_var_out


# ── Grid proximity helpers ────────────────────────────────────────────────────

def find_prox(TGM_):
    """
    Return the closest stellar-grid point to *TGM_* (Teff, logg, [M/H]).

    Parameters
    ----------
    TGM_ : array_like, shape (≥3,)

    Returns
    -------
    TGM_prox : ndarray, shape (3,)
    """
    T = float(data_lib.temp_grid[np.argmin(np.abs(data_lib.temp_grid - TGM_[0]))])
    G = float(data_lib.logg_grid[np.argmin(np.abs(data_lib.logg_grid - TGM_[1]))])
    M = float(data_lib.met_grid[np.argmin(np.abs(data_lib.met_grid - TGM_[2]))])

    # Avoid the empty corner of the grid (hot + low gravity)
    if T > 5600 and G < 1.4:
        if (T - 5600) / 200.0 > (1.4 - G) / 0.4:
            G = 1.4
        else:
            T = 5600.0

    return np.array([T, G, M], dtype=DP)


def make_gridS(TGM_prox):
    """
    Build a 5-point local sub-grid centred on the closest grid node.

    Parameters
    ----------
    TGM_prox : array_like, shape (3,) — closest grid point

    Returns
    -------
    temp_gridS, logg_gridS, met_gridS : ndarrays, shape (5,)
    """
    def _subgrid(arr, val):
        pos = int(np.argmin(np.abs(arr - val)))
        inf = max(0, pos - 2)
        sup = min(len(arr), pos + 3)
        return arr[inf:sup].copy()

    return (_subgrid(data_lib.temp_gridL, TGM_prox[0]),
            _subgrid(data_lib.logg_gridL, TGM_prox[1]),
            _subgrid(data_lib.met_gridL, TGM_prox[2]))


def find_proxS(TGM_, temp_gridS, logg_gridS, met_gridS, TGM_prox):
    """
    Find proximity flags for the current TGM estimate relative to the
    local sub-grid and the large grid.

    Returns
    -------
    TGM_prox  : ndarray, shape (3,) — updated closest large-grid point
    flag_lim  : bool — True if TGM_ is on/beyond the large-grid boundary
    flag_limS : bool — True if TGM_ is outside the local sub-grid
    flag_move : bool — True if TGM_ moved more than 1 sub-grid step
    """
    TGM_prox = np.asarray(TGM_prox, dtype=DP)
    flag_lim = flag_limS = flag_move = False

    tL = data_lib.temp_gridL
    gL = data_lib.logg_gridL
    mL = data_lib.met_gridL

    Td = float(temp_gridS[np.argmin(np.abs(temp_gridS - TGM_[0]))])
    Gd = float(logg_gridS[np.argmin(np.abs(logg_gridS - TGM_[1]))])
    Md = float(met_gridS[np.argmin(np.abs(met_gridS - TGM_[2]))])

    TGM_prox[:] = find_prox([Td, Gd, Md])

    # Check large-grid boundary
    if ((Td == tL[0] and TGM_[0] <= tL[0]) or
            (Td == tL[-1] and TGM_[0] >= tL[-1]) or
            (Gd == gL[0] and TGM_[1] <= gL[0]) or
            (Gd == gL[-1] and TGM_[1] >= gL[-1]) or
            (Md == mL[0] and TGM_[2] <= mL[0]) or
            (Md == mL[-1] and TGM_[2] >= mL[-1])):
        flag_lim = True
        return TGM_prox, flag_lim, flag_limS, flag_move

    # Check sub-grid boundary
    if (TGM_[0] <= temp_gridS[0] or TGM_[0] >= temp_gridS[-1] or
            TGM_[1] <= logg_gridS[0] or TGM_[1] >= logg_gridS[-1] or
            TGM_[2] <= met_gridS[0] or TGM_[2] >= met_gridS[-1]):
        flag_limS = True
        return TGM_prox, flag_lim, flag_limS, flag_move

    # Check movement > 1 step
    if (TGM_[0] <= temp_gridS[1] or TGM_[0] >= temp_gridS[-2] or
            TGM_[1] <= logg_gridS[1] or TGM_[1] >= logg_gridS[-2] or
            TGM_[2] <= met_gridS[1] or TGM_[2] >= met_gridS[-2]):
        flag_move = True

    return TGM_prox, flag_lim, flag_limS, flag_move


# ── Abundance array allocation ────────────────────────────────────────────────

def alloc_ABD(ele_in, dim_in, n_ele_symb):
    """
    Allocate and initialise the element-abundance arrays in *share*.

    Parameters
    ----------
    ele_in      : ndarray, shape (dim_in,) — all element IDs in the line list
    dim_in      : int
    n_ele_symb  : int — total number of element symbols

    Returns
    -------
    dim_ele : int — number of distinct measured elements
    """
    import error as err

    mask_ele = np.zeros(n_ele_symb + 1, dtype=bool)
    for i in range(2, n_ele_symb + 1):
        if np.any(np.abs(ele_in - i) < 0.3):
            mask_ele[i] = True

    if not mask_ele[26]:
        err.stop_msg("no iron lines!")

    dim_ele = int(np.count_nonzero(mask_ele))
    share.dim_ele = dim_ele

    share.ABD = np.zeros(dim_ele, dtype=DP)
    share.ABD_old = np.full(dim_ele, -0.1, dtype=DP)
    share.up_ABD = np.zeros(dim_ele, dtype=DP)
    share.lo_ABD = np.zeros(dim_ele, dtype=DP)
    share.residABD = np.zeros(dim_ele, dtype=DP)
    share.ABD_mask = np.ones(dim_ele, dtype=bool)
    share.alpha_mask = np.zeros(dim_ele, dtype=bool)
    share.write_ABD_mask = np.zeros(dim_ele, dtype=bool)
    share.X_abd = np.zeros((70, dim_ele), dtype=DP)

    # Build ele2meas: Fe first, then all others
    ele2meas = np.zeros(dim_ele, dtype=np.int16)
    ele2meas[0] = 26
    j = 1
    for i in range(2, n_ele_symb + 1):
        if mask_ele[i] and i != 26:
            ele2meas[j] = i
            j += 1
    share.ele2meas = ele2meas

    # If ele2write not yet set, mirror ele2meas
    if share.ele2write is None:
        share.ele2write = ele2meas.copy()

    global ABD_mask, write_ABD_mask
    ABD_mask = share.ABD_mask
    write_ABD_mask = share.write_ABD_mask

    return dim_ele


# ── [α/Fe] computation ────────────────────────────────────────────────────────

def compute_aFe(ABD, ele2meas, ABD_mask_in):
    """
    Compute the mean [α/Fe] from C, O, Mg, Si, Ca, Ti abundances.

    Parameters
    ----------
    ABD       : ndarray, shape (N,)
    ele2meas  : ndarray, shape (N,)
    ABD_mask_in : ndarray, bool, shape (N,)

    Returns
    -------
    alpha_mean : float
    """
    alpha_elem = np.array([6, 8, 12, 14, 20, 22])
    alpha_mask = np.zeros(len(ele2meas), dtype=bool)
    for ae in alpha_elem:
        alpha_mask |= (ele2meas == ae)
    alpha_mask &= ABD_mask_in
    share.alpha_mask = alpha_mask

    if np.count_nonzero(alpha_mask) > 0:
        _, _, _ = avg_std_(ABD, alpha_mask)
        alpha_mean, _, _ = avg_std_(ABD, alpha_mask)
    else:
        alpha_mean = 0.0
    return alpha_mean


# ── Line selection & weight setting ──────────────────────────────────────────

def select_lines(fun, coeff, pars, flag):
    """
    Select spectral lines for measurement and assign pixel weights.

    Parameters
    ----------
    fun   : callable — ew_poly4 or ew_poly6
    coeff : ndarray  — polynomial coefficient array
    pars  : ndarray, shape (4,) — current stellar parameters
    flag  : bool — if True, update ABD_mask
    """
    import space_pars as sp
    import read_sp_ll as rsl
    from func_poly import poly4_transform

    wave_ll = share.wave_ll
    ele_ll = share.ele_ll
    dim_ll = share.dim_ll
    dim_rej = share.dim_rej
    wave_rej = share.wave_rej
    rad_rej = share.rad_rej
    w_sp = share.w_sp
    ew_arr = share.ew
    TGM = share.TGM
    ABD = share.ABD
    ele2meas = share.ele2meas
    sn_var = share.sn_var
    dimsp = share.dimsp

    add_pix = rsl.add_pix_sp
    c = 1.0 + TGM[4] / 299792.0
    sigma = TGM[3]

    share.select_ll_mask = np.zeros(dim_ll, dtype=bool)
    select2 = np.zeros(dim_ll, dtype=bool)
    share.weights = np.ones(dimsp, dtype=DP)

    # Compute EWs
    X = poly4_transform(np.array([TGM[0], TGM[1], TGM[2], 1e-4], dtype=DP),
                        ABD, share.dim_ele, share.X_abd)
    for i in range(dim_ll):
        ew_arr[i] = fun(X, ele_ll[i], coeff[:, i], ele2meas, share.X_abd)

    # Select lines with sufficient EW relative to local S/N
    for i in range(dim_ll):
        neigh = (np.abs(wave_ll - wave_ll[i]) < sigma)
        ew_neigh = float(np.sum(ew_arr[neigh]))
        line_pos = int(np.argmin(np.abs(wave_ll[i] - w_sp)))
        local_sn = sn_var[line_pos]

        if ew_neigh / (2.5 * sigma) > 1.0 / local_sn:
            nearby = (np.abs(ele_ll - wave_ll[i]) < 3 * sigma) & (ew_arr > 0.1 / local_sn)
            share.select_ll_mask |= nearby
            select2 |= nearby

    # Set low weights around rejected lines
    for i in range(dim_rej):
        rad = max(3.0 * sigma, rad_rej[i])
        share.select_ll_mask |= ((wave_ll >= (wave_rej[i] - rad_rej[i])) &
                                  (wave_ll <= (wave_rej[i] + rad_rej[i])))
        share.weights[np.abs(w_sp - wave_rej[i] * c) <= rad] = 1e-6

    # Padding regions
    share.weights[:add_pix] = 1e-6
    share.weights[dimsp - add_pix:] = 1e-6

    # Excluded wavelength intervals
    if sp.N_w_int > 1:
        for i in range(sp.N_w_int - 1):
            share.weights[(w_sp > sp.w_sup[i]) & (w_sp < sp.w_inf[i + 1])] = 1e-6

    # Cosmic rays
    share.weights[share.cosmic_mask] = 1e-6

    if flag:
        share.ABD_mask[:] = False
        for i in range(len(ele2meas)):
            mask_ew = (np.abs(ele_ll - ele2meas[i]) < 0.3) & select2
            if (np.count_nonzero(mask_ew) > 0 and
                    ABD[i] < 0.7 and ABD[i] > -0.5):
                share.ABD_mask[i] = True
            else:
                share.ABD_mask[i] = False
                ABD[i] = 0.0

        share.write_ABD_mask[:] = False
        for i in range(len(share.ele2write)):
            share.write_ABD_mask |= (ele2meas == share.ele2write[i])


# ── Parameter normalisation helpers ──────────────────────────────────────────

def normalize_pars(tgm, tgm_mask):
    """
    Normalise TGM: Teff/5000, logg/3, then pack by *tgm_mask*.

    Returns
    -------
    tgmx : ndarray, shape (n_free,)
    """
    tgm_local = tgm.copy().astype(DP)
    tgm_local[0] /= 5000.0
    tgm_local[1] /= 3.0
    return tgm_local[tgm_mask]


def denormalize_pars(tgmx, tgm_mask):
    """
    Unpack *tgmx* into the full TGM vector and de-normalise.

    Returns
    -------
    tgm_local : ndarray, shape (6,)
    """
    tgm_local = share.TGM.copy().astype(DP)
    tgm_local[tgm_mask] = tgmx
    if tgm_mask[0]:
        tgm_local[0] *= 5000.0
    if tgm_mask[1]:
        tgm_local[1] *= 3.0
    return tgm_local
