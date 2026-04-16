"""
make_model.py — SP_Ace spectrum model construction
Translated from make_model.f95 (Corrado Boeche, 2016, GPLv3)
"""

import numpy as np
import share
import data_lib
from num_type import DP
from func_poly import (poly6_transform, ew_poly6,
                       poly4_transform, ew_poly4, voigt_logg)

_SQRT2LN2 = np.float64(1.17741002252)  # sqrt(2*ln2)


def _clamp_params(params):
    """Apply physical parameter limits in-place."""
    params[0] = np.clip(params[0], 3600., 7400.)
    params[1] = np.clip(params[1], 0.1, 5.5)
    params[2] = np.clip(params[2], -2.6, 0.5)
    params[3] = np.clip(abs(params[3]), 0.01, 5.0)


def gammaL(i_line, x0, x1, x2, x3):
    """
    Compute the Lorentzian broadening parameter γ_L.

    Parameters
    ----------
    i_line : int   — 0-based line index
    x0     : float — equivalent width
    x1     : float — Gaussian FWHM (gamG)
    x2     : float — log g
    x3     : float — Teff

    Returns
    -------
    gaml : float  (>= 0.01)
    """
    x = np.array([1.0, x0, x1, x2, x3], dtype=DP)
    count = 0
    transf = np.empty(15, dtype=DP)
    for i in range(5):
        for j in range(i, 5):
            transf[count] = x[i] * x[j]
            count += 1

    # All line types use gamL_coeff in this version (see commented code)
    result = float(np.dot(transf, data_lib.gamL_coeff))
    return max(result, 0.01)


def _build_voigt_model(model, wave_ll, wave_center_ll, select_ll_mask,
                       ew_arr, gamg, c_doppler, st_pars):
    """
    Subtract Voigt absorption profiles from *model* for all selected lines.
    Updates *model* in place.
    """
    w_sp = share.w_sp
    dimsp = share.dimsp

    for i in range(share.dim_ll):
        if not select_ll_mask[i]:
            continue
        ew_val = ew_arr[i]
        w_center = int(wave_center_ll[i])
        gaml = gammaL(i, ew_val, gamg, st_pars[1], st_pars[0])
        wc = wave_ll[i] * c_doppler

        # Blue wing
        for j in range(w_center - 1, -1, -1):
            vv = voigt_logg(w_sp[j], wc, gamg, ew_val, gaml)
            if vv > 1e-4:
                model[j] -= vv
            else:
                break

        # Red wing
        for j in range(w_center, dimsp):
            vv = voigt_logg(w_sp[j], wc, gamg, ew_val, gaml)
            if vv > 1e-4:
                model[j] -= vv
            else:
                break


def make_model_TGM_quick(model, params):
    """
    Build the quick (6th-degree polynomial) spectral model.

    Parameters
    ----------
    model  : ndarray, shape (dimsp,) — modified in place, initialised to 1
    params : ndarray, shape (≥5,)   — [Teff, logg, [M/H], sigma, RV, ...]
    """
    params = np.array(params, dtype=DP)
    _clamp_params(params)

    st_pars = np.array([params[0], params[1], params[2], 1e-4], dtype=DP)
    sig = abs(params[3])
    gamg = 2.0 * sig * _SQRT2LN2
    c = 1.0 + params[4] / 299792.0

    model[:] = 1.0
    X = poly6_transform(st_pars)

    select = share.select_ll_mask
    wave_ll = share.wave_ll
    wave_center_ll = share.wave_center_ll
    coeff = share.coeff_4deg_quick
    ew_arr = share.ew

    for i in range(share.dim_ll):
        if select[i]:
            ew_arr[i] = ew_poly6(X, coeff[:, i])

    _build_voigt_model(model, wave_ll, wave_center_ll, select,
                       ew_arr, gamg, c, st_pars)


def make_model_TGM(model, params):
    """
    Build the full (4th-degree polynomial) spectral model using element
    abundances from *share.ABD*.

    Parameters
    ----------
    model  : ndarray, shape (dimsp,)
    params : ndarray, shape (≥5,)
    """
    params = np.array(params, dtype=DP)
    _clamp_params(params)

    st_pars = np.array([params[0], params[1], params[2], 1e-4], dtype=DP)
    sig = abs(params[3])
    gamg = 2.0 * sig * _SQRT2LN2
    c = 1.0 + params[4] / 299792.0

    model[:] = 1.0
    X = poly4_transform(st_pars, share.ABD, share.dim_ele, share.X_abd)

    select = share.select_ll_mask
    wave_ll = share.wave_ll
    wave_center_ll = share.wave_center_ll
    coeff = share.coeff_4deg
    ew_arr = share.ew

    for i in range(share.dim_ll):
        if select[i]:
            ew_arr[i] = ew_poly4(X, share.ele_ll[i], coeff[:, i],
                                 share.ele2meas, share.X_abd)

    _build_voigt_model(model, wave_ll, wave_center_ll, select,
                       ew_arr, gamg, c, st_pars)


def make_model_ABD(model, params):
    """
    Build a spectral model with fixed TGM and variable abundances.

    Parameters
    ----------
    model  : ndarray, shape (dimsp,)
    params : ndarray, shape (dim_ele,) — element abundances (clamped)
    """
    params = np.clip(params, -0.6, 0.8)
    TGM = share.TGM

    st_pars = np.array([TGM[0], TGM[1], TGM[2], 1e-4], dtype=DP)
    sig = abs(TGM[3])
    gamg = 2.0 * sig * _SQRT2LN2
    c = 1.0 + TGM[4] / 299792.0

    model[:] = 1.0
    X = poly4_transform(st_pars, params, share.dim_ele, share.X_abd)

    select = share.select_ll_mask
    wave_ll = share.wave_ll
    wave_center_ll = share.wave_center_ll
    coeff = share.coeff_4deg
    ew_arr = share.ew

    for i in range(share.dim_ll):
        if select[i]:
            ew_arr[i] = ew_poly4(X, share.ele_ll[i], coeff[:, i],
                                 share.ele2meas, share.X_abd)

    _build_voigt_model(model, wave_ll, wave_center_ll, select,
                       ew_arr, gamg, c, st_pars)


def make_model_ABDerr(model, parTGM, parABD):
    """
    Build model for uncertainty estimation with separate TGM and ABD vectors.

    Parameters
    ----------
    model  : ndarray, shape (dimsp,)
    parTGM : ndarray, shape (≥5,)
    parABD : ndarray, shape (dim_ele,)
    """
    parABD = np.clip(parABD, -0.6, 0.8)

    st_pars = np.array([parTGM[0], parTGM[1], parTGM[2], 1e-4], dtype=DP)
    sig = abs(parTGM[3])
    gamg = 2.0 * sig * _SQRT2LN2
    c = 1.0 + parTGM[4] / 299792.0

    model[:] = 1.0
    X = poly4_transform(st_pars, parABD, share.dim_ele, share.X_abd)

    select = share.select_ll_mask
    wave_ll = share.wave_ll
    wave_center_ll = share.wave_center_ll
    coeff = share.coeff_4deg
    ew_arr = share.ew

    for i in range(share.dim_ll):
        if select[i]:
            ew_arr[i] = ew_poly4(X, share.ele_ll[i], coeff[:, i],
                                 share.ele2meas, share.X_abd)

    _build_voigt_model(model, wave_ll, wave_center_ll, select,
                       ew_arr, gamg, c, st_pars)
