"""
fit_cont.py — SP_Ace continuum fitting
Translated from fit_cont.f95 (Corrado Boeche, 2016, GPLv3)
"""

import numpy as np
from num_type import DP
from stats import avg_std_


def fit_cont(f_sp, f_model, weights, rad_pix):
    """
    Estimate the pseudo-continuum correction *cont* for each pixel.

    For every pixel i that has weight > 1e-3 a local window of ±rad_pix[i]
    pixels is used to compute:
      - the average flux and residual (f_sp - f_model)
      - a 2σ lower clip and 4σ upper clip
      - the continuum offset to apply

    If rad_pix[0] >= dimsp the same scalar correction is applied to the
    whole spectrum.

    Parameters
    ----------
    f_sp    : ndarray, shape (N,) — observed flux
    f_model : ndarray, shape (N,) — model flux
    weights : ndarray, shape (N,) — pixel weights
    rad_pix : ndarray, shape (N,) — half-window radius in pixels

    Returns
    -------
    cont : ndarray, shape (N,)
    """
    f_sp = np.asarray(f_sp, dtype=DP)
    f_model = np.asarray(f_model, dtype=DP)
    weights = np.asarray(weights, dtype=DP)
    rad_pix = np.asarray(rad_pix, dtype=np.int32)

    dimsp = len(f_sp)
    cont = np.ones(dimsp, dtype=DP)
    resid = f_sp - f_model

    if rad_pix[0] < dimsp:
        dimsp_loop = dimsp
    else:
        dimsp_loop = 1

    for i in range(dimsp_loop):
        weig = weights[i] if dimsp_loop > 1 else 1.0

        if weig > 1e-3:
            mask = np.zeros(dimsp, dtype=bool)
            iinf = max(0, i - int(rad_pix[i]))
            isup = min(dimsp, i + int(rad_pix[i]) + 1)

            mask[iinf:isup] = weights[iinf:isup] > 1e-3

            if np.any(mask[iinf:isup]):
                avg_f, var, sig = avg_std_(f_sp[iinf:isup],
                                           mask[iinf:isup])
                avg_r, var, sig_r = avg_std_(resid[iinf:isup],
                                             mask[iinf:isup])
                avg_mm, var, sig_mm = avg_std_(f_model[iinf:isup],
                                               mask[iinf:isup])
                # 2σ lower clip
                mask[iinf:isup] &= f_sp[iinf:isup] > (avg_f - 2.0 * sig_r)
                # 4σ upper clip
                mask[iinf:isup] &= f_sp[iinf:isup] < (1.0 + 4.0 * sig_r)

            if np.count_nonzero(mask) > 1:
                avg_f, var, sig = avg_std_(f_sp[iinf:isup], mask[iinf:isup])
                avg_m, var, sig = avg_std_(f_model[iinf:isup], mask[iinf:isup])
                avg_r, var, sig_r = avg_std_(resid[iinf:isup], mask[iinf:isup])

                # also need avg_mm for the full local region
                avg_mm_full, _, _ = avg_std_(
                    f_model[iinf:isup],
                    weights[iinf:isup] > 1e-3
                )
            else:
                avg_f = 1.0
                avg_m = 1.0
                avg_mm_full = 1.0
                sig_r = 0.0

            set_val = 4.0 * (1.0 - np.exp((1.0 - avg_mm_full) ** 3))
            cont[i] = 1.0 + (avg_f - avg_m - set_val)
        else:
            cont[i] = 1.0

    if rad_pix[0] >= dimsp:
        cont[:] = cont[0]

    return cont
