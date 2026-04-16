"""
read_sp_ll.py — SP_Ace spectrum & line-list reader
Translated from read_sp_ll.f95 (Corrado Boeche, 2016, GPLv3)
"""

import numpy as np
import share
import error
import data_lib
import space_pars
from num_type import DP, I4B

# ── Module-level masks (allocated by read_files) ───────────────────────────────
wave_mask = None
llist_mask = None
llist_mask_w = None
llist_mask_ll = None

add_pix_sp: int = 0
add_pix_sp_l: int = 0
add_pix_sp_u: int = 0

MAX_SIZE_SP = 128_000


def read_files(name1, name2, name4):
    """
    Read the observed spectrum (*name1*), absorption-line list (*name2*),
    and optionally a user rejection list (*name4*).  Populates all
    relevant arrays in the *share* module.

    Parameters
    ----------
    name1 : str — observed-spectrum file
    name2 : str — line-list file
    name4 : str — user rejection-list file (may be empty if flag_rej=False)
    """
    global wave_mask, llist_mask, llist_mask_w, llist_mask_ll
    global add_pix_sp, add_pix_sp_l, add_pix_sp_u

    # ── Read spectrum ──────────────────────────────────────────────────────
    print(f"spectrum: {name1}")
    try:
        raw = np.loadtxt(name1, max_rows=MAX_SIZE_SP)
    except OSError:
        error.error_msg(9, "I cannot open the spectrum! (maybe wrong file?)")

    if raw.ndim == 1:
        raw = raw.reshape(1, -1)

    wave_raw = raw[:, 0]
    flux_raw = raw[:, 1]

    # Mask out the telluric band 6860–8400 Å
    tel = (wave_raw >= 6860.) & (wave_raw <= 8400.)
    wave_raw = wave_raw[~tel]
    flux_raw = flux_raw[~tel]

    # Keep only pixels inside the user-supplied wavelength intervals
    keep = np.zeros(len(wave_raw), dtype=bool)
    for i in range(space_pars.N_w_int):
        lo = space_pars.w_inf[i]
        hi = space_pars.w_sup[i]
        keep |= (wave_raw >= lo) & (wave_raw <= hi)
        # Add 4 Å padding on each side (set flux=1)
        pad_lo = (wave_raw > (lo - 4.0)) & (wave_raw < lo)
        pad_hi = (wave_raw > hi) & (wave_raw < (hi + 4.0))
        keep |= pad_lo | pad_hi

    wave = wave_raw[keep]
    flux = flux_raw[keep]
    # Restore flux=1 in padding regions
    for i in range(space_pars.N_w_int):
        lo = space_pars.w_inf[i]
        hi = space_pars.w_sup[i]
        pad = ((wave > (lo - 4.0)) & (wave < lo)) | ((wave > hi) & (wave < (hi + 4.0)))
        flux[pad] = 1.0

    sp_dim_trim = len(wave)
    if sp_dim_trim < 1:
        error.stop_msg("No readable spectrum or wrong wave_lims.")

    add_pix_sp_l = max(5, int(round(4.0 / (wave[1] - wave[0]))))
    add_pix_sp_u = max(5, int(round(4.0 / (wave[-1] - wave[-2]))))
    add_pix_sp = max(add_pix_sp_l, add_pix_sp_u, 5)

    # ── Read line list ─────────────────────────────────────────────────────
    try:
        ll_raw = np.loadtxt(name2)
    except OSError:
        error.stop_msg("I cannot open the line list")

    if ll_raw.ndim == 1:
        ll_raw = ll_raw.reshape(1, -1)

    w_ll  = ll_raw[:, 0]
    e_ll  = ll_raw[:, 1]
    ex_ll = ll_raw[:, 2] if ll_raw.shape[1] > 2 else np.zeros(len(w_ll))

    # ── Determine wavelength bounds from both spectrum and line list ────────
    wave_lbound = max(space_pars.w_inf[0], w_ll[0], wave[0])
    wave_ubound = min(space_pars.w_sup[space_pars.N_w_int - 1], w_ll[-1])

    sp_mask = (wave > wave_lbound) & (wave < wave_ubound)
    dimsp = int(np.count_nonzero(sp_mask))
    share.dimsp = dimsp

    # ── Build the padded spectrum arrays ──────────────────────────────────
    pad = add_pix_sp
    total = dimsp + 2 * pad

    w_sp = np.empty(total, dtype=DP)
    f_sp = np.ones(total, dtype=DP)

    w_sp[pad: pad + dimsp] = wave[sp_mask]

    # Extrapolate wavelength grid at both ends
    dw_l = w_sp[pad + 1] - w_sp[pad]
    dw_u = w_sp[pad + dimsp - 1] - w_sp[pad + dimsp - 2]
    for i in range(1, pad + 1):
        w_sp[pad - i] = w_sp[pad] - dw_l * i
        w_sp[pad + dimsp + i - 1] = w_sp[pad + dimsp - 1] + dw_u * i

    f_sp[pad: pad + dimsp] = flux[sp_mask]

    share.w_sp = w_sp
    share.f_sp = f_sp
    share.f_sp_norm = np.ones(total, dtype=DP)
    share.f_model = np.ones(total, dtype=DP)
    share.cont = np.ones(total, dtype=DP)
    share.cont0 = np.ones(total, dtype=DP)
    share.weights = np.ones(total, dtype=DP)
    share.cosmic_mask = np.zeros(total, dtype=bool)
    share.sn_var = np.full(total, space_pars.sn_ratio, dtype=DP)

    # Update dimsp to the full padded size
    share.dimsp = total

    # rad_pix for continuum normalisation
    rad_pix = np.empty(total, dtype=np.int32)
    if share.norm_rad < (w_sp[-1] - w_sp[0]):
        nr = max(share.norm_rad, 5.0)
        dw = np.diff(w_sp)
        dw = np.append(dw, dw[-1])
        rad_pix = np.maximum(1, np.round(nr / dw).astype(np.int32))
    else:
        rad_pix[:] = total
    share.rad_pix = rad_pix

    # Fix rad_pix at interval boundaries when multiple intervals used
    if space_pars.N_w_int > 1:
        for i in range(1, space_pars.N_w_int):
            c = int(np.argmin(np.abs(w_sp - (space_pars.w_inf[i] - 4.0))))
            if c + 1 < total:
                rad_pix[c] = rad_pix[c + 1]

    # ── Select lines within wavelength range ───────────────────────────────
    lm_w = np.zeros(len(w_ll), dtype=bool)
    lm_l = np.zeros(len(w_ll), dtype=bool)
    for i in range(space_pars.N_w_int):
        lm_w |= (w_ll > space_pars.w_inf[i]) & (w_ll < space_pars.w_sup[i])
        lm_l |= (w_ll > wave[0]) & (w_ll < wave[-1])

    llist_mask = lm_w & lm_l
    dim_ll = int(np.count_nonzero(llist_mask))
    share.dim_ll = dim_ll

    wave_ll = w_ll[llist_mask]
    ele_ll  = e_ll[llist_mask]
    Ex_inf  = ex_ll[llist_mask]

    # Dispersion & center pixel per line
    disp = np.empty(dim_ll, dtype=DP)
    wave_center_ll = np.empty(dim_ll, dtype=np.int32)
    for i in range(dim_ll):
        c = int(np.argmin(np.abs(w_sp - wave_ll[i])))
        wave_center_ll[i] = c
        if c + 1 < total:
            disp[i] = w_sp[c + 1] - w_sp[c]
        else:
            disp[i] = w_sp[c] - w_sp[c - 1]

    # Map molecules to dummy atomic numbers
    mapping = {106.: 95., 107.: 96., 112.: 97., 114.: 98.,
               606.: 99., 607.: 100., 814.: 101.}
    for src, dst in mapping.items():
        ele_ll[np.abs(ele_ll - src) < 0.3] = dst

    # If alpha-element flag: merge Mg/Si/Ca/Ti → 94, others → Fe (26)
    if space_pars.flag_alpha:
        alpha = (np.abs(ele_ll - 12) < 0.3) | (np.abs(ele_ll - 14) < 0.3) | \
                (np.abs(ele_ll - 20) < 0.3) | (np.abs(ele_ll - 22) < 0.3)
        ele_ll[alpha] = 94.0
        other = (np.abs(ele_ll - 6) > 0.3) & (np.abs(ele_ll - 7) > 0.3) & \
                (np.abs(ele_ll - 8) > 0.3) & (np.abs(ele_ll - 94) > 0.3) & \
                (ele_ll > 94.0)
        ele_ll[other] = 26.0

    # Flag special lines
    flag_lines = np.zeros(dim_ll, dtype=np.int16)
    # H lines (Hα, Hβ — element=1)
    h_mask = (ele_ll == 1.0) & (wave_ll < 6563.)
    flag_lines[h_mask] = 1
    # NaI doublet
    na_mask = (ele_ll == 11.0) & (wave_ll > 5889.) & (wave_ll < 5896.)
    flag_lines[na_mask] = 2

    share.wave_ll = wave_ll
    share.ele_ll = ele_ll
    share.Ex_inf = Ex_inf
    share.ew = np.zeros(dim_ll, dtype=DP)
    share.disp = disp
    share.wave_center_ll = wave_center_ll
    share.flag_lines = flag_lines
    share.coeff_4deg = np.zeros((70, dim_ll), dtype=DP)
    share.coeff_4deg_quick = np.zeros((84, dim_ll), dtype=DP)

    # ── Build rejection list ───────────────────────────────────────────────
    rej_wave = np.concatenate([data_lib.w_rej_op, data_lib.w_rej_nlte,
                               data_lib.w_rej_unknown, data_lib.w_rej_bad])
    rej_rad  = np.concatenate([data_lib.r_rej_op, data_lib.r_rej_nlte,
                               data_lib.r_rej_unknown, data_lib.r_rej_bad])

    if space_pars.flag_rej and name4:
        try:
            user_rej = np.loadtxt(name4)
            if user_rej.ndim == 1:
                user_rej = user_rej.reshape(1, -1)
            rej_wave = np.append(rej_wave, user_rej[:, 0])
            rej_rad  = np.append(rej_rad,  user_rej[:, 1])
        except OSError:
            error.stop_msg("I cannot open the reject line list")

    # Sort rejection list
    order = np.argsort(rej_wave)
    rej_wave = rej_wave[order]
    rej_rad  = rej_rad[order]

    # Select rejection entries that overlap the spectrum
    rej_mask = np.zeros(len(rej_wave), dtype=bool)
    for i in range(space_pars.N_w_int):
        lo = max(space_pars.w_inf[i], wave[0])
        hi = min(space_pars.w_sup[i], wave[-1])
        rej_mask |= ((rej_wave + rej_rad) >= lo) & ((rej_wave - rej_rad) <= hi)

    share.wave_rej = rej_wave[rej_mask]
    share.rad_rej  = rej_rad[rej_mask]
    share.dim_rej  = int(np.count_nonzero(rej_mask))


def _sorting(x, x1):
    """Sort both arrays x and x1 by the values in x (ascending)."""
    order = np.argsort(x)
    return x[order].copy(), x1[order].copy()
