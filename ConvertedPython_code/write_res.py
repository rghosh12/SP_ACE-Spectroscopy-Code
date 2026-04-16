"""
write_res.py — SP_Ace result writer
Translated from write_res.f95 (Corrado Boeche, 2016, GPLv3)
"""

import numpy as np
import share
import data_lib
import space_pars
from num_type import DP
import uncertains2


def write_res(flag, chisq, conv):
    """
    Write stellar parameters, element abundances, equivalent widths, and
    the model spectrum to output files.

    Parameters
    ----------
    flag  : bool  — True if uncertainty estimation succeeded
    chisq : float — reduced chi-squared of the best fit
    conv  : int   — convergence flag
    """
    TGM = share.TGM
    ABD = share.ABD
    ele2meas = share.ele2meas
    ele2write = share.ele2write
    up_TGM = uncertains2.up_TGM
    lo_TGM = uncertains2.lo_TGM
    up_ABD = share.up_ABD
    lo_ABD = share.lo_ABD
    ABD_mask = share.ABD_mask
    select_ll_mask = share.select_ll_mask.copy()
    dim_rej = share.dim_rej
    wave_rej = share.wave_rej
    rad_rej = share.rad_rej
    sigma = TGM[3]
    null_val = space_pars.null_val
    flag_alpha = space_pars.flag_alpha

    # Remove rejected lines from selection mask before counting
    for i in range(dim_rej):
        rad = max(3.0 * sigma, rad_rej[i])
        mask = (share.wave_ll >= (wave_rej[i] - rad_rej[i])) & \
               (share.wave_ll <= (wave_rej[i] + rad_rej[i]))
        select_ll_mask[mask] = False

    null_str = {"NaN": "  NaN", "null": " null", "-9.99": "-9.99"}.get(null_val, " null")

    # ── Format TGM values ──────────────────────────────────────────────────
    def fmt_teff(v):
        return f"{int(round(v)):5d}"

    def fmt_logg(v):
        return f"{v:5.2f}"

    def fmt_met(v):
        return f"{v:5.2f}"

    def fmt_fwhm(v):
        return f"{v:6.2f}"

    def fmt_rv(v):
        return f"{v:7.1f}"

    def fmt_abd(v):
        return f"{v:5.2f}"

    Teff_s = fmt_teff(TGM[0])
    logg_s = fmt_logg(TGM[1])
    met_s  = fmt_met(TGM[2])
    fwhm_s = fmt_fwhm(abs(TGM[3]) * 2.35)
    rv_s   = fmt_rv(TGM[4])

    if flag:
        T_lo = f"{int(round(lo_TGM[0])):5d}" if lo_TGM[0] >= data_lib.temp_gridL[0] else null_str
        T_hi = f"{int(round(up_TGM[0])):5d}" if up_TGM[0] <= data_lib.temp_gridL[-1] else null_str
        G_lo = f"{lo_TGM[1]:5.2f}"  if lo_TGM[1] >= data_lib.logg_gridL[0] else null_str
        G_hi = f"{up_TGM[1]:5.2f}"  if up_TGM[1] <= data_lib.logg_gridL[-1] else null_str
        M_lo = f"{lo_TGM[2]:5.2f}"  if lo_TGM[2] >= data_lib.met_gridL[0]  else null_str
        M_hi = f"{up_TGM[2]:5.2f}"  if up_TGM[2] <= data_lib.met_gridL[-1]  else null_str
    else:
        T_lo = T_hi = G_lo = G_hi = M_lo = M_hi = null_str

    # ABD formatted values
    n_abd = len(ABD) if ABD is not None else 0
    ABDc   = [null_str] * n_abd
    loABDc = [null_str] * n_abd
    hiABDc = [null_str] * n_abd

    if ABD is not None and ABD_mask is not None:
        for j in range(n_abd):
            if -0.5 <= ABD[j] <= 0.7 and ABD_mask[j]:
                ABDc[j]   = fmt_abd(ABD[j]    + TGM[2])
                loABDc[j] = fmt_abd(lo_ABD[j] + TGM[2]) if lo_ABD[j] >= -0.5 else null_str
                hiABDc[j] = fmt_abd(up_ABD[j] + TGM[2]) if up_ABD[j] <=  0.7 else null_str

    # ── Constraints for reduced chi-sq ────────────────────────────────────
    n_tgm_free = int(np.sum(share.TGM_mask)) if share.TGM_mask is not None else 3
    n_abd_free = int(np.sum(ABD_mask)) if ABD_mask is not None else 0
    constraints = n_tgm_free + n_abd_free
    n_weig = share.n_weig if share.n_weig else 1
    chi_red = chisq / max(n_weig - constraints, 1)

    # ── Build header and result line ──────────────────────────────────────
    hdr = ""
    line = ""

    hdr  += "  conv";  line += f"   {int(conv):3d}"
    hdr  += "        RV";  line += f"   {rv_s}"
    hdr  += "    FWHM";    line += f"  {fwhm_s}"
    hdr  += "     S/N";    line += f"   {share.sn:5.1f}"
    hdr  += "    chisq";   line += f"   {chi_red:6.2f}"

    hdr  += "  Teff   T_l   T_h"
    line += f" {Teff_s} {T_lo} {T_hi}"
    hdr  += "   logg   L_l   L_h"
    line += f"  {logg_s} {G_lo} {G_hi}"
    hdr  += "    MH  MH_l  MH_h"
    line += f" {met_s} {M_lo} {M_hi}"

    def _write_element(ele_int, j_idx):
        nonlocal hdr, line
        sym = data_lib.ELE_symb[int(ele_int) - 1]
        tag = f"{sym}_l  {sym}_h {sym}_N"
        mask_el = (share.ele_ll.astype(int) == int(ele2meas[j_idx])) & select_ll_mask
        N_lin = int(np.count_nonzero(mask_el))
        hdr  += f"   {sym}  {tag}"
        line += f" {ABDc[j_idx]:6s} {loABDc[j_idx]:6s} {hiABDc[j_idx]:6s} {N_lin:4d}"

    if flag_alpha:
        # Write "metals" (Fe proxy)
        mask_fe = (share.ele_ll.astype(int) == 26) & select_ll_mask
        sym93 = data_lib.ELE_symb[92]
        sym94 = data_lib.ELE_symb[93]
        N_fe = int(np.count_nonzero(mask_fe))
        j_fe = int(np.argmin(np.abs(ele2meas - 26)))
        hdr  += f"   {sym93}  {sym93}_l  {sym93}_h {sym93}_N"
        line += (f" {ABDc[j_fe]:6s} {loABDc[j_fe]:6s} {hiABDc[j_fe]:6s} {N_fe:4d}"
                 if N_fe > 0 else f" {null_str:5s} {null_str:6s} {null_str:6s} {0:4d}")
        # Write alpha
        mask_al = (share.ele_ll.astype(int) == 94) & select_ll_mask
        N_al = int(np.count_nonzero(mask_al))
        j_al = int(np.argmin(np.abs(ele2meas - 94)))
        hdr  += f"   {sym94}  {sym94}_l  {sym94}_h {sym94}_N"
        line += (f" {ABDc[j_al]:6s} {loABDc[j_al]:6s} {hiABDc[j_al]:6s} {N_al:4d}"
                 if N_al > 0 else f" {null_str:6s} {null_str:6s} {null_str:6s} {0:4d}")
    else:
        for k, ew_k in enumerate(ele2write):
            if np.any(ele2meas == ew_k):
                j_idx = int(np.argmin(np.abs(ele2meas - ew_k)))
                _write_element(ew_k, j_idx)
            else:
                sym = data_lib.ELE_symb[int(ew_k) - 1]
                tag = f"{sym}_l  {sym}_h {sym}_N"
                hdr  += f"   {sym}  {tag}"
                line += f" {null_str:6s} {null_str:6s} {null_str:6s} {0:4d}"

    # ── Write _TGM_ABD.dat ────────────────────────────────────────────────
    base = share.space_params_file
    if base.endswith(".par"):
        base = base[:-4]

    with open(base + "_TGM_ABD.dat", "w") as f:
        f.write(hdr.strip() + "\n")
        f.write(line.strip() + "\n")

    # ── Write _ew_meas.dat ────────────────────────────────────────────────
    with open(base + "_ew_meas.dat", "w") as f:
        for i in range(share.dim_ll):
            if select_ll_mask[i]:
                f.write(f"{share.wave_ll[i]:8.3f}  {share.ele_ll[i]:5.1f}  "
                        f"{share.ew[i] * 1000.:6.1f}\n")

    # ── Write _model.dat ──────────────────────────────────────────────────
    with open(base + "_model.dat", "w") as f:
        dimsp = share.dimsp
        for ii in range(dimsp):
            f.write(f"{share.w_sp[ii]:9.3f} {share.f_sp[ii]:8.5f} "
                    f"{share.f_sp_norm[ii]:8.5f} {share.f_model[ii]:8.5f} "
                    f"{share.cont[ii]:8.5f} {share.weights[ii]:4.2f} "
                    f"{int(round(share.sn_var[ii])):4d}\n")
