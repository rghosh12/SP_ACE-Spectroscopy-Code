"""
space.py — SP_Ace main program
Translated from space.f95 (Corrado Boeche, 2020, GPLv3)

SP_Ace derives stellar parameters (Teff, log g, [M/H]) and element
abundances from optical stellar spectra under LTE assumptions.

Usage
-----
    python space.py space.par
"""

import sys
import numpy as np
from num_type import DP

# ── Imports (order mirrors Fortran USE statements) ────────────────────────────
import share
import space_pars as sp
import read_sp_ll
import read_GCOG
import utils
import fit_cont as fc
import make_model as mm
import minimize as mini
import uncertains2 as unc
import write_res as wr
import error as err
from func_poly import ew_poly4


def update_gridS_ll(TGM_prox, temp_gridS, logg_gridS, met_gridS):
    """
    Update the local sub-grid and re-load the GCOG coefficients for the
    grid point closest to the current TGM estimate.
    """
    import read_GCOG as rg
    new_prox = utils.find_prox(share.TGM[:3])
    temp_gridS, logg_gridS, met_gridS = utils.make_gridS(new_prox)
    rg.load_GCOG_4deg(new_prox[0], new_prox[1], new_prox[2])
    return new_prox, temp_gridS, logg_gridS, met_gridS


def main(par_file):
    # ── Read parameter file ────────────────────────────────────────────────
    share.space_params_file = par_file
    sp.read_space_pars(par_file)

    # ── Read spectrum and line list ────────────────────────────────────────
    read_sp_ll.read_files(share.obs_sp_file, sp.llist,
                          sp.llist_rej if sp.flag_rej else "")

    # ── Allocate abundance arrays ──────────────────────────────────────────
    dim_ele = utils.alloc_ABD(share.ele_ll, share.dim_ll, share.n_ele_symb)

    # ── Initialise TGM ────────────────────────────────────────────────────
    share.TGM[0] = 5000.0   # Teff (K)
    share.TGM[1] = 2.5      # log g
    share.TGM[2] = 0.0      # [M/H]
    share.TGM[3] = sp.fwhm / 2.35  # sigma (Å)
    share.TGM[4] = sp.rv_ini        # RV (km/s)
    share.TGM[5] = 1.0              # normalisation factor

    # Apply forced parameters
    if sp.TGM_force[0] > 0:
        share.TGM[0] = sp.TGM_force[0]
        share.TGM_mask[0] = False
    if sp.TGM_force[1] > 0:
        share.TGM[1] = sp.TGM_force[1]
        share.TGM_mask[1] = False

    # ── Noise model ───────────────────────────────────────────────────────
    share.sig_noise = np.ones(share.dimsp, dtype=DP) / sp.sn_ratio

    # ── Load quick (6-degree) GCOG ─────────────────────────────────────────
    read_GCOG.load_GCOG_4deg_quick()

    # ── Initialise cosmic-ray mask ────────────────────────────────────────
    share.cosmic_mask[:] = False

    # ── Initial line selection (quick model) ──────────────────────────────
    utils.select_lines(ew_poly4, share.coeff_4deg,
                       np.array([share.TGM[0], share.TGM[1],
                                 share.TGM[2], share.TGM[3]], dtype=DP),
                       flag=True)

    # ── Find closest grid point & build local sub-grid ────────────────────
    TGM_prox = utils.find_prox(share.TGM[:3])
    temp_gridS, logg_gridS, met_gridS = utils.make_gridS(TGM_prox)
    read_GCOG.load_GCOG_4deg(TGM_prox[0], TGM_prox[1], TGM_prox[2])

    # ── NORMALISATION LOOP (max 30 iterations) ────────────────────────────
    print("Starting normalisation loop...")
    for norm_iter in range(30):
        model = np.ones(share.dimsp, dtype=DP)
        mm.make_model_TGM_quick(model, share.TGM.copy())
        share.f_model = model.copy()

        share.cont = fc.fit_cont(share.f_sp, share.f_model,
                                 share.weights, share.rad_pix)
        share.f_sp_norm = share.f_sp / share.cont

        # Estimate S/N if not provided
        if sp.sn_flag:
            share.sn = utils.new_sn(share.sn)
            if share.sn > 0:
                share.sig_noise[:] = 1.0 / share.sn
        else:
            share.sn = sp.sn_ratio
            share.sig_noise[:] = 1.0 / share.sn

        # Update local S/N per pixel
        share.sn_var = utils.find_sn_var(share.f_sp_norm, share.f_model,
                                          share.sn_var)

        # Quick TGM minimisation
        tgmx = utils.normalize_pars(share.TGM, share.TGM_mask)
        dim_var = len(tgmx)

        def chi_q(x):
            return mini.chi_TGM_Q(share.dimsp, dim_var, x,
                                   share.f_sp, share.cont, share.sig_noise,
                                   share.TGM_mask)

        tgmx, _, info = mini.lmdif1(chi_q, share.dimsp, dim_var, tgmx)
        share.TGM = utils.denormalize_pars(tgmx, share.TGM_mask)

        utils.select_lines(ew_poly4, share.coeff_4deg,
                           share.TGM[:4], flag=False)
        utils.cosmic_rej()

        if info in (1, 2, 3):
            break

    # ── OUTER LOOP: full TGM fitting ──────────────────────────────────────
    print("Starting TGM fitting...")
    max_outer = 20
    flag_lim = flag_limS = flag_move = False
    info_tgm = 0  # last LM exit code from outer TGM loop (Fortran lmdif1 info)

    for outer in range(max_outer):
        # Load GCOG for current grid point
        TGM_prox, flag_lim, flag_limS, flag_move = utils.find_proxS(
            share.TGM[:3], temp_gridS, logg_gridS, met_gridS, TGM_prox
        )

        if flag_move:
            TGM_prox, temp_gridS, logg_gridS, met_gridS = \
                update_gridS_ll(TGM_prox, temp_gridS, logg_gridS, met_gridS)

        utils.select_lines(ew_poly4, share.coeff_4deg, share.TGM[:4], flag=True)

        tgmx = utils.normalize_pars(share.TGM, share.TGM_mask)
        dim_var = len(tgmx)

        def chi_tgm(x):
            return mini.chi_TGM(share.dimsp, dim_var, x,
                                 share.f_sp_norm, share.sig_noise,
                                 share.TGM_mask)

        tgmx, _, info_tgm = mini.lmdif1(chi_tgm, share.dimsp, dim_var, tgmx)
        share.TGM = utils.denormalize_pars(tgmx, share.TGM_mask)

        # ── INNER LOOP: abundance fitting ──────────────────────────────
        if sp.flag_ABD_loop:
            for inner in range(10):
                packed_abd = share.ABD[share.ABD_mask]
                dim_abd = len(packed_abd)
                if dim_abd == 0:
                    break

                def chi_abd(x):
                    return mini.chi_ABD(share.dimsp, dim_abd, x,
                                        share.f_sp_norm, share.sig_noise,
                                        share.ABD_mask)

                packed_abd, _, info_abd = mini.lmdif1(
                    chi_abd, share.dimsp, dim_abd, packed_abd)
                share.ABD[share.ABD_mask] = packed_abd

                if np.max(np.abs(share.ABD - share.ABD_old)) < 0.01:
                    break
                share.ABD_old = share.ABD.copy()

        # Continuum update
        model = np.ones(share.dimsp, dtype=DP)
        mm.make_model_TGM(model, share.TGM.copy())
        share.f_model = model.copy()
        share.cont = fc.fit_cont(share.f_sp, share.f_model,
                                  share.weights, share.rad_pix)
        share.f_sp_norm = share.f_sp / share.cont

        if sp.sn_flag:
            share.sn = utils.new_sn(share.sn)
            if share.sn > 0:
                share.sig_noise[:] = 1.0 / share.sn

        share.sn_var = utils.find_sn_var(share.f_sp_norm, share.f_model,
                                          share.sn_var)

        if flag_lim or info_tgm in (1, 2, 3):
            break

    # ── Count weighted pixels ─────────────────────────────────────────────
    share.n_weig = int(np.count_nonzero(share.weights > 1e-3))

    # ── Compute final chi-squared ─────────────────────────────────────────
    model_final = np.ones(share.dimsp, dtype=DP)
    mm.make_model_TGM(model_final, share.TGM.copy())
    share.f_model = model_final.copy()
    chisq = float(np.sum(
        ((share.f_sp_norm - share.f_model) / share.sig_noise) ** 2
    ))

    # ── Uncertainty estimation ────────────────────────────────────────────
    flag_err = False
    if sp.error_est:
        try:
            unc.TGM_errors(chisq)
            flag_err = True
        except Exception as e:
            print(f"Warning: uncertainty estimation failed: {e}")
            unc.TGM_errors_null()
    else:
        unc.TGM_errors_null()

    # ── Write results ─────────────────────────────────────────────────────
    conv = int(info_tgm)
    wr.write_res(flag_err, chisq, conv)
    print("SP_Ace finished successfully.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python space.py <parameter_file>")
        sys.exit(1)
    main(sys.argv[1])
