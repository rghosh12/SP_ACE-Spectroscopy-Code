"""
uncertains2.py — SP_Ace uncertainty estimation
Translated from uncertains2.f95 (Corrado Boeche, 2016, GPLv3)

Uses scipy.linalg.lstsq instead of the Fortran LAPACK DGELSS wrapper.
"""

import numpy as np
from scipy.linalg import lstsq
import share
from num_type import DP
from utils import ABD_mask, write_ABD_mask

# ── Module-level state ────────────────────────────────────────────────────────
up_TGM = np.zeros(3, dtype=DP)
lo_TGM = np.zeros(3, dtype=DP)

GRADO = 2
DIM_COEFF = 10          # for grado=2
PROB = np.array([3.53, 6.25, 8.02], dtype=DP)  # 68.3%, 90%, 95.4%

_matrixA = np.zeros((27, 10), dtype=DP)
_vecB = np.zeros((27, 1), dtype=DP)
_TGM_test = np.zeros((27, 3), dtype=DP)
_chisq_target = 0.0
_pos = 0        # 0-based (Fortran 1-based → subtract 1)
_n_eq = 27
_count_point = 0
_flag_TGM = True
_flag_write = False


def normalize3D(tgm):
    tgm = np.asarray(tgm, dtype=DP).copy()
    tgm[0] /= 5000.0
    tgm[1] /= 3.0
    return tgm


def denormalize3D(tgm):
    tgm = np.asarray(tgm, dtype=DP).copy()
    tgm[0] *= 5000.0
    tgm[1] *= 3.0
    return tgm


def eval_poly(tgm, coeff):
    """Evaluate the 3D quadratic polynomial at *tgm* using *coeff*."""
    tgm = np.asarray(tgm, dtype=DP)
    val = 0.0
    count = 0
    for i1 in range(GRADO + 1):
        a = tgm[0] ** i1
        for i2 in range(GRADO - i1 + 1):
            b = tgm[1] ** i2
            for i3 in range(GRADO - i1 - i2 + 1):
                val += coeff[count] * a * b * (tgm[2] ** i3)
                count += 1
    return val


def _build_matrix_row(tgm):
    """Return the 10-element row vector for the quadratic basis at *tgm*."""
    tgm = np.asarray(tgm, dtype=DP)
    row = np.empty(DIM_COEFF, dtype=DP)
    count = 0
    for i1 in range(GRADO + 1):
        a = tgm[0] ** i1
        for i2 in range(GRADO - i1 + 1):
            b = tgm[1] ** i2
            for i3 in range(GRADO - i1 - i2 + 1):
                row[count] = a * b * (tgm[2] ** i3)
                count += 1
    return row


def chisq_tgm_e(x):
    """χ² evaluated at TGM point *x* (normalised, shape (3,))."""
    from make_model import make_model_TGM
    xTGM = share.TGM.copy().astype(DP)
    xTGM[:3] = denormalize3D(x)
    model = np.ones(share.dimsp, dtype=DP)
    make_model_TGM(model, xTGM)
    return float(np.sum(((share.f_sp_norm - model) / share.sig_noise) ** 2))


def chisq_abd_e(x):
    """χ² evaluated at ABD point *x* = [Teff/5000, logg/3, ABD[pos]] (norm.)."""
    from make_model import make_model_ABDerr
    parTGM = share.TGM.copy().astype(DP)
    parTGM[:2] = denormalize3D(x[:2])
    parABD = share.ABD.copy()
    parABD[_pos] = x[2]
    model = np.ones(share.dimsp, dtype=DP)
    make_model_ABDerr(model, parTGM, parABD)
    return float(np.sum(((share.f_sp_norm - model) / share.sig_noise) ** 2))


def _chisq(tgm_norm):
    if _flag_TGM:
        return chisq_tgm_e(tgm_norm)
    return chisq_abd_e(tgm_norm)


def prepare_matrix(m_start, m_end):
    global _matrixA, _vecB
    for i in range(m_start, m_end + 1):
        _matrixA[i, :] = _build_matrix_row(_TGM_test[i, :])
        _vecB[i, 0] = _chisq(_TGM_test[i, :])


def coeff_poly_find(TGM_ini, step):
    """Fit the chi^2 surface with a 3D quadratic polynomial."""
    global _n_eq, _matrixA, _vecB, _TGM_test
    _n_eq = 27
    define_grid(TGM_ini, step)
    prepare_matrix(0, _n_eq - 1)
    A_sub = _matrixA[:_n_eq, :]
    b_sub = _vecB[:_n_eq, :]
    coeff, _, _, _ = lstsq(A_sub, b_sub, cond=1e-18)
    return coeff[:DIM_COEFF, 0]


def define_grid(TGM_ini, step):
    """Place 27 test points around *TGM_ini* with spacing *step*."""
    global _TGM_test, _n_eq
    from data_lib import temp_gridL, logg_gridL, met_gridL
    TGM_ = denormalize3D(TGM_ini)
    step_ = denormalize3D(step)

    # Clamp to grid
    TGM_[0] = np.clip(TGM_[0], temp_gridL[0], temp_gridL[-1])
    TGM_[1] = np.clip(TGM_[1], logg_gridL[0], logg_gridL[-1])
    if _flag_TGM:
        TGM_[2] = np.clip(TGM_[2], met_gridL[0], met_gridL[-1])
    else:
        TGM_[2] = np.clip(TGM_[2], -0.4, 0.8)

    TGM_n = normalize3D(TGM_)
    step_n = normalize3D(step_)
    count = 0
    for k in [-1, 0, 1]:
        for j in [-1, 0, 1]:
            for i in [-1, 0, 1]:
                if not (k == 0 and j == 0 and i == 0):
                    pt = TGM_n + step_n * np.array([i, j, k], dtype=DP)
                    # Find chi-sq border along this direction
                    pt_border, chi, found = find_chisq_border(TGM_n, pt)
                    _TGM_test[count, :] = pt_border if found else TGM_n
                    count += 1
    _TGM_test[count, :] = TGM_n
    _n_eq = 27


def stick_one_point(TGM_loc, coeff):
    """Add a new chi-sq evaluation point and refit the polynomial."""
    global _count_point, _matrixA, _vecB, _TGM_test
    _count_point += 1
    _TGM_test[_count_point, :] = TGM_loc
    _matrixA[_count_point, :] = _build_matrix_row(TGM_loc)
    _vecB[_count_point, 0] = _chisq(TGM_loc)
    A_sub = _matrixA[:_n_eq, :]
    b_sub = _vecB[:_n_eq, :]
    coeff_new, _, _, _ = lstsq(A_sub, b_sub, cond=1e-60)
    return coeff_new[:DIM_COEFF, 0]


def check_up_limits(step_loc, TGM_loc):
    from data_lib import temp_gridL, logg_gridL, met_gridL
    TGM_d = denormalize3D(TGM_loc)
    flag = False
    if TGM_d[0] > temp_gridL[-1]:
        TGM_d[0] = temp_gridL[-1] - step_loc[0]
        flag = True
    if TGM_d[1] > logg_gridL[-1]:
        TGM_d[1] = logg_gridL[-1] - step_loc[1]
        flag = True
    if _flag_TGM:
        if TGM_d[2] > met_gridL[-1]:
            TGM_d[2] = met_gridL[-1] - step_loc[2]
            flag = True
    else:
        if TGM_d[2] > 0.8:
            TGM_d[2] = 0.8 - step_loc[2]
            flag = True
    return normalize3D(TGM_d), flag


def check_lo_limits(step_loc, TGM_loc):
    from data_lib import temp_gridL, logg_gridL, met_gridL
    TGM_d = denormalize3D(TGM_loc)
    flag = False
    if TGM_d[0] < temp_gridL[0]:
        TGM_d[0] = temp_gridL[0] + step_loc[0]
        flag = True
    if TGM_d[1] < logg_gridL[0]:
        TGM_d[1] = logg_gridL[0] + step_loc[1]
        flag = True
    if _flag_TGM:
        if TGM_d[2] < met_gridL[0]:
            TGM_d[2] = met_gridL[0] + step_loc[2]
            flag = True
    else:
        if TGM_d[2] < -0.4:
            TGM_d[2] = -0.4 + step_loc[2]
            flag = True
    return normalize3D(TGM_d), flag


def find_extremes(coeff, chisq_best, flag_up_low):
    """
    Analytically solve the quadratic chi-sq surface for the extremes.

    Returns
    -------
    TGM_loc  : ndarray, shape (3,)
    flag_D1  : bool — True if the ellipsoid fit did not converge
    """
    # Map coeff to quadratic form: a*x^2 + b*y^2 + c*z^2 + d*xy + e*xz +
    #                              f*yz + g*x + h*y + i*z + j = chi
    c10 = coeff[9]  # Teff^2
    c9  = coeff[8]  # Teff*logg
    c8  = coeff[7]  # Teff*met
    c7  = coeff[6]  # Teff
    c6  = coeff[5]  # logg^2
    c5  = coeff[4]  # logg*met
    c4  = coeff[3]  # logg
    c3  = coeff[2]  # met^2
    c2  = coeff[1]  # met
    c1  = coeff[0]  # const

    # Build per-pos mapping
    mats = {
        0: (c10, c9, c8, c7, c6, c9, c5, c4, c3, c5, c8, c2),
        1: (c6,  c9, c5, c4, c10, c9, c8, c7, c3, c8, c5, c2),
        2: (c3,  c5, c8, c2, c6,  c5, c9, c4, c10, c9, c8, c7),
    }

    pos_local = _pos if _flag_TGM else 2
    am, d, e, g, bm, _d2, f, h, cm, _f2, _e2, i_c = mats[pos_local]

    j_c = c1 - (chisq_best + PROB[0])

    A1 = -cm * d * e * f + cm**2 * d**2 + am * cm * f**2 + bm * cm * e**2 - 4 * am * bm * cm**2
    B1 = (-cm * d * f * i_c - cm * e * f * h + 2 * cm**2 * d * h + 2 * bm * cm * e * i_c
          + cm * f**2 * g - 4 * bm * cm**2 * g)
    C1 = -cm * f * h * i_c + cm**2 * h**2 + cm * f**2 * j_c + bm * cm * i_c**2 - 4 * bm * cm**2 * j_c
    D1 = B1**2 - 4.0 * A1 * C1

    cond1 = d**2 - 4 * am * bm
    cond2 = e**2 - 4 * am * cm
    cond3 = f**2 - 4 * bm * cm

    if D1 > 0.0 and cond1 < 0 and cond2 < 0 and cond3 < 0:
        x_m = (-B1 + np.sqrt(D1)) / (2 * A1)
        x_p = (-B1 - np.sqrt(D1)) / (2 * A1)
        y_m = -(2 * e * f * x_m + 2 * f * i_c - 4 * cm * d * x_m - 4 * cm * h) / (2 * (f**2 - 4 * bm * cm))
        y_p = -(2 * e * f * x_p + 2 * f * i_c - 4 * cm * d * x_p - 4 * cm * h) / (2 * (f**2 - 4 * bm * cm))
        z_m = -(e * x_m + f * y_m + i_c) / (2 * cm)
        z_p = -(e * x_p + f * y_p + i_c) / (2 * cm)

        order = {0: [0, 1, 2], 1: [1, 0, 2], 2: [2, 1, 0]}
        idx = order[pos_local]
        up_pt = np.array([x_p, y_p, z_p], dtype=DP)
        lo_pt = np.array([x_m, y_m, z_m], dtype=DP)
        pts = [up_pt[idx], lo_pt[idx]]
        TGM_loc = pts[0] if flag_up_low else pts[1]
        return TGM_loc, False
    return np.zeros(3, dtype=DP), True


def est_errors_p(coeff, chisq_best, chisq_tol, step, TGM_p):
    """Estimate upper error on parameter at position *_pos*."""
    global _count_point
    coeff_loc = coeff.copy()
    TGM_loc = TGM_p.copy()
    pos_loc = _pos if _flag_TGM else 2
    flag_solved = False

    for _ in range(20):
        TGM_loc, flag_D1 = find_extremes(coeff_loc, chisq_best, True)
        if flag_D1:
            coeff_loc = stick_one_point(TGM_loc, coeff_loc)
            continue
        if TGM_loc[pos_loc] <= TGM_p[pos_loc]:
            break
        TGM_loc, flag_out = check_up_limits(step, TGM_loc)
        chisq = _chisq(TGM_loc)
        chisq_poly = eval_poly(TGM_loc, coeff_loc)

        if chisq < _chisq_target - chisq_tol:
            if chisq_poly < _chisq_target - chisq_tol and flag_out:
                break
            coeff_loc = stick_one_point(TGM_loc, coeff_loc)
        elif chisq > _chisq_target + chisq_tol:
            coeff_loc = stick_one_point(TGM_loc, coeff_loc)
        else:
            flag_solved = True
            break

    if flag_solved:
        TGM_loc = denormalize3D(TGM_loc)
        if _flag_TGM:
            up_TGM[_pos] = TGM_loc[_pos]
        else:
            share.up_ABD[_pos] = TGM_loc[2]
    else:
        if _flag_TGM:
            up_TGM[_pos] = 10000.0
        else:
            share.up_ABD[_pos] = 10000.0


def est_errors_m(coeff, chisq_best, chisq_tol, step, TGM_m):
    """Estimate lower error on parameter at position *_pos*."""
    global _count_point
    coeff_loc = coeff.copy()
    TGM_loc = TGM_m.copy()
    pos_loc = _pos if _flag_TGM else 2
    flag_solved = False

    for _ in range(20):
        TGM_loc, flag_D1 = find_extremes(coeff_loc, chisq_best, False)
        if flag_D1:
            coeff_loc = stick_one_point(TGM_loc, coeff_loc)
            continue
        TGM_loc, flag_out = check_lo_limits(step, TGM_loc)
        if TGM_loc[pos_loc] >= TGM_m[pos_loc]:
            break
        chisq = _chisq(TGM_loc)
        chisq_poly = eval_poly(TGM_loc, coeff_loc)

        if chisq < _chisq_target - chisq_tol:
            if chisq_poly < _chisq_target - chisq_tol and flag_out:
                break
            coeff_loc = stick_one_point(TGM_loc, coeff_loc)
        elif chisq > _chisq_target + chisq_tol:
            coeff_loc = stick_one_point(TGM_loc, coeff_loc)
        else:
            flag_solved = True
            break

    if flag_solved:
        TGM_loc = denormalize3D(TGM_loc)
        if _flag_TGM:
            lo_TGM[_pos] = TGM_loc[_pos]
        else:
            share.lo_ABD[_pos] = TGM_loc[2]
    else:
        if _flag_TGM:
            lo_TGM[_pos] = -10000.0
        else:
            share.lo_ABD[_pos] = -10000.0


def TGM_errors(chisq_best):
    """
    Estimate 1σ uncertainties on Teff, logg, [M/H], and all abundances.

    Parameters
    ----------
    chisq_best : float — best-fit chi-squared value
    """
    global _chisq_target, _flag_TGM, _pos, _count_point

    _chisq_target = chisq_best + PROB[0]

    step_raw = np.array([5.0, 0.02, 0.01], dtype=DP)
    step = normalize3D(step_raw)
    chisq_tol = 0.2

    # ── TGM errors ────────────────────────────────────────────────────────
    _flag_TGM = True
    TGM_temp = normalize3D(share.TGM[:3])

    for pos in range(3):
        _pos = pos
        _count_point = 0
        coeff = coeff_poly_find(TGM_temp, step)
        est_errors_p(coeff, chisq_best, chisq_tol, step, TGM_temp)

    TGM_temp = normalize3D(share.TGM[:3])
    for pos in range(3):
        _pos = pos
        _count_point = 0
        coeff = coeff_poly_find(TGM_temp, -step)
        est_errors_m(coeff, chisq_best, chisq_tol, step, TGM_temp)

    # ── Abundance errors ──────────────────────────────────────────────────
    _flag_TGM = False
    step_abd = np.array([0.001, 0.03, 0.01], dtype=DP)
    ABD = share.ABD
    wm = share.write_ABD_mask if share.write_ABD_mask is not None else np.ones(len(ABD), dtype=bool)
    am = share.ABD_mask if share.ABD_mask is not None else np.ones(len(ABD), dtype=bool)

    for i in range(len(ABD)):
        _pos = i
        if wm[i] and am[i] and -0.599 < ABD[i] < 0.799:
            TGM_temp2 = np.array([share.TGM[0], share.TGM[1], ABD[i]], dtype=DP)
            TGM_temp2 = normalize3D(TGM_temp2)
            _count_point = 0
            coeff = coeff_poly_find(TGM_temp2, step_abd)
            est_errors_p(coeff, chisq_best, chisq_tol, step_abd, TGM_temp2)
            coeff = coeff_poly_find(TGM_temp2, -step_abd)
            _count_point = 0
            est_errors_m(coeff, chisq_best, chisq_tol, step_abd, TGM_temp2)
        else:
            share.up_ABD[i] = 10.0
            share.lo_ABD[i] = -10.0


def TGM_errors_null():
    """Set all uncertainties to null/placeholder values."""
    for i in range(3):
        up_TGM[i] = 9999.0
        lo_TGM[i] = -9999.0
    ABD = share.ABD
    wm = share.write_ABD_mask
    if ABD is not None and wm is not None:
        for i in range(len(ABD)):
            if wm[i]:
                share.up_ABD[i] = 9.99
                share.lo_ABD[i] = -9.99


def find_chisq_border(point1, point2):
    """
    Locate the chi-sq = chisq_target boundary along the direction
    from *point1* toward *point2*.

    Returns
    -------
    tgm_out   : ndarray, shape (3,)
    chi       : float
    flag_found: bool
    """
    v_comp = point2 - point1
    t = np.array([0.0, 1.0, 2.0], dtype=DP)
    vec = _set_t_vec(point1, v_comp, t)
    chi = -1.0
    tgm_loc = point1.copy()

    for _ in range(10):
        t_p, t_m, flag_solved = _find_zeros(t, vec)
        if not flag_solved:
            t += np.array([0.0, 1.0, 2.0])
            if np.any(t > 100):
                return point1 + v_comp * t[1], chi, True
            vec = _set_t_vec(point1, v_comp, t)
            continue

        tgm_loc = point1 + v_comp * t_p
        chi = _chisq(tgm_loc)
        if abs(chi - _chisq_target) < 0.05:
            break
        idx = int(np.argmax(np.abs(vec - chi)))
        vec[idx] = chi
        t[idx] = t_p

    return tgm_loc, chi, True


def _set_t_vec(point1, v_comp, t):
    vec = np.empty(3, dtype=DP)
    for i in range(3):
        vec[i] = _chisq(point1 + v_comp * t[i])
    return vec


def _find_zeros(t, vec):
    """Fit a parabola through (t, vec) and find where it equals chisq_target."""
    mat = np.column_stack([t**2, t, np.ones(3)])
    det = np.linalg.det(mat)
    if abs(det) < 1e-30:
        return t[1], t[1], False

    def _det_swap(col, vals):
        m = mat.copy()
        m[:, col] = vals
        return np.linalg.det(m)

    a = _det_swap(0, vec) / det
    b = _det_swap(1, vec) / det
    c = _det_swap(2, vec) / det - _chisq_target
    arg = b**2 - 4 * a * c
    if arg < 0.0:
        return t[1], t[1], False

    t1_p = (-b + np.sqrt(arg)) / (2 * a)
    t1_m = (-b - np.sqrt(arg)) / (2 * a)
    return max(t1_p, t1_m), min(t1_p, t1_m), True
