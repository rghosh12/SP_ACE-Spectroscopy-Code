"""
minimize.py — SP_Ace χ² computation & Levenberg-Marquardt minimisation
Translated from minimize.f95 (Corrado Boeche, 2016, GPLv3)

The core lmdif/lmdif1/lmpar/qrfac/qrsolv routines are faithful ports of the
MINPACK routines by Jorge More, Burton Garbow, Kenneth Hillstrom
(Argonne National Laboratory, 1980), as adapted to Fortran 90 by John Burkardt
and further modified by Corrado Boeche.
"""

import numpy as np
import share
import error
from num_type import DP
from utils import denormalize_pars


# ── χ² evaluation routines ────────────────────────────────────────────────────

def chi_TGM_Q(dimsp, dim_var, tgmx, f_sp, cont, sig_noise, TGM_mask):
    """
    Discrepancy vector for the quick (6th-degree) TGM model.

    Returns
    -------
    discrep : ndarray, shape (dimsp,)
    """
    from make_model import make_model_TGM_quick
    tgm_local = denormalize_pars(tgmx, TGM_mask)
    sp_norm = f_sp / tgm_local[5] / cont
    model = np.ones(dimsp, dtype=DP)
    make_model_TGM_quick(model, tgm_local)
    return (sp_norm - model) / sig_noise


def chi_TGM(dimsp, dim_var, tgmx, f_sp_norm, sig_noise, TGM_mask):
    """
    Discrepancy vector for the full (4th-degree) TGM model.

    Returns
    -------
    discrep : ndarray, shape (dimsp,)
    """
    from make_model import make_model_TGM
    tgm_local = denormalize_pars(tgmx, TGM_mask)
    model = np.ones(dimsp, dtype=DP)
    make_model_TGM(model, tgm_local)
    return (f_sp_norm - model) / sig_noise


def chi_ABD(dimsp, dim_var, abdx, f_sp_norm, sig_noise, ABD_mask):
    """
    Discrepancy vector for the abundance model, with out-of-range penalty.

    Returns
    -------
    discrep : ndarray, shape (dimsp,)
    """
    from make_model import make_model_ABD
    abd = share.ABD
    abd_in = np.where(ABD_mask, abdx, abd).copy()
    mask_internal = (abd_in > 0.8) | (abd_in < -0.6)
    penal = np.abs(abd_in)
    model = np.ones(dimsp, dtype=DP)
    make_model_ABD(model, abd_in)
    discrep = (f_sp_norm - model) / sig_noise

    if np.count_nonzero(mask_internal) > 0:
        penal[mask_internal] = np.abs(penal[mask_internal] - 2.0) / 1e6
        penalty = np.sum(penal[mask_internal]) / dimsp
        discrep += np.sign(discrep) * penalty

    return discrep


# ── MINPACK routines (Levenberg-Marquardt) ────────────────────────────────────

def enorm(x):
    """Euclidean norm of vector x."""
    return float(np.sqrt(np.sum(np.asarray(x, dtype=DP) ** 2)))


def qrsolv(r, ipvt, diag, qtb):
    """
    Solve the system  [R; sqrt(par)*D]*x = [Q'b; 0]  in the least-squares
    sense, where R is upper-triangular.

    Parameters
    ----------
    r    : ndarray, shape (n, n) — modified in place
    ipvt : ndarray, shape (n,)  — permutation
    diag : ndarray, shape (n,)  — diagonal of D
    qtb  : ndarray, shape (n,)  — first n elements of Q'b

    Returns
    -------
    x     : ndarray, shape (n,)
    sdiag : ndarray, shape (n,)
    """
    n = len(qtb)
    r = r.copy()
    wa = qtb.copy()
    x = np.array([r[j, j] for j in range(n)], dtype=DP)
    sdiag = np.zeros(n, dtype=DP)

    for j in range(n):
        r[j:n, j] = r[j, j:n].copy()
        x[j] = r[j, j]
    wa[:] = qtb

    for j in range(n):
        l = ipvt[j]
        if diag[l] != 0.0:
            sdiag[j:] = 0.0
            sdiag[j] = diag[l]
            qtbpj = 0.0
            for k in range(j, n):
                if sdiag[k] != 0.0:
                    if abs(r[k, k]) < abs(sdiag[k]):
                        cotan = r[k, k] / sdiag[k]
                        s = 0.5 / np.sqrt(0.25 + 0.25 * cotan ** 2)
                        c = s * cotan
                    else:
                        t = sdiag[k] / r[k, k]
                        c = 0.5 / np.sqrt(0.25 + 0.25 * t ** 2)
                        s = c * t
                    r[k, k] = c * r[k, k] + s * sdiag[k]
                    temp = c * wa[k] + s * qtbpj
                    qtbpj = -s * wa[k] + c * qtbpj
                    wa[k] = temp
                    for i in range(k + 1, n):
                        temp = c * r[i, k] + s * sdiag[i]
                        sdiag[i] = -s * r[i, k] + c * sdiag[i]
                        r[i, k] = temp
        sdiag[j] = r[j, j]
        r[j, j] = x[j]

    nsing = n
    for j in range(n):
        if sdiag[j] == 0.0 and nsing == n:
            nsing = j
        if nsing < n:
            wa[j] = 0.0

    for j in range(nsing - 1, -1, -1):
        s2 = np.dot(wa[j + 1:nsing], r[j + 1:nsing, j])
        wa[j] = (wa[j] - s2) / sdiag[j]

    x_out = np.empty(n, dtype=DP)
    for j in range(n):
        x_out[ipvt[j]] = wa[j]
    return x_out, sdiag


def qrfac(a, pivot=True):
    """
    QR factorisation of a with optional column pivoting.

    Returns
    -------
    a      : ndarray — factored form (upper triangle = R, lower = Householder)
    ipvt   : ndarray, shape (n,) — column permutation
    rdiag  : ndarray, shape (n,)
    acnorm : ndarray, shape (n,)
    """
    m, n = a.shape
    a = a.copy()
    acnorm = np.array([enorm(a[:, j]) for j in range(n)], dtype=DP)
    rdiag = acnorm.copy()
    wa = acnorm.copy()
    ipvt = np.arange(n, dtype=np.int32) if pivot else np.arange(n, dtype=np.int32)
    epsmch = np.finfo(DP).eps
    minmn = min(m, n)

    for j in range(minmn):
        if pivot:
            kmax = j + int(np.argmax(rdiag[j:]))
            if kmax != j:
                a[:, [j, kmax]] = a[:, [kmax, j]]
                rdiag[kmax] = rdiag[j]
                wa[kmax] = wa[j]
                ipvt[[j, kmax]] = ipvt[[kmax, j]]

        ajnorm = enorm(a[j:, j])
        if ajnorm != 0.0:
            if a[j, j] < 0.0:
                ajnorm = -ajnorm
            a[j:, j] /= ajnorm
            a[j, j] += 1.0
            for k in range(j + 1, n):
                temp = np.dot(a[j:, j], a[j:, k]) / a[j, j]
                a[j:, k] -= temp * a[j:, j]
                if pivot and rdiag[k] != 0.0:
                    t = a[j, k] / rdiag[k]
                    rdiag[k] *= np.sqrt(max(0.0, 1.0 - t ** 2))
                    if 0.05 * (rdiag[k] / wa[k]) ** 2 <= epsmch:
                        rdiag[k] = enorm(a[j + 1:, k]) if j + 1 < m else 0.0
                        wa[k] = rdiag[k]
        rdiag[j] = -ajnorm

    return a, ipvt, rdiag, acnorm


def fdjac2(fcn, m, n, x, fvec, epsfcn=0.0):
    """
    Forward-difference Jacobian approximation.

    Returns
    -------
    fjac : ndarray, shape (m, n)
    """
    epsmch = np.finfo(DP).eps
    eps = np.sqrt(max(epsfcn, epsmch))
    fjac = np.empty((m, n), dtype=DP)
    x = x.copy()
    for j in range(n):
        temp = x[j]
        h = eps * abs(temp) if temp != 0.0 else eps
        x[j] = temp + h
        wa = fcn(x)
        x[j] = temp
        fjac[:, j] = (wa - fvec) / h
    return fjac


def lmpar(r, ipvt, diag, qtb, delta, par):
    """
    Compute the Levenberg-Marquardt parameter.

    Returns
    -------
    par    : float
    x      : ndarray, shape (n,)
    sdiag  : ndarray, shape (n,)
    """
    n = len(qtb)
    r = r.copy()
    dwarf = np.finfo(DP).tiny

    nsing = n
    wa1 = qtb.copy()
    for j in range(n):
        if r[j, j] == 0.0 and nsing == n:
            nsing = j
        if nsing < n:
            wa1[j] = 0.0

    for k in range(nsing - 1, -1, -1):
        wa1[k] /= r[k, k]
        wa1[:k] -= r[:k, k] * wa1[k]

    x = np.empty(n, dtype=DP)
    for j in range(n):
        x[ipvt[j]] = wa1[j]

    wa2 = diag * x
    dxnorm = enorm(wa2)
    fp = dxnorm - delta
    if fp <= 0.1 * delta:
        return 0.0, x, np.zeros(n, dtype=DP)

    parl = 0.0
    if n <= nsing:
        wa1_p = np.array([diag[ipvt[j]] * wa2[ipvt[j]] / dxnorm for j in range(n)])
        for j in range(n):
            wa1_p[j] = (wa1_p[j] - np.dot(wa1_p[:j], r[:j, j])) / r[j, j]
        temp = enorm(wa1_p)
        parl = ((fp / delta) / temp) / temp

    wa1_g = np.array([np.dot(qtb[:j + 1], r[:j + 1, j]) / diag[ipvt[j]] for j in range(n)])
    gnorm = enorm(wa1_g)
    paru = gnorm / delta if gnorm / delta != 0.0 else dwarf / min(delta, 0.1)

    par = np.clip(par, parl, paru)
    if par == 0.0:
        par = gnorm / dxnorm

    for _ in range(10):
        if par == 0.0:
            par = max(dwarf, 0.001 * paru)
        wa1_s = np.sqrt(par) * diag
        x, sdiag = qrsolv(r, ipvt, wa1_s, qtb)
        wa2 = diag * x
        dxnorm = enorm(wa2)
        fp_old = fp
        fp = dxnorm - delta
        if abs(fp) <= 0.1 * delta:
            break
        if parl == 0.0 and fp <= fp_old < 0.0:
            break
        wa1_c = np.array([diag[ipvt[j]] * wa2[ipvt[j]] / dxnorm for j in range(n)])
        for j in range(n):
            wa1_c[j] /= sdiag[j]
            wa1_c[j + 1:] -= r[j + 1:, j] * wa1_c[j]
        temp = enorm(wa1_c)
        parc = ((fp / delta) / temp) / temp
        if fp > 0.0:
            parl = max(parl, par)
        else:
            paru = min(paru, par)
        par = max(parl, par + parc)

    return par, x, sdiag


def lmdif(fcn, m, n, x, ftol=1.49012e-8, xtol=1.49012e-8, gtol=0.0,
          maxfev=None, epsfcn=0.0, mode=1, factor=100.0, nprint=0):
    """
    Levenberg-Marquardt minimiser (core).

    Parameters
    ----------
    fcn    : callable(x) → ndarray(m,)
    m, n   : int
    x      : ndarray, shape (n,)  — modified in place

    Returns
    -------
    x    : ndarray, shape (n,)
    fvec : ndarray, shape (m,)
    info : int
    nfev : int
    """
    if maxfev is None:
        maxfev = 200 * (n + 1)

    x = np.array(x, dtype=DP)
    epsmch = np.finfo(DP).eps
    info = 0
    nfev = 0
    xnorm = 0.0
    delta = 0.0
    par = 0.0
    diag = np.zeros(n, dtype=DP)

    fvec = fcn(x)
    nfev = 1
    fnorm = enorm(fvec)
    iter_ = 1

    while True:
        fjac = fdjac2(fcn, m, n, x, fvec, epsfcn)
        nfev += n

        fjac, ipvt, wa1, wa2 = qrfac(fjac, pivot=True)

        if iter_ == 1:
            if mode != 2:
                diag[:] = wa2
                diag[wa2 == 0.0] = 1.0
            wa3 = diag * x
            xnorm = enorm(wa3)
            delta = factor * xnorm if xnorm != 0.0 else factor

        wa4 = fvec.copy()
        for j in range(n):
            if fjac[j, j] != 0.0:
                s = -np.dot(wa4[j:], fjac[j:, j]) / fjac[j, j]
                wa4[j:] += fjac[j:, j] * s
            fjac[j, j] = wa1[j]

        qtf = wa4[:n].copy()

        gnorm = 0.0
        if fnorm != 0.0:
            for j in range(n):
                l = ipvt[j]
                if wa2[l] != 0.0:
                    s2 = sum(fjac[i, j] * (qtf[i] / fnorm) for i in range(j + 1))
                    gnorm = max(gnorm, abs(s2 / wa2[l]))

        if gnorm <= gtol:
            info = 4
            break

        if mode != 2:
            diag = np.maximum(diag, wa2)

        while True:
            par, wa1_lm, _ = lmpar(fjac[:n, :], ipvt, diag, qtf, delta, par)
            wa1_neg = -wa1_lm
            wa2_new = x + wa1_neg
            wa3_new = diag * wa1_neg
            pnorm = enorm(wa3_new)

            if iter_ == 1:
                delta = min(delta, pnorm)

            wa4_new = fcn(wa2_new)
            nfev += 1
            fnorm1 = enorm(wa4_new)

            actred = (1.0 - (fnorm1 / fnorm) ** 2) if 0.1 * fnorm1 < fnorm else -1.0

            wa3_p = np.zeros(n, dtype=DP)
            for j in range(n):
                wa3_p[:j + 1] += fjac[:j + 1, j] * wa1_neg[ipvt[j]]
            temp1 = enorm(wa3_p) / fnorm
            temp2 = (np.sqrt(par) * pnorm) / fnorm
            prered = temp1 ** 2 + temp2 ** 2 / 0.5
            dirder = -(temp1 ** 2 + temp2 ** 2)

            ratio = actred / prered if prered != 0.0 else 0.0

            if ratio <= 0.25:
                temp = 0.5 if actred >= 0.0 else 0.5 * dirder / (dirder + 0.5 * actred)
                if 0.1 * fnorm1 >= fnorm or temp < 0.1:
                    temp = 0.1
                delta = temp * min(delta, pnorm / 0.1)
                par /= temp
            else:
                if par == 0.0 or ratio >= 0.75:
                    delta = 2.0 * pnorm
                    par *= 0.5

            if ratio >= 1e-4:
                x = wa2_new.copy()
                wa2 = diag * x
                fvec = wa4_new.copy()
                xnorm = enorm(wa2)
                fnorm = fnorm1
                iter_ += 1

            if abs(actred) <= ftol and prered <= ftol and 0.5 * ratio <= 1.0:
                info = 1
            if delta <= xtol * xnorm:
                info = 2
            if abs(actred) <= ftol and prered <= ftol and 0.5 * ratio <= 1.0 and info == 2:
                info = 3
            if info != 0:
                break
            if nfev >= maxfev:
                info = 5
            if (abs(actred) <= epsmch and prered <= epsmch and 0.5 * ratio <= 1.0):
                info = 6
            if delta <= epsmch * xnorm:
                info = 7
            if gnorm <= epsmch:
                info = 8
            if info != 0:
                break
            if ratio >= 1e-4:
                break

        if info != 0:
            break

    return x, fvec, info, nfev


def lmdif1(fcn, m, n, x, tol=1.49012e-8):
    """
    Simplified interface to the Levenberg-Marquardt minimiser.

    Parameters
    ----------
    fcn : callable(x) → ndarray(m,)
    m,n : int
    x   : ndarray, shape (n,)  — initial estimate; overwritten
    tol : float

    Returns
    -------
    x    : ndarray, shape (n,)
    fvec : ndarray, shape (m,)
    info : int
    """
    x, fvec, info, _ = lmdif(fcn, m, n, x, ftol=tol, xtol=tol)
    if info == 8:
        info = 4
    return x, fvec, info
