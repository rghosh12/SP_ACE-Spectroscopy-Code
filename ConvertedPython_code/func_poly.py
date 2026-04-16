"""
func_poly.py — SP_Ace polynomial & Voigt profile functions
Translated from func_poly.f95 (Corrado Boeche, 2016, GPLv3)
"""

import numpy as np
from num_type import DP

# ── Constants ─────────────────────────────────────────────────────────────────
_PI = np.float64(3.14159265358979)
_SQRTLN2 = np.float64(0.832554611)
_SQRTPI = np.float64(1.772453851)
_PIxSQRTPi = np.float64(5.568327996831707)       # pi * sqrtpi
_SQRTPIxSQRTLN2 = np.float64(1.4756646266356057) # sqrtpi * sqrtln2
_SQRTPI2 = np.float64(0.28209479177387814)        # 0.5 / sqrtpi
_SQRTLN2x2 = np.float64(1.6651092223153954)       # 2 * sqrtln2

_A = np.array([-1.2150, -1.3509, -1.2150, -1.3509], dtype=DP)
_B = np.array([ 1.2359,  0.3786, -1.2359, -0.3786], dtype=DP)
_C = np.array([-0.3085,  0.5906, -0.3085,  0.5906], dtype=DP)
_D = np.array([ 0.0210, -1.1858, -0.0210,  1.1858], dtype=DP)


# ── 4th-degree polynomial (70-term) ──────────────────────────────────────────

def ew_poly4(X, ele, coeff, ele2meas, X_abd):
    """
    Compute the equivalent width from a 4th-degree polynomial model.

    Parameters
    ----------
    X       : ndarray, shape (70,)  — polynomial basis vector
    ele     : float                 — atomic number of the line
    coeff   : ndarray, shape (70,)  — polynomial coefficients
    ele2meas: ndarray               — array of element atomic numbers
    X_abd   : ndarray, shape (70, dim_ele) — abundance transform matrix

    Returns
    -------
    ew : float  (clamped to >= 1e-6)
    """
    ele_pos = int(np.argmin(np.abs(ele - ele2meas)))
    X_pars = X * X_abd[:, ele_pos]
    result = float(np.dot(X_pars, coeff))
    return max(result, 1e-6)


def poly4_transform(pars, ABD, dim_ele, X_abd_out):
    """
    Build the 70-element polynomial basis vector X for a 4th-degree expansion
    over (Teff, logg, [M/H], micro-turbulence) and update *X_abd_out*.

    Parameters
    ----------
    pars      : array_like, shape (4,)  — [Teff, logg, [M/H], xi] (normalised)
    ABD       : array_like, shape (dim_ele,) — element abundances
    dim_ele   : int
    X_abd_out : ndarray, shape (70, dim_ele) — updated in place

    Returns
    -------
    X : ndarray, shape (70,)
    """
    pars = np.asarray(pars, dtype=DP)
    ABD = np.asarray(ABD, dtype=DP)
    grado = 4
    X = np.zeros(70, dtype=DP)
    X_exp = np.zeros(70, dtype=np.int16)

    count = 0
    for i1 in range(grado + 1):
        a = pars[0] ** i1
        for i2 in range(grado - i1 + 1):
            b = a * pars[1] ** i2
            for i3 in range(grado - i2 - i1 + 1):
                c = b * pars[2] ** i3
                for i4 in range(grado - i3 - i2 - i1 + 1):
                    X[count] = c
                    X_exp[count] = i4
                    count += 1

    for i1 in range(dim_ele):
        X_abd_out[:, i1] = ABD[i1] ** X_exp

    return X


# ── 6th-degree polynomial (84-term) ──────────────────────────────────────────

def ew_poly6(X, coeff):
    """
    Compute the equivalent width from a 6th-degree polynomial model.

    Parameters
    ----------
    X     : ndarray, shape (84,)
    coeff : ndarray, shape (84,)

    Returns
    -------
    ew : float  (clamped to >= 1e-6)
    """
    result = float(np.dot(X, coeff))
    return max(result, 1e-6)


def poly6_transform(pars):
    """
    Build the 84-element polynomial basis vector X for a 6th-degree expansion.

    Parameters
    ----------
    pars : array_like, shape (4,)  — [Teff, logg, [M/H], xi]

    Returns
    -------
    X : ndarray, shape (84,)
    """
    pars = np.asarray(pars, dtype=DP)
    grado = 6
    X = np.zeros(84, dtype=DP)
    count = 0
    for i1 in range(grado + 1):
        a = pars[0] ** i1
        for i2 in range(grado - i1 + 1):
            b = a * pars[1] ** i2
            for i3 in range(grado - i2 - i1 + 1):
                X[count] = b * pars[2] ** i3
                count += 1
    return X


# ── Voigt profile ─────────────────────────────────────────────────────────────

def voigt_logg(w, mu, gamG, ew_val, gamL):
    """
    Evaluate the Voigt absorption profile at wavelength *w*.

    Uses the pseudo-Voigt approximation of Thompson, Cox & Hastings.

    Parameters
    ----------
    w     : float  — wavelength (Å)
    mu    : float  — line centre (Å)
    gamG  : float  — Gaussian FWHM parameter
    ew_val: float  — equivalent width (Å)
    gamL  : float  — Lorentzian HWHM

    Returns
    -------
    absorption : float  (>= 1e-6)
    """
    sigmaL = gamL * _SQRTPI2
    aL = ew_val / (sigmaL * _PIxSQRTPi)
    X = (w - mu) * _SQRTLN2x2 / gamG
    Y = gamL * _SQRTLN2 / gamG

    V = (_C * (Y - _A) + _D * (X - _B)) / ((Y - _A) ** 2 + (X - _B) ** 2)
    result = float(np.sum(V)) * (gamL * aL * _SQRTPIxSQRTLN2 / gamG)
    return max(result, 1e-6)
