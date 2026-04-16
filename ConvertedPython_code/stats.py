"""
stats.py — SP_Ace statistics routines
Translated from stats.f95 (Corrado Boeche, 2016, GPLv3)
"""

import numpy as np
from num_type import DP


def avg_std_mask(dat, mask):
    """
    Compute mean, variance, and standard deviation of *dat* over *mask*.

    Parameters
    ----------
    dat  : array_like, shape (N,)
    mask : boolean array_like, shape (N,)

    Returns
    -------
    mean : float
    var  : float
    sdev : float
    """
    dat = np.asarray(dat, dtype=DP)
    mask = np.asarray(mask, dtype=bool)
    n = int(np.count_nonzero(mask))
    if n == 0:
        return 0.0, 0.0, 0.0
    mean = np.sum(dat[mask]) / n
    discr = dat - mean
    var = np.sum((discr**2)[mask]) / n
    sdev = np.sqrt(var)
    return mean, var, sdev


def avg_std(dat):
    """
    Compute mean, variance, and standard deviation of the full array *dat*.

    Parameters
    ----------
    dat : array_like, shape (N,)

    Returns
    -------
    mean : float
    var  : float
    sdev : float
    """
    dat = np.asarray(dat, dtype=DP)
    n = dat.size
    mean = np.sum(dat) / n
    discr = dat - mean
    var = np.sum(discr**2) / n
    sdev = np.sqrt(var)
    return mean, var, sdev


def avg_std_(dat, mask=None):
    """
    Dispatcher matching the Fortran INTERFACE avg_std_:
    - If *mask* is provided, calls avg_std_mask(dat, mask)
    - Otherwise calls avg_std(dat)
    """
    if mask is not None:
        return avg_std_mask(dat, mask)
    return avg_std(dat)
