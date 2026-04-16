"""
read_GCOG.py — SP_Ace GCOG (Gaussian Curve Of Growth) library reader
Translated from read_GCOG.f95 (Corrado Boeche, 2016, GPLv3)
"""

import numpy as np
import share
import error
import space_pars
import read_sp_ll
from num_type import DP


def load_GCOG_4deg(temp, logg, met):
    """
    Load the 4th-degree GCOG coefficient file for a given stellar grid point.

    File-name convention (from Fortran):
        {GCOGlib}/{met_label}{temp_label}{logg_label}
    where
        met_label  = 'am{|met*10|:02d}-' (negative) or 'ap{met*10:02d}-'
        temp_label = 't{temp:04d}-'
        logg_label = 'g{logg*10:02d}-4degpoly.dat'

    Parameters
    ----------
    temp : float — effective temperature (K)
    logg : float — log surface gravity
    met  : float — metallicity [M/H]
    """
    logg_label = f"g{int(round(logg * 10)):02d}-4degpoly.dat"
    temp_label = f"t{int(round(temp)):04d}-"
    if met < 0:
        met_label = f"am{int(round(-met * 10)):02d}-"
    else:
        met_label = f"ap{int(round(met * 10)):02d}-"

    gcog_dir = space_pars.GCOGlib.rstrip("/")
    file_GCOG = f"{gcog_dir}/{met_label}{temp_label}{logg_label}"

    try:
        fh = open(file_GCOG, "r")
    except OSError:
        error.error_msg(10, "I cannot open the GCOG library!")
        return

    llist_mask = read_sp_ll.llist_mask
    dim_ll = share.dim_ll

    i = 0  # index into selected lines
    for k, selected in enumerate(llist_mask):
        line = fh.readline()
        if not line:
            break
        if selected:
            tokens = line.split()
            if len(tokens) >= 70:
                share.coeff_4deg[:, i] = [float(t) for t in tokens[:70]]
            i += 1

    fh.close()


def load_GCOG_4deg_quick():
    """
    Load the 6th-degree 'quick' GCOG coefficient file
    (``space_6degpoly.dat``) into *share.coeff_4deg_quick*.
    """
    gcog_dir = space_pars.GCOGlib.rstrip("/")
    file_GCOG = f"{gcog_dir}/space_6degpoly.dat"

    try:
        fh = open(file_GCOG, "r")
    except OSError:
        error.error_msg(11, "I cannot open the space_6degpoly.dat file!")
        return

    llist_mask = read_sp_ll.llist_mask

    i = 0
    for k, selected in enumerate(llist_mask):
        line = fh.readline()
        if not line:
            break
        if selected:
            tokens = line.split()
            if len(tokens) >= 84:
                share.coeff_4deg_quick[:, i] = [float(t) for t in tokens[:84]]
            i += 1

    fh.close()
