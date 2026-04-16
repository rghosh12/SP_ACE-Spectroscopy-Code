"""
error.py — SP_Ace error-handling routines
Translated from error.f95 (Corrado Boeche, 2016, GPLv3)
"""

import sys
import numpy as np
import share
import data_lib


def error_msg(conv, text):
    """
    Write null results to the output file and exit.

    Parameters
    ----------
    conv : int   — convergence flag
    text : str   — error description
    """
    import space_pars  # imported here to avoid circular imports at module load

    print(f"{text}  SP_Ace exit with no results!")

    null_str = {
        "NaN":   "  NaN",
        "null":  " null",
        "-9.99": "-9.99",
    }.get(space_pars.null_val, " null")

    header = ""
    line = ""

    # convergence
    header += "  conv"
    line += f"   {int(conv):2d}"

    # RV
    header += "      RV"
    line += f"   {null_str}"

    # FWHM
    header += "    FWHM"
    line += f"   {null_str}"

    # S/N
    header += "     S/N"
    line += f"   {null_str}"

    # chisq
    header += "    chisq"
    line += f"    {null_str}"

    # Teff
    header += "  Teff   T_l   T_h"
    line += f" {null_str} {null_str} {null_str}"

    # logg
    header += "   logg   L_l   L_h"
    line += f"  {null_str} {null_str} {null_str}"

    # [M/H]
    header += "    MH  MH_l  MH_h"
    line += f" {null_str} {null_str} {null_str}"

    # elements
    ele2write = share.ele2write if share.ele2write is not None else []
    for ele_int in ele2write:
        sym = data_lib.ELE_symb[int(ele_int) - 1]  # 1-indexed in Fortran
        string = f"{sym}_l  {sym}_h {sym}_N"
        header += f"   {sym}  {string}"
        line += f" {null_str:6s} {null_str:6s} {null_str:6s} {0:4d}"

    file_out = space_pars.space_params_file
    if file_out.endswith(".par"):
        file_out = file_out[:-4]
    file_out += "_TGM_ABD.dat"

    with open(file_out, "w") as f:
        f.write(header.strip() + "\n")
        f.write(line.strip() + "\n")

    print("SP_Ace exits with no solution.")
    sys.exit(1)


def stop_msg(text):
    """
    Print an error message, write it to space_msg.txt, and exit.

    Parameters
    ----------
    text : str
    """
    print(f"{text}, SP_Ace stops!")
    line = f"{share.obs_sp_file} {text}, SP_Ace stops!"
    with open("space_msg.txt", "w") as f:
        f.write(line + "\n")
    sys.exit(1)
