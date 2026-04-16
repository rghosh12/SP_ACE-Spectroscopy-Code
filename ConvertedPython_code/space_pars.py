"""
space_pars.py — SP_Ace parameter-file reader
Translated from space_pars.f95 (Corrado Boeche, 2016, GPLv3)
"""

import sys
import numpy as np
import share
import data_lib

# ── Module-level parameters ────────────────────────────────────────────────────
GCOGlib: str = ""
llist: str = ""
llist_rej: str = ""
null_val: str = "null"
fwhm: float = 0.0
w_inf = np.zeros(10)
w_sup = np.zeros(10)
sn_ratio: float = 100.0
rv_ini: float = 0.0
TGM_force = np.array([-99., -99.])
N_w_int: int = 0
error_est: bool = False
flag_rej: bool = False
sn_flag: bool = True
flag_norm: bool = True
flag_alpha: bool = False
flag_ABD_loop: bool = True
flag_salaris_MH: bool = True


def stop_local(text):
    print(f"{text}, SP_Ace stops!")
    with open("space_msg.txt", "w") as f:
        f.write(f"### {share.obs_sp_file}\n")
        f.write(f"    {text}, SP_Ace stops!\n")
    sys.exit(1)


def read_space_pars(file_pars):
    """
    Read the SP_Ace parameter file ('space.par') and populate
    module-level variables and share.* globals.

    Parameters
    ----------
    file_pars : str — path to the parameter file
    """
    global GCOGlib, llist, llist_rej, null_val, fwhm, w_inf, w_sup
    global sn_ratio, rv_ini, TGM_force, N_w_int
    global error_est, flag_rej, sn_flag, flag_norm
    global flag_alpha, flag_ABD_loop, flag_salaris_MH

    n_ele_symb = len(data_lib.ELE_symb)
    share.n_ele_symb = n_ele_symb

    # Allowed elements to output (besides Fe=26)
    ele_allowed = np.array(
        [6, 11, 12, 13, 14, 20, 21, 22, 23, 24, 25, 27, 28, 29,
         30, 39, 40, 56, 57, 58, 60], dtype=np.int16
    )

    # Defaults
    sn_ratio = 100.0
    rv_ini = 0.0
    sn_flag = True
    TGM_force = np.array([-99., -99.])
    error_est = False
    flag_rej = False
    flag_alpha = False
    flag_ABD_loop = True
    flag_salaris_MH = True
    flag_norm = True
    null_val = "null"
    share.norm_rad = 30.0
    share.sigma = 0.4
    n_ele = 0
    flag_ele2write = False
    flag_keys = [False, False, False, False]
    mask_2write = np.zeros(len(ele_allowed), dtype=bool)
    ele = np.zeros(101, dtype=np.int16)

    try:
        fh = open(file_pars, "r")
    except OSError:
        print(f"SPACE cannot open {file_pars} STOP!")
        sys.exit(1)

    for line in fh:
        line = line.rstrip("\n")
        if not line.strip() or line.strip().startswith("#"):
            continue
        parts = line.split(None, 1)
        keyword = parts[0]
        rest = parts[1].strip() if len(parts) > 1 else ""

        if keyword == "obs_sp_file":
            flag_keys[0] = True
            share.obs_sp_file = rest.strip()
        elif keyword == "sn_ratio":
            try:
                v = float(rest)
                if v >= 1:
                    sn_ratio = v
                    sn_flag = False
                else:
                    print("sn_ratio wrong! SP_Ace estimates it itself!")
            except ValueError:
                print("sn_ratio wrong! SP_Ace estimates it itself!")
        elif keyword == "GCOGlib":
            flag_keys[1] = True
            GCOGlib = rest.strip()
            llist = GCOGlib.rstrip("/") + "/linelist.dat"
        elif keyword == "llist_rej":
            llist_rej = rest.strip()
            flag_rej = True
        elif keyword == "fwhm":
            flag_keys[2] = True
            fwhm = float(rest)
            share.sigma = fwhm / 2.35
        elif keyword == "wave_lims":
            flag_keys[3] = True
            tokens = rest.split()
            if len(tokens) % 2 != 0:
                stop_local("wave_lims must have an even number of values")
            pairs = [(float(tokens[i]), float(tokens[i + 1]))
                     for i in range(0, len(tokens), 2)]
            N_w_int = len(pairs)
            for i, (lo, hi) in enumerate(pairs):
                if lo > hi:
                    stop_local("wrong wavelength limits!")
                w_inf[i] = lo
                w_sup[i] = hi
        elif keyword == "RV_ini":
            try:
                v = float(rest)
                if abs(v) <= 9999:
                    rv_ini = v
                else:
                    print("Rv beyond the limit of +-9999 km/sec! SP_Ace starts from Rv=0!")
                    rv_ini = 0.0
            except ValueError:
                rv_ini = 0.0
        elif keyword == "ele2write":
            flag_ele2write = True
            tokens = rest.split()
            n_ele = 0
            for tok in tokens:
                try:
                    v = int(tok)
                    ele[n_ele] = v
                    n_ele += 1
                except ValueError:
                    pass
        elif keyword == "T_force":
            v = float(rest)
            if v < 3600 or v > 7400:
                stop_local("T_force must be 3600<T_force<7400")
            TGM_force[0] = v
        elif keyword == "G_force":
            v = float(rest)
            if v < 0.2 or v > 5.0:
                stop_local("G_force must be 0.2<G_force<5.4")
            TGM_force[1] = v
        elif keyword == "error_est":
            error_est = True
        elif keyword == "no_norm":
            flag_norm = False
        elif keyword == "alpha":
            flag_alpha = True
        elif keyword == "ABD_loop":
            flag_ABD_loop = True
        elif keyword == "Salaris_MH":
            flag_ABD_loop = True
            flag_salaris_MH = True
        elif keyword == "norm_rad":
            share.norm_rad = float(rest)
        elif keyword == "null_value":
            v = rest.strip()
            if v not in ("-9.99", "null", "NaN"):
                stop_local("null_value is wrong!")
            null_val = v
        else:
            stop_local(f"There is an unrecognized keyword: {keyword}")

    fh.close()

    if not all(flag_keys):
        stop_local("one necessary keyword is missing,")

    # Build ele2write
    if flag_alpha:
        share.ele2write = np.array([26, 94], dtype=np.int16)
    elif flag_ele2write:
        for i in range(n_ele):
            mask_2write |= (ele[i] == ele_allowed)
        selected = ele_allowed[mask_2write]
        share.ele2write = np.concatenate([[26], selected]).astype(np.int16)
    else:
        share.ele2write = np.concatenate([[26], ele_allowed]).astype(np.int16)
