"""
Parse ASCII outputs written by ``write_res`` / Fortran ``write_res.f95``.

Output basenames are derived from the parameter file stem, e.g.
``myrun.par`` → ``myrun_TGM_ABD.dat``, ``myrun_model.dat``, ``myrun_ew_meas.dat``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


def output_stem(par_path: str | Path) -> str:
    """Stem used by SP_Ace for output files (``space.par`` → ``space``)."""
    p = Path(par_path)
    name = p.name
    if name.endswith(".par"):
        name = name[:-4]
    return str(p.with_name(name))


def model_dat_path(par_path: str | Path) -> Path:
    return Path(f"{output_stem(par_path)}_model.dat")


def tgm_abd_path(par_path: str | Path) -> Path:
    return Path(f"{output_stem(par_path)}_TGM_ABD.dat")


def ew_meas_path(par_path: str | Path) -> Path:
    return Path(f"{output_stem(par_path)}_ew_meas.dat")


@dataclass
class ModelColumns:
    wave: np.ndarray
    f_sp: np.ndarray
    f_sp_norm: np.ndarray
    f_model: np.ndarray
    cont: np.ndarray
    weights: np.ndarray
    sn_var: np.ndarray


def load_model_dat(path: str | Path) -> ModelColumns:
    """
    Load ``*_model.dat``: wavelength, raw flux, normalised flux, model,
    continuum, weights, S/N per pixel (last column is ``sn_var`` rounded in writer).
    """
    path = Path(path)
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 7:
        raise ValueError(f"Expected ≥7 columns in {path}, got {data.shape[1]}")
    return ModelColumns(
        wave=data[:, 0],
        f_sp=data[:, 1],
        f_sp_norm=data[:, 2],
        f_model=data[:, 3],
        cont=data[:, 4],
        weights=data[:, 5],
        sn_var=data[:, 6],
    )


def load_ew_meas(path: str | Path) -> np.ndarray:
    """Load ``*_ew_meas.dat`` as float array: wavelength, element id, EW (mÅ)."""
    path = Path(path)
    return np.loadtxt(path, comments="#")


def read_tgm_abd_table(path: str | Path) -> tuple[list[str], list[str]]:
    """
    Read ``*_TGM_ABD.dat``: first line = header tokens, second = value tokens.

    Returns (header_tokens, value_tokens). Length mismatch is preserved so
    callers can detect format drift vs Fortran.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").strip().splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path}: expected at least 2 lines")
    hdr = lines[0].split()
    vals = lines[1].split()
    return hdr, vals


def tgm_abd_as_dict(path: str | Path) -> dict[str, Any]:
    """
    Pair header and value tokens by position (best-effort).

    If lengths differ, returns ``{"_raw_header": ..., "_raw_values": ...}``.
    """
    hdr, vals = read_tgm_abd_table(path)
    if len(hdr) == len(vals):
        return dict(zip(hdr, vals))
    return {"_raw_header": hdr, "_raw_values": vals, "_mismatch": True}
