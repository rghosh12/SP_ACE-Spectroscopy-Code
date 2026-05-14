"""
Multi-panel visualization for a **completed Fortran SP_Ace run** (``*_model.dat``,
``*_TGM_ABD.dat``, ``*_ew_meas.dat``).

Use this at the wrapper stage so plots reflect the **reference** Fortran pipeline
before trusting the Python port.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .parse_outputs import (
    ModelColumns,
    ew_meas_path,
    load_ew_meas,
    load_model_dat,
    model_dat_path,
    tgm_abd_as_dict,
    tgm_abd_path,
)


@dataclass
class FortranRunOutputs:
    """Paths and arrays for one run identified by ``space.par`` stem."""

    par_path: Path
    model_path: Path
    tgm_path: Path
    ew_path: Path
    model: ModelColumns
    tgm: dict[str, Any]
    ew: np.ndarray | None


def load_fortran_run_outputs(par_path: str | Path) -> FortranRunOutputs:
    """
    Load ``model.dat``, ``TGM_ABD.dat``, and ``ew_meas.dat`` (if present).

    ``par_path`` must be the same path convention SP_Ace used when writing
    outputs (stem matches ``write_res`` / Fortran).
    """
    par_path = Path(par_path).resolve()
    mp = model_dat_path(par_path)
    tp = tgm_abd_path(par_path)
    ep = ew_meas_path(par_path)
    if not mp.is_file():
        raise FileNotFoundError(f"Missing model file: {mp}")
    if not tp.is_file():
        raise FileNotFoundError(f"Missing TGM/ABD file: {tp}")
    model = load_model_dat(mp)
    tgm = tgm_abd_as_dict(tp)
    ew = load_ew_meas(ep) if ep.is_file() else None
    return FortranRunOutputs(par_path, mp, tp, ep, model, tgm, ew)


def _tgm_title_lines(tgm: dict[str, Any]) -> str:
    if tgm.get("_mismatch"):
        return "TGM_ABD: header/value length mismatch — check parse drift"
    parts = []
    for key in ("conv", "RV", "FWHM", "S/N", "chisq", "Teff", "logg", "MH"):
        if key in tgm:
            parts.append(f"{key}={tgm[key]}")
    return "  |  ".join(parts) if parts else "TGM_ABD (parsed)"


def plot_run_dashboard(
    par_path: str | Path,
    out_png: str | Path,
    *,
    weight_threshold: float = 0.01,
    title: str | None = None,
    dpi: int = 150,
) -> Path:
    """
    One figure, **four panels** (2×2):

    1. **Observed vs model** — ``f_sp_norm`` and ``f_model`` vs wavelength (Fortran reference).
    2. **Residuals** — ``f_sp_norm - f_model`` (absolute) and ``(f_sp_norm - f_model) / f_model`` (relative, axis on right).
    3. **Mask / weights** — pixel ``weights`` from ``*_model.dat``; shaded where ``weight <= weight_threshold``.
    4. **Line EW contributions** — ``*_ew_meas.dat`` stem plot (λ vs mÅ); skipped with a note if missing.

    Suptitle pulls **conv**, **chisq**, **S/N** from ``*_TGM_ABD.dat`` when keys match ``write_res`` headers.
    """
    import matplotlib.pyplot as plt

    run = load_fortran_run_outputs(par_path)
    m = run.model
    wave = m.wave
    obs = m.f_sp_norm
    mod = m.f_model
    wts = m.weights
    low_w = wts <= float(weight_threshold)

    with np.errstate(divide="ignore", invalid="ignore"):
        resid_abs = obs - mod
        resid_rel = resid_abs / mod

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    ax00, ax01 = axes[0]
    ax10, ax11 = axes[1]

    # (0,0) spectrum
    ax00.plot(wave, obs, lw=0.55, label="obs (f_sp_norm)", alpha=0.9)
    ax00.plot(wave, mod, lw=0.55, label="model (f_model)", alpha=0.9)
    ax00.set_ylabel("flux (norm.)")
    ax00.legend(loc="upper right", fontsize=8)
    ax00.grid(True, alpha=0.25)
    ax00.set_title("Observed vs best-fit model (Fortran outputs)")

    # (0,1) residuals
    ax01.plot(wave, resid_abs, lw=0.5, color="C0", label="obs − model")
    ax01.axhline(0.0, color="k", lw=0.4)
    ax01.set_ylabel("Δflux", color="C0")
    ax01.tick_params(axis="y", labelcolor="C0")
    ax01r = ax01.twinx()
    ax01r.plot(wave, resid_rel, lw=0.45, color="C1", alpha=0.75, label="(obs−model)/model")
    ax01r.set_ylabel("(obs − model) / model", color="C1")
    ax01r.tick_params(axis="y", labelcolor="C1")
    ax01.set_title("Residuals")
    ax01.grid(True, alpha=0.2)

    # (1,0) weights / mask
    ax10.fill_between(wave, 0, 1, where=low_w, color="red", alpha=0.12, label=f"weight ≤ {weight_threshold}")
    ax10.plot(wave, wts, lw=0.45, color="C2", label="pixel weight")
    ax10.set_ylim(-0.05, 1.05)
    ax10.set_xlabel("wavelength (Å)")
    ax10.set_ylabel("weight")
    ax10.legend(loc="lower right", fontsize=7)
    ax10.set_title("Mask proxy (weights from model.dat)")
    ax10.grid(True, alpha=0.25)

    # (1,1) EW lines
    ew = run.ew
    if ew is not None and ew.size > 0:
        if ew.ndim == 1:
            ew = ew.reshape(1, -1)
        wl = ew[:, 0]
        ewm = ew[:, 2]
        ax11.stem(wl, ewm, linefmt="C3-", markerfmt="C3o", basefmt=" ")
        ax11.set_xlabel("wavelength (Å)")
        ax11.set_ylabel("EW (mÅ)")
        ax11.set_title("Measured EWs (selected lines)")
        ax11.grid(True, alpha=0.25)
    else:
        ax11.text(0.5, 0.5, "No ew_meas.dat (or empty)", ha="center", va="center", transform=ax11.transAxes)
        ax11.set_axis_off()

    st = title or _tgm_title_lines(run.tgm)
    fig.suptitle(st, fontsize=10)

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    return out_png
