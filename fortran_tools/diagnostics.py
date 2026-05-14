"""
Diagnostic plots for SP_Ace ``*_model.dat`` outputs (observed vs model, residuals).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .parse_outputs import load_model_dat


def plot_model_diagnostic(
    model_dat: str | Path,
    out_png: str | Path,
    title: str | None = None,
) -> Path:
    """
    Save a two-panel figure: spectrum + model, and normalised residuals.

    Requires matplotlib.
    """
    import matplotlib.pyplot as plt

    model_dat = Path(model_dat)
    out_png = Path(out_png)
    m = load_model_dat(model_dat)

    wave = m.wave
    obs = m.f_sp_norm
    mod = m.f_model
    sig = np.where(m.sn_var > 0, m.sn_var.astype(float), np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        resid = (obs - mod) / mod
        resid_sigma = (obs - mod) * sig / mod

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, constrained_layout=True)
    axes[0].plot(wave, obs, lw=0.6, label="f_sp_norm", alpha=0.85)
    axes[0].plot(wave, mod, lw=0.6, label="f_model", alpha=0.85)
    axes[0].set_ylabel("flux (norm.)")
    axes[0].legend(loc="upper right", fontsize=8)
    axes[0].grid(True, alpha=0.3)
    if title:
        axes[0].set_title(title)

    axes[1].plot(wave, resid, lw=0.5, color="C2", label="(obs-mod)/mod")
    axes[1].plot(wave, resid_sigma, lw=0.5, color="C3", alpha=0.7, label="~resid×(S/N)")
    axes[1].axhline(0.0, color="k", lw=0.5)
    axes[1].set_xlabel("wavelength (Å)")
    axes[1].set_ylabel("residual")
    axes[1].legend(loc="upper right", fontsize=8)
    axes[1].grid(True, alpha=0.3)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return out_png
