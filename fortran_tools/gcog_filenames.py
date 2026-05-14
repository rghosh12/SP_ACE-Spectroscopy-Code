"""
GCOG file naming — must stay identical to ``read_GCOG.f95`` / ``read_GCOG.py``.

Used for parity tests and tooling without importing the full SP_Ace Python port.
"""

from __future__ import annotations


def gcog_4deg_filename(temp: float, logg: float, met: float) -> str:
    """
    Return the filename *segment* ``{met_label}{temp_label}{logg_label}``,
    matching ``load_GCOG_4deg`` in ``ConvertedPython_code/read_GCOG.py``.
    """
    logg_label = f"g{int(round(logg * 10)):02d}-4degpoly.dat"
    temp_label = f"t{int(round(temp)):04d}-"
    if met < 0:
        met_label = f"am{int(round(-met * 10)):02d}-"
    else:
        met_label = f"ap{int(round(met * 10)):02d}-"
    return f"{met_label}{temp_label}{logg_label}"


def gcog_4deg_path(gcoglib: str, temp: float, logg: float, met: float) -> str:
    """Full path to the 4th-degree GCOG coefficient file for one grid node."""
    root = gcoglib.rstrip("/")
    return f"{root}/{gcog_4deg_filename(temp, logg, met)}"


def gcog_quick_poly_path(gcoglib: str) -> str:
    """Path to ``space_6degpoly.dat`` (quick coefficients)."""
    return f"{gcoglib.rstrip('/')}/space_6degpoly.dat"
