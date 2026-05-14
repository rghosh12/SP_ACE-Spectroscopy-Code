#!/usr/bin/env python3
"""Print a short text summary from SP_Ace ``*_TGM_ABD.dat`` and ``*_model.dat``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fortran_tools.parse_outputs import (
    load_model_dat,
    model_dat_path,
    read_tgm_abd_table,
    tgm_abd_path,
)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Summarise SP_Ace outputs for a given space.par path (same stem)."
    )
    ap.add_argument("par_file", help="Path to parameter file (e.g. myrun.par)")
    args = ap.parse_args()
    par = Path(args.par_file)

    tgm_p = tgm_abd_path(par)
    mod_p = model_dat_path(par)
    if not tgm_p.is_file():
        print(f"Missing {tgm_p}", file=sys.stderr)
        return 1
    hdr, vals = read_tgm_abd_table(tgm_p)
    print("=== TGM_ABD ===")
    for h, v in zip(hdr, vals):
        print(f"  {h}: {v}")
    if len(hdr) != len(vals):
        print(f"  (warning: {len(hdr)} header tokens vs {len(vals)} value tokens)", file=sys.stderr)

    if mod_p.is_file():
        m = load_model_dat(mod_p)
        resid = (m.f_sp_norm - m.f_model) / m.f_model
        print("\n=== model.dat ===")
        print(f"  pixels: {len(m.wave)}")
        print(f"  wavelength range: {m.wave[0]:.2f} – {m.wave[-1]:.2f} Å")
        print(f"  median |(f_norm - f_mod)/f_mod|: {float(np.median(np.abs(resid))):.5f}")
    else:
        print(f"\n(no {mod_p} — model summary skipped)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
