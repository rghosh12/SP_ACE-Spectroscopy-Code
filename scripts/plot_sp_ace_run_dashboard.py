#!/usr/bin/env python3
"""
Build a multi-panel dashboard PNG from a completed Fortran SP_Ace run (same stem as ``space.par``).

Example::

    python scripts/plot_sp_ace_run_dashboard.py C:\\work\\myrun.par -o C:\\work\\myrun_dashboard.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fortran_tools.visualization import plot_run_dashboard


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Plot Fortran SP_Ace run: obs vs model, residuals, weights, EWs."
    )
    ap.add_argument("par_file", help="Path to space.par (outputs share the same stem)")
    ap.add_argument("-o", "--output", default="sp_ace_dashboard.png", help="Output PNG path")
    ap.add_argument(
        "--weight-threshold",
        type=float,
        default=0.01,
        help="Highlight pixels with model.dat weight at or below this value",
    )
    ap.add_argument("--title", default=None)
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()
    out = plot_run_dashboard(
        args.par_file,
        args.output,
        weight_threshold=args.weight_threshold,
        title=args.title,
        dpi=args.dpi,
    )
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
