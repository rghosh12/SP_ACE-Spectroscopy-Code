#!/usr/bin/env python3
"""Plot diagnostics from ``*_model.dat`` produced by SP_Ace."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fortran_tools.diagnostics import plot_model_diagnostic


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot SP_Ace model diagnostic PNG.")
    ap.add_argument("model_dat", help="Path to *_model.dat")
    ap.add_argument("-o", "--output", default="sp_ace_diagnostic.png", help="Output PNG path")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()
    plot_model_diagnostic(args.model_dat, args.output, title=args.title)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
