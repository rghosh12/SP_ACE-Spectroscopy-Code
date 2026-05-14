#!/usr/bin/env python3
"""
Plot SP_Ace ``*_model.dat`` diagnostics (roadmap name; same as ``plot_sp_ace_diagnostics.py``).

See GitHub issue #2.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def main() -> int:
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location(
        "_plot_sp_ace_diagnostics", here / "plot_sp_ace_diagnostics.py"
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load plot_sp_ace_diagnostics.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return int(mod.main())


if __name__ == "__main__":
    raise SystemExit(main())
