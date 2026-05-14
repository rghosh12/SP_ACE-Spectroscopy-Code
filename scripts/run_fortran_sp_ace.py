#!/usr/bin/env python3
"""
Run the Fortran SP_Ace ``space`` program with a parameter file.

Examples
--------
  set SP_ACE_EXE=C:\\path\\to\\space.exe
  python scripts/run_fortran_sp_ace.py path\\to\\space.par

  python scripts/run_fortran_sp_ace.py --write-par work\\run.par \\
      --obs-sp ..\\fast_spectrum.txt --gcoglib C:\\GCOG\\lib --fwhm 2.0 \\
      --wave-lims 4500 6800
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fortran_tools.launcher import find_sp_ace_executable, run_sp_ace
from fortran_tools.space_par_build import write_minimal_space_par


def main() -> int:
    p = argparse.ArgumentParser(description="Run Fortran SP_Ace (PROGRAM space).")
    p.add_argument("par_file", nargs="?", help="Parameter file (space.par)")
    p.add_argument("--write-par", metavar="PATH", help="Write minimal space.par to PATH then run it")
    p.add_argument("--obs-sp", help="obs_sp_file path for --write-par")
    p.add_argument("--gcoglib", help="GCOG library directory for --write-par")
    p.add_argument("--fwhm", type=float, help="FWHM (Å) for --write-par")
    p.add_argument("--wave-lims", nargs="+", type=float, metavar="L", help="Pairs: lo hi lo hi ...")
    p.add_argument("--exe", help="Path to Fortran ``space`` binary (else SP_ACE_EXE)")
    p.add_argument("--cwd", help="Working directory (default: par file directory)")
    p.add_argument("--timeout", type=float, default=None, help="Subprocess timeout (s)")
    args = p.parse_args()

    par: Path | None = None
    if args.write_par:
        if not all([args.obs_sp, args.gcoglib, args.fwhm is not None, args.wave_lims]):
            p.error("--write-par requires --obs-sp, --gcoglib, --fwhm, and --wave-lims")
        lims = list(zip(args.wave_lims[::2], args.wave_lims[1::2]))
        par = write_minimal_space_par(
            args.write_par,
            obs_sp_file=args.obs_sp,
            gcoglib=args.gcoglib,
            fwhm=args.fwhm,
            wave_lims=lims,
        )
    elif args.par_file:
        par = Path(args.par_file)
    else:
        p.error("Provide par_file or --write-par ...")

    exe = find_sp_ace_executable(args.exe)
    if exe is None:
        print(
            "No executable found. Build Fortran-95_code (see Makefile) and set SP_ACE_EXE.",
            file=sys.stderr,
        )
        return 2

    cwd = Path(args.cwd) if args.cwd else None
    res = run_sp_ace(par, executable=args.exe, cwd=cwd, timeout=args.timeout)
    sys.stdout.write(res.stdout)
    sys.stderr.write(res.stderr)
    return int(res.returncode != 0)


if __name__ == "__main__":
    raise SystemExit(main())
