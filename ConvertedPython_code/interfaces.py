"""
interfaces.py — SP_Ace interface declarations
Translated from interfaces.f95 (Corrado Boeche, 2016, GPLv3)

In Fortran, this module existed solely to provide explicit interfaces for
external subroutines so the compiler could type-check call sites.  Python
does not need explicit interface declarations — all function signatures are
resolved at runtime.  This file is therefore a no-op stub kept for
structural completeness.

The two subroutines declared here are implemented in:
  - fit_cont.py   → fit_cont()
  - write_res.py  → write_res()
"""
