"""
Fortran SP_Ace orchestration, output I/O, and diagnostics.

This package does not reimplement the physics; it wraps the compiled
``space`` program, parses its text outputs, and plots diagnostics. Python
parity with Fortran lives in ``ConvertedPython_code/`` and must be
validated separately (see ``parity/README.md``).
"""

from .launcher import RunResult, find_sp_ace_executable, run_sp_ace
from .space_par_build import write_minimal_space_par

__all__ = [
    "RunResult",
    "find_sp_ace_executable",
    "run_sp_ace",
    "write_minimal_space_par",
]
