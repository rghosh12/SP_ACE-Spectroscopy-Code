"""
Run the compiled Fortran ``PROGRAM space`` binary with a parameter file.

Set ``SP_ACE_EXE`` to the full path of the executable, or pass ``executable=``
explicitly. The Fortran program reads ``getarg(1)`` as the ``space.par``
path (default ``space.par`` if empty).
"""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class RunResult:
    returncode: int
    stdout: str
    stderr: str
    par_path: Path
    cwd: Path


def find_sp_ace_executable(explicit: Optional[str | Path] = None) -> Optional[Path]:
    """
    Resolve the SP_Ace Fortran executable.

    Search order:
    1. ``explicit`` argument
    2. Environment variable ``SP_ACE_EXE``
    3. ``./space`` / ``./space.exe`` under ``Fortran-95_code`` (if present)
    """
    if explicit is not None:
        p = Path(explicit).expanduser()
        return p if p.is_file() else None

    env = os.environ.get("SP_ACE_EXE")
    if env:
        p = Path(env).expanduser()
        return p if p.is_file() else None

    repo = Path(__file__).resolve().parents[1] / "Fortran-95_code"
    for name in ("space.exe", "space"):
        cand = repo / name
        if cand.is_file():
            return cand
    return None


def run_sp_ace(
    par_path: str | Path,
    *,
    executable: Optional[str | Path] = None,
    cwd: Optional[str | Path] = None,
    timeout: Optional[float] = None,
    env: Optional[dict[str, str]] = None,
) -> RunResult:
    """
    Invoke the Fortran SP_Ace binary.

    Parameters
    ----------
    par_path
        Path to ``space.par`` (or equivalent). Passed as argv[1] to Fortran.
    executable
        Override path to the ``space`` binary; otherwise ``find_sp_ace_executable``.
    cwd
        Working directory for the child process. Defaults to the parameter
        file's parent so relative paths inside ``space.par`` resolve like Fortran.
    timeout
        Optional subprocess timeout in seconds.
    """
    par_path = Path(par_path).resolve()
    exe = find_sp_ace_executable(executable)
    if exe is None:
        raise FileNotFoundError(
            "SP_Ace Fortran executable not found. Build Fortran-95_code (see Makefile) "
            "and set SP_ACE_EXE to the binary path, or pass executable=."
        )

    work = Path(cwd).resolve() if cwd is not None else par_path.parent
    # Fortran uses getarg(1); use name when cwd matches par directory
    arg1 = par_path.name if par_path.parent == work else str(par_path)

    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)

    proc = subprocess.run(
        [str(exe), arg1],
        cwd=str(work),
        capture_output=True,
        text=True,
        timeout=timeout,
        env=merged_env,
        check=False,
    )
    return RunResult(
        returncode=proc.returncode,
        stdout=proc.stdout or "",
        stderr=proc.stderr or "",
        par_path=par_path,
        cwd=work,
    )
