# Parity work: Python vs Fortran

This folder documents **line-by-line** and **numerical** parity between:

- ``Fortran-95_code/`` (reference)
- ``ConvertedPython_code/`` (port)
- ``fortran_tools/`` (I/O helpers that must match Fortran file conventions only)

## GCOG coefficient files

1. **Filename convention** — Implemented in ``fortran_tools/gcog_filenames.py`` and
   unit-tested against the logic in ``ConvertedPython_code/read_GCOG.py``.
2. **Row layout** — Each selected line in ``linelist.dat`` corresponds to one row
   in the GCOG ``*4degpoly.dat`` file: **70** floats for the 4th-degree load
   (**84** for ``space_6degpoly.dat``). Any Python reader must preserve **token
   order** and **masking** from ``read_GCOG.f95``.

## Suggested workflow

1. Run Fortran on a fixed ``space.par`` with a small wavelength interval.
2. Dump intermediate arrays from Fortran (temporary ``print`` or binary dump)
   *or* compare final ``*_model.dat`` / ``*_TGM_ABD.dat`` after the Python port
   reaches the same stage.
3. Extend ``tests/`` with golden files once you have permission to redistribute
   a minimal GCOG snippet (copyright may restrict redistributing full libraries).

## Optimizer

Do **not** replace the Fortran LM / χ² loop with SciPy until Stages 3–4 in the
README roadmap are satisfied. ``fortran_tools/optimizer_placeholder.py`` exists
only as an explicit guardrail.
