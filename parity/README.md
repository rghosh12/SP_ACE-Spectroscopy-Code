# Parity work: Python vs Fortran

This folder documents **line-by-line** and **numerical** parity between:

- ``Fortran-95_code/`` (reference)
- ``ConvertedPython_code/`` (port)
- ``fortran_tools/`` (I/O helpers that must match Fortran file conventions only)

## GCOG coefficient files

1. **Filename convention** — Implemented in ``fortran_tools/gcog_filenames.py`` and
   unit-tested against the logic in ``ConvertedPython_code/read_GCOG.py``.
2. **Row layout** — See **[`GCOG_LIBRARY_FORMAT.md`](GCOG_LIBRARY_FORMAT.md)** for file paths, 70 vs 84 floats, and in-memory targets. Structural tests live in ``tests/test_gcog_filenames.py`` and ``tests/test_gcog_line_tokens.py``.

## Model spectrum parity

See **[`MODEL_PARITY_PLAN.md`](MODEL_PARITY_PLAN.md)** for the Fortran-vs-Python model comparison plan ([GitHub #4](https://github.com/rghosh12/SP_ACE-Spectroscopy-Code/issues/4)).

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
