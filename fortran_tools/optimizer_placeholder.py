"""
Reserved for a future Python-side optimizer (e.g. scipy.optimize.least_squares).

**Do not implement or wire this in until:**
  - GCOG ingestion matches Fortran (see ``parity/README.md``),
  - synthetic spectra from Python match Fortran,
  - and regression tests on reference stars pass.

The authoritative fitting loop remains the Fortran code in ``Fortran-95_code/``.
"""


class OptimizerNotValidatedError(RuntimeError):
    pass


def run_scipy_least_squares_fit(*_args, **_kwargs):
    """Placeholder — raises until parity stages are complete."""
    raise OptimizerNotValidatedError(
        "Python optimizer is intentionally disabled until Fortran parity "
        "is demonstrated for GCOG load, model spectrum, and benchmarks."
    )
