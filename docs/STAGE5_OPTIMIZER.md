# Stage 5: Python optimizer ([issue #5](https://github.com/rghosh12/SP_ACE-Spectroscopy-Code/issues/5))

**Blocked by:** [issue #3](https://github.com/rghosh12/SP_ACE-Spectroscopy-Code/issues/3) (GCOG parity) and [issue #4](https://github.com/rghosh12/SP_ACE-Spectroscopy-Code/issues/4) (model spectrum parity).

## Policy

- The **Fortran** `minimize.f95` / LM loop remains authoritative until Stages 3–4 pass automated checks.
- `fortran_tools/optimizer_placeholder.py` raises **`OptimizerNotValidatedError`** if a SciPy-based path is invoked early — this is intentional.

## When Stages 3–4 are satisfied

1. Wrap the same χ² residuals the Fortran LM sees (vector-valued for `least_squares`) or match the scalar objective exactly.
2. Compare **final** `TGM` / `ABD` to Fortran on reference cases ([issue #7](https://github.com/rghosh12/SP_ACE-Spectroscopy-Code/issues/7)).
3. Document runtime and any acceptable numerical drift.
