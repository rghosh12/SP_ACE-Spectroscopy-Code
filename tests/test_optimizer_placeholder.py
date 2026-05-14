import pytest

from fortran_tools.optimizer_placeholder import (
    OptimizerNotValidatedError,
    run_scipy_least_squares_fit,
)


def test_optimizer_placeholder_raises():
    with pytest.raises(OptimizerNotValidatedError):
        run_scipy_least_squares_fit()
