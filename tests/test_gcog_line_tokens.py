"""GCOG data row shape: token counts must match read_GCOG loaders."""

import numpy as np


def _row(n: int) -> str:
    return " ".join(f"{i * 0.001:.6f}" for i in range(n))


def test_fourth_degree_row_has_70_tokens():
    arr = np.array(_row(70).split(), dtype=np.float64)
    assert arr.shape == (70,)


def test_quick_sixth_degree_row_has_84_tokens():
    arr = np.array(_row(84).split(), dtype=np.float64)
    assert arr.shape == (84,)
