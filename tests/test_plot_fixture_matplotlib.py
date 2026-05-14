"""Plot diagnostics on bundled fixture (issue #2)."""

from pathlib import Path

import pytest

pytest.importorskip("matplotlib")

from fortran_tools.diagnostics import plot_model_diagnostic


def test_plot_model_fixture_writes_png(tmp_path):
    fix = Path(__file__).resolve().parent / "fixtures" / "example_model_snippet.dat"
    out = tmp_path / "diag.png"
    plot_model_diagnostic(fix, out, title="fixture")
    assert out.is_file()
    assert out.stat().st_size > 400
