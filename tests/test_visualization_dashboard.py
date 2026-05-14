"""Fortran-run dashboard (multi-panel) from fixtures."""

import shutil
from pathlib import Path

import pytest

pytest.importorskip("matplotlib")

from fortran_tools.visualization import plot_run_dashboard


def test_dashboard_from_fixture_run(tmp_path):
    par = tmp_path / "fixturerun.par"
    par.write_text("# dummy par\n", encoding="utf-8")
    stem = tmp_path / "fixturerun"
    shutil.copy(
        Path(__file__).parent / "fixtures" / "example_model_snippet.dat",
        stem.parent / f"{stem.name}_model.dat",
    )
    shutil.copy(
        Path(__file__).parent / "fixtures" / "minimal_TGM_ABD.dat",
        stem.parent / f"{stem.name}_TGM_ABD.dat",
    )
    shutil.copy(
        Path(__file__).parent / "fixtures" / "minimal_ew_meas.dat",
        stem.parent / f"{stem.name}_ew_meas.dat",
    )
    out = tmp_path / "dash.png"
    plot_run_dashboard(par, out)
    assert out.is_file() and out.stat().st_size > 2000
