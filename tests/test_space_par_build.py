"""Minimal ``space.par`` writer."""

import tempfile
from pathlib import Path

from fortran_tools.space_par_build import write_minimal_space_par


def test_write_minimal_space_par_required_keys():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "x.par"
        write_minimal_space_par(
            p,
            obs_sp_file="obs.dat",
            gcoglib="/g/",
            fwhm=2.0,
            wave_lims=[(4500.0, 6800.0)],
        )
        text = p.read_text(encoding="utf-8")
        assert "obs_sp_file" in text
        assert "GCOGlib" in text
        assert "fwhm" in text
        assert "wave_lims" in text
        assert "obs.dat" in text
