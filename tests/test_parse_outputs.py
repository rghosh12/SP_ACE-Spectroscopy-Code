"""Output parsers for ``*_model.dat`` / ``*_TGM_ABD.dat``."""

import tempfile
from pathlib import Path

import pytest

from fortran_tools.parse_outputs import (
    load_model_dat,
    model_dat_path,
    output_stem,
    read_tgm_abd_table,
    tgm_abd_as_dict,
)


def test_output_stem():
    assert output_stem(Path("/tmp/run.par")).replace("\\", "/").endswith("/tmp/run")


def test_model_dat_path():
    p = model_dat_path("/tmp/run.par")
    assert str(p).replace("\\", "/").endswith("run_model.dat")


def test_load_model_dat_roundtrip():
    txt = " 4500.050  0.94800  0.91400  0.92000  1.00000  1.00   100\n"
    with tempfile.TemporaryDirectory() as d:
        fp = Path(d) / "m.dat"
        fp.write_text(txt, encoding="utf-8")
        m = load_model_dat(fp)
        assert m.wave[0] == pytest.approx(4500.05)
        assert m.f_sp_norm[0] == pytest.approx(0.914)


def test_read_tgm_abd_table():
    with tempfile.TemporaryDirectory() as d:
        fp = Path(d) / "t.dat"
        fp.write_text("  conv   Teff\n   1  5000\n", encoding="utf-8")
        h, v = read_tgm_abd_table(fp)
        assert h == ["conv", "Teff"]
        assert v == ["1", "5000"]
        dct = tgm_abd_as_dict(fp)
        assert dct["conv"] == "1" and dct["Teff"] == "5000"
