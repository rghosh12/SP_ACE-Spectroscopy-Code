"""GCOG path strings must match ``read_GCOG.py`` / ``read_GCOG.f95``."""

from fortran_tools.gcog_filenames import gcog_4deg_filename, gcog_4deg_path


def test_gcog_4deg_filename_solar_metallicity():
    assert gcog_4deg_filename(5750.0, 4.5, 0.0) == "ap00-t5750-g45-4degpoly.dat"


def test_gcog_4deg_filename_negative_metallicity():
    assert gcog_4deg_filename(5750.0, 4.5, -0.5) == "am05-t5750-g45-4degpoly.dat"


def test_gcog_4deg_path_joins_directory():
    p = gcog_4deg_path("/data/GCOG", 5000.0, 2.5, 0.1)
    assert p.startswith("/data/GCOG/")
    assert p.endswith("4degpoly.dat")
