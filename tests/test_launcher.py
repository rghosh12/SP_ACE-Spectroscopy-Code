"""Fortran launcher — executable discovery and failure modes."""

import os

import pytest

from fortran_tools.launcher import find_sp_ace_executable, run_sp_ace


def test_run_sp_ace_raises_without_executable(monkeypatch, tmp_path):
    monkeypatch.delenv("SP_ACE_EXE", raising=False)
    par = tmp_path / "space.par"
    par.write_text("dummy", encoding="utf-8")
    with pytest.raises(FileNotFoundError, match="SP_Ace Fortran executable"):
        run_sp_ace(par)


def test_find_executable_explicit_missing(tmp_path):
    assert find_sp_ace_executable(tmp_path / "nonexistent.exe") is None


def test_find_executable_from_env(tmp_path, monkeypatch):
    fake = tmp_path / "space.exe"
    fake.write_bytes(b"")
    monkeypatch.setenv("SP_ACE_EXE", str(fake))
    assert find_sp_ace_executable() == fake
