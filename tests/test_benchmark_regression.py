"""Reference-star regression placeholder (issue #7)."""

import os

import pytest


@pytest.mark.benchmark
def test_benchmark_cases_when_enabled():
    if os.environ.get("BENCHMARK_DATA") != "1":
        pytest.skip("Set BENCHMARK_DATA=1 and add benchmarks/cases/... (see benchmarks/REFERENCE_STAR_TEMPLATE.md)")
    cases = os.environ.get("BENCHMARK_CASES_DIR", "benchmarks/cases")
    assert os.path.isdir(cases), f"Expected directory {cases}"
    assert any(os.scandir(cases)), f"No cases under {cases}"
