#!/usr/bin/env python3
"""Regression tests for rates, paper-style tables, and draft figures."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "study_analysis", ROOT / "scripts" / "analyze_study.py")
assert SPEC and SPEC.loader
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def row(case: str, method: str, resolution: int, error: float,
        order: int = 1) -> dict[str, str]:
    return {
        "case": case, "kind": "convergence", "method": method,
        "order": str(order), "resolution": str(resolution), "ranks": "1",
        "scale": "1", "runtime_seconds": "0.1", "profile_prefix": "",
        "dofs": str((order + 1) * resolution), "status": "ok",
        "l1_error": str(error), "l2_error": str(1.2 * error),
        "linf_error": str(2.0 * error),
        "cell_average_l2_error": str(0.2 * error),
        "downwind_l2_error": str(0.3 * error),
        "filter_applications": "0", "active_elements": "0",
    }


def expect_raises(function, exception: type[Exception]) -> None:
    try:
        function()
    except exception:
        return
    raise RuntimeError(f"expected {exception.__name__}")


def main() -> None:
    raw = [row("synthetic", "dg", 10, 0.04),
           row("synthetic", "dg", 20, 0.01)]
    summary = ANALYSIS.summarize(raw)
    finest = max(summary, key=lambda item: int(item["resolution"]))
    rate = float(finest["observed_rate"])
    if abs(rate - 2.0) > 1e-12:
        raise RuntimeError(f"expected second order, obtained {rate}")

    table_raw = []
    for method in ANALYSIS.HEADLINE:
        for resolution, error in ((16, 1.0e-2), (32, 2.5e-3)):
            table_raw.append(row("advection_1d", method, resolution, error))
            table_raw.append(row("advection_superconvergence", method,
                                 resolution, error * 0.1))
    table_summary = ANALYSIS.summarize(table_raw)

    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        figure_dir = root / "figures"
        report_dir = root / "report"
        ANALYSIS.convergence_figures(summary, figure_dir, report_dir)
        for suffix in ("pdf", "png"):
            figure = figure_dir / f"convergence_synthetic_p1.{suffix}"
            if not figure.exists() or figure.stat().st_size == 0:
                raise RuntimeError(f"missing synthetic {suffix} figure")

        table_dir = root / "tables"
        ANALYSIS.write_advection_tables(table_summary, table_dir)
        tables = [table_dir / f"advection_1d_{method}.tex"
                  for method in ANALYSIS.HEADLINE]
        contents = [path.read_text() for path in tables]
        required = (r"\multirow{2}{*}{1}", "1.00e-02", "2.00", " & --",
                    r"$e_5$ & rate")
        for content in contents:
            for token in required:
                if token not in content:
                    raise RuntimeError(f"paper table is missing {token!r}")
            if "resizebox" in content:
                raise RuntimeError("paper table must not use resizebox")
        row_lines = [[line for line in content.splitlines()
                      if line.startswith((r"\multirow", " &"))]
                     for content in contents]
        if not all(len(lines) == len(row_lines[0]) for lines in row_lines):
            raise RuntimeError("headline method tables use different row sets")

        sensor = root / "sensor.pdf"
        x = np.linspace(0.05, 0.95, 10)
        indicator = np.asarray([0.01, 0.03, 0.02, 0.15, 0.8,
                                0.9, 0.2, 0.04, 0.02, 0.01])
        ANALYSIS.kxrcf_sensor_figure(x, indicator, 0.25,
                                    indicator > 0.25, sensor)
        if not sensor.exists() or sensor.stat().st_size == 0:
            raise RuntimeError("missing spatial KXRCF style fixture")

    if ANALYSIS.mesh_label(32, 1) != "32":
        raise RuntimeError("incorrect 1D mesh label")
    if ANALYSIS.mesh_label(16, 2) != r"$16\times 16$":
        raise RuntimeError("incorrect 2D mesh label")
    expect_raises(lambda: ANALYSIS.report_asset_root("quick", "final"),
                  ValueError)

    incomplete = table_summary[:-1]
    expect_raises(lambda: ANALYSIS.assert_identical_method_rows(
        incomplete, "advection_superconvergence"), ValueError)
    print("Study-analysis style and guard tests passed.")


if __name__ == "__main__":
    main()
