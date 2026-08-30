#!/usr/bin/env python3
"""Validate manifest expansion and the draft execution/report guards."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RUNNER = load("study_runner", ROOT / "scripts" / "run_study.py")
FINAL_CHECK = load("final_check", ROOT / "scripts" / "check_final_report.py")


def main() -> None:
    document, profile = RUNNER.load_manifest("full")
    if profile.get("status") != "draft":
        raise RuntimeError("full profile should remain draft in this phase")
    runs = RUNNER.expand_cases(document, profile, set())
    if not runs:
        raise RuntimeError("full design expanded to no rows")

    with tempfile.TemporaryDirectory() as temporary:
        plan = Path(temporary) / "full.csv"
        seconds = RUNNER.write_study_plan(document, profile, runs, plan)
        with plan.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        if len(rows) != len(runs):
            raise RuntimeError("dry-run plan lost expanded rows")
        if any(row["case_status"] != "proposed" for row in rows):
            raise RuntimeError("candidate cases must remain proposed")
        required = ("scientific_question", "expected_contribution", "metrics",
                    "reference_check", "intended_artifact", "estimated_seconds")
        if any(not row[field] for row in rows for field in required):
            raise RuntimeError("expanded plan contains incomplete design metadata")
        target = float(profile["target_runtime_hours"])
        hours = seconds / 3600.0
        if not 0.75 * target <= hours <= 1.25 * target:
            raise RuntimeError(f"estimate {hours:.2f} h is not near {target:g} h")

    blocked = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "run_study.py"),
         "--profile", "full", "--skip-build"], cwd=ROOT, text=True,
        capture_output=True, check=False)
    if blocked.returncode != 3 or "Refusing to execute" not in blocked.stderr:
        raise RuntimeError("draft full-profile execution was not refused")
    if not FINAL_CHECK.validation_errors():
        raise RuntimeError("final report guard unexpectedly accepted draft state")
    print(f"Study-design guard tests passed ({len(runs)} planned rows).")


if __name__ == "__main__":
    main()
