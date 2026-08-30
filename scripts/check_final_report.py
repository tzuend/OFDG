#!/usr/bin/env python3
"""Refuse a final report build until approval and generated findings exist."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "experiments" / "study_manifest.json"
REQUIRED = (
    ROOT / "report" / "generated" / "final" / "results.tex",
    ROOT / "report" / "generated" / "final" / "discussion.tex",
    ROOT / "report" / "generated" / "final" / "conclusion.tex",
)


def validation_errors() -> list[str]:
    document = json.loads(MANIFEST.read_text())
    errors = []
    status = document["profiles"]["full"].get("status", "unspecified")
    if status != "approved":
        errors.append(f"full study status is '{status}', not 'approved'")
    for path in REQUIRED:
        if not path.is_file() or not path.read_text().strip():
            errors.append(f"missing generated final section: {path.relative_to(ROOT)}")
    return errors


def main() -> int:
    errors = validation_errors()
    if errors:
        print("Final report build refused:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        print("Review and approve the manifest, run the final study, and generate "
              "the final findings first.", file=sys.stderr)
        return 3
    print("Final report approval and generated-section checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
