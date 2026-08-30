#!/usr/bin/env python3
"""Backward-compatible entry point for the manifest-driven study runner."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compatibility wrapper; prefer scripts/run_study.py.")
    parser.add_argument("--matrix", choices=("quick", "full", "hpc"),
                        default="quick")
    args, remainder = parser.parse_known_args()
    print("run_measurements.py is deprecated; forwarding to run_study.py",
          file=sys.stderr)
    command = [sys.executable, str(ROOT / "scripts" / "run_study.py"),
               "--profile", args.matrix, *remainder]
    return subprocess.run(command, cwd=ROOT, check=False).returncode


if __name__ == "__main__":
    sys.exit(main())
