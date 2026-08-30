#!/usr/bin/env python3
"""Process-level checks for contracts that deliberately terminate via MFEM."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def run(binary: Path, mode: str, mesh: Path | None = None) -> subprocess.CompletedProcess[str]:
    command = [str(binary), mode]
    if mesh is not None:
        command.append(str(mesh))
    return subprocess.run(command, capture_output=True, text=True, timeout=20, check=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("--mfem-data", type=Path, required=True)
    arguments = parser.parse_args()

    supported = run(arguments.binary, "supported")
    if supported.returncode != 0:
        raise RuntimeError(f"supported mesh was rejected:\n{supported.stderr}")

    rejected = [
        ("mixed", arguments.mfem_data / "square-mixed.mesh", "multiple element geometries"),
        ("curved", arguments.mfem_data / "disc-nurbs.mesh", "curvilinear or non-affine"),
        ("variable-order", None, "variable polynomial order"),
        ("modified", None, "changed after construction"),
    ]
    for mode, mesh, expected in rejected:
        result = run(arguments.binary, mode, mesh)
        output = result.stdout + result.stderr
        if result.returncode == 0:
            raise RuntimeError(f"unsupported {mode} case unexpectedly succeeded")
        if expected not in output:
            raise RuntimeError(
                f"unsupported {mode} case did not report {expected!r}:\n{output}"
            )

    print("Supported-space contract tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
