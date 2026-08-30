#!/usr/bin/env python3
"""Bounded end-to-end checks for the maintained MPI example drivers."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
import os
from pathlib import Path
import shlex
import signal
import subprocess
import sys
from typing import Sequence


RELATIVE_TOLERANCE = 1.0e-10
ABSOLUTE_TOLERANCE = 1.0e-12
CONSERVATION_TOLERANCE = 1.0e-11
TIMEOUT_SECONDS = 30.0


@dataclass(frozen=True)
class Case:
    name: str
    program: str
    arguments: tuple[str, ...]
    exact_fields: tuple[str, ...]
    numeric_fields: tuple[str, ...]


COMMON_EXACT_FIELDS = (
    "method", "cadence", "order", "steps", "dofs",
    "filter_applications", "active_elements",
)

CASES = (
    Case(
        "advection-dg", "advection",
        ("-d", "1", "-n", "8", "-o", "1", "-tf", "0.001",
         "-dt", "0.001", "-method", "dg", "-no-vis"),
        COMMON_EXACT_FIELDS + ("dimension", "elements"),
        ("time", "l1_error", "l2_error", "linf_error",
         "cell_average_l2_error", "downwind_l2_error",
         "conservation_drift", "minimum", "maximum"),
    ),
    Case(
        "burgers-dg", "burgers",
        ("-p", "1", "-n", "8", "-o", "1", "-tf", "0.001",
         "-dt", "0.001", "-method", "dg", "-no-vis"),
        COMMON_EXACT_FIELDS + ("problem", "elements"),
        ("time", "l1_error", "l2_error", "linf_error",
         "conservation_drift", "min", "max"),
    ),
    Case(
        "euler-dg", "euler",
        ("-p", "4", "-r", "0", "-n", "8", "-o", "1",
         "-tf", "0.001", "-dt", "0.001", "-method", "dg",
         "-no-positivity", "-no-vis"),
        COMMON_EXACT_FIELDS + (
            "problem", "positivity_applications", "limited_elements",
            "rejected_steps",
        ),
        ("time", "l1_error", "l2_error", "linf_error",
         "state_l2_error", "conservation_drift", "min_density",
         "min_pressure"),
    ),
    Case(
        "advection-ofdg-kxrcf", "advection",
        ("-d", "1", "-n", "8", "-o", "1", "-tf", "0.001",
         "-dt", "0.001", "-method", "ofdg-kxrcf", "-no-vis"),
        COMMON_EXACT_FIELDS + ("dimension", "elements"),
        ("time", "l1_error", "l2_error", "linf_error",
         "cell_average_l2_error", "downwind_l2_error",
         "conservation_drift", "minimum", "maximum"),
    ),
    Case(
        "euler-oedg-positivity", "euler",
        ("-p", "4", "-r", "0", "-n", "8", "-o", "1",
         "-tf", "0.001", "-dt", "0.001", "-method", "oedg",
         "-positivity", "-no-vis"),
        COMMON_EXACT_FIELDS + (
            "problem", "positivity_applications", "limited_elements",
            "rejected_steps",
        ),
        ("time", "l1_error", "l2_error", "linf_error",
         "state_l2_error", "conservation_drift", "min_density",
         "min_pressure"),
    ),
)


def terminate_process_group(process: subprocess.Popen[str]) -> None:
    """Terminate an MPI launch and all ranks without leaving CI processes."""
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=2.0)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def run_case(command: Sequence[str], timeout: float) -> str:
    environment = os.environ.copy()
    if hasattr(os, "geteuid") and os.geteuid() == 0:
        environment.setdefault("OMPI_ALLOW_RUN_AS_ROOT", "1")
        environment.setdefault("OMPI_ALLOW_RUN_AS_ROOT_CONFIRM", "1")

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=environment,
        start_new_session=True,
    )
    try:
        output, _ = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired as error:
        terminate_process_group(process)
        partial = error.output or ""
        raise RuntimeError(
            f"command timed out after {timeout:g} seconds:\n"
            f"{' '.join(command)}\n{partial}"
        ) from error

    if process.returncode != 0:
        raise RuntimeError(
            f"command exited with {process.returncode}:\n"
            f"{' '.join(command)}\n{output}"
        )
    return output


def parse_summary(output: str, case: Case, ranks: int) -> dict[str, str]:
    lines = [line.strip() for line in output.splitlines()
             if line.startswith("method=")]
    if len(lines) != 1:
        raise RuntimeError(
            f"{case.name} with {ranks} rank(s) produced {len(lines)} "
            f"root summaries instead of one:\n{output}"
        )
    fields: dict[str, str] = {}
    for token in lines[0].split():
        if "=" in token:
            key, value = token.split("=", 1)
            fields[key] = value
    missing = [field for field in case.exact_fields + case.numeric_fields
               if field not in fields]
    if missing:
        raise RuntimeError(
            f"{case.name} summary is missing {', '.join(missing)}:\n{lines[0]}"
        )
    return fields


def compare_summaries(case: Case, serial: dict[str, str],
                      parallel: dict[str, str]) -> None:
    for field in case.exact_fields:
        if serial[field] != parallel[field]:
            raise RuntimeError(
                f"{case.name}: {field} differs between one and two ranks: "
                f"{serial[field]} != {parallel[field]}"
            )

    for field in case.numeric_fields:
        serial_value = float(serial[field])
        parallel_value = float(parallel[field])
        if not math.isfinite(serial_value) or not math.isfinite(parallel_value):
            raise RuntimeError(f"{case.name}: non-finite {field}")
        if not math.isclose(serial_value, parallel_value,
                            rel_tol=RELATIVE_TOLERANCE,
                            abs_tol=ABSOLUTE_TOLERANCE):
            raise RuntimeError(
                f"{case.name}: {field} differs between one and two ranks: "
                f"{serial_value:.17g} != {parallel_value:.17g}"
            )

    if abs(float(parallel["conservation_drift"])) > CONSERVATION_TOLERANCE:
        raise RuntimeError(
            f"{case.name}: parallel conservation drift exceeds "
            f"{CONSERVATION_TOLERANCE:g}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--example-dir", type=Path,
                        default=Path("build/release/examples"))
    parser.add_argument("--mpiexec", default=os.environ.get("MPIEXEC", "mpirun"))
    parser.add_argument("--timeout", type=float, default=TIMEOUT_SECONDS)
    arguments = parser.parse_args()
    mpi_command = shlex.split(arguments.mpiexec)

    for case in CASES:
        executable = arguments.example_dir / case.program
        if not executable.is_file():
            raise RuntimeError(f"missing example executable: {executable}")
        summaries = {}
        for ranks in (1, 2):
            command = [*mpi_command, "-np", str(ranks), str(executable),
                       *case.arguments]
            output = run_case(command, arguments.timeout)
            summaries[ranks] = parse_summary(output, case, ranks)
        compare_summaries(case, summaries[1], summaries[2])

    print(f"Parallel example tests passed ({len(CASES)} cases, 1/2 ranks).")


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, ValueError) as error:
        print(f"Parallel example tests failed: {error}", file=sys.stderr)
        raise SystemExit(1)
