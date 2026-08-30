#!/usr/bin/env python3
"""Run the manifest-defined DG/OFDG/OEDG numerical study."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import itertools
import json
import os
from pathlib import Path
import platform
import re
import shlex
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "experiments" / "study_manifest.json"
BIN = ROOT / "build" / "release" / "examples"
REFERENCE_BIN = ROOT / "build" / "release" / "benchmarks" / "euler_reference_weno"

FIELDS = [
    "case", "kind", "program", "method", "cadence", "order", "resolution",
    "resolution_flag", "dofs", "ranks", "repeat", "scale", "cfl", "dt",
    "rk_method", "time", "steps", "l1_error", "l2_error", "linf_error",
    "state_l2_error", "cell_average_l2_error", "downwind_l2_error",
    "reference_l1_error", "conservation_drift",
    "total_variation", "minimum", "maximum", "min_density", "min_pressure",
    "overshoot", "undershoot", "filter_applications", "active_elements",
    "positivity_applications", "limited_elements", "rejected_steps",
    "failure_time", "runtime_seconds", "wall_time_seconds", "status",
    "exit_code", "failure",
    "profile_prefix", "reference", "log", "command", "measured_at_utc",
    "commit", "compiler", "mfem", "host", "manifest_profile"
]

PLAN_FIELDS = [
    "case", "case_status", "profile_status", "scientific_question",
    "expected_contribution", "equation", "domain", "boundary_conditions",
    "method", "order", "resolution", "resolution_flag", "ranks", "scale",
    "time_integrator", "metrics", "reference_check", "priority",
    "intended_artifact", "estimated_seconds", "command"
]


def load_manifest(profile: str) -> tuple[dict, dict]:
    document = json.loads(MANIFEST.read_text())
    return document, document["profiles"][profile]


def git_value(*args: str) -> str:
    result = subprocess.run(["git", *args], cwd=ROOT, text=True,
                            capture_output=True, check=False)
    return result.stdout.strip() if result.returncode == 0 else "unknown"


def version_line(command: list[str]) -> str:
    result = subprocess.run(command, cwd=ROOT, text=True,
                            capture_output=True, check=False)
    output = (result.stdout or result.stderr).splitlines()
    return output[0].strip() if output else "unknown"


def expand_cases(document: dict, profile: dict, only: set[str]) -> list[dict]:
    methods = list(document["methods"]) + list(document.get("ablation_methods", []))
    expanded: list[dict] = []
    for case in profile["cases"]:
        if only and case["name"] not in only:
            continue
        case_methods = case.get("methods", methods)
        ranks = case.get("ranks", [1])
        scales = case.get("scales", [1.0])
        if case.get("zip_order_resolution"):
            order_resolutions = zip(case["orders"], case["resolutions"])
        else:
            order_resolutions = itertools.product(case["orders"],
                                                   case["resolutions"])
        pairs = list(order_resolutions)
        repeats = profile["timing_repeats"] if case["kind"] == "timing" else 1
        for (order, resolution), method, rank_count, scale, repeat in itertools.product(
                pairs, case_methods, ranks, scales, range(1, repeats + 1)):
            expanded.append({**case, "order": order, "resolution": resolution,
                             "method": method, "rank_count": rank_count,
                             "scale": scale, "repeat": repeat})
    return expanded


def ensure_references(profile: dict, output_dir: Path,
                      required: set[str]) -> dict[str, str]:
    paths: dict[str, str] = {}
    reference_dir = output_dir / "references"
    reference_dir.mkdir(parents=True, exist_ok=True)
    for specification in profile.get("references", []):
        if specification["name"] not in required:
            continue
        destination = reference_dir / f"{specification['name']}.csv"
        paths[specification["name"]] = str(destination)
        if destination.exists():
            continue
        command = [str(REFERENCE_BIN), "--problem", str(specification["problem"]),
                   "--cells", str(specification["cells"]), "--final-time",
                   str(specification["final_time"]), "--scale",
                   str(specification.get("scale", 1.0)), "--output", str(destination)]
        subprocess.run(command, cwd=ROOT, check=True)
    return paths


def command_for(run: dict, profile_prefix: Path) -> list[str]:
    command = [str(BIN / run["program"]), *map(str, run.get("args", [])),
               run["resolution_flag"], str(run["resolution"]),
               "-o", str(run["order"]), "-method", run["method"],
               "-cadence", "auto"]
    if run["program"] == "euler":
        command.extend(["-scale", str(run["scale"])])
    if run["kind"] in {"profile", "failure"}:
        command.extend(["-profile", str(profile_prefix)])
    if run["rank_count"] > 1:
        command = ["mpirun", "--oversubscribe", "-np",
                   str(run["rank_count"]), *command]
    return command


def parse_metrics(output: str) -> dict[str, str]:
    line = ""
    for candidate in output.splitlines():
        if "runtime_seconds=" in candidate:
            line = candidate
    return dict(re.findall(r"([A-Za-z_][A-Za-z0-9_]*)=([^\s]+)", line))


def requested_value(arguments: list[str], flag: str, default: str = "") -> str:
    try:
        return str(arguments[arguments.index(flag) + 1])
    except (ValueError, IndexError):
        return default


def run_one(index: int, total: int, run: dict, output_dir: Path,
            references: dict[str, str], metadata: dict[str, str]) -> dict[str, str]:
    raw_dir = output_dir / "raw"
    profile_dir = output_dir / "profiles"
    raw_dir.mkdir(parents=True, exist_ok=True)
    profile_dir.mkdir(parents=True, exist_ok=True)
    stem = (f"{index:04d}_{run['name']}_{run['method']}_p{run['order']}_"
            f"r{run['resolution']}_np{run['rank_count']}_s{run['scale']}_"
            f"rep{run['repeat']}")
    profile_prefix = profile_dir / stem
    command = command_for(run, profile_prefix)
    print(f"[{index}/{total}] {run['name']} {run['method']} "
          f"p={run['order']} resolution={run['resolution']} np={run['rank_count']}",
          flush=True)
    started = time.perf_counter()
    result = subprocess.run(command, cwd=raw_dir, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            check=False)
    wall_time = time.perf_counter() - started
    log = raw_dir / f"{stem}.log"
    log.write_text(result.stdout)
    metrics = parse_metrics(result.stdout)
    failure_match = re.search(r"(Nonphysical Euler state[^\n]*|"
                              r"Euler step remained inadmissible[^\n]*)",
                              result.stdout)
    failure_time_match = re.search(r"\bt=([0-9.eE+-]+)", result.stdout)
    base_args = list(map(str, run.get("args", [])))
    row = {field: "" for field in FIELDS}
    row.update(metadata)
    row.update({
        "case": run["name"], "kind": run["kind"],
        "program": run["program"], "method": run["method"],
        "cadence": metrics.get("cadence", ""), "order": str(run["order"]),
        "resolution": str(run["resolution"]),
        "resolution_flag": run["resolution_flag"],
        "ranks": str(run["rank_count"]), "repeat": str(run["repeat"]),
        "scale": str(run["scale"]), "cfl": requested_value(base_args, "-c"),
        "dt": requested_value(base_args, "-dt", "cfl"),
        "rk_method": requested_value(base_args, "-s", "4"),
        "wall_time_seconds": f"{wall_time:.9g}",
        "status": "ok" if result.returncode == 0 else "failed",
        "exit_code": str(result.returncode),
        "failure": failure_match.group(1) if failure_match else "",
        "failure_time": (failure_time_match.group(1)
                         if failure_time_match else ""),
        "profile_prefix": str(profile_prefix),
        "reference": references.get(run.get("reference", ""), ""),
        "log": str(log.relative_to(output_dir)),
        "command": shlex.join(command),
    })
    for field in FIELDS:
        if field in metrics:
            row[field] = metrics[field]
    if not row["minimum"]:
        row["minimum"] = metrics.get("min", "")
    if not row["maximum"]:
        row["maximum"] = metrics.get("max", "")
    return row


def run_warmup(run: dict, output_dir: Path) -> None:
    """Perform the single unrecorded warm-up required by timing cases."""
    warmup_prefix = output_dir / "profiles" / "timing_warmup"
    command = command_for(run, warmup_prefix)
    print(f"[warm-up] {run['name']} {run['method']} p={run['order']} "
          f"resolution={run['resolution']} np={run['rank_count']}", flush=True)
    result = subprocess.run(command, cwd=output_dir, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            check=False)
    if result.returncode != 0:
        raise RuntimeError("timing warm-up failed:\n" + result.stdout)


def write_rows(destination: Path, rows: list[dict[str, str]]) -> None:
    """Atomically checkpoint all completed rows after every run."""
    temporary = destination.with_suffix(".csv.tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(destination)


def write_study_plan(document: dict, profile: dict, runs: list[dict],
                     destination: Path) -> float:
    """Write the expanded proposed matrix without executing any solver."""
    design = document.get("study_design", {})
    destination.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    total_seconds = 0.0
    safety_factor = float(profile.get("runtime_estimate_safety_factor", 1.0))
    for run in runs:
        if run["name"] not in design:
            raise ValueError(f"missing study_design metadata for {run['name']}")
        metadata = design[run["name"]]
        estimated = float(metadata["estimated_seconds_per_run"]) * safety_factor
        total_seconds += estimated
        planned_profile = ROOT / "measurements" / "study" / "planned" / run["name"]
        rows.append({
            "case": run["name"],
            "case_status": metadata["status"],
            "profile_status": profile.get("status", "unspecified"),
            "scientific_question": metadata["question"],
            "expected_contribution": metadata["contribution"],
            "equation": metadata["equation"],
            "domain": metadata["domain"],
            "boundary_conditions": metadata["boundary_conditions"],
            "method": run["method"],
            "order": str(run["order"]),
            "resolution": str(run["resolution"]),
            "resolution_flag": run["resolution_flag"],
            "ranks": str(run["rank_count"]),
            "scale": str(run["scale"]),
            "time_integrator": metadata["time_integrator"],
            "metrics": "; ".join(metadata["metrics"]),
            "reference_check": metadata["reference_check"],
            "priority": metadata["priority"],
            "intended_artifact": metadata["artifact"],
            "estimated_seconds": f"{estimated:.6g}",
            "command": shlex.join(command_for(run, planned_profile)),
        })
    with destination.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, PLAN_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return total_seconds


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("quick", "full", "hpc"),
                        default="quick")
    parser.add_argument("--only", action="append", default=[],
                        help="Run only a named case; may be repeated.")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--dry-run", action="store_true",
                        help="expand and validate the design without running solvers")
    parser.add_argument("--plan-output", type=Path,
                        help="CSV destination used with --dry-run")
    args = parser.parse_args()

    document, profile = load_manifest(args.profile)
    output_dir = (args.output or ROOT / "measurements" / "study" / args.profile).resolve()
    runs = expand_cases(document, profile, set(args.only))

    if args.dry_run:
        destination = (args.plan_output or
                       ROOT / "measurements" / "study" / "plans" /
                       f"{args.profile}.csv").resolve()
        total_seconds = write_study_plan(document, profile, runs, destination)
        target = profile.get("target_runtime_hours")
        target_text = f"; target={target} h" if target is not None else ""
        print(f"Wrote {destination} ({len(runs)} proposed rows, "
              f"estimated {total_seconds / 3600:.2f} h{target_text}).")
        return 0

    profile_status = profile.get("status", "unspecified")
    if profile_status != "approved":
        print(f"Refusing to execute profile '{args.profile}': manifest status is "
              f"'{profile_status}', not 'approved'. Use --dry-run to inspect it.",
              file=sys.stderr)
        return 3

    output_dir.mkdir(parents=True, exist_ok=True)
    if not args.skip_build:
        subprocess.run(["make", "release", "reference-solver"], cwd=ROOT,
                       check=True)

    required_references = {run["reference"] for run in runs
                           if run.get("reference")}
    references = ensure_references(profile, output_dir, required_references)
    metadata = {
        "measured_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "commit": git_value("rev-parse", "HEAD"),
        "compiler": version_line(["mpicxx", "--version"]),
        "mfem": git_value("-C", "../mfem", "rev-parse", "HEAD"),
        "host": platform.node(), "manifest_profile": args.profile,
    }
    destination = output_dir / "results.csv"
    rows = []
    warmed: set[tuple] = set()
    for i, run in enumerate(runs, 1):
        warmup_key = (run["name"], run["method"], run["order"],
                      run["resolution"], run["rank_count"], run["scale"])
        if run["kind"] == "timing" and warmup_key not in warmed:
            run_warmup(run, output_dir)
            warmed.add(warmup_key)
        rows.append(run_one(i, len(runs), run, output_dir,
                            references, metadata))
        write_rows(destination, rows)
    print(f"Wrote {destination} ({len(rows)} rows)")
    return 0 if all(row["status"] == "ok" or row["kind"] == "failure"
                    for row in rows) else 2


if __name__ == "__main__":
    sys.exit(main())
