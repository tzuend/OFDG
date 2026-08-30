#!/usr/bin/env python3
"""Validate the independent WENO solver on smooth and exact Sod solutions."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def run(executable: Path, problem: int, cells: int, output: Path) -> tuple[float, list[dict[str, str]]]:
    completed = subprocess.run(
        [str(executable), "--problem", str(problem), "--cells", str(cells),
         "--output", str(output)], text=True, capture_output=True, check=True)
    match = re.search(r"density_l1_error=([^\s]+)", completed.stdout)
    smooth_error = float(match.group(1)) if match else math.nan
    with output.open(newline="") as stream:
        return smooth_error, list(csv.DictReader(stream))


def sod_density(x: float, time: float = 0.2) -> float:
    gamma = 1.4
    p_star = 0.30313017805064685
    u_star = 0.9274526200489499
    rho_left, pressure_left = 1.0, 1.0
    rho_right, pressure_right = 0.125, 0.1
    sound_left = math.sqrt(gamma * pressure_left / rho_left)
    sound_right = math.sqrt(gamma * pressure_right / rho_right)
    xi = x / time
    fan_head = -sound_left
    fan_tail = u_star - sound_left * (p_star / pressure_left) ** (
        (gamma - 1.0) / (2.0 * gamma))
    shock = sound_right * math.sqrt(
        (gamma + 1.0) / (2.0 * gamma) * p_star / pressure_right
        + (gamma - 1.0) / (2.0 * gamma))
    rho_star_left = rho_left * (p_star / pressure_left) ** (1.0 / gamma)
    ratio = p_star / pressure_right
    rho_star_right = rho_right * (
        (ratio + (gamma - 1.0) / (gamma + 1.0)) /
        ((gamma - 1.0) / (gamma + 1.0) * ratio + 1.0))
    if xi <= fan_head:
        return rho_left
    if xi < fan_tail:
        sound = 2.0 / (gamma + 1.0) * (
            sound_left - 0.5 * (gamma - 1.0) * xi)
        return rho_left * (sound / sound_left) ** (2.0 / (gamma - 1.0))
    if xi < u_star:
        return rho_star_left
    if xi < shock:
        return rho_star_right
    return rho_right


def sod_error(rows: list[dict[str, str]]) -> float:
    if len(rows) < 2:
        return math.inf
    dx = float(rows[1]["x"]) - float(rows[0]["x"])
    return sum(abs(float(row["density"]) - sod_density(float(row["x"])))
               for row in rows) * dx


def main() -> None:
    executable = Path(sys.argv[1])
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        smooth_40, _ = run(executable, 4, 40, directory / "smooth40.csv")
        smooth_80, _ = run(executable, 4, 80, directory / "smooth80.csv")
        if not (smooth_80 < smooth_40 / 4.0):
            raise RuntimeError(
                f"smooth WENO validation did not converge: {smooth_40}, {smooth_80}")
        _, sod_100 = run(executable, 6, 100, directory / "sod100.csv")
        _, sod_200 = run(executable, 6, 200, directory / "sod200.csv")
        error_100, error_200 = sod_error(sod_100), sod_error(sod_200)
        if not (error_200 < error_100):
            raise RuntimeError(
                f"Sod reference validation did not improve: {error_100}, {error_200}")
    print("Independent WENO smooth and exact-Riemann tests passed.")


if __name__ == "__main__":
    main()
