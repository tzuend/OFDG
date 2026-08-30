#!/usr/bin/env python3
"""Summarize study data and generate deliberately separated report assets.

Quick-profile outputs are presentation prototypes only.  The path-selection
guard in :func:`report_asset_root` prevents them from entering the directory
reserved for approved final evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import statistics
import sys
from collections import defaultdict


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "experiments" / "study_manifest.json"
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "measurements" / ".matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HEADLINE = ("dg", "ofdg-kxrcf", "oedg")
COLORS = {"dg": "#222222", "ofdg-kxrcf": "#0072B2",
          "oedg": "#D55E00", "ofdg": "#009E73"}
LABELS = {"dg": "RKDG", "ofdg-kxrcf": "adapted OFDG--KXRCF",
          "oedg": "2024 OEDG", "ofdg": "OFDG (all cells)"}
LINESTYLES = {"dg": "-", "ofdg-kxrcf": "--", "oedg": "-."}
MARKERS = {"dg": "o", "ofdg-kxrcf": "s", "oedg": "^"}


def number(value: str | float | int | None) -> float:
    try:
        result = float(value)
        return result if math.isfinite(result) else math.nan
    except (TypeError, ValueError):
        return math.nan


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def rebase_recorded_paths(rows: list[dict[str, str]], data_root: Path) -> None:
    """Make a retained study checkpoint portable across checkout locations.

    Raw rows retain the paths emitted by the original run for provenance. If
    such a path is unavailable, use a same-named profile or reference beside
    the input CSV. Rows from a live run continue to use their original paths.
    """
    for row in rows:
        prefix = row.get("profile_prefix", "")
        if prefix and not profile_files(prefix):
            local_prefix = data_root / "profiles" / Path(prefix).name
            if profile_files(str(local_prefix)):
                row["profile_prefix"] = str(local_prefix)

        reference = row.get("reference", "")
        if reference and not Path(reference).exists():
            local_reference = data_root / "references" / Path(reference).name
            if local_reference.exists():
                row["reference"] = str(local_reference)


def profile_files(prefix: str) -> list[Path]:
    if not prefix:
        return []
    base = Path(prefix)
    return sorted(base.parent.glob(base.name + ".rank*.csv"))


def load_profile(prefix: str) -> dict[str, np.ndarray]:
    records: list[dict[str, str]] = []
    for path in profile_files(prefix):
        with path.open(newline="") as stream:
            records.extend(csv.DictReader(stream))
    if not records:
        return {}
    result: dict[str, np.ndarray] = {}
    for key in records[0]:
        if key == "sample":
            result[key] = np.asarray([record[key] for record in records],
                                     dtype=object)
        else:
            result[key] = np.asarray([number(record[key]) for record in records])
    return result


def reference_error(row: dict[str, str], profile: dict[str, np.ndarray]) -> float:
    reference_path = row.get("reference", "")
    if not reference_path or not profile or not Path(reference_path).exists():
        return math.nan
    with Path(reference_path).open(newline="") as stream:
        reference = list(csv.DictReader(stream))
    mask = profile["sample"] == "center"
    x = profile["x"][mask]
    density_key = "density" if "density" in profile else "value"
    density = profile[density_key][mask] / number(row.get("scale", "1"))
    ordering = np.argsort(x)
    x = x[ordering]
    density = density[ordering]
    ref_x = np.asarray([float(item["x"]) for item in reference])
    ref_density = np.asarray([float(item["density"]) for item in reference])
    interpolated = np.interp(x, ref_x, ref_density)
    return float(np.trapezoid(np.abs(density - interpolated), x) /
                 max(x[-1] - x[0], np.finfo(float).eps))


def summarize(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Combine repetitions and compute rates without discarding raw metrics."""
    groups: dict[tuple, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (row["case"], row["kind"], row["method"], row["order"],
               row["resolution"], row["ranks"], row["scale"])
        groups[key].append(row)

    summaries: list[dict[str, str]] = []
    for members in groups.values():
        representative = members[0]
        runtimes = [number(item.get("runtime_seconds")) for item in members]
        runtimes = [value for value in runtimes if math.isfinite(value)]
        item = dict(representative)
        item["runtime_median"] = (f"{statistics.median(runtimes):.12g}"
                                  if runtimes else "nan")
        item["runtime_spread"] = (f"{statistics.pstdev(runtimes):.12g}"
                                  if len(runtimes) > 1 else "0")
        profile = load_profile(representative.get("profile_prefix", ""))
        item["reference_l1_error"] = f"{reference_error(representative, profile):.12g}"
        if profile:
            mask = profile["sample"] == "center"
            key = "density" if "density" in profile else "value"
            values = profile[key][mask]
            ordering = np.argsort(profile["x"][mask])
            values = values[ordering]
            item["total_variation"] = f"{np.abs(np.diff(values)).sum():.12g}"
            item["profile_minimum"] = f"{values.min():.12g}"
            item["profile_maximum"] = f"{values.max():.12g}"
        else:
            item["total_variation"] = "nan"
            item["profile_minimum"] = "nan"
            item["profile_maximum"] = "nan"
        applications = number(item.get("filter_applications"))
        active = number(item.get("active_elements"))
        item["active_per_application"] = (
            f"{active / applications:.12g}"
            if applications > 0 and math.isfinite(active) else "nan")
        for metric in ("l1_error", "l2_error", "linf_error",
                       "cell_average_l2_error", "downwind_l2_error"):
            item[f"observed_{metric}_rate"] = "nan"
        item["observed_rate"] = "nan"  # compatibility alias for L2
        summaries.append(item)

    rate_groups: dict[tuple, list[dict[str, str]]] = defaultdict(list)
    for item in summaries:
        if item["kind"] == "convergence" and item["status"] == "ok":
            rate_groups[(item["case"], item["method"], item["order"],
                         item["ranks"], item["scale"])].append(item)
    for members in rate_groups.values():
        members.sort(key=lambda item: number(item["resolution"]))
        previous = {metric: math.nan for metric in (
            "l1_error", "l2_error", "linf_error",
            "cell_average_l2_error", "downwind_l2_error")}
        for item in members:
            for metric, old_error in previous.items():
                error = number(item.get(metric))
                if math.isfinite(old_error) and old_error > 0 and error > 0:
                    item[f"observed_{metric}_rate"] = (
                        f"{math.log(old_error / error, 2):.6g}")
                previous[metric] = error
            item["observed_rate"] = item["observed_l2_error_rate"]
    return summaries


def report_asset_root(profile: str, mode: str) -> Path | None:
    """Resolve a report destination while enforcing draft/final separation."""
    if mode == "none":
        return None
    if mode == "auto":
        mode = "draft" if profile == "quick" else "none"
    if mode == "none":
        return None
    if mode == "final":
        manifest = json.loads(MANIFEST.read_text())
        status = manifest["profiles"][profile].get("status")
        if profile != "full" or status != "approved":
            raise ValueError("final report assets require the approved full profile")
    if profile == "quick" and mode != "draft":
        raise ValueError("quick data may only be written to draft report assets")
    return ROOT / "report" / "generated" / mode


def save_figure(fig: plt.Figure, name: str, figure_dir: Path,
                report_figure_dir: Path | None) -> None:
    destinations = [figure_dir]
    if report_figure_dir is not None:
        destinations.append(report_figure_dir)
    for directory in destinations:
        directory.mkdir(parents=True, exist_ok=True)
        fig.savefig(directory / f"{name}.pdf", bbox_inches="tight")
        fig.savefig(directory / f"{name}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def convergence_figures(summary: list[dict[str, str]], figure_dir: Path,
                        report_figure_dir: Path | None) -> None:
    """One case and one polynomial degree per axis, with three methods max."""
    cases = sorted({item["case"] for item in summary
                    if item["kind"] == "convergence"
                    and item["case"] != "advection_superconvergence"})
    for case in cases:
        orders = sorted({int(item["order"]) for item in summary
                         if item["case"] == case and item["status"] == "ok"})
        for order in orders:
            fig, ax = plt.subplots(figsize=(5.9, 4.0))
            all_points: list[tuple[float, float]] = []
            for method in HEADLINE:
                values = sorted((item for item in summary
                                 if item["case"] == case
                                 and item["method"] == method
                                 and int(item["order"]) == order
                                 and item["status"] == "ok"),
                                key=lambda item: number(item["dofs"]))
                points = [(number(item["dofs"]), number(item["l2_error"]))
                          for item in values]
                points = [(x, y) for x, y in points
                          if x > 0 and y > 0 and math.isfinite(x + y)]
                if not points:
                    continue
                all_points.extend(points)
                ax.loglog([point[0] for point in points],
                          [point[1] for point in points],
                          color=COLORS[method], marker=MARKERS[method],
                          linestyle=LINESTYLES[method], label=LABELS[method])
            if not all_points:
                plt.close(fig)
                continue
            dimension = 2 if any(token in case for token in ("2d", "vortex")) else 1
            x0 = min(point[0] for point in all_points)
            x1 = max(point[0] for point in all_points)
            if x1 > x0:
                y0 = max(point[1] for point in all_points)
                exponent = -(order + 1) / dimension
                guide_x = np.asarray([x0, x1])
                guide_y = 1.35 * y0 * (guide_x / x0) ** exponent
                ax.loglog(guide_x, guide_y, color="#777777", linewidth=1.0,
                          linestyle=":", label=rf"reference $h^{{{order + 1}}}$")
            ax.set_xlabel("global degrees of freedom")
            ax.set_ylabel(r"$L^2$ error")
            ax.set_title(f"{case.replace('_', ' ')}, $P^{order}$")
            ax.grid(True, which="both", alpha=0.22)
            ax.legend(fontsize=8)
            save_figure(fig, f"convergence_{case}_p{order}", figure_dir,
                        report_figure_dir)


def _selected_profile(summary: list[dict[str, str]], case: str, method: str,
                      order: int, scale: float = 1.0) -> tuple[dict, dict]:
    item = next((row for row in summary if row["case"] == case
                 and row["method"] == method and int(row["order"]) == order
                 and row["status"] == "ok"
                 and math.isclose(number(row["scale"]), scale)), None)
    return (item or {}), load_profile(item.get("profile_prefix", "")) if item else {}


def lax_profile_figure(summary: list[dict[str, str]], figure_dir: Path,
                       report_figure_dir: Path | None, order: int = 2) -> None:
    fig, ax = plt.subplots(figsize=(6.3, 4.0))
    plotted = False
    reference_path = ""
    for method in HEADLINE:
        item, profile = _selected_profile(summary, "lax_scaled", method, order)
        if not profile:
            continue
        mask = profile["sample"] == "center"
        ordering = np.argsort(profile["x"][mask])
        ax.plot(profile["x"][mask][ordering], profile["density"][mask][ordering],
                color=COLORS[method], linestyle=LINESTYLES[method],
                label=LABELS[method])
        plotted = True
        reference_path = reference_path or item.get("reference", "")
    if reference_path and Path(reference_path).exists():
        with Path(reference_path).open(newline="") as stream:
            reference = list(csv.DictReader(stream))
        ax.plot([number(row["x"]) for row in reference],
                [number(row["density"]) for row in reference],
                color="#777777", linewidth=1.0, label="WENO5 reference")
    if not plotted:
        plt.close(fig)
        return
    ax.set_xlabel("x")
    ax.set_ylabel(r"density $\rho$")
    ax.set_title(rf"Scaled Lax problem, $\lambda=1$, $P^{order}$")
    ax.grid(True, alpha=0.22)
    ax.legend(fontsize=8)
    save_figure(fig, f"profile_lax_scaled_p{order}", figure_dir,
                report_figure_dir)


def scale_invariance_figure(summary: list[dict[str, str]], figure_dir: Path,
                            report_figure_dir: Path | None,
                            order: int = 2) -> None:
    fig, axes = plt.subplots(1, len(HEADLINE), figsize=(10.4, 3.3),
                             sharex=True, sharey=True)
    plotted = False
    for ax, method in zip(axes, HEADLINE):
        rows = sorted((row for row in summary if row["case"] == "lax_scaled"
                       and row["method"] == method
                       and int(row["order"]) == order
                       and row["status"] == "ok"),
                      key=lambda row: number(row["scale"]))
        for row in rows:
            profile = load_profile(row.get("profile_prefix", ""))
            if not profile:
                continue
            mask = profile["sample"] == "center"
            ordering = np.argsort(profile["x"][mask])
            scale = number(row["scale"])
            ax.plot(profile["x"][mask][ordering],
                    profile["density"][mask][ordering] / scale,
                    label=rf"$\lambda={scale:g}$")
            plotted = True
        ax.set_title(LABELS[method])
        ax.set_xlabel("x")
        ax.grid(True, alpha=0.22)
    axes[0].set_ylabel(r"normalized density $\rho/\lambda$")
    axes[-1].legend(fontsize=8)
    if plotted:
        save_figure(fig, f"scale_invariance_p{order}", figure_dir,
                    report_figure_dir)
    else:
        plt.close(fig)


def riemann_2d_figures(summary: list[dict[str, str]], figure_dir: Path,
                       report_figure_dir: Path | None, order: int = 1) -> None:
    profiles: dict[str, dict[str, np.ndarray]] = {}
    all_density: list[float] = []
    for method in HEADLINE:
        _, profile = _selected_profile(summary, "riemann_2d", method, order)
        if not profile:
            continue
        profiles[method] = profile
        mask = profile["sample"] == "center"
        all_density.extend(profile["density"][mask])
    if not profiles:
        return
    levels = np.linspace(min(all_density), max(all_density), 25)
    fig, axes = plt.subplots(1, len(HEADLINE), figsize=(10.6, 3.25),
                             sharex=True, sharey=True)
    contour = None
    for ax, method in zip(axes, HEADLINE):
        profile = profiles.get(method)
        if profile is None:
            ax.set_axis_off()
            continue
        mask = profile["sample"] == "center"
        contour = ax.tricontourf(profile["x"][mask], profile["y"][mask],
                                 profile["density"][mask], levels=levels,
                                 cmap="viridis", extend="both")
        ax.set_title(LABELS[method])
        ax.set_xlabel("x")
        ax.set_aspect("equal")
    axes[0].set_ylabel("y")
    if contour is not None:
        fig.colorbar(contour, ax=axes.tolist(), shrink=0.84, label="density")
    save_figure(fig, f"contour_riemann_2d_p{order}", figure_dir,
                report_figure_dir)

    fig, ax = plt.subplots(figsize=(6.3, 3.8))
    for method in HEADLINE:
        profile = profiles.get(method)
        if profile is None:
            continue
        center_mask = profile["sample"] == "center"
        distinct_y = np.unique(np.round(profile["y"][center_mask], 13))
        selected_y = distinct_y[np.argmin(np.abs(distinct_y - 0.5))]
        tolerance = max(1e-12, 1e-10 * max(1.0, abs(selected_y)))
        mask = center_mask & (np.abs(profile["y"] - selected_y) <= tolerance)
        ordering = np.argsort(profile["x"][mask])
        ax.plot(profile["x"][mask][ordering], profile["density"][mask][ordering],
                color=COLORS[method], linestyle=LINESTYLES[method],
                label=LABELS[method])
    ax.set_xlabel("x")
    ax.set_ylabel("density")
    ax.set_title(rf"2D Riemann horizontal cut, $P^{order}$")
    ax.grid(True, alpha=0.22)
    ax.legend(fontsize=8)
    save_figure(fig, f"linecut_riemann_2d_p{order}", figure_dir,
                report_figure_dir)


def kxrcf_sensor_figure(x: np.ndarray, indicator: np.ndarray, threshold: float,
                        active: np.ndarray, destination: Path) -> None:
    """Create the approved spatial sensor layout; used by style tests for now."""
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.plot(x, indicator, color=COLORS["ofdg-kxrcf"], label="KXRCF indicator")
    ax.axhline(threshold, color="#555555", linestyle="--", label="threshold")
    ax.scatter(x[active], indicator[active], color="#D55E00", marker="o",
               s=22, label="active cell")
    ax.set_xlabel("cell center")
    ax.set_ylabel("indicator")
    ax.grid(True, alpha=0.22)
    ax.legend(fontsize=8)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, bbox_inches="tight")
    plt.close(fig)


def mesh_label(resolution: str | int, dimension: int = 1) -> str:
    return str(resolution) if dimension == 1 else rf"${resolution}\times {resolution}$"


def _row_key(row: dict[str, str]) -> tuple[int, float]:
    return int(row["order"]), number(row["resolution"])


def assert_identical_method_rows(summary: list[dict[str, str]], case: str,
                                 methods: tuple[str, ...] = HEADLINE) -> None:
    row_sets = {}
    for method in methods:
        row_sets[method] = {_row_key(row) for row in summary
                            if row["case"] == case and row["method"] == method
                            and row["kind"] == "convergence"
                            and row["status"] == "ok"}
    baseline = row_sets[methods[0]]
    differing = {method: rows for method, rows in row_sets.items()
                 if rows != baseline}
    if differing:
        raise ValueError(f"method tables do not have identical rows: {differing}")


def _metric_rate(values: list[float]) -> list[float]:
    rates = [math.nan]
    for previous, current in zip(values, values[1:]):
        rates.append(math.log(previous / current, 2)
                     if previous > 0 and current > 0 else math.nan)
    return rates


def _format_error(value: float) -> str:
    return f"{value:.2e}" if math.isfinite(value) else "--"


def _format_rate(value: float) -> str:
    return f"{value:.2f}" if math.isfinite(value) else "--"


def write_advection_table(summary: list[dict[str, str]], method: str,
                          destination: Path, draft: bool = True) -> None:
    """Write the reference-paper e1--e5 table for one method."""
    base = {(row["order"], row["resolution"]): row for row in summary
            if row["case"] == "advection_1d" and row["method"] == method
            and row["status"] == "ok"}
    super_rows = {(row["order"], row["resolution"]): row for row in summary
                  if row["case"] == "advection_superconvergence"
                  and row["method"] == method and row["status"] == "ok"}
    if set(base) != set(super_rows):
        raise ValueError(f"e1--e3 and e4--e5 rows differ for {method}")
    groups: dict[int, list[tuple[dict[str, str], dict[str, str]]]] = defaultdict(list)
    for key in sorted(base, key=lambda item: (int(item[0]), number(item[1]))):
        groups[int(key[0])].append((base[key], super_rows[key]))

    destination.parent.mkdir(parents=True, exist_ok=True)
    caption_prefix = "Preliminary development run: " if draft else ""
    lines = [
        "% Generated by scripts/analyze_study.py; do not edit.",
        "\\begin{table}[htbp]", "\\centering", "\\small",
        "\\setlength{\\tabcolsep}{3.0pt}",
        f"\\caption{{{caption_prefix}errors and convergence rates for "
        f"{LABELS[method]} on smooth one-dimensional advection.}}",
        f"\\label{{tab:advection-{method}}}",
        "\\begin{tabular}{c r rr rr rr rr rr}", "\\toprule",
        "$k$ & $N_x$ & $e_1$ & rate & $e_2$ & rate & $e_3$ & rate & "
        "$e_4$ & rate & $e_5$ & rate \\\\", "\\midrule",
    ]
    for group_index, (order, pairs) in enumerate(sorted(groups.items())):
        metric_values = {
            "e1": [number(base_row.get("l1_error")) for base_row, _ in pairs],
            "e2": [number(base_row.get("l2_error")) for base_row, _ in pairs],
            "e3": [number(base_row.get("linf_error")) for base_row, _ in pairs],
            "e4": [number(super_row.get("cell_average_l2_error"))
                   for _, super_row in pairs],
            "e5": [number(super_row.get("downwind_l2_error"))
                   for _, super_row in pairs],
        }
        rates = {metric: _metric_rate(values)
                 for metric, values in metric_values.items()}
        for row_index, (base_row, _) in enumerate(pairs):
            degree = (f"\\multirow{{{len(pairs)}}}{{*}}{{{order}}}"
                      if row_index == 0 else "")
            cells = [degree, base_row["resolution"]]
            for metric in ("e1", "e2", "e3", "e4", "e5"):
                cells.extend((_format_error(metric_values[metric][row_index]),
                              _format_rate(rates[metric][row_index])))
            lines.append(" & ".join(cells) + r" \\")
        if group_index + 1 < len(groups):
            lines.append("\\midrule")
    lines.extend(("\\bottomrule", "\\end{tabular}", "\\end{table}", ""))
    destination.write_text("\n".join(lines))


def write_advection_tables(summary: list[dict[str, str]], table_dir: Path,
                           draft: bool = True) -> None:
    assert_identical_method_rows(summary, "advection_1d")
    assert_identical_method_rows(summary, "advection_superconvergence")
    for method in HEADLINE:
        write_advection_table(summary, method,
                              table_dir / f"advection_1d_{method}.tex", draft)


def write_standard_convergence_table(summary: list[dict[str, str]], case: str,
                                     method: str, destination: Path,
                                     dimension: int = 1,
                                     draft: bool = True) -> None:
    """Write scalar or primary-density L1/L2/Linf pairs for one method."""
    rows = sorted((row for row in summary if row["case"] == case
                   and row["method"] == method and row["kind"] == "convergence"
                   and row["status"] == "ok"), key=_row_key)
    if not rows:
        return
    groups: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[int(row["order"])].append(row)
    destination.parent.mkdir(parents=True, exist_ok=True)
    prefix = "Preliminary development run: " if draft else ""
    equation_word = "density " if case.startswith("euler") else ""
    resolution_header = "$N_x$" if dimension == 1 else "$N_x\\times N_y$"
    lines = [
        "% Generated by scripts/analyze_study.py; do not edit.",
        "\\begin{table}[htbp]", "\\centering", "\\small",
        "\\setlength{\\tabcolsep}{4.0pt}",
        f"\\caption{{{prefix}{equation_word}errors and convergence rates for "
        f"{LABELS[method]} on {case.replace('_', ' ')}.}}",
        f"\\label{{tab:{case}-{method}}}",
        "\\begin{tabular}{c r rr rr rr}", "\\toprule",
        f"$k$ & {resolution_header} & $L^1$ error & rate & $L^2$ error & "
        "rate & $L^\\infty$ error & rate \\\\", "\\midrule",
    ]
    metrics = ("l1_error", "l2_error", "linf_error")
    for group_index, (order, group_rows) in enumerate(sorted(groups.items())):
        values = {metric: [number(row.get(metric)) for row in group_rows]
                  for metric in metrics}
        rates = {metric: _metric_rate(metric_values)
                 for metric, metric_values in values.items()}
        for row_index, row in enumerate(group_rows):
            degree = (f"\\multirow{{{len(group_rows)}}}{{*}}{{{order}}}"
                      if row_index == 0 else "")
            cells = [degree, mesh_label(row["resolution"], dimension)]
            for metric in metrics:
                cells.extend((_format_error(values[metric][row_index]),
                              _format_rate(rates[metric][row_index])))
            lines.append(" & ".join(cells) + r" \\")
        if group_index + 1 < len(groups):
            lines.append("\\midrule")
    lines.extend(("\\bottomrule", "\\end{tabular}", "\\end{table}", ""))
    destination.write_text("\n".join(lines))


def write_standard_tables(summary: list[dict[str, str]], table_dir: Path,
                          draft: bool = True) -> None:
    cases = (("burgers_smooth", 1), ("euler_smooth_1d", 1),
             ("advection_2d_triangular", 2))
    for case, dimension in cases:
        if not any(row["case"] == case for row in summary):
            continue
        assert_identical_method_rows(summary, case)
        for method in HEADLINE:
            write_standard_convergence_table(
                summary, case, method, table_dir / f"{case}_{method}.tex",
                dimension=dimension, draft=draft)


def write_summary(destination: Path, summary: list[dict[str, str]]) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not summary:
        raise ValueError("study input contains no rows")
    with destination.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("quick", "full", "hpc"),
                        default="quick")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--report-mode", choices=("auto", "draft", "final", "none"),
                        default="auto")
    args = parser.parse_args()

    study_dir = ROOT / "measurements" / "study" / args.profile
    source = args.input or study_dir / "results.csv"
    rows = read_rows(source)
    rebase_recorded_paths(rows, source.resolve().parent)
    summary = summarize(rows)
    destination = study_dir / "summary.csv"
    write_summary(destination, summary)

    report_root = report_asset_root(args.profile, args.report_mode)
    report_figures = report_root / "figures" if report_root else None
    report_tables = report_root / "tables" if report_root else None
    figure_dir = study_dir / "figures"
    convergence_figures(summary, figure_dir, report_figures)
    lax_profile_figure(summary, figure_dir, report_figures)
    scale_invariance_figure(summary, figure_dir, report_figures)
    riemann_2d_figures(summary, figure_dir, report_figures)
    if report_tables is not None:
        write_advection_tables(summary, report_tables,
                               draft=args.report_mode != "final")
        write_standard_tables(summary, report_tables,
                              draft=args.report_mode != "final")
    print(f"Wrote {destination} and separated publication-style assets")
    return 0


if __name__ == "__main__":
    sys.exit(main())
