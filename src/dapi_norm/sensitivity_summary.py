from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np


SENSITIVITY_LONG_COLUMNS = [
    "run_id",
    "image_id",
    "method",
    "foreground_method",
    "per_nucleus",
    "region_fraction",
    "unseeded_foreground_fraction",
    "qc_status",
    "qc_flags",
]

SENSITIVITY_RUN_SUMMARY_COLUMNS = [
    "run_id",
    "n_images",
    "n_reviewable",
    "n_manual_review",
    "n_rejected",
    "expected_order_preserved",
    "challenge_all_rejected",
    "challenge_all_zero_per_nucleus",
    "manual_review_present",
    "hard_reject_present",
    "max_challenge_per_nucleus",
]

SENSITIVITY_IMAGE_SUMMARY_COLUMNS = [
    "image_id",
    "n_runs",
    "min_per_nucleus",
    "max_per_nucleus",
    "mean_per_nucleus",
    "coefficient_of_variation",
    "all_runs_rejected",
    "any_manual_review",
    "any_qc_reject",
]


def build_sensitivity_rows(
    *,
    table_specs: list[tuple[str, Path]],
    positions: list[str] | None = None,
) -> list[dict[str, Any]]:
    selected_positions = None
    if positions is not None:
        selected_positions = {_normalize_position(position) for position in positions}
    rows: list[dict[str, Any]] = []
    for run_id, path in table_specs:
        for source_row in _read_seeded_csv(path):
            image_id = _normalize_position(source_row.get("image_id", ""))
            if selected_positions is not None and image_id not in selected_positions:
                continue
            rows.append(
                {
                    "run_id": run_id,
                    "image_id": image_id,
                    "method": source_row.get("method", ""),
                    "foreground_method": source_row.get("foreground_method", ""),
                    "per_nucleus": _float_or_nan(
                        source_row.get(
                            "seeded_region_intensity_per_DAPI_positive_nucleus"
                        )
                    ),
                    "region_fraction": _float_or_nan(
                        source_row.get("seeded_region_fraction")
                    ),
                    "unseeded_foreground_fraction": _float_or_nan(
                        source_row.get("unseeded_foreground_fraction")
                    ),
                    "qc_status": source_row.get("qc_status", ""),
                    "qc_flags": source_row.get("qc_flags", ""),
                }
            )
    return sorted(rows, key=lambda row: (str(row["run_id"]), _image_sort_key(str(row["image_id"]))))


def summarize_runs(
    rows: list[dict[str, Any]],
    *,
    ordered_positions: list[str],
    challenge_positions: list[str],
) -> list[dict[str, Any]]:
    grouped = _group_by(rows, key="run_id")
    ordered = [_normalize_position(position) for position in ordered_positions]
    challenges = [_normalize_position(position) for position in challenge_positions]
    summary: list[dict[str, Any]] = []
    for run_id in sorted(grouped):
        run_rows = grouped[run_id]
        by_image = {str(row["image_id"]): row for row in run_rows}
        challenge_rows = [by_image[position] for position in challenges if position in by_image]
        challenge_values = [
            _float_or_nan(row.get("per_nucleus"))
            for row in challenge_rows
            if np.isfinite(_float_or_nan(row.get("per_nucleus")))
        ]
        summary.append(
            {
                "run_id": run_id,
                "n_images": len(run_rows),
                "n_reviewable": sum(
                    str(row.get("qc_status", "")) == "reviewable_not_validated"
                    for row in run_rows
                ),
                "n_manual_review": sum(
                    str(row.get("qc_status", "")) == "needs_manual_review"
                    for row in run_rows
                ),
                "n_rejected": sum(
                    str(row.get("qc_status", "")) == "reject_qc_failure"
                    for row in run_rows
                ),
                "expected_order_preserved": _strict_descending_order(by_image, ordered),
                "challenge_all_rejected": bool(challenges)
                and len(challenge_rows) == len(challenges)
                and all(
                    str(row.get("qc_status", "")) == "reject_qc_failure"
                    for row in challenge_rows
                ),
                "challenge_all_zero_per_nucleus": bool(challenges)
                and len(challenge_rows) == len(challenges)
                and all(
                    _is_zero_or_near_zero(row.get("per_nucleus"))
                    for row in challenge_rows
                ),
                "manual_review_present": any(
                    str(row.get("qc_status", "")) == "needs_manual_review"
                    for row in run_rows
                ),
                "hard_reject_present": any(
                    str(row.get("qc_status", "")) == "reject_qc_failure"
                    for row in run_rows
                ),
                "max_challenge_per_nucleus": max(challenge_values)
                if challenge_values
                else float("nan"),
            }
        )
    return summary


def summarize_image_stability(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped = _group_by(rows, key="image_id")
    summary: list[dict[str, Any]] = []
    for image_id in sorted(grouped, key=_image_sort_key):
        image_rows = grouped[image_id]
        values = [
            _float_or_nan(row.get("per_nucleus"))
            for row in image_rows
            if np.isfinite(_float_or_nan(row.get("per_nucleus")))
        ]
        mean = float(np.mean(values)) if values else float("nan")
        std = float(np.std(values, ddof=0)) if values else float("nan")
        coefficient_of_variation = std / mean if values and mean != 0 else float("nan")
        qc_statuses = [str(row.get("qc_status", "")) for row in image_rows]
        summary.append(
            {
                "image_id": image_id,
                "n_runs": len(image_rows),
                "min_per_nucleus": min(values) if values else float("nan"),
                "max_per_nucleus": max(values) if values else float("nan"),
                "mean_per_nucleus": mean,
                "coefficient_of_variation": coefficient_of_variation,
                "all_runs_rejected": bool(qc_statuses)
                and all(status == "reject_qc_failure" for status in qc_statuses),
                "any_manual_review": any(status == "needs_manual_review" for status in qc_statuses),
                "any_qc_reject": any(status == "reject_qc_failure" for status in qc_statuses),
            }
        )
    return summary


def write_sensitivity_outputs(
    rows: list[dict[str, Any]],
    *,
    output_dir: Path,
    ordered_positions: list[str],
    challenge_positions: list[str],
) -> dict[str, Path]:
    if not rows:
        raise ValueError("Cannot write sensitivity outputs for an empty row set")
    output_dir.mkdir(parents=True, exist_ok=True)
    long_csv = output_dir / "sensitivity_long.csv"
    run_summary_csv = output_dir / "sensitivity_run_summary.csv"
    image_summary_csv = output_dir / "sensitivity_image_summary.csv"
    report = output_dir / "sensitivity_report.md"
    plot = output_dir / "sensitivity_per_nucleus_by_method.png"

    run_summary = summarize_runs(
        rows,
        ordered_positions=ordered_positions,
        challenge_positions=challenge_positions,
    )
    image_summary = summarize_image_stability(rows)
    _write_csv(long_csv, rows, SENSITIVITY_LONG_COLUMNS)
    _write_csv(run_summary_csv, run_summary, SENSITIVITY_RUN_SUMMARY_COLUMNS)
    _write_csv(image_summary_csv, image_summary, SENSITIVITY_IMAGE_SUMMARY_COLUMNS)
    _write_report(
        report,
        rows=rows,
        run_summary=run_summary,
        image_summary=image_summary,
        ordered_positions=ordered_positions,
        challenge_positions=challenge_positions,
    )
    _write_plot(plot, rows, challenge_positions=challenge_positions)
    return {
        "long_csv": long_csv,
        "run_summary_csv": run_summary_csv,
        "image_summary_csv": image_summary_csv,
        "report": report,
        "plot": plot,
    }


def _read_seeded_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])
        required = {
            "image_id",
            "seeded_region_intensity_per_DAPI_positive_nucleus",
            "seeded_region_fraction",
            "unseeded_foreground_fraction",
            "qc_status",
        }
        missing = sorted(required - fieldnames)
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
        return list(reader)


def _write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: _format_csv_value(row.get(column, "")) for column in columns})


def _write_report(
    path: Path,
    *,
    rows: list[dict[str, Any]],
    run_summary: list[dict[str, Any]],
    image_summary: list[dict[str, Any]],
    ordered_positions: list[str],
    challenge_positions: list[str],
) -> None:
    lines = [
        "# Region-Restricted Method Sensitivity",
        "",
        "Status: exploratory robustness diagnostic; this is not manual validation.",
        "",
        "Purpose: compare whether region-restricted aSMA-per-nucleus values are stable enough to justify manual review, and whether challenge fields are rejected instead of producing plausible-looking numbers.",
        "",
        f"- source rows: {len(rows)}",
        f"- methods/runs compared: {len(run_summary)}",
        f"- ordering check: {' > '.join(_normalize_position(pos) for pos in ordered_positions) or 'not set'}",
        f"- challenge fields expected to reject: {', '.join(_normalize_position(pos) for pos in challenge_positions) or 'not set'}",
        "",
        "Run-level summary:",
        "",
        "| Run | n | Reviewable | Manual review | Rejected | Expected order | Challenge rejected | Challenge zero | Max challenge per nucleus |",
        "| --- | ---: | ---: | ---: | ---: | --- | --- | --- | ---: |",
    ]
    for row in run_summary:
        lines.append(
            "| {run_id} | {n_images} | {n_reviewable} | {n_manual_review} | {n_rejected} | "
            "{expected_order_preserved} | {challenge_all_rejected} | "
            "{challenge_all_zero_per_nucleus} | {max_challenge_per_nucleus} |".format(
                **_format_report_row(row)
            )
        )
    lines.extend(
        [
            "",
            "Image-level stability:",
            "",
            "| Image | n runs | Min | Max | Mean | CV | All rejected | Any manual review |",
            "| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |",
        ]
    )
    for row in image_summary:
        lines.append(
            "| {image_id} | {n_runs} | {min_per_nucleus} | {max_per_nucleus} | "
            "{mean_per_nucleus} | {coefficient_of_variation} | {all_runs_rejected} | "
            "{any_manual_review} |".format(**_format_report_row(row))
        )
    lines.extend(
        [
            "",
            "Interpretation guardrails:",
            "",
            "- Preserved ordering is a robustness signal only; it does not validate masks or cell boundaries.",
            "- Challenge-field rejection is required before considering a method reviewable for this dataset.",
            "- Precision, recall, F1, and IoU are intentionally absent because no manual ground truth masks were supplied.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_plot(
    path: Path,
    rows: list[dict[str, Any]],
    *,
    challenge_positions: list[str],
) -> None:
    image_ids = sorted({str(row["image_id"]) for row in rows}, key=_image_sort_key)
    run_ids = sorted({str(row["run_id"]) for row in rows})
    by_key = {(str(row["run_id"]), str(row["image_id"])): row for row in rows}
    challenge_set = {_normalize_position(position) for position in challenge_positions}
    expected_ids = [image_id for image_id in image_ids if image_id not in challenge_set]
    challenge_ids = [image_id for image_id in image_ids if image_id in challenge_set]
    panel_defs = []
    if expected_ids:
        panel_defs.append(("Expected ordering fields", expected_ids))
    if challenge_ids:
        panel_defs.append(("Challenge fields", challenge_ids))
    if not panel_defs:
        panel_defs.append(("Fields", image_ids))

    width = min(0.8 / max(len(run_ids), 1), 0.24)
    width_ratios = [max(1, len(ids)) for _title, ids in panel_defs]
    fig, axes = plt.subplots(
        1,
        len(panel_defs),
        figsize=(max(9.0, 1.15 * len(image_ids) + 2.0), 5.4),
        constrained_layout=True,
        width_ratios=width_ratios,
        squeeze=False,
    )
    palette = _palette(len(run_ids))
    method_handles: list[Any] = []
    method_labels: list[str] = []
    for ax, (title, panel_ids) in zip(axes.ravel(), panel_defs, strict=False):
        x = np.arange(len(panel_ids))
        y_values = [
            _float_or_nan(row.get("per_nucleus"))
            for image_id in panel_ids
            for run_id in run_ids
            for row in [by_key.get((run_id, image_id))]
            if row is not None and np.isfinite(_float_or_nan(row.get("per_nucleus")))
        ]
        y_max = max(y_values) if y_values else 1.0
        ax.set_ylim(0, y_max * 1.25 if y_max > 0 else 1.0)
        for run_index, run_id in enumerate(run_ids):
            values: list[float] = []
            colors: list[str] = []
            edgecolors: list[str] = []
            hatches: list[str] = []
            for image_id in panel_ids:
                row = by_key.get((run_id, image_id))
                values.append(_float_or_nan(row.get("per_nucleus")) if row else float("nan"))
                colors.append(palette[run_index])
                qc_status = str(row.get("qc_status", "")) if row else ""
                edgecolors.append(_qc_edge_color(qc_status))
                hatches.append(_qc_hatch(qc_status))
            bars = ax.bar(
                x + (run_index - (len(run_ids) - 1) / 2) * width,
                values,
                width=width,
                label=run_id,
                color=colors,
                edgecolor=edgecolors,
                linewidth=1.1,
                zorder=2,
            )
            for bar, hatch in zip(bars, hatches, strict=False):
                bar.set_hatch(hatch)
        ax.set_xticks(x)
        ax.set_xticklabels(panel_ids, rotation=35, ha="right")
        ax.set_title(title)
        ax.grid(axis="y", color="#dee2e6", linewidth=0.8, zorder=0)
        if not method_handles:
            method_handles, method_labels = ax.get_legend_handles_labels()
    axes.ravel()[0].set_ylabel("CH2/aSMA intensity per DAPI-positive nucleus")
    fig.suptitle("Region-restricted sensitivity by method\nExploratory; no method accepted")
    qc_handles = [
        Patch(facecolor="white", edgecolor="#212529", label="Reviewable, not validated"),
        Patch(facecolor="white", edgecolor="#f08c00", hatch="///", label="Manual review"),
        Patch(facecolor="white", edgecolor="#c92a2a", hatch="xxx", label="QC reject"),
    ]
    fig.legend(
        handles=[*method_handles, *qc_handles],
        labels=[*method_labels, *[handle.get_label() for handle in qc_handles]],
        loc="center left",
        ncols=1,
        bbox_to_anchor=(1.0, 0.5),
        fontsize=8,
    )
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def _strict_descending_order(by_image: dict[str, dict[str, Any]], positions: list[str]) -> bool:
    if len(positions) < 2:
        return False
    values: list[float] = []
    for position in positions:
        row = by_image.get(position)
        if row is None:
            return False
        value = _float_or_nan(row.get("per_nucleus"))
        if not np.isfinite(value):
            return False
        values.append(value)
    return all(left > right for left, right in zip(values, values[1:], strict=False))


def _group_by(rows: list[dict[str, Any]], *, key: str) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row[key])].append(row)
    return dict(grouped)


def _normalize_position(position: str) -> str:
    return position.strip().upper().replace(" ", "")


def _image_sort_key(image_id: str) -> tuple[int, str]:
    digits = "".join(char for char in image_id if char.isdigit())
    return (int(digits) if digits else 10**9, image_id)


def _float_or_nan(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _format_csv_value(value: Any) -> Any:
    if isinstance(value, float) and not np.isfinite(value):
        return ""
    return value


def _format_report_row(row: dict[str, Any]) -> dict[str, Any]:
    formatted = dict(row)
    for key in [
        "max_challenge_per_nucleus",
        "min_per_nucleus",
        "max_per_nucleus",
        "mean_per_nucleus",
        "coefficient_of_variation",
    ]:
        if key in formatted:
            formatted[key] = _format_number(formatted[key])
    for key in [
        "expected_order_preserved",
        "challenge_all_rejected",
        "challenge_all_zero_per_nucleus",
        "all_runs_rejected",
        "any_manual_review",
    ]:
        if key in formatted:
            formatted[key] = "yes" if formatted[key] else "no"
    return formatted


def _format_number(value: Any) -> str:
    numeric = _float_or_nan(value)
    if not np.isfinite(numeric):
        return "NA"
    if numeric == 0:
        return "0"
    if abs(numeric) < 0.01:
        return f"{numeric:.3g}"
    return f"{numeric:.3e}"


def _is_zero_or_near_zero(value: Any) -> bool:
    numeric = _float_or_nan(value)
    return bool(np.isfinite(numeric) and abs(numeric) < 1e-9)


def _palette(count: int) -> list[str]:
    base = ["#1c7ed6", "#0ca678", "#9c36b5", "#f76707", "#1098ad", "#e03131"]
    if count <= len(base):
        return base[:count]
    cmap = plt.get_cmap("tab20")
    return [cmap(index % cmap.N) for index in range(count)]


def _qc_edge_color(qc_status: str) -> str:
    if qc_status == "reject_qc_failure":
        return "#c92a2a"
    if qc_status == "needs_manual_review":
        return "#f08c00"
    return "#212529"


def _qc_hatch(qc_status: str) -> str:
    if qc_status == "reject_qc_failure":
        return "xxx"
    if qc_status == "needs_manual_review":
        return "///"
    return ""
