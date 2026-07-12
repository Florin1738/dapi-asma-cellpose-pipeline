from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np


TRIAGE_COLUMNS = [
    "image_id",
    "whole_field_raw_per_nucleus",
    "whole_field_raw_integrated",
    "dapi_nucleus_count",
    "seeded_watershed_per_nucleus",
    "seeded_watershed_region_fraction",
    "seeded_watershed_unseeded_foreground_fraction",
    "seeded_watershed_qc_status",
    "seeded_watershed_qc_flags",
    "seeded_random_walker_per_nucleus",
    "seeded_random_walker_region_fraction",
    "seeded_random_walker_unseeded_foreground_fraction",
    "seeded_random_walker_qc_status",
    "seeded_random_walker_qc_flags",
    "seeded_propagation_per_nucleus",
    "seeded_propagation_region_fraction",
    "seeded_propagation_unseeded_foreground_fraction",
    "seeded_propagation_qc_status",
    "seeded_propagation_qc_flags",
    "cellpose_candidate_per_nucleus",
    "cellpose_candidate_region_fraction",
    "cellpose_excluded_corrected_fraction",
    "cellpose_qc_status",
    "cellpose_qc_flags",
    "missing_sources",
    "any_region_restricted_qc_reject",
    "accepted_region_restricted_method",
    "triage_status",
]


def triage_region_restricted_status(
    qc_statuses: Iterable[str],
    *,
    manual_validation_available: bool,
) -> str:
    statuses = {status for status in qc_statuses if status}
    if not statuses:
        return "region_restricted_sources_missing"
    if statuses == {"reject_qc_failure"}:
        return "all_region_restricted_methods_rejected"
    if "reject_qc_failure" in statuses:
        return "mixed_region_restricted_qc_rejection_not_validated"
    if "needs_manual_review" in statuses:
        return "manual_review_flagged_not_validated"
    if manual_validation_available:
        return "manual_validation_available_review_required"
    return "not_validated_manual_validation_required"


def build_method_triage_rows(
    *,
    pi_metrics_path: Path | None,
    seeded_watershed_path: Path | None,
    seeded_random_walker_path: Path | None,
    seeded_propagation_path: Path | None = None,
    cellpose_paths: list[Path] | None,
    manual_validation_available: bool = False,
    positions: list[str] | None = None,
) -> list[dict[str, Any]]:
    pi_rows = _read_optional_csv_by_key(pi_metrics_path, key_column="loc")
    watershed_rows = _read_optional_csv_by_key(seeded_watershed_path, key_column="image_id")
    random_walker_rows = _read_optional_csv_by_key(
        seeded_random_walker_path,
        key_column="image_id",
    )
    propagation_rows = _read_optional_csv_by_key(seeded_propagation_path, key_column="image_id")
    cellpose_rows: dict[str, dict[str, str]] = {}
    for path in cellpose_paths or []:
        cellpose_rows.update(_read_optional_csv_by_key(path, key_column="image_id"))

    image_ids = _select_image_ids(
        pi_rows=pi_rows,
        watershed_rows=watershed_rows,
        random_walker_rows=random_walker_rows,
        propagation_rows=propagation_rows,
        cellpose_rows=cellpose_rows,
        positions=positions,
    )
    expected_sources = _expected_sources(
        pi_metrics_path=pi_metrics_path,
        seeded_watershed_path=seeded_watershed_path,
        seeded_random_walker_path=seeded_random_walker_path,
        seeded_propagation_path=seeded_propagation_path,
        cellpose_paths=cellpose_paths,
    )
    rows: list[dict[str, Any]] = []
    for image_id in image_ids:
        pi = pi_rows.get(image_id, {})
        watershed = watershed_rows.get(image_id, {})
        random_walker = random_walker_rows.get(image_id, {})
        propagation = propagation_rows.get(image_id, {})
        cellpose = cellpose_rows.get(image_id, {})
        qc_statuses = [
            watershed.get("qc_status", ""),
            random_walker.get("qc_status", ""),
            propagation.get("qc_status", ""),
            cellpose.get("qc_status", ""),
        ]
        triage_status = triage_region_restricted_status(
            qc_statuses,
            manual_validation_available=manual_validation_available,
        )
        row = {
            "image_id": image_id,
            "whole_field_raw_per_nucleus": _float_or_nan(pi.get("raw_per_nuc")),
            "whole_field_raw_integrated": _float_or_nan(pi.get("raw")),
            "dapi_nucleus_count": _int_or_none(pi.get("nuclei")),
            "seeded_watershed_per_nucleus": _float_or_nan(
                watershed.get("seeded_region_intensity_per_DAPI_positive_nucleus")
            ),
            "seeded_watershed_region_fraction": _float_or_nan(
                watershed.get("seeded_region_fraction")
            ),
            "seeded_watershed_unseeded_foreground_fraction": _float_or_nan(
                watershed.get("unseeded_foreground_fraction")
            ),
            "seeded_watershed_qc_status": watershed.get("qc_status", ""),
            "seeded_watershed_qc_flags": watershed.get("qc_flags", ""),
            "seeded_random_walker_per_nucleus": _float_or_nan(
                random_walker.get("seeded_region_intensity_per_DAPI_positive_nucleus")
            ),
            "seeded_random_walker_region_fraction": _float_or_nan(
                random_walker.get("seeded_region_fraction")
            ),
            "seeded_random_walker_unseeded_foreground_fraction": _float_or_nan(
                random_walker.get("unseeded_foreground_fraction")
            ),
            "seeded_random_walker_qc_status": random_walker.get("qc_status", ""),
            "seeded_random_walker_qc_flags": random_walker.get("qc_flags", ""),
            "seeded_propagation_per_nucleus": _float_or_nan(
                propagation.get("seeded_region_intensity_per_DAPI_positive_nucleus")
            ),
            "seeded_propagation_region_fraction": _float_or_nan(
                propagation.get("seeded_region_fraction")
            ),
            "seeded_propagation_unseeded_foreground_fraction": _float_or_nan(
                propagation.get("unseeded_foreground_fraction")
            ),
            "seeded_propagation_qc_status": propagation.get("qc_status", ""),
            "seeded_propagation_qc_flags": propagation.get("qc_flags", ""),
            "cellpose_candidate_per_nucleus": _float_or_nan(
                cellpose.get("target_integrated_intensity_per_DAPI_positive_nucleus")
            ),
            "cellpose_candidate_region_fraction": _float_or_nan(
                cellpose.get("candidate_region_fraction")
            ),
            "cellpose_excluded_corrected_fraction": _float_or_nan(
                cellpose.get("excluded_region_background_corrected_fraction")
            ),
            "cellpose_qc_status": cellpose.get("qc_status", ""),
            "cellpose_qc_flags": cellpose.get("qc_flags", ""),
            "missing_sources": _missing_sources(
                image_id=image_id,
                pi_rows=pi_rows,
                watershed_rows=watershed_rows,
                random_walker_rows=random_walker_rows,
                propagation_rows=propagation_rows,
                cellpose_rows=cellpose_rows,
                expected_sources=expected_sources,
            ),
            "any_region_restricted_qc_reject": any(
                status == "reject_qc_failure" for status in qc_statuses
            ),
            "accepted_region_restricted_method": False,
            "triage_status": triage_status,
        }
        rows.append(row)
    return rows


def write_method_triage_outputs(rows: list[dict[str, Any]], output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "method_triage_summary.csv"
    report_path = output_dir / "method_triage_report.md"
    plot_path = output_dir / "method_triage_per_nucleus_comparison.png"
    _write_csv(csv_path, rows)
    _write_report(report_path, rows)
    _write_plot(plot_path, rows)
    return {"csv": csv_path, "report": report_path, "plot": plot_path}


def _read_csv_by_key(path: Path, *, key_column: str) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if key_column not in (reader.fieldnames or []):
            raise ValueError(f"{path} is missing key column {key_column!r}")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            rows[row[key_column].strip().upper()] = row
        return rows


def _read_optional_csv_by_key(
    path: Path | None,
    *,
    key_column: str,
) -> dict[str, dict[str, str]]:
    if path is None or not path.exists():
        return {}
    return _read_csv_by_key(path, key_column=key_column)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRIAGE_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {column: _format_csv_value(row.get(column, "")) for column in TRIAGE_COLUMNS}
            )


def _write_report(path: Path, rows: list[dict[str, Any]]) -> None:
    total = len(rows)
    rejected = sum(bool(row["any_region_restricted_qc_reject"]) for row in rows)
    accepted = sum(bool(row["accepted_region_restricted_method"]) for row in rows)
    lines = [
        "# Method Triage Report",
        "",
        "Status: exploratory comparison; no manual ground truth supplied.",
        "",
        "No region-restricted method is accepted as validated whole-cell segmentation.",
        "",
        f"- images compared: {total}",
        f"- images with at least one region-restricted QC rejection: {rejected}",
        f"- accepted region-restricted methods: {accepted}",
        "",
        "| Image | Whole-field raw/nucleus | Seeded Otsu | Random walker | Propagation | Cellpose candidate | Any QC reject | Missing sources | Triage |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            "| {image_id} | {whole_field_raw_per_nucleus} | "
            "{seeded_watershed_per_nucleus} | {seeded_random_walker_per_nucleus} | "
            "{seeded_propagation_per_nucleus} | {cellpose_candidate_per_nucleus} | "
            "{any_region_restricted_qc_reject} | {missing_sources} | {triage_status} |".format(
                **_format_report_row(row),
            )
        )
    lines.extend(
        [
            "",
            "Triage status keys:",
            "",
            "- `not_validated_manual_validation_required`: available region-restricted methods are not QC-rejected, but no manual validation supports acceptance.",
            "- `manual_review_flagged_not_validated`: at least one available region-restricted method has a manual-review QC flag but none are hard-rejected.",
            "- `mixed_region_restricted_qc_rejection_not_validated`: at least one available region-restricted method is QC-rejected and at least one is not.",
            "- `all_region_restricted_methods_rejected`: every available region-restricted method is QC-rejected.",
            "- `region_restricted_sources_missing`: no region-restricted source method rows were available for that image.",
            "- `manual_validation_available_review_required`: manual validation was declared available, but this triage still does not auto-accept a method.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("Cannot plot an empty triage table")
    labels = [str(row["image_id"]) for row in rows]
    series = [
        ("Whole-field raw/nucleus", "whole_field_raw_per_nucleus", "#4c6ef5"),
        ("Seeded Otsu", "seeded_watershed_per_nucleus", "#15aabf"),
        ("Random walker", "seeded_random_walker_per_nucleus", "#12b886"),
        ("Propagation", "seeded_propagation_per_nucleus", "#9c36b5"),
        ("Cellpose candidate", "cellpose_candidate_per_nucleus", "#f76707"),
    ]
    x = np.arange(len(rows))
    width = 0.16
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(rows)), 5), constrained_layout=True)
    finite_values: list[float] = []
    for row in rows:
        for _label, key, _color in series:
            value = _float_or_nan(row.get(key))
            if np.isfinite(value):
                finite_values.append(value)
    y_max = max(finite_values) if finite_values else 1.0
    ax.set_ylim(0, y_max * 1.25)
    for index, row in enumerate(rows):
        style = _triage_plot_style(str(row.get("triage_status", "")))
        if style["background"] is not None:
            ax.axvspan(
                index - 0.5,
                index + 0.5,
                color=style["background"],
                alpha=0.35,
                zorder=0,
            )
        if style["annotation"]:
            ax.text(
                index,
                y_max * 1.08,
                style["annotation"],
                ha="center",
                va="bottom",
                fontsize=8,
                color=style["text_color"],
            )
    for index, (label, key, color) in enumerate(series):
        values = [_float_or_nan(row.get(key)) for row in rows]
        ax.bar(
            x + (index - 2.0) * width,
            values,
            width=width,
            label=label,
            color=color,
            edgecolor="white",
            linewidth=0.4,
            zorder=2,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("CH2/aSMA intensity per DAPI-positive nucleus")
    ax.set_title(
        "Per-nucleus CH2/aSMA by method\n"
        "Region-restricted methods exploratory; none accepted"
    )
    method_handles, method_labels = ax.get_legend_handles_labels()
    status_handles = [
        Patch(facecolor="#fff3bf", alpha=0.35, label="No QC reject, not validated"),
        Patch(facecolor="#ffe8cc", alpha=0.35, label="Manual review flagged"),
        Patch(facecolor="#ffc9c9", alpha=0.35, label="Mixed QC rejection"),
        Patch(facecolor="#ffa8a8", alpha=0.35, label="All available methods rejected"),
    ]
    ax.legend(handles=[*method_handles, *status_handles], labels=[*method_labels, *[handle.get_label() for handle in status_handles]], loc="best")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _format_report_row(row: dict[str, Any]) -> dict[str, Any]:
    formatted = dict(row)
    for key in [
        "whole_field_raw_per_nucleus",
        "seeded_watershed_per_nucleus",
        "seeded_random_walker_per_nucleus",
        "seeded_propagation_per_nucleus",
        "cellpose_candidate_per_nucleus",
    ]:
        formatted[key] = _format_report_value(row.get(key))
    formatted["missing_sources"] = row.get("missing_sources") or "none"
    formatted["any_region_restricted_qc_reject"] = (
        "yes" if row.get("any_region_restricted_qc_reject") else "no"
    )
    return formatted


def _triage_plot_style(triage_status: str) -> dict[str, str | None]:
    if triage_status == "not_validated_manual_validation_required":
        return {
            "background": "#fff3bf",
            "annotation": "not\nvalidated",
            "text_color": "#5c4b00",
        }
    if triage_status == "manual_review_flagged_not_validated":
        return {
            "background": "#ffe8cc",
            "annotation": "manual\nreview",
            "text_color": "#d9480f",
        }
    if triage_status == "mixed_region_restricted_qc_rejection_not_validated":
        return {
            "background": "#ffc9c9",
            "annotation": "mixed\nQC reject",
            "text_color": "#c92a2a",
        }
    if triage_status == "all_region_restricted_methods_rejected":
        return {
            "background": "#ffa8a8",
            "annotation": "all\nrejected",
            "text_color": "#a61e4d",
        }
    if triage_status == "region_restricted_sources_missing":
        return {
            "background": "#dee2e6",
            "annotation": "sources\nmissing",
            "text_color": "#495057",
        }
    if triage_status == "manual_validation_available_review_required":
        return {
            "background": "#d0ebff",
            "annotation": "manual\nreview",
            "text_color": "#1864ab",
        }
    return {"background": None, "annotation": None, "text_color": "#212529"}


def _select_image_ids(
    *,
    pi_rows: dict[str, dict[str, str]],
    watershed_rows: dict[str, dict[str, str]],
    random_walker_rows: dict[str, dict[str, str]],
    propagation_rows: dict[str, dict[str, str]],
    cellpose_rows: dict[str, dict[str, str]],
    positions: list[str] | None,
) -> list[str]:
    if positions is not None:
        return [_normalize_position(position) for position in positions]
    return sorted(
        set(pi_rows)
        | set(watershed_rows)
        | set(random_walker_rows)
        | set(propagation_rows)
        | set(cellpose_rows),
        key=_image_sort_key,
    )


def _missing_sources(
    *,
    image_id: str,
    pi_rows: dict[str, dict[str, str]],
    watershed_rows: dict[str, dict[str, str]],
    random_walker_rows: dict[str, dict[str, str]],
    propagation_rows: dict[str, dict[str, str]],
    cellpose_rows: dict[str, dict[str, str]],
    expected_sources: set[str],
) -> str:
    missing: list[str] = []
    if "whole_field" in expected_sources and image_id not in pi_rows:
        missing.append("whole_field")
    if "seeded_watershed" in expected_sources and image_id not in watershed_rows:
        missing.append("seeded_watershed")
    if "seeded_random_walker" in expected_sources and image_id not in random_walker_rows:
        missing.append("seeded_random_walker")
    if "seeded_propagation" in expected_sources and image_id not in propagation_rows:
        missing.append("seeded_propagation")
    if "cellpose" in expected_sources and image_id not in cellpose_rows:
        missing.append("cellpose")
    return ";".join(missing)


def _expected_sources(
    *,
    pi_metrics_path: Path | None,
    seeded_watershed_path: Path | None,
    seeded_random_walker_path: Path | None,
    seeded_propagation_path: Path | None,
    cellpose_paths: list[Path] | None,
) -> set[str]:
    expected: set[str] = set()
    if pi_metrics_path is not None:
        expected.add("whole_field")
    if seeded_watershed_path is not None:
        expected.add("seeded_watershed")
    if seeded_random_walker_path is not None:
        expected.add("seeded_random_walker")
    if seeded_propagation_path is not None:
        expected.add("seeded_propagation")
    if cellpose_paths:
        expected.add("cellpose")
    return expected


def _normalize_position(position: str) -> str:
    return position.strip().upper().replace(" ", "")


def _format_report_value(value: Any) -> str:
    numeric = _float_or_nan(value)
    if not np.isfinite(numeric):
        return "NA"
    return f"{numeric:.3e}"


def _format_csv_value(value: Any) -> Any:
    if isinstance(value, float) and not np.isfinite(value):
        return ""
    return value


def _float_or_nan(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _int_or_none(value: Any) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _image_sort_key(image_id: str) -> tuple[int, str]:
    digits = "".join(char for char in image_id if char.isdigit())
    return (int(digits) if digits else 10**9, image_id)
