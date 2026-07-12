from __future__ import annotations

from collections import Counter
import csv
from html import escape
import os
from pathlib import Path
from typing import Any


REVIEW_COLUMNS = [
    "image_id",
    "source_id",
    "qc_status",
    "qc_flags",
    "dapi_positive_nucleus_count",
    "cellpose_object_count",
    "candidate_region_fraction",
    "dapi_nuclei_centroid_coverage_fraction",
    "excluded_region_background_corrected_fraction",
    "target_integrated_intensity_per_DAPI_positive_nucleus",
    "qc_panel_path",
    "excluded_signal_check_path",
]

STATUS_ORDER = {
    "reject_qc_failure": 0,
    "needs_manual_review": 1,
    "reviewable_not_validated": 2,
}


def write_cellpose_review_report(
    *,
    metrics_csv: Path,
    output_dir: Path,
    title: str = "Cellpose CH2+CH4 Candidate Region Visual Review",
) -> dict[str, Path]:
    metrics_path = Path(metrics_csv)
    output_path = Path(output_dir)
    rows = _read_metrics(metrics_path)
    if not rows:
        raise ValueError(f"No Cellpose metrics rows found in {metrics_path}")
    prepared = [_prepare_row(row, metrics_path=metrics_path, output_dir=output_path) for row in rows]
    prepared.sort(key=_sort_key)

    output_path.mkdir(parents=True, exist_ok=True)
    html_path = output_path / "index.html"
    summary_path = output_path / "README.md"
    html_path.write_text(_render_html(prepared, title=title), encoding="utf-8")
    summary_path.write_text(_render_summary(prepared, title=title), encoding="utf-8")
    return {"html": html_path, "summary": summary_path}


def _read_metrics(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Cellpose metrics CSV missing: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    missing = [column for column in REVIEW_COLUMNS if column not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"Cellpose metrics CSV is missing required columns: {', '.join(missing)}")
    return rows


def _prepare_row(
    row: dict[str, str],
    *,
    metrics_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    image_id = _image_id(row)
    panel_path = _resolve_existing_path(row["qc_panel_path"], base_path=metrics_path)
    if not panel_path.exists():
        raise FileNotFoundError(f"QC panel missing for {image_id}: {panel_path}")
    excluded_path = _resolve_existing_path(row["excluded_signal_check_path"], base_path=metrics_path)
    if not excluded_path.exists():
        raise FileNotFoundError(f"excluded-signal panel missing for {image_id}: {excluded_path}")
    return {
        "image_id": image_id,
        "source_id": row.get("source_id", image_id),
        "qc_status": row["qc_status"],
        "qc_flags": row["qc_flags"],
        "dapi_positive_nucleus_count": _int_or_zero(row["dapi_positive_nucleus_count"]),
        "cellpose_object_count": _int_or_zero(row["cellpose_object_count"]),
        "candidate_region_fraction": _float_or_nan(row["candidate_region_fraction"]),
        "dapi_nuclei_centroid_coverage_fraction": _float_or_nan(
            row["dapi_nuclei_centroid_coverage_fraction"]
        ),
        "excluded_region_background_corrected_fraction": _float_or_nan(
            row["excluded_region_background_corrected_fraction"]
        ),
        "target_integrated_intensity_per_DAPI_positive_nucleus": _float_or_nan(
            row["target_integrated_intensity_per_DAPI_positive_nucleus"]
        ),
        "qc_panel_href": _relative_href(panel_path, output_dir),
        "excluded_signal_check_href": _relative_href(excluded_path, output_dir),
    }


def _render_html(rows: list[dict[str, Any]], *, title: str) -> str:
    status_counts = Counter(row["qc_status"] for row in rows)
    cards = "\n".join(_render_card(row) for row in rows)
    status_summary = "\n".join(
        f"<li><strong>{escape(status)}</strong>: {status_counts[status]}</li>"
        for status in sorted(status_counts, key=lambda item: STATUS_ORDER.get(item, 99))
    )
    return (
        "<!doctype html>\n"
        "<html lang=\"en\">\n"
        "<head>\n"
        "  <meta charset=\"utf-8\">\n"
        "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        f"  <title>{escape(title)}</title>\n"
        "  <style>\n"
        "    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 24px; color: #17202a; background: #f6f8fa; }\n"
        "    header, section.card { max-width: 1180px; margin: 0 auto 18px auto; }\n"
        "    h1 { font-size: 28px; margin: 0 0 8px 0; }\n"
        "    .caveat { background: #fff4cc; border: 1px solid #e1b100; padding: 12px 14px; border-radius: 6px; }\n"
        "    .counts { display: flex; gap: 18px; flex-wrap: wrap; padding-left: 20px; }\n"
        "    section.card { background: white; border: 1px solid #d8dee4; border-radius: 8px; padding: 14px; }\n"
        "    .topline { display: flex; align-items: baseline; justify-content: space-between; gap: 12px; border-bottom: 1px solid #eaeef2; padding-bottom: 8px; margin-bottom: 12px; }\n"
        "    .status { font-weight: 700; }\n"
        "    .status.reject_qc_failure { color: #b42318; }\n"
        "    .status.needs_manual_review { color: #a15c00; }\n"
        "    .status.reviewable_not_validated { color: #176a3a; }\n"
        "    .metrics { display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 8px 14px; font-size: 14px; margin-bottom: 12px; }\n"
        "    .metric span { display: block; color: #57606a; font-size: 12px; }\n"
        "    .flags { font-size: 13px; color: #424a53; margin-bottom: 12px; }\n"
        "    .images { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 12px; }\n"
        "    figure { margin: 0; }\n"
        "    figcaption { font-size: 13px; color: #57606a; margin-bottom: 4px; }\n"
        "    img { width: 100%; height: auto; border: 1px solid #d8dee4; background: #fff; }\n"
        "  </style>\n"
        "</head>\n"
        "<body>\n"
        "  <header>\n"
        f"    <h1>{escape(title)}</h1>\n"
        "    <p class=\"caveat\"><strong>Qualitative QC only.</strong> These panels help decide whether the Cellpose CH2+CH4 candidate regions are visually plausible. No precision, recall, F1, or IoU is reported here because real manual/reference masks are still required for quantitative validation.</p>\n"
        f"    <p>Fields reviewed: {len(rows)}</p>\n"
        f"    <ul class=\"counts\">{status_summary}</ul>\n"
        "  </header>\n"
        f"{cards}\n"
        "</body>\n"
        "</html>\n"
    )


def _render_card(row: dict[str, Any]) -> str:
    status = str(row["qc_status"])
    return (
        f"<section class=\"card\" id=\"{escape(row['image_id'])}\">\n"
        "  <div class=\"topline\">\n"
        f"    <h2>{escape(row['image_id'])}</h2>\n"
        f"    <span class=\"status {escape(status)}\">{escape(status)}</span>\n"
        "  </div>\n"
        "  <div class=\"metrics\">\n"
        f"    <div class=\"metric\"><span>DAPI-positive nuclei</span>{row['dapi_positive_nucleus_count']}</div>\n"
        f"    <div class=\"metric\"><span>Cellpose objects</span>{row['cellpose_object_count']}</div>\n"
        f"    <div class=\"metric\"><span>CH2 in candidate region / DAPI-positive nucleus</span>{_format_scientific(row['target_integrated_intensity_per_DAPI_positive_nucleus'])}</div>\n"
        f"    <div class=\"metric\"><span>Candidate region area</span>{_format_percent(row['candidate_region_fraction'])}</div>\n"
        f"    <div class=\"metric\"><span>DAPI centroid coverage</span>{_format_percent(row['dapi_nuclei_centroid_coverage_fraction'])}</div>\n"
        f"    <div class=\"metric\"><span>CH2 outside candidate region</span>{_format_percent(row['excluded_region_background_corrected_fraction'])}</div>\n"
        "  </div>\n"
        f"  <div class=\"flags\"><strong>Flags:</strong> {escape(str(row['qc_flags']).replace(';', '; '))}</div>\n"
        "  <div class=\"images\">\n"
        "    <figure>\n"
        "      <figcaption>Full QC panel: DAPI nuclei, raw CH2, Cellpose candidate region, combined overlay</figcaption>\n"
        f"      <a href=\"{escape(row['qc_panel_href'])}\"><img src=\"{escape(row['qc_panel_href'])}\" alt=\"{escape(row['image_id'])} Cellpose QC panel\"></a>\n"
        "    </figure>\n"
        "    <figure>\n"
        "      <figcaption>Excluded-signal check: magenta marks high displayed CH2 outside the candidate region</figcaption>\n"
        f"      <a href=\"{escape(row['excluded_signal_check_href'])}\"><img src=\"{escape(row['excluded_signal_check_href'])}\" alt=\"{escape(row['image_id'])} excluded signal check\"></a>\n"
        "    </figure>\n"
        "  </div>\n"
        "</section>\n"
    )


def _render_summary(rows: list[dict[str, Any]], *, title: str) -> str:
    status_counts = Counter(row["qc_status"] for row in rows)
    lines = [
        f"# {title}",
        "",
        "Status: qualitative visual QC report; not quantitative validation.",
        "",
        "No precision, recall, F1, or IoU is reported here because completed manual/reference masks are required first.",
        "",
        "## QC Status Counts",
        "",
    ]
    for status in sorted(status_counts, key=lambda item: STATUS_ORDER.get(item, 99)):
        lines.append(f"- {status}: {status_counts[status]}")
    lines.extend(
        [
            "",
            "## Review File",
            "",
            "- `index.html` links each full Cellpose QC panel and excluded-signal check.",
            "",
        ]
    )
    return "\n".join(lines)


def _sort_key(row: dict[str, Any]) -> tuple[int, float, str]:
    status_rank = STATUS_ORDER.get(str(row["qc_status"]), 99)
    excluded = row["excluded_region_background_corrected_fraction"]
    return status_rank, -float(excluded), str(row["image_id"])


def _resolve_existing_path(value: str, *, base_path: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    candidates = [
        path,
        base_path.parent / path,
        base_path.parent.parent / path,
        Path.cwd() / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return Path.cwd() / path


def _relative_href(path: Path, output_dir: Path) -> str:
    return Path(os.path.relpath(path.resolve(), output_dir.resolve())).as_posix()


def _image_id(row: dict[str, str]) -> str:
    return row.get("image_id", row.get("source_id", "")).strip().upper().replace(" ", "")


def _int_or_zero(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def _float_or_nan(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _format_scientific(value: float) -> str:
    if value != value:
        return "NA"
    return f"{value:.3e}"


def _format_percent(value: float) -> str:
    if value != value:
        return "NA"
    return f"{100.0 * value:.1f}%"
