from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np
import tifffile
import yaml

from dapi_norm.image_arrays import read_primary_intensity_plane

CELLPOSE_PLATE_SUMMARY_COLUMNS = [
    "plate",
    "source_id",
    "location",
    "target_channel_id",
    "dapi_channel_id",
    "whole_field_target_integrated_raw",
    "whole_field_target_integrated_raw_per_DAPI_positive_nucleus",
    "whole_field_ch2_integrated_raw",
    "whole_field_ch2_integrated_raw_per_DAPI_positive_nucleus",
    "dapi_positive_nucleus_count",
    "cellpose_masked_target_integrated_raw",
    "cellpose_masked_target_integrated_raw_per_DAPI_positive_nucleus",
    "cellpose_masked_target_integrated_background_corrected",
    "cellpose_masked_target_integrated_background_corrected_per_DAPI_positive_nucleus",
    "cellpose_masked_ch2_integrated_raw",
    "cellpose_masked_ch2_integrated_raw_per_DAPI_positive_nucleus",
    "cellpose_masked_ch2_integrated_background_corrected",
    "cellpose_masked_ch2_integrated_background_corrected_per_DAPI_positive_nucleus",
    "dapi_anchored_cellpose_target_integrated_raw",
    "dapi_anchored_cellpose_target_integrated_raw_per_DAPI_positive_nucleus",
    "dapi_anchored_cellpose_target_integrated_background_corrected",
    "dapi_anchored_cellpose_target_integrated_background_corrected_per_DAPI_positive_nucleus",
    "dapi_anchored_cellpose_ch2_integrated_raw",
    "dapi_anchored_cellpose_ch2_integrated_raw_per_DAPI_positive_nucleus",
    "dapi_anchored_cellpose_ch2_integrated_background_corrected",
    "dapi_anchored_cellpose_ch2_integrated_background_corrected_per_DAPI_positive_nucleus",
    "cellpose_masked_area_px",
    "dapi_anchored_cellpose_masked_area_px",
    "cellpose_masked_area_per_DAPI_positive_nucleus",
    "dapi_anchored_cellpose_masked_area_per_DAPI_positive_nucleus",
    "no_dapi_cellpose_object_count_excluded_in_anchored_variant",
    "cellpose_object_count",
    "qc_status",
    "qc_flags",
    "source_warnings",
    "cellpose_mask_path",
    "dapi_nuclei_mask_path",
    "source_qc_panel_path",
    "source_excluded_signal_check_path",
    "target_path",
    "dapi_path",
    "ch2_path",
    "ch4_path",
]


def build_cellpose_plate_summary(runs_root: Path | str) -> list[dict[str, Any]]:
    root = Path(runs_root)
    rows: list[dict[str, Any]] = []
    for metrics_path in sorted(root.glob("*/*/summaries/cellpose_cell_region_image_metrics.csv")):
        run_dir = metrics_path.parents[1]
        plate = _plate_name(metrics_path.parents[2].name)
        run_name = run_dir.name
        config = _read_yaml(run_dir / "logs" / "config_resolved.yaml")
        image_lookup = {
            str(item.get("source_id") or item.get("location")): item
            for item in config.get("image_inputs", [])
        }
        for metric_row in _read_csv(metrics_path):
            image_id = metric_row["image_id"]
            source_id = _source_id(run_name, metric_row)
            image_record = image_lookup.get(metric_row.get("source_id", "")) or image_lookup.get(image_id)
            if image_record is None:
                raise ValueError(f"No image record for {source_id} in {run_dir / 'logs' / 'config_resolved.yaml'}")
            target_path = _resolve_reference(
                str(image_record.get("target_path", image_record["ch2_path"])),
                run_dir,
            )
            dapi_path = _resolve_reference(
                str(image_record.get("dapi_path", image_record.get("ch4_path", ""))),
                run_dir,
            )
            rows.append(
                _summary_row(
                    plate,
                    source_id,
                    image_id,
                    metric_row,
                    target_path,
                    dapi_path,
                    run_dir,
                    target_channel_id=str(
                        image_record.get(
                            "target_channel_id",
                            metric_row.get("target_channel_id", "CH2"),
                        )
                    ),
                    dapi_channel_id=str(
                        image_record.get("dapi_channel_id", metric_row.get("dapi_channel_id", "CH4"))
                    ),
                )
            )
    return sorted(rows, key=_summary_sort_key)


def write_cellpose_plate_summary_csv(path: Path | str, rows: list[dict[str, Any]]) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CELLPOSE_PLATE_SUMMARY_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return output_path


def write_cellpose_plate_summary_markdown(path: Path | str, rows: list[dict[str, Any]]) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "| "
        + " | ".join(CELLPOSE_PLATE_SUMMARY_COLUMNS)
        + " |",
        "| " + " | ".join(["---"] * len(CELLPOSE_PLATE_SUMMARY_COLUMNS)) + " |",
    ]
    for row in rows:
        values = [_format_markdown_value(row.get(column, "")) for column in CELLPOSE_PLATE_SUMMARY_COLUMNS]
        lines.append("| " + " | ".join(values) + " |")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def _summary_row(
    plate: str,
    source_id: str,
    location: str,
    metric_row: dict[str, str],
    ch2_path: Path,
    ch4_path: Path,
    run_dir: Path,
    target_channel_id: str = "CH2",
    dapi_channel_id: str = "CH4",
) -> dict[str, Any]:
    nuclei = _to_int(metric_row["dapi_positive_nucleus_count"])
    whole_field = _whole_field_raw_ch2(ch2_path)
    candidate_area = _to_int(metric_row["candidate_region_area_px"])
    anchored_area = _to_int(metric_row.get("dapi_anchored_candidate_region_area_px", "0"))
    cellpose_raw = _to_float(metric_row["target_integrated_raw_in_cellpose_region"])
    cellpose_background_corrected = _to_float(
        metric_row.get(
            "target_integrated_background_corrected_in_cellpose_region",
            metric_row["target_integrated_raw_in_cellpose_region"],
        )
    )
    anchored_raw = _to_float(metric_row.get("dapi_anchored_target_integrated_raw", "nan"))
    anchored_background_corrected = _to_float(
        metric_row.get("dapi_anchored_target_integrated_background_corrected", "nan")
    )
    return {
        "plate": plate,
        "source_id": source_id,
        "location": location,
        "target_channel_id": target_channel_id,
        "dapi_channel_id": dapi_channel_id,
        "whole_field_target_integrated_raw": whole_field,
        "whole_field_target_integrated_raw_per_DAPI_positive_nucleus": _per_nucleus(
            whole_field, nuclei
        ),
        "whole_field_ch2_integrated_raw": whole_field,
        "whole_field_ch2_integrated_raw_per_DAPI_positive_nucleus": _per_nucleus(
            whole_field, nuclei
        ),
        "dapi_positive_nucleus_count": nuclei,
        "cellpose_masked_target_integrated_raw": cellpose_raw,
        "cellpose_masked_target_integrated_raw_per_DAPI_positive_nucleus": _per_nucleus(
            cellpose_raw, nuclei
        ),
        "cellpose_masked_target_integrated_background_corrected": cellpose_background_corrected,
        "cellpose_masked_target_integrated_background_corrected_per_DAPI_positive_nucleus": _per_nucleus(
            cellpose_background_corrected, nuclei
        ),
        "cellpose_masked_ch2_integrated_raw": cellpose_raw,
        "cellpose_masked_ch2_integrated_raw_per_DAPI_positive_nucleus": _per_nucleus(
            cellpose_raw, nuclei
        ),
        "dapi_anchored_cellpose_target_integrated_raw": anchored_raw,
        "dapi_anchored_cellpose_target_integrated_raw_per_DAPI_positive_nucleus": _per_nucleus(
            anchored_raw, nuclei
        ),
        "dapi_anchored_cellpose_target_integrated_background_corrected": anchored_background_corrected,
        "dapi_anchored_cellpose_target_integrated_background_corrected_per_DAPI_positive_nucleus": _per_nucleus(
            anchored_background_corrected, nuclei
        ),
        "cellpose_masked_ch2_integrated_background_corrected": cellpose_background_corrected,
        "cellpose_masked_ch2_integrated_background_corrected_per_DAPI_positive_nucleus": _per_nucleus(
            cellpose_background_corrected, nuclei
        ),
        "dapi_anchored_cellpose_ch2_integrated_raw": anchored_raw,
        "dapi_anchored_cellpose_ch2_integrated_raw_per_DAPI_positive_nucleus": _per_nucleus(
            anchored_raw, nuclei
        ),
        "dapi_anchored_cellpose_ch2_integrated_background_corrected": anchored_background_corrected,
        "dapi_anchored_cellpose_ch2_integrated_background_corrected_per_DAPI_positive_nucleus": _per_nucleus(
            anchored_background_corrected, nuclei
        ),
        "cellpose_masked_area_px": candidate_area,
        "dapi_anchored_cellpose_masked_area_px": anchored_area,
        "cellpose_masked_area_per_DAPI_positive_nucleus": (
            candidate_area / nuclei if nuclei else float("nan")
        ),
        "dapi_anchored_cellpose_masked_area_per_DAPI_positive_nucleus": _to_float(
            metric_row.get("dapi_anchored_positive_area_per_DAPI_positive_nucleus", "nan")
        ),
        "no_dapi_cellpose_object_count_excluded_in_anchored_variant": _to_int(
            metric_row.get("dapi_anchored_excluded_no_dapi_object_count", "0")
        ),
        "cellpose_object_count": _to_int(metric_row["cellpose_object_count"]),
        "qc_status": metric_row.get("qc_status", ""),
        "qc_flags": metric_row.get("qc_flags", ""),
        "source_warnings": metric_row.get("warnings", ""),
        "cellpose_mask_path": str(_resolve_reference(metric_row.get("mask_path", ""), run_dir)),
        "dapi_nuclei_mask_path": str(_resolve_reference(metric_row.get("nuclei_mask_path", ""), run_dir)),
        "source_qc_panel_path": _resolve_optional_reference(metric_row.get("qc_panel_path", ""), run_dir),
        "source_excluded_signal_check_path": _resolve_optional_reference(
            metric_row.get("excluded_signal_check_path", ""), run_dir
        ),
        "target_path": str(ch2_path),
        "dapi_path": str(ch4_path),
        "ch2_path": str(ch2_path),
        "ch4_path": str(ch4_path),
    }


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _read_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _source_id(run_name: str, metric_row: dict[str, str]) -> str:
    raw_source_id = metric_row.get("source_id") or metric_row["image_id"]
    if "/" in raw_source_id:
        return raw_source_id
    return f"{run_name}/{raw_source_id}"


def _plate_name(value: str) -> str:
    return value.replace("_", " ").title()


def _summary_sort_key(row: dict[str, Any]) -> tuple[int, str, int]:
    plate_number = 999
    parts = str(row["plate"]).split()
    if parts and parts[-1].isdigit():
        plate_number = int(parts[-1])
    return (plate_number, str(row["source_id"]).rsplit("/", 1)[0], _xy_number(str(row["location"])))


def _xy_number(value: str) -> int:
    digits = "".join(character for character in value if character.isdigit())
    return int(digits) if digits else 99999


def _resolve_reference(value: str, base_dir: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    for candidate in [path, base_dir / path, base_dir.parent / path, base_dir.parent.parent / path]:
        if candidate.exists():
            return candidate
    return path


def _resolve_optional_reference(value: str, base_dir: Path) -> str:
    if not value:
        return ""
    return str(_resolve_reference(value, base_dir))


def _whole_field_raw_ch2(path: Path) -> float:
    try:
        image, _pages = read_primary_intensity_plane(path)
    except Exception:
        image = np.asarray(tifffile.imread(path))
    return float(np.sum(np.asarray(image, dtype=np.float64)))


def _to_float(value: str | int | float) -> float:
    return float(value)


def _to_int(value: str | int | float) -> int:
    return int(float(value))


def _per_nucleus(value: float, nuclei: int) -> float:
    return value / nuclei if nuclei else float("nan")


def _format_markdown_value(value: Any) -> str:
    if isinstance(value, float):
        if np.isnan(value):
            return "nan"
        if abs(value) >= 1_000_000:
            return f"{value:.3e}"
        return f"{value:.6g}"
    return str(value)
