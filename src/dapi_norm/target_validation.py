from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np
import tifffile

from dapi_norm.image_arrays import read_primary_intensity_plane

REQUIRED_PLOTS = [
    "normalized_intensity_by_well.png",
    "target_integrated_vs_nucleus_count.png",
]


def validate_target_outputs(output_dir: Path | str) -> dict[str, Any]:
    output_path = Path(output_dir)
    summary_path = output_path / "summaries" / "image_level_summary.csv"
    well_path = output_path / "summaries" / "well_level_summary.csv"
    config_path = output_path / "logs" / "config_resolved.yaml"
    run_log_path = output_path / "logs" / "run_log.txt"
    contact_sheet_path = output_path / "qc_contact_sheet.png"
    for required_path in [summary_path, well_path, config_path, run_log_path, contact_sheet_path]:
        if not required_path.exists():
            raise FileNotFoundError(
                "Required target per-DAPI-positive-nucleus artifact missing: "
                f"{required_path}"
            )

    rows = _read_csv(summary_path)
    well_rows = _read_csv(well_path)
    _validate_well_rows(rows, well_rows)
    if not rows:
        raise ValueError(f"No rows found in {summary_path}")

    for plot_name in REQUIRED_PLOTS:
        plot_path = output_path / "plots" / plot_name
        if not plot_path.exists():
            raise FileNotFoundError(f"Required plot missing: {plot_path}")

    for row in rows:
        _validate_row(row, output_path)

    return {
        "summary_rows": len(rows),
        "plots_exist": True,
        "qc_overlays_exist": True,
        "formulas_match": True,
    }


def _validate_row(row: dict[str, str], output_path: Path) -> None:
    image_id = row["image_id"]
    target_path = _resolve_reference(row["input_path"], output_path)
    mask_path = _resolve_reference(row["mask_path"], output_path)
    qc_path = _resolve_reference(row["qc_overlay_path"], output_path)
    for path in [target_path, mask_path, qc_path]:
        if not path.exists():
            raise FileNotFoundError(f"{image_id} referenced artifact missing: {path}")

    target, _page_count = read_primary_intensity_plane(target_path)
    mask = tifffile.imread(mask_path)
    mask_label_count = int(np.count_nonzero(np.unique(mask)))
    raw_count = int(row["raw_nucleus_count"])
    filtered_count = int(row["filtered_nucleus_count"])
    if mask_label_count != raw_count or mask_label_count != filtered_count:
        raise ValueError(
            f"{image_id} mask label count {mask_label_count} does not match "
            f"raw_nucleus_count={raw_count} and filtered_nucleus_count={filtered_count}"
        )

    target_float = target.astype(np.float64)
    background_percentile = _background_percentile_from_method(row["background_method"])
    background_value = float(np.percentile(target_float, background_percentile))
    corrected_sum = float(np.sum(np.clip(target_float - background_value, 0, None)))
    raw_sum = float(np.sum(target_float))
    mean_raw = float(np.mean(target_float))
    nucleus_count = filtered_count
    per_nucleus_endpoint = corrected_sum / nucleus_count if nucleus_count > 0 else float("nan")

    _assert_close(image_id, "background_value_per_px", background_value, row["background_value_per_px"])
    _assert_close(image_id, "target_integrated_raw", raw_sum, row["target_integrated_raw"])
    _assert_close(image_id, "target_mean_raw", mean_raw, row["target_mean_raw"])
    _assert_close(
        image_id,
        "target_integrated_background_corrected",
        corrected_sum,
        row["target_integrated_background_corrected"],
    )
    _assert_close(
        image_id,
        "target_integrated_intensity_per_DAPI_positive_nucleus",
        per_nucleus_endpoint,
        row["target_integrated_intensity_per_DAPI_positive_nucleus"],
    )


def _assert_close(image_id: str, label: str, expected: float, observed_text: str) -> None:
    observed = float(observed_text)
    if not np.isclose(observed, expected, rtol=1e-9, atol=1e-6):
        raise ValueError(f"{image_id} {label} mismatch: expected {expected}, observed {observed}")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _validate_well_rows(
    image_rows: list[dict[str, str]], well_rows: list[dict[str, str]]
) -> None:
    if len(image_rows) != len(well_rows):
        raise ValueError(
            f"image_level_summary rows ({len(image_rows)}) do not match well_level_summary rows "
            f"({len(well_rows)})"
        )
    image_by_id = {row["image_id"]: row for row in image_rows}
    well_by_id = {row["image_id"]: row for row in well_rows}
    if image_by_id.keys() != well_by_id.keys():
        raise ValueError("well_level_summary image IDs do not match image_level_summary")
    for image_id, image_row in image_by_id.items():
        if image_row != well_by_id[image_id]:
            raise ValueError(f"well_level_summary row diverges from image_level_summary for {image_id}")


def _background_percentile_from_method(method: str) -> float:
    prefix = "percentile_"
    if not method.startswith(prefix):
        raise ValueError(f"Unsupported background_method: {method}")
    return float(method.removeprefix(prefix).replace("p", "."))


def _resolve_reference(value: str, output_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    candidates = [
        path,
        output_path / path,
        output_path.parent / path,
        output_path.parent.parent / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[-1]
