from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np
from skimage.measure import regionprops

SUMMARY_COLUMNS = [
    "image_id",
    "input_path",
    "backend",
    "model_name",
    "channel_id",
    "candidate_stain",
    "channel_identity_confirmed",
    "nucleus_count",
    "mask_path",
    "qc_montage_path",
    "warnings",
]

PER_NUCLEUS_COLUMNS = [
    "image_id",
    "input_path",
    "nucleus_id",
    "x_centroid",
    "y_centroid",
    "area_px",
    "bbox_min_row",
    "bbox_min_col",
    "bbox_max_row",
    "bbox_max_col",
    "touches_border",
    "kept_after_filtering",
]


def summarize_labeled_mask(
    *,
    image_id: str,
    input_path: Path,
    mask: np.ndarray,
    backend: str,
    model_name: str,
    channel_id: str,
    candidate_stain: str,
    mask_path: Path | None = None,
    qc_montage_path: Path | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    label_mask = np.asarray(mask)
    if label_mask.ndim != 2:
        raise ValueError("Expected a 2-D labeled mask")

    input_path_text = str(input_path)
    nucleus_rows: list[dict[str, Any]] = []
    for region in regionprops(label_mask):
        min_row, min_col, max_row, max_col = region.bbox
        y_centroid, x_centroid = region.centroid
        nucleus_rows.append(
            {
                "image_id": image_id,
                "input_path": input_path_text,
                "nucleus_id": int(region.label),
                "x_centroid": float(x_centroid),
                "y_centroid": float(y_centroid),
                "area_px": int(region.area),
                "bbox_min_row": int(min_row),
                "bbox_min_col": int(min_col),
                "bbox_max_row": int(max_row),
                "bbox_max_col": int(max_col),
                "touches_border": _touches_border(region.bbox, label_mask.shape),
                "kept_after_filtering": True,
            }
        )

    summary = {
        "image_id": image_id,
        "input_path": input_path_text,
        "backend": backend,
        "model_name": model_name,
        "channel_id": channel_id,
        "candidate_stain": candidate_stain,
        "channel_identity_confirmed": False,
        "nucleus_count": len(nucleus_rows),
        "mask_path": "" if mask_path is None else str(mask_path),
        "qc_montage_path": "" if qc_montage_path is None else str(qc_montage_path),
        "warnings": "channel_identity_unconfirmed",
    }
    return summary, nucleus_rows


def _touches_border(bbox: tuple[int, int, int, int], shape: tuple[int, int]) -> bool:
    min_row, min_col, max_row, max_col = bbox
    height, width = shape
    return min_row == 0 or min_col == 0 or max_row == height or max_col == width


def write_nuclei_count_tables(
    *,
    output_dir: Path,
    summary_rows: list[dict[str, Any]],
    nucleus_rows: list[dict[str, Any]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "nucleus_counts.csv", SUMMARY_COLUMNS, summary_rows)
    _write_csv(output_dir / "per_nucleus_locations.csv", PER_NUCLEUS_COLUMNS, nucleus_rows)


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
