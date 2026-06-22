from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Any

import numpy as np
import tifffile
import yaml


def validate_count_outputs(output_dir: Path | str) -> dict[str, Any]:
    output_path = Path(output_dir)
    summary_path = output_path / "summaries" / "nucleus_counts.csv"
    per_nucleus_path = output_path / "summaries" / "per_nucleus_locations.csv"
    config_path = output_path / "logs" / "config_resolved.yaml"
    run_log_path = output_path / "logs" / "run_log.txt"
    contact_sheet_path = output_path / "qc_contact_sheet.png"

    for required_path in [
        summary_path,
        per_nucleus_path,
        config_path,
        run_log_path,
        contact_sheet_path,
    ]:
        if not required_path.exists():
            raise FileNotFoundError(f"Required output artifact is missing: {required_path}")

    summary_rows = _read_csv(summary_path)
    per_nucleus_rows = _read_csv(per_nucleus_path)
    config = _read_config(config_path)
    if not summary_rows:
        raise ValueError(f"No image rows found in {summary_path}")
    image_shapes = _image_shapes_by_id(config)

    per_image_row_counts: dict[str, int] = {}
    for row in per_nucleus_rows:
        image_id = row["image_id"]
        per_image_row_counts[image_id] = per_image_row_counts.get(image_id, 0) + 1

    total_count = 0
    for row in summary_rows:
        image_id = row["image_id"]
        expected_count = int(row["nucleus_count"])
        total_count += expected_count
        per_rows = per_image_row_counts.get(image_id, 0)
        if per_rows != expected_count:
            raise ValueError(
                f"{image_id} has {per_rows} per-nucleus rows but nucleus_counts.csv reports "
                f"{expected_count}"
            )

        mask_path = _resolve_output_reference(row["mask_path"], output_path)
        qc_path = _resolve_output_reference(row["qc_montage_path"], output_path)
        if not mask_path.exists():
            raise FileNotFoundError(f"Mask listed for {image_id} is missing: {mask_path}")
        if not qc_path.exists():
            raise FileNotFoundError(f"QC montage listed for {image_id} is missing: {qc_path}")
        mask = tifffile.imread(mask_path)
        expected_shape = image_shapes.get(image_id)
        if expected_shape is None:
            raise ValueError(f"No image shape metadata found for {image_id} in {config_path}")
        if tuple(mask.shape) != expected_shape:
            raise ValueError(
                f"{image_id} mask shape {tuple(mask.shape)} does not match extracted image shape "
                f"{expected_shape}"
            )
        label_count = int(np.count_nonzero(np.unique(mask)))
        if label_count != expected_count:
            raise ValueError(
                f"{image_id} mask has {label_count} nonzero labels but CSV reports "
                f"{expected_count}"
            )

    if sum(per_image_row_counts.values()) != total_count:
        raise ValueError(
            "Total per-nucleus rows do not match total nucleus count: "
            f"{sum(per_image_row_counts.values())} != {total_count}"
        )

    _validate_config_model_hash(config, output_path)
    return {
        "summary_rows": len(summary_rows),
        "total_nucleus_count": total_count,
        "per_nucleus_rows": len(per_nucleus_rows),
        "mask_counts_match_csv": True,
    }


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_config(config_path: Path) -> dict[str, Any]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ValueError(f"Expected mapping config in {config_path}")
    return config


def _image_shapes_by_id(config: dict[str, Any]) -> dict[str, tuple[int, int]]:
    records = config.get("image_inputs")
    if not records:
        raise ValueError("No image shape metadata found in config_resolved.yaml")

    shapes: dict[str, tuple[int, int]] = {}
    for record in records:
        if not isinstance(record, dict):
            continue
        image_id = record.get("image_id")
        shape = record.get("extracted_plane_shape")
        if isinstance(image_id, str) and _is_2d_shape(shape):
            shapes[image_id] = (int(shape[0]), int(shape[1]))
    if not shapes:
        raise ValueError("No image shape metadata found in config_resolved.yaml")
    return shapes


def _is_2d_shape(value: Any) -> bool:
    return (
        isinstance(value, list | tuple)
        and len(value) == 2
        and all(isinstance(item, int) for item in value)
    )


def _validate_config_model_hash(config: dict[str, Any], output_path: Path) -> None:
    model = config.get("model", {}) if isinstance(config, dict) else {}
    model_path = model.get("path")
    expected_sha256 = model.get("sha256")
    if not model_path or not expected_sha256:
        return

    path = _resolve_output_reference(model_path, output_path)
    if not path.exists():
        raise FileNotFoundError(f"Model path recorded in config is missing: {path}")
    actual_sha256 = _sha256(path)
    if actual_sha256 != expected_sha256:
        raise ValueError(
            f"Model SHA256 mismatch for {path}: expected {expected_sha256}, got {actual_sha256}"
        )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_output_reference(value: str, output_path: Path) -> Path:
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
