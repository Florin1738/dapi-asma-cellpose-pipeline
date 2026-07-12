from __future__ import annotations

import csv
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import tifffile

from dapi_norm.segmentation_validation import evaluate_instance_mask_iou


VALID_COMMIT_STATUSES = {"auto", "complete_non_empty", "confirmed_empty"}
STATUS_COLUMNS = [
    "image_id",
    "manual_reference_mask_path",
    "annotation_panel_path",
    "status",
    "labeler",
    "completed_date",
    "notes",
]


def commit_manual_reference_mask(
    *,
    package_dir: Path,
    image_id: str,
    labels_path: Path,
    status: str = "auto",
    labeler: str = "",
    completed_date: str | None = None,
    notes: str = "",
    npz_key: str = "manual_reference_labels",
) -> dict[str, Any]:
    prepared = _prepare_manual_reference_mask_commit(
        package_dir=package_dir,
        image_id=image_id,
        labels_path=labels_path,
        status=status,
        npz_key=npz_key,
    )

    tifffile.imwrite(prepared["reference_path"], prepared["labels_u32"], photometric="minisblack")
    _update_status_rows(
        prepared["status_rows"],
        image_id=prepared["image_id"],
        status=prepared["status"],
        labeler=labeler,
        completed_date=completed_date or _today_iso(),
        notes=notes,
    )
    _write_csv(prepared["status_path"], prepared["status_rows"])
    return _commit_summary(prepared)


def validate_manual_reference_mask_commit(
    *,
    package_dir: Path,
    image_id: str,
    labels_path: Path,
    status: str = "auto",
    npz_key: str = "manual_reference_labels",
) -> dict[str, Any]:
    prepared = _prepare_manual_reference_mask_commit(
        package_dir=package_dir,
        image_id=image_id,
        labels_path=labels_path,
        status=status,
        npz_key=npz_key,
    )
    return _commit_summary(prepared)


def _prepare_manual_reference_mask_commit(
    *,
    package_dir: Path,
    image_id: str,
    labels_path: Path,
    status: str,
    npz_key: str,
) -> dict[str, Any]:
    package_path = Path(package_dir)
    normalized_image_id = _normalize_image_id(image_id)
    requested_status = status.strip().lower()
    if requested_status not in VALID_COMMIT_STATUSES:
        raise ValueError(
            f"status must be one of {', '.join(sorted(VALID_COMMIT_STATUSES))}; "
            f"got {status!r}"
        )

    manifest_rows = _read_csv(package_path / "manual_validation_manifest.csv")
    status_rows = _read_csv(package_path / "manual_labeling_status.csv")
    manifest_row = _single_row_for_image(manifest_rows, normalized_image_id, "manifest")
    status_row = _single_row_for_image(status_rows, normalized_image_id, "status")
    reference_path = _resolve_package_path(
        status_row["manual_reference_mask_path"],
        package_dir=package_path,
    )
    manifest_reference_path = _resolve_package_path(
        manifest_row["manual_reference_mask_path"],
        package_dir=package_path,
    )
    if reference_path.resolve() != manifest_reference_path.resolve():
        raise ValueError(
            f"{normalized_image_id} status/manifest reference path mismatch: "
            f"status={reference_path.resolve()} manifest={manifest_reference_path.resolve()}"
        )
    if not reference_path.exists():
        raise FileNotFoundError(f"reference mask does not exist: {reference_path}")

    existing_reference = np.asarray(tifffile.imread(reference_path))
    labels = _read_label_image(labels_path, npz_key=npz_key)
    if labels.shape != existing_reference.shape:
        raise ValueError(
            f"{normalized_image_id} shape mismatch: "
            f"edited labels={labels.shape} reference={existing_reference.shape}"
        )
    labels_u32 = _validate_reference_labels(normalized_image_id, labels)
    positive_label_count = _positive_label_count(labels_u32)
    final_status = _final_status(requested_status, positive_label_count)

    return {
        "image_id": normalized_image_id,
        "status": final_status,
        "positive_label_count": positive_label_count,
        "foreground_area_px": int(np.count_nonzero(labels_u32 > 0)),
        "reference_path": reference_path,
        "labels_u32": labels_u32,
        "status_rows": status_rows,
        "status_path": package_path / "manual_labeling_status.csv",
    }


def _commit_summary(prepared: dict[str, Any]) -> dict[str, Any]:
    return {
        "image_id": prepared["image_id"],
        "status": prepared["status"],
        "positive_label_count": prepared["positive_label_count"],
        "foreground_area_px": prepared["foreground_area_px"],
        "manual_reference_mask_path": prepared["reference_path"],
        "status_path": prepared["status_path"],
    }


def _read_label_image(path: Path, *, npz_key: str) -> np.ndarray:
    labels_path = Path(path)
    if not labels_path.exists():
        raise FileNotFoundError(f"edited label image does not exist: {labels_path}")
    if labels_path.suffix.lower() == ".npz":
        with np.load(labels_path) as bundle:
            if npz_key not in bundle.files:
                raise ValueError(f"{labels_path} does not contain NPZ key {npz_key!r}")
            return np.asarray(bundle[npz_key])
    return np.asarray(tifffile.imread(labels_path))


def _validate_reference_labels(image_id: str, labels: np.ndarray) -> np.ndarray:
    array = np.asarray(labels)
    if array.ndim != 2:
        raise ValueError(f"{image_id} reference labels must be a 2-D image; got {array.shape}")
    if np.issubdtype(array.dtype, np.bool_):
        raise ValueError(f"{image_id} reference labels must use an integer dtype, not bool")
    if not np.issubdtype(array.dtype, np.integer):
        raise ValueError(f"{image_id} reference labels must use an integer dtype; got {array.dtype}")
    if np.any(array < 0):
        raise ValueError(f"{image_id} reference labels contain negative values")
    if array.size and int(np.max(array)) > np.iinfo(np.uint32).max:
        raise ValueError(f"{image_id} reference labels contain a value that exceeds uint32")
    array_u32 = array.astype(np.uint32, copy=False)
    try:
        evaluate_instance_mask_iou(
            image_id=image_id,
            candidate_mask=np.zeros(array_u32.shape, dtype=np.uint32),
            reference_mask=array_u32,
            iou_threshold=0.5,
        )
    except ValueError as exc:
        raise ValueError(f"{image_id} invalid instance labels: {exc}") from exc
    return array_u32


def _final_status(requested_status: str, positive_label_count: int) -> str:
    if requested_status == "auto":
        if positive_label_count == 0:
            raise ValueError("empty labels require status='confirmed_empty'")
        return "complete_non_empty"
    if requested_status == "complete_non_empty" and positive_label_count == 0:
        raise ValueError("status='complete_non_empty' requires at least one positive label")
    if requested_status == "confirmed_empty" and positive_label_count > 0:
        raise ValueError("status='confirmed_empty' requires an empty label image")
    return requested_status


def _update_status_rows(
    rows: list[dict[str, str]],
    *,
    image_id: str,
    status: str,
    labeler: str,
    completed_date: str,
    notes: str,
) -> None:
    for row in rows:
        if _normalize_image_id(row["image_id"]) != image_id:
            continue
        row["image_id"] = image_id
        row["status"] = status
        row["labeler"] = labeler
        row["completed_date"] = completed_date
        row["notes"] = notes
        return
    raise KeyError(f"No status row found for {image_id}")


def _positive_label_count(labels: np.ndarray) -> int:
    positive = np.unique(labels)
    return int(np.count_nonzero(positive > 0))


def _single_row_for_image(
    rows: list[dict[str, str]],
    image_id: str,
    source_name: str,
) -> dict[str, str]:
    matches = [row for row in rows if _normalize_image_id(row["image_id"]) == image_id]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {source_name} row for {image_id}, found {len(matches)}")
    return matches[0]


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"required CSV does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    columns = list(rows[0]) if rows else STATUS_COLUMNS
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _resolve_package_path(value: str, *, package_dir: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    candidates = [
        path,
        package_dir / path,
        package_dir.parent / path,
        package_dir.parent.parent / path,
        Path.cwd() / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return Path.cwd() / path


def _normalize_image_id(value: str) -> str:
    return value.strip().upper().replace(" ", "")


def _today_iso() -> str:
    return datetime.now(UTC).date().isoformat()
