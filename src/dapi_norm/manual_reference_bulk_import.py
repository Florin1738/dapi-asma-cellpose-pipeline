from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np
import tifffile

from dapi_norm.manual_reference_commit import (
    commit_manual_reference_mask,
    validate_manual_reference_mask_commit,
)


BULK_IMPORT_COLUMNS = [
    "image_id",
    "action",
    "status_after",
    "positive_label_count",
    "foreground_area_px",
    "editable_reference_labels_path",
    "notes",
]


def import_raw_annotation_labels(
    *,
    package_dir: Path,
    raw_export_manifest_path: Path,
    output_dir: Path,
    labeler: str,
    completed_date: str | None = None,
    notes: str = "",
    confirm_empty_ids: set[str] | None = None,
    require_all_decisions: bool = False,
) -> dict[str, Path]:
    package_path = Path(package_dir)
    manifest_path = Path(raw_export_manifest_path)
    output_path = Path(output_dir)
    confirm_empty = {_normalize_image_id(value) for value in (confirm_empty_ids or set())}
    raw_rows = _read_csv(manifest_path)
    _validate_raw_manifest_rows(raw_rows, package_dir=package_path)
    decisions = [
        _decision_for_raw_row(row, package_dir=package_path, confirm_empty_ids=confirm_empty)
        for row in raw_rows
    ]
    skipped = [row["image_id"] for row in decisions if row["action"] == "skipped_empty_unconfirmed"]
    if require_all_decisions and skipped:
        raise ValueError("unconfirmed empty labels for image IDs: " + ", ".join(skipped))
    _validate_committable_decisions(decisions, package_dir=package_path)

    summary_rows: list[dict[str, Any]] = []
    for decision in decisions:
        image_id = decision["image_id"]
        labels_path = decision["labels_path"]
        action = decision["action"]
        if action == "skipped_empty_unconfirmed":
            summary_rows.append(
                {
                    "image_id": image_id,
                    "action": action,
                    "status_after": "not_started",
                    "positive_label_count": 0,
                    "foreground_area_px": 0,
                    "editable_reference_labels_path": labels_path,
                    "notes": "empty label image not imported without explicit confirmation",
                }
            )
            continue

        result = commit_manual_reference_mask(
            package_dir=package_path,
            image_id=image_id,
            labels_path=labels_path,
            status="confirmed_empty" if action == "committed_confirmed_empty" else "auto",
            labeler=labeler,
            completed_date=completed_date,
            notes=notes,
        )
        summary_rows.append(
            {
                "image_id": image_id,
                "action": action,
                "status_after": result["status"],
                "positive_label_count": result["positive_label_count"],
                "foreground_area_px": result["foreground_area_px"],
                "editable_reference_labels_path": labels_path,
                "notes": notes,
            }
        )

    summary_path = output_path / "bulk_import_summary.csv"
    _write_csv(summary_path, summary_rows, BULK_IMPORT_COLUMNS)
    return {"summary": summary_path}


def _validate_committable_decisions(
    decisions: list[dict[str, Any]],
    *,
    package_dir: Path,
) -> None:
    for decision in decisions:
        action = decision["action"]
        if action == "skipped_empty_unconfirmed":
            continue
        validate_manual_reference_mask_commit(
            package_dir=package_dir,
            image_id=decision["image_id"],
            labels_path=decision["labels_path"],
            status="confirmed_empty" if action == "committed_confirmed_empty" else "auto",
        )


def _validate_raw_manifest_rows(raw_rows: list[dict[str, str]], *, package_dir: Path) -> None:
    image_ids = [_normalize_image_id(row.get("image_id", "")) for row in raw_rows]
    _require_unique("raw annotation export manifest", image_ids)
    package_manifest_rows = _read_csv(package_dir / "manual_validation_manifest.csv")
    package_status_rows = _read_csv(package_dir / "manual_labeling_status.csv")

    for row, image_id in zip(raw_rows, image_ids, strict=True):
        package_manifest_row = _single_row_for_image(
            package_manifest_rows,
            image_id,
            "manual_validation_manifest.csv",
        )
        package_status_row = _single_row_for_image(
            package_status_rows,
            image_id,
            "manual_labeling_status.csv",
        )
        package_reference_path = _resolve_path(
            package_manifest_row["manual_reference_mask_path"],
            package_dir=package_dir,
        )
        status_reference_path = _resolve_path(
            package_status_row["manual_reference_mask_path"],
            package_dir=package_dir,
        )
        if package_reference_path.resolve() != status_reference_path.resolve():
            raise ValueError(
                f"{image_id} package manifest/status reference path mismatch: "
                f"manifest={package_reference_path.resolve()} "
                f"status={status_reference_path.resolve()}"
            )

        raw_reference_value = row.get("authoritative_reference_mask_path", "").strip()
        if not raw_reference_value:
            raise ValueError(f"{image_id} raw export manifest is missing authoritative_reference_mask_path")
        raw_reference_path = _resolve_path(raw_reference_value, package_dir=package_dir)
        if raw_reference_path.resolve() != package_reference_path.resolve():
            raise ValueError(
                f"{image_id} raw export authoritative reference path mismatch: "
                f"raw_export={raw_reference_path.resolve()} "
                f"package={package_reference_path.resolve()}"
            )

        candidate_included = row.get("candidate_mask_path_included", "").strip().lower()
        if candidate_included and candidate_included not in {"false", "0", "no"}:
            raise ValueError(
                f"{image_id} raw export manifest is not raw-only; "
                "candidate_mask_path_included must be false"
            )


def _decision_for_raw_row(
    row: dict[str, str],
    *,
    package_dir: Path,
    confirm_empty_ids: set[str],
) -> dict[str, Any]:
    image_id = _normalize_image_id(row["image_id"])
    labels_path = _resolve_path(row["editable_reference_labels_path"], package_dir=package_dir)
    if not labels_path.exists():
        raise FileNotFoundError(f"editable reference labels do not exist for {image_id}: {labels_path}")
    labels = np.asarray(tifffile.imread(labels_path))
    positive_label_count = int(np.count_nonzero(np.unique(labels) > 0))
    if positive_label_count:
        action = "committed_non_empty"
    elif image_id in confirm_empty_ids:
        action = "committed_confirmed_empty"
    else:
        action = "skipped_empty_unconfirmed"
    return {
        "image_id": image_id,
        "labels_path": labels_path,
        "action": action,
        "positive_label_count": positive_label_count,
    }


def _resolve_path(value: str, *, package_dir: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    candidates = [
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


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"required CSV does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _single_row_for_image(
    rows: list[dict[str, str]],
    image_id: str,
    source_name: str,
) -> dict[str, str]:
    matches = [row for row in rows if _normalize_image_id(row.get("image_id", "")) == image_id]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {source_name} row for {image_id}, found {len(matches)}")
    return matches[0]


def _require_unique(source: str, image_ids: list[str]) -> None:
    duplicates = sorted({image_id for image_id in image_ids if image_ids.count(image_id) > 1})
    if duplicates:
        raise ValueError(f"duplicate image_id rows in {source}: " + ", ".join(duplicates))


def _write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
