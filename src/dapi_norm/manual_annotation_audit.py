from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image, ImageDraw, ImageFont
import tifffile

from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.segmentation_validation import evaluate_instance_mask_iou


AUDIT_COLUMNS = [
    "image_id",
    "status",
    "labeler",
    "completed_date",
    "manual_reference_mask_path",
    "annotation_panel_path",
    "reference_mask_exists",
    "annotation_panel_exists",
    "mask_shape",
    "mask_dtype",
    "mask_state",
    "positive_label_count",
    "foreground_area_px",
    "status_mask_consistent",
    "package_has_positive_reference",
    "validation_ready_image",
    "blocking_reasons",
]

COMPLETE_STATUSES = {"complete_non_empty", "confirmed_empty"}


def run_manual_annotation_audit(*, package_dir: Path, output_dir: Path) -> dict[str, Path]:
    manifest_path = package_dir / "manual_validation_manifest.csv"
    status_path = package_dir / "manual_labeling_status.csv"
    status_rows = _read_status_rows(status_path)
    manifest_rows = _read_manifest_rows(manifest_path)
    _validate_manifest_status_sets(manifest_rows, status_rows)
    status_by_image = {row["image_id"]: row for row in status_rows}
    audit_rows = [
        _audit_one_manifest_row(
            row,
            status_row=status_by_image[row["image_id"]],
            package_dir=package_dir,
        )
        for row in manifest_rows
    ]
    _apply_package_level_readiness(audit_rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "manual_annotation_audit.csv"
    report_path = output_dir / "manual_annotation_audit_report.md"
    contact_sheet_path = output_dir / "manual_annotation_status_contact_sheet.png"
    _write_csv(csv_path, audit_rows, AUDIT_COLUMNS)
    _write_report(report_path, audit_rows)
    _write_contact_sheet(contact_sheet_path, audit_rows)
    return {
        "csv": csv_path,
        "report": report_path,
        "contact_sheet": contact_sheet_path,
    }


def _audit_one_manifest_row(
    manifest_row: dict[str, str],
    *,
    status_row: dict[str, str],
    package_dir: Path,
) -> dict[str, Any]:
    image_id = manifest_row["image_id"]
    status = status_row["status"].strip().lower()
    reference_mask_path = _resolve_package_path(
        status_row["manual_reference_mask_path"],
        package_dir=package_dir,
    )
    annotation_panel_path = _resolve_package_path(
        status_row["annotation_panel_path"],
        package_dir=package_dir,
    )
    manifest_reference = _resolve_package_path(
        manifest_row["manual_reference_mask_path"],
        package_dir=package_dir,
    )
    manifest_annotation = _resolve_package_path(
        manifest_row["annotation_panel_path"],
        package_dir=package_dir,
    )
    reference_exists = reference_mask_path.exists()
    annotation_exists = annotation_panel_path.exists()
    blocking: list[str] = []
    mask_shape = ""
    mask_dtype = ""
    positive_label_count = 0
    foreground_area_px = 0
    mask_state = "missing"

    if not reference_exists:
        blocking.append("reference_mask_missing")
    else:
        mask = np.asarray(tifffile.imread(reference_mask_path))
        mask_shape = "x".join(str(part) for part in mask.shape)
        mask_dtype = str(mask.dtype)
        positive_labels = np.unique(mask)
        positive_labels = positive_labels[positive_labels > 0]
        positive_label_count = int(len(positive_labels))
        foreground_area_px = int(np.count_nonzero(mask > 0))
        mask_state = "non_empty" if positive_label_count else "empty"
        if mask.ndim != 2:
            blocking.append("reference_mask_not_2d")
        if not np.issubdtype(mask.dtype, np.integer):
            blocking.append("reference_mask_not_integer")
        if mask.ndim == 2 and (
            np.issubdtype(mask.dtype, np.integer) or np.issubdtype(mask.dtype, np.bool_)
        ):
            try:
                evaluate_instance_mask_iou(
                    image_id=image_id,
                    candidate_mask=np.zeros(mask.shape, dtype=np.uint32),
                    reference_mask=mask,
                    iou_threshold=0.5,
                )
            except ValueError:
                blocking.append("reference_mask_invalid_instance_labels")
        expected_shape = _expected_shape_from_manifest(
            manifest_row,
            package_dir=package_dir,
        )
        if expected_shape is not None and mask.shape != expected_shape:
            blocking.append("reference_mask_shape_mismatch")
    if not annotation_exists:
        blocking.append("annotation_panel_missing")
    if status not in COMPLETE_STATUSES:
        blocking.append("status_not_complete")
    if status == "complete_non_empty" and positive_label_count == 0:
        blocking.append("complete_non_empty_mask_empty")
    if status == "confirmed_empty" and positive_label_count > 0:
        blocking.append("confirmed_empty_mask_non_empty")
    if reference_mask_path.resolve() != manifest_reference.resolve():
        blocking.append("status_manifest_reference_path_mismatch")
    if annotation_panel_path.resolve() != manifest_annotation.resolve():
        blocking.append("status_manifest_annotation_panel_path_mismatch")

    status_mask_consistent = not any(
        reason
        in {
            "complete_non_empty_mask_empty",
            "confirmed_empty_mask_non_empty",
            "status_manifest_reference_path_mismatch",
            "status_manifest_annotation_panel_path_mismatch",
            "reference_mask_invalid_instance_labels",
            "reference_mask_shape_mismatch",
        }
        for reason in blocking
    )
    validation_ready = not blocking
    return {
        "image_id": image_id,
        "status": status,
        "labeler": status_row.get("labeler", ""),
        "completed_date": status_row.get("completed_date", ""),
        "manual_reference_mask_path": str(reference_mask_path),
        "annotation_panel_path": str(annotation_panel_path),
        "reference_mask_exists": reference_exists,
        "annotation_panel_exists": annotation_exists,
        "mask_shape": mask_shape,
        "mask_dtype": mask_dtype,
        "mask_state": mask_state,
        "positive_label_count": positive_label_count,
        "foreground_area_px": foreground_area_px,
        "status_mask_consistent": status_mask_consistent,
        "package_has_positive_reference": "",
        "validation_ready_image": validation_ready,
        "blocking_reasons": ";".join(blocking),
    }


def _validate_manifest_status_sets(
    manifest_rows: list[dict[str, str]],
    status_rows: list[dict[str, str]],
) -> None:
    manifest_ids = [row["image_id"] for row in manifest_rows]
    status_ids = [row["image_id"] for row in status_rows]
    duplicate_manifest = sorted({image_id for image_id in manifest_ids if manifest_ids.count(image_id) > 1})
    duplicate_status = sorted({image_id for image_id in status_ids if status_ids.count(image_id) > 1})
    if duplicate_manifest:
        raise ValueError("duplicate manifest rows: " + ", ".join(duplicate_manifest))
    if duplicate_status:
        raise ValueError("duplicate status rows: " + ", ".join(duplicate_status))
    missing_status = sorted(set(manifest_ids) - set(status_ids))
    extra_status = sorted(set(status_ids) - set(manifest_ids))
    if missing_status:
        raise ValueError("missing status rows for manifest images: " + ", ".join(missing_status))
    if extra_status:
        raise ValueError("status rows not present in manifest: " + ", ".join(extra_status))


def _apply_package_level_readiness(rows: list[dict[str, Any]]) -> None:
    package_has_positive_reference = any(int(row["positive_label_count"]) > 0 for row in rows)
    for row in rows:
        row["package_has_positive_reference"] = package_has_positive_reference
        blocking = [reason for reason in str(row["blocking_reasons"]).split(";") if reason]
        if not package_has_positive_reference:
            blocking.append("package_all_reference_masks_empty")
        row["blocking_reasons"] = ";".join(blocking)
        row["validation_ready_image"] = not blocking


def _write_report(path: Path, rows: list[dict[str, Any]]) -> None:
    status_counts = Counter(row["status"] for row in rows)
    ready_count = sum(bool(row["validation_ready_image"]) for row in rows)
    positive_reference_count = sum(int(row["positive_label_count"]) > 0 for row in rows)
    overall_ready = ready_count == len(rows) and bool(rows)
    lines = [
        "# Manual Annotation Audit",
        "",
        f"Overall validation-ready: `{overall_ready}`",
        "",
        "## Summary",
        "",
        f"- images: {len(rows)}",
        f"- validation-ready images: {ready_count}",
        f"- blocked images: {len(rows) - ready_count}",
        f"- positive-reference images: {positive_reference_count}",
        "",
        "Status counts:",
        "",
    ]
    for status, count in sorted(status_counts.items()):
        lines.append(f"- {status}: {count}")
    lines.extend(
        [
            "",
            "## Image Status",
            "",
            "| Image | Status | Ready | Labels | Area px | Blocking reasons |",
            "| --- | --- | --- | ---: | ---: | --- |",
        ]
    )
    for row in rows:
        lines.append(
            "| {image_id} | {status} | {validation_ready_image} | "
            "{positive_label_count} | {foreground_area_px} | {blocking_reasons} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "Interpretation: this audit checks whether manual/reference masks are complete enough to run quantitative validation. It does not validate any candidate segmentation method by itself.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_contact_sheet(path: Path, rows: list[dict[str, Any]]) -> None:
    tiles = [_make_status_tile(row) for row in rows]
    if not tiles:
        Image.new("RGB", (480, 120), "white").save(path)
        return
    tile_width = max(tile.width for tile in tiles)
    tile_height = max(tile.height for tile in tiles)
    columns = min(3, len(tiles))
    rows_n = int(np.ceil(len(tiles) / columns))
    sheet = Image.new("RGB", (tile_width * columns, tile_height * rows_n), "white")
    for index, tile in enumerate(tiles):
        x = (index % columns) * tile_width
        y = (index // columns) * tile_height
        sheet.paste(tile, (x, y))
    sheet.save(path)


def _make_status_tile(row: dict[str, Any]) -> Image.Image:
    panel_path = Path(str(row["annotation_panel_path"]))
    try:
        image = Image.open(panel_path).convert("RGB")
    except (FileNotFoundError, OSError):
        image = Image.new("RGB", (420, 260), "#343a40")
    image.thumbnail((420, 220))
    ready = bool(row["validation_ready_image"])
    color = "#2b8a3e" if ready else "#c92a2a"
    tile = Image.new("RGB", (460, 310), "white")
    tile.paste(image, ((460 - image.width) // 2, 42))
    draw = ImageDraw.Draw(tile)
    draw.rectangle([0, 0, 459, 309], outline=color, width=8)
    font = ImageFont.load_default()
    header = f"{row['image_id']} | {row['status']} | ready={ready}"
    draw.text((14, 12), header, fill=color, font=font)
    footer = (
        f"labels={row['positive_label_count']} area={row['foreground_area_px']} "
        f"{row['blocking_reasons']}"
    )
    draw.text((14, 282), footer[:90], fill="#212529", font=font)
    return tile


def _read_status_rows(path: Path) -> list[dict[str, str]]:
    rows = _read_csv(path)
    required = {"image_id", "manual_reference_mask_path", "annotation_panel_path", "status"}
    _require_columns(path, rows, required)
    return [_normalize_row(row) for row in rows]


def _read_manifest_rows(path: Path) -> list[dict[str, str]]:
    rows = _read_csv(path)
    required = {"image_id", "manual_reference_mask_path", "annotation_panel_path"}
    _require_columns(path, rows, required)
    return [_normalize_row(row) for row in rows]


def _expected_shape_from_manifest(
    manifest_row: dict[str, str],
    *,
    package_dir: Path,
) -> tuple[int, ...] | None:
    for column in ["candidate_mask_path", "nuclei_mask_path"]:
        path_text = manifest_row.get(column, "").strip()
        if not path_text:
            continue
        path = _resolve_package_path(path_text, package_dir=package_dir)
        if not path.exists():
            continue
        return tuple(np.asarray(tifffile.imread(path)).shape)
    ch2_path = manifest_row.get("ch2_path", "").strip()
    if ch2_path:
        path = _resolve_package_path(ch2_path, package_dir=package_dir)
        if path.exists():
            image, _ = read_primary_intensity_plane(path)
            return tuple(image.shape)
    return None


def _resolve_package_path(value: str, *, package_dir: Path) -> Path:
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
    return package_dir / path


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Required CSV missing: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _require_columns(path: Path, rows: list[dict[str, str]], required: set[str]) -> None:
    observed = set(rows[0]) if rows else set()
    missing = sorted(required - observed)
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def _normalize_row(row: dict[str, str]) -> dict[str, str]:
    normalized = dict(row)
    normalized["image_id"] = normalized["image_id"].strip().upper().replace(" ", "")
    normalized["status"] = normalized.get("status", "").strip().lower()
    return normalized


def _write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
