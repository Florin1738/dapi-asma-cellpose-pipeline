from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np
import tifffile

from dapi_norm.image_arrays import read_primary_intensity_plane


HANDOFF_COLUMNS = [
    "image_id",
    "validation_task",
    "layer_bundle_path",
    "manual_reference_mask_path",
    "annotation_panel_path",
    "guide_panel_path",
    "ch2_path",
    "ch4_path",
    "candidate_mask_path",
    "nuclei_mask_path",
    "status",
    "bundle_shape",
    "instructions",
]


def prepare_manual_annotation_handoff(*, package_dir: Path, output_dir: Path) -> dict[str, Path]:
    manifest_path = package_dir / "manual_validation_manifest.csv"
    status_path = package_dir / "manual_labeling_status.csv"
    manifest_rows = _read_manifest(manifest_path)
    status_rows = _read_status(status_path)
    _validate_manifest_status_match(manifest_rows, status_rows)
    status_by_image = {row["image_id"]: row for row in status_rows}

    layer_dir = output_dir / "layers_npz"
    layer_dir.mkdir(parents=True, exist_ok=True)
    handoff_rows: list[dict[str, Any]] = []
    for row in manifest_rows:
        image_id = row["image_id"]
        status = status_by_image[image_id]
        package_paths = _resolve_row_paths(row, status, package_dir=package_dir)
        ch2, _ = read_primary_intensity_plane(package_paths["ch2_path"])
        ch4, _ = read_primary_intensity_plane(package_paths["ch4_path"])
        candidate_labels = np.asarray(
            tifffile.imread(package_paths["candidate_mask_path"]),
            dtype=np.uint32,
        )
        nuclei_labels = np.asarray(
            tifffile.imread(package_paths["nuclei_mask_path"]),
            dtype=np.uint32,
        )
        manual_reference_path = package_paths["manual_reference_mask_path"]
        manual_reference_labels = np.asarray(tifffile.imread(manual_reference_path), dtype=np.uint32)
        shapes = {
            "ch2": ch2.shape,
            "ch4": ch4.shape,
            "candidate": candidate_labels.shape,
            "nuclei": nuclei_labels.shape,
            "manual": manual_reference_labels.shape,
        }
        if len(set(shapes.values())) != 1:
            raise ValueError(f"{image_id} layer shape mismatch: {shapes}")
        bundle_path = layer_dir / f"{image_id}_annotation_layers.npz"
        np.savez_compressed(
            bundle_path,
            ch2=ch2,
            ch4=ch4,
            candidate_labels=candidate_labels,
            nuclei_labels=nuclei_labels,
            manual_reference_labels=manual_reference_labels,
        )
        handoff_rows.append(
            {
                "image_id": image_id,
                "validation_task": row.get("validation_task", ""),
                "layer_bundle_path": str(bundle_path),
                "manual_reference_mask_path": str(manual_reference_path.resolve()),
                "annotation_panel_path": str(package_paths["annotation_panel_path"].resolve()),
                "guide_panel_path": str(package_paths["guide_panel_path"].resolve()),
                "ch2_path": str(package_paths["ch2_path"].resolve()),
                "ch4_path": str(package_paths["ch4_path"].resolve()),
                "candidate_mask_path": str(package_paths["candidate_mask_path"].resolve()),
                "nuclei_mask_path": str(package_paths["nuclei_mask_path"].resolve()),
                "status": status["status"],
                "bundle_shape": "x".join(str(part) for part in ch2.shape),
                "instructions": "edit manual_reference_labels only; do not edit candidate_labels",
            }
        )

    manifest_out = output_dir / "annotation_handoff_manifest.csv"
    readme_out = output_dir / "README.md"
    _write_csv(manifest_out, handoff_rows, HANDOFF_COLUMNS)
    _write_readme(readme_out, package_dir=package_dir, image_count=len(handoff_rows))
    return {
        "manifest": manifest_out,
        "readme": readme_out,
        "layer_dir": layer_dir,
    }


def _resolve_row_paths(
    row: dict[str, str],
    status: dict[str, str],
    *,
    package_dir: Path,
) -> dict[str, Path]:
    paths = {
        "ch2_path": _resolve_package_path(row["ch2_path"], package_dir=package_dir),
        "ch4_path": _resolve_package_path(row["ch4_path"], package_dir=package_dir),
        "candidate_mask_path": _resolve_package_path(row["candidate_mask_path"], package_dir=package_dir),
        "nuclei_mask_path": _resolve_package_path(row["nuclei_mask_path"], package_dir=package_dir),
        "manual_reference_mask_path": _resolve_package_path(row["manual_reference_mask_path"], package_dir=package_dir),
        "annotation_panel_path": _resolve_package_path(row["annotation_panel_path"], package_dir=package_dir),
        "guide_panel_path": _resolve_package_path(row["guide_panel_path"], package_dir=package_dir),
    }
    _require_status_path_match_manifest(
        row["image_id"],
        "manual_reference_mask_path",
        paths["manual_reference_mask_path"],
        _resolve_package_path(status["manual_reference_mask_path"], package_dir=package_dir),
    )
    _require_status_path_match_manifest(
        row["image_id"],
        "annotation_panel_path",
        paths["annotation_panel_path"],
        _resolve_package_path(status["annotation_panel_path"], package_dir=package_dir),
    )
    return paths


def _require_status_path_match_manifest(
    image_id: str,
    column: str,
    manifest_path: Path,
    status_path: Path,
) -> None:
    if manifest_path.resolve() != status_path.resolve():
        raise ValueError(
            f"{image_id} status/manifest path mismatch for {column}: "
            f"manifest={manifest_path.resolve()} status={status_path.resolve()}"
        )


def _write_readme(path: Path, *, package_dir: Path, image_count: int) -> None:
    lines = [
        "# Manual Annotation Handoff",
        "",
        f"Source package: `{package_dir}`",
        f"Images bundled: `{image_count}`",
        "",
        "Each `layers_npz/*_annotation_layers.npz` contains:",
        "",
        "- `ch2`: raw CH2/aSMA image",
        "- `ch4`: raw CH4/DAPI image",
        "- `candidate_labels`: automated candidate segmentation for comparison only",
        "- `nuclei_labels`: DAPI nucleus labels",
        "- `manual_reference_labels`: current editable manual/reference labels",
        "",
        "Use these bundles in napari or another label-editing tool. Do not edit candidate_labels or nuclei_labels. Edit only manual_reference_labels as an integer instance-label image, export it to a scratch TIFF or NPZ, then commit it back with `scripts/commit_manual_reference_mask.py`. Do not manually overwrite `reference_masks_to_fill/` or hand-edit `manual_labeling_status.csv`.",
        "",
        "Recommended safe commit command after exporting an edited label TIFF or NPZ:",
        "",
        "```bash",
        ".venv/bin/python scripts/commit_manual_reference_mask.py \\",
        "  --package manual_validation/package \\",
        "  --image-id XY22 \\",
        "  --labels path/to/edited_manual_reference_labels.tif \\",
        "  --labeler YOUR_INITIALS \\",
        '  --notes "brief annotation note"',
        "```",
        "",
        "For an intentionally empty field, pass `--status confirmed_empty`; the default `auto` status refuses to commit an empty mask so an accidental blank export is not silently treated as reviewed. Boolean exports are rejected. Integer label images must use distinct IDs for distinct objects; a single connected object may use label `1`.",
        "",
        "After committing edited labels, rerun the annotation audit, and only then run quantitative validation.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _validate_manifest_status_match(
    manifest_rows: list[dict[str, str]],
    status_rows: list[dict[str, str]],
) -> None:
    manifest_ids = [row["image_id"] for row in manifest_rows]
    status_ids = [row["image_id"] for row in status_rows]
    _require_unique_image_ids("manual_validation_manifest.csv", manifest_ids)
    _require_unique_image_ids("manual_labeling_status.csv", status_ids)
    missing = sorted(set(manifest_ids) - set(status_ids))
    extra = sorted(set(status_ids) - set(manifest_ids))
    if missing:
        raise ValueError("missing status rows for manifest images: " + ", ".join(missing))
    if extra:
        raise ValueError("status rows not present in manifest: " + ", ".join(extra))


def _require_unique_image_ids(source_name: str, image_ids: list[str]) -> None:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for image_id in image_ids:
        if image_id in seen:
            duplicates.add(image_id)
        seen.add(image_id)
    if duplicates:
        raise ValueError(
            f"duplicate image_id rows in {source_name}: " + ", ".join(sorted(duplicates))
        )


def _read_manifest(path: Path) -> list[dict[str, str]]:
    rows = _read_csv(path)
    required = {
        "image_id",
        "validation_task",
        "ch2_path",
        "ch4_path",
        "candidate_mask_path",
        "nuclei_mask_path",
        "manual_reference_mask_path",
        "annotation_panel_path",
        "guide_panel_path",
    }
    _require_columns(path, rows, required)
    return [_normalize_row(row) for row in rows]


def _read_status(path: Path) -> list[dict[str, str]]:
    rows = _read_csv(path)
    required = {"image_id", "manual_reference_mask_path", "annotation_panel_path", "status"}
    _require_columns(path, rows, required)
    return [_normalize_row(row) for row in rows]


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


def _write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
