from __future__ import annotations

import csv
import shlex
from pathlib import Path
from typing import Any

import numpy as np
import tifffile

from dapi_norm.image_arrays import read_primary_intensity_plane


RAW_EXPORT_COLUMNS = [
    "image_id",
    "status",
    "ch2_export_path",
    "ch4_export_path",
    "editable_reference_labels_path",
    "authoritative_reference_mask_path",
    "commit_command",
    "confirm_empty_command",
    "candidate_mask_path_included",
]


def prepare_raw_annotation_export(*, package_dir: Path, output_dir: Path) -> dict[str, Path]:
    package_path = Path(package_dir)
    output_path = Path(output_dir)
    _refuse_existing_candidate_files(output_path)
    manifest_rows = _read_csv(package_path / "manual_validation_manifest.csv")
    status_rows = _read_csv(package_path / "manual_labeling_status.csv")
    _validate_manifest_status_match(manifest_rows, status_rows)
    status_by_image = {_image_id(row): row for row in status_rows}

    prepared_rows: list[dict[str, Any]] = []
    for row in manifest_rows:
        image_id = _image_id(row)
        _require_status_path_match_manifest(
            image_id,
            "manual_reference_mask_path",
            _resolve_package_path(row["manual_reference_mask_path"], package_dir=package_path),
            _resolve_package_path(status_by_image[image_id]["manual_reference_mask_path"], package_dir=package_path),
        )
        authoritative_reference = _resolve_package_path(
            row["manual_reference_mask_path"],
            package_dir=package_path,
        )
        ch2, _ = read_primary_intensity_plane(
            _resolve_package_path(row["ch2_path"], package_dir=package_path)
        )
        ch4, _ = read_primary_intensity_plane(
            _resolve_package_path(row["ch4_path"], package_dir=package_path)
        )
        reference = np.asarray(tifffile.imread(authoritative_reference), dtype=np.uint32)
        if ch2.shape != ch4.shape or ch2.shape != reference.shape:
            raise ValueError(
                f"{image_id} raw export shape mismatch: "
                f"ch2={ch2.shape}, ch4={ch4.shape}, reference={reference.shape}"
            )
        prepared_rows.append(
            {
                "image_id": image_id,
                "status": status_by_image[image_id]["status"],
                "authoritative_reference": authoritative_reference,
                "ch2": ch2,
                "ch4": ch4,
                "reference": reference,
            }
        )

    export_rows: list[dict[str, Any]] = []
    for prepared in prepared_rows:
        image_id = prepared["image_id"]
        image_dir = output_path / image_id
        image_dir.mkdir(parents=True, exist_ok=True)
        ch2_export = image_dir / f"{image_id}_CH2_raw.tif"
        ch4_export = image_dir / f"{image_id}_CH4_raw.tif"
        editable_reference = image_dir / f"{image_id}_editable_manual_reference_labels.tif"
        tifffile.imwrite(ch2_export, prepared["ch2"], photometric="minisblack")
        tifffile.imwrite(ch4_export, prepared["ch4"], photometric="minisblack")
        tifffile.imwrite(editable_reference, prepared["reference"], photometric="minisblack")
        export_rows.append(
            {
                "image_id": image_id,
                "status": prepared["status"],
                "ch2_export_path": str(ch2_export),
                "ch4_export_path": str(ch4_export),
                "editable_reference_labels_path": str(editable_reference),
                "authoritative_reference_mask_path": str(prepared["authoritative_reference"].resolve()),
                "commit_command": _commit_command(package_path, image_id, editable_reference),
                "confirm_empty_command": _confirm_empty_command(package_path, image_id, editable_reference),
                "candidate_mask_path_included": False,
            }
        )
    manifest_out = output_path / "raw_annotation_export_manifest.csv"
    readme_out = output_path / "README.md"
    _write_csv(manifest_out, export_rows, RAW_EXPORT_COLUMNS)
    example_image_id = str(export_rows[0]["image_id"]) if export_rows else "XY##"
    _write_readme(
        readme_out,
        package_dir=package_path,
        image_count=len(export_rows),
        example_image_id=example_image_id,
    )
    return {"manifest": manifest_out, "readme": readme_out}


def _write_readme(path: Path, *, package_dir: Path, image_count: int, example_image_id: str) -> None:
    lines = [
        "# Raw-Only Manual Annotation Export",
        "",
        f"Source package: `{package_dir}`",
        f"Images exported: `{image_count}`",
        "",
        "Each image folder contains:",
        "",
        "- `*_CH2_raw.tif`: raw CH2/aSMA image",
        "- `*_CH4_raw.tif`: raw CH4/DAPI image",
        "- `*_editable_manual_reference_labels.tif`: scratch editable copy of the current manual/reference labels",
        "",
        "Automated candidate masks are intentionally not exported in this folder. Use this export for drawing manual/reference labels without seeing the candidate segmentation.",
        "",
        "After editing the scratch label TIFF, replace `YOUR_INITIALS` and commit it with the `scripts/commit_manual_reference_mask.py` command template in `raw_annotation_export_manifest.csv`. Do not manually overwrite reference_masks_to_fill/ or hand-edit manual_labeling_status.csv.",
        f"The per-image commands include explicit IDs, for example `--image-id {example_image_id}`.",
        "",
        "For intentionally empty fields, use the `confirm_empty_command` from the manifest.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _commit_command(package_dir: Path, image_id: str, editable_reference: Path) -> str:
    return shlex.join(
        [
            ".venv/bin/python",
            "scripts/commit_manual_reference_mask.py",
            "--package",
            str(package_dir),
            "--image-id",
            image_id,
            "--labels",
            str(editable_reference),
            "--labeler",
            "YOUR_INITIALS",
            "--notes",
            "manual/reference labels drawn from raw-only export",
        ]
    )


def _confirm_empty_command(package_dir: Path, image_id: str, editable_reference: Path) -> str:
    return shlex.join(
        [
            ".venv/bin/python",
            "scripts/commit_manual_reference_mask.py",
            "--package",
            str(package_dir),
            "--image-id",
            image_id,
            "--labels",
            str(editable_reference),
            "--status",
            "confirmed_empty",
            "--labeler",
            "YOUR_INITIALS",
            "--notes",
            "reviewed raw-only export; no traceable aSMA-associated region",
        ]
    )


def _refuse_existing_candidate_files(output_path: Path) -> None:
    if not output_path.exists():
        return
    candidate_files = [
        path
        for path in output_path.rglob("*")
        if path.is_file() and any("candidate" in part.lower() for part in path.relative_to(output_path).parts)
    ]
    if candidate_files:
        examples = ", ".join(str(path) for path in candidate_files[:5])
        raise ValueError(
            "raw-only export output already contains candidate-looking files; "
            f"remove or choose a clean output directory before export: {examples}"
        )


def _validate_manifest_status_match(
    manifest_rows: list[dict[str, str]],
    status_rows: list[dict[str, str]],
) -> None:
    manifest_ids = [_image_id(row) for row in manifest_rows]
    status_ids = [_image_id(row) for row in status_rows]
    missing = sorted(set(manifest_ids) - set(status_ids))
    extra = sorted(set(status_ids) - set(manifest_ids))
    if missing:
        raise ValueError("missing status rows for manifest images: " + ", ".join(missing))
    if extra:
        raise ValueError("status rows not present in manifest: " + ", ".join(extra))
    _require_unique("manual_validation_manifest.csv", manifest_ids)
    _require_unique("manual_labeling_status.csv", status_ids)


def _require_unique(source: str, image_ids: list[str]) -> None:
    duplicates = sorted({image_id for image_id in image_ids if image_ids.count(image_id) > 1})
    if duplicates:
        raise ValueError(f"duplicate image_id rows in {source}: " + ", ".join(duplicates))


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


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"required CSV does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _resolve_package_path(value: str, *, package_dir: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    candidates = [path, package_dir / path, package_dir.parent / path, Path.cwd() / path]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return package_dir / path


def _image_id(row: dict[str, str]) -> str:
    return row.get("image_id", "").strip().upper().replace(" ", "")
