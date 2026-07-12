from __future__ import annotations

import csv
import os
from pathlib import Path

import numpy as np
from PIL import Image
import pytest
import tifffile

from dapi_norm.manual_annotation_audit import run_manual_annotation_audit


def test_run_manual_annotation_audit_reports_ready_and_blocked_images(tmp_path: Path):
    package_dir = _write_annotation_package(tmp_path)

    outputs = run_manual_annotation_audit(
        package_dir=package_dir,
        output_dir=tmp_path / "audit",
    )

    assert outputs["csv"].exists()
    assert outputs["report"].exists()
    assert outputs["contact_sheet"].exists()
    rows = _read_rows(outputs["csv"])
    by_image = {row["image_id"]: row for row in rows}
    assert by_image["XY01"]["validation_ready_image"] == "True"
    assert by_image["XY01"]["positive_label_count"] == "1"
    assert by_image["XY02"]["validation_ready_image"] == "True"
    assert by_image["XY02"]["mask_state"] == "empty"
    assert by_image["XY03"]["validation_ready_image"] == "False"
    assert "status_not_complete" in by_image["XY03"]["blocking_reasons"]
    report = outputs["report"].read_text(encoding="utf-8")
    assert "Overall validation-ready: `False`" in report
    assert "not_started: 1" in report


def test_run_manual_annotation_audit_flags_status_mask_mismatch(tmp_path: Path):
    package_dir = _write_annotation_package(
        tmp_path,
        xy01_status="complete_non_empty",
        xy01_non_empty=False,
    )

    outputs = run_manual_annotation_audit(
        package_dir=package_dir,
        output_dir=tmp_path / "audit",
    )

    rows = _read_rows(outputs["csv"])
    xy01 = {row["image_id"]: row for row in rows}["XY01"]
    assert xy01["validation_ready_image"] == "False"
    assert "complete_non_empty_mask_empty" in xy01["blocking_reasons"]


def test_run_manual_annotation_audit_requires_manifest_status_exact_match(tmp_path: Path):
    package_dir = _write_annotation_package(tmp_path)
    status_path = package_dir / "manual_labeling_status.csv"
    rows = _read_rows(status_path)
    _write_csv(status_path, [row for row in rows if row["image_id"] != "XY03"])

    with pytest.raises(ValueError, match="missing status rows"):
        run_manual_annotation_audit(
            package_dir=package_dir,
            output_dir=tmp_path / "audit_missing",
        )

    rows.append(
        {
            **rows[0],
            "image_id": "XY99",
        }
    )
    _write_csv(status_path, rows)
    with pytest.raises(ValueError, match="status rows not present in manifest"):
        run_manual_annotation_audit(
            package_dir=package_dir,
            output_dir=tmp_path / "audit_extra",
        )


def test_run_manual_annotation_audit_blocks_all_empty_package(tmp_path: Path):
    package_dir = _write_annotation_package(
        tmp_path,
        xy01_status="confirmed_empty",
        xy01_non_empty=False,
        xy03_status="confirmed_empty",
    )

    outputs = run_manual_annotation_audit(
        package_dir=package_dir,
        output_dir=tmp_path / "audit",
    )

    rows = _read_rows(outputs["csv"])
    assert all(row["validation_ready_image"] == "False" for row in rows)
    assert all("package_all_reference_masks_empty" in row["blocking_reasons"] for row in rows)
    report = outputs["report"].read_text(encoding="utf-8")
    assert "positive-reference images: 0" in report


def test_run_manual_annotation_audit_flags_validator_incompatible_masks(tmp_path: Path):
    package_dir = _write_annotation_package(tmp_path)
    xy01_path = package_dir / "reference_masks_to_fill" / "XY01_manual_reference_labels.tif"
    disconnected_single_label = np.zeros((12, 14), dtype=np.uint32)
    disconnected_single_label[1:3, 1:3] = 1
    disconnected_single_label[8:10, 8:10] = 1
    tifffile.imwrite(xy01_path, disconnected_single_label)

    outputs = run_manual_annotation_audit(
        package_dir=package_dir,
        output_dir=tmp_path / "audit",
    )

    xy01 = {row["image_id"]: row for row in _read_rows(outputs["csv"])}["XY01"]
    assert xy01["validation_ready_image"] == "False"
    assert "reference_mask_invalid_instance_labels" in xy01["blocking_reasons"]


def test_run_manual_annotation_audit_flags_annotation_panel_path_mismatch(tmp_path: Path):
    package_dir = _write_annotation_package(tmp_path)
    status_path = package_dir / "manual_labeling_status.csv"
    rows = _read_rows(status_path)
    rows[0]["annotation_panel_path"] = str((package_dir / "wrong_panel.png").resolve())
    _write_csv(status_path, rows)

    outputs = run_manual_annotation_audit(
        package_dir=package_dir,
        output_dir=tmp_path / "audit",
    )

    xy01 = {row["image_id"]: row for row in _read_rows(outputs["csv"])}["XY01"]
    assert xy01["validation_ready_image"] == "False"
    assert "status_manifest_annotation_panel_path_mismatch" in xy01["blocking_reasons"]


def test_run_manual_annotation_audit_flags_shape_mismatch_against_manifest_candidate(
    tmp_path: Path,
):
    package_dir = _write_annotation_package(tmp_path)
    xy01_path = package_dir / "reference_masks_to_fill" / "XY01_manual_reference_labels.tif"
    wrong_shape = np.zeros((5, 5), dtype=np.uint32)
    wrong_shape[1:3, 1:3] = 1
    tifffile.imwrite(xy01_path, wrong_shape)

    outputs = run_manual_annotation_audit(
        package_dir=package_dir,
        output_dir=tmp_path / "audit",
    )

    xy01 = {row["image_id"]: row for row in _read_rows(outputs["csv"])}["XY01"]
    assert xy01["validation_ready_image"] == "False"
    assert "reference_mask_shape_mismatch" in xy01["blocking_reasons"]


def test_run_manual_annotation_audit_resolves_package_relative_paths_from_other_cwd(
    tmp_path: Path,
):
    package_dir = _write_annotation_package(
        tmp_path,
        manifest_paths="package_relative",
    )
    outside_cwd = tmp_path / "outside"
    outside_cwd.mkdir()
    previous_cwd = Path.cwd()
    try:
        os.chdir(outside_cwd)
        outputs = run_manual_annotation_audit(
            package_dir=package_dir.resolve(),
            output_dir=tmp_path / "audit",
        )
    finally:
        os.chdir(previous_cwd)

    rows = _read_rows(outputs["csv"])
    by_image = {row["image_id"]: row for row in rows}
    assert by_image["XY01"]["validation_ready_image"] == "True"
    assert "status_manifest_reference_path_mismatch" not in by_image["XY01"]["blocking_reasons"]
    assert "status_manifest_annotation_panel_path_mismatch" not in by_image["XY01"]["blocking_reasons"]


def test_run_manual_annotation_audit_resolves_repo_relative_paths_from_other_cwd(
    tmp_path: Path,
):
    package_dir = _write_annotation_package(
        tmp_path,
        manifest_paths="repo_relative",
    )
    outside_cwd = tmp_path / "outside"
    outside_cwd.mkdir()
    previous_cwd = Path.cwd()
    try:
        os.chdir(outside_cwd)
        outputs = run_manual_annotation_audit(
            package_dir=package_dir.resolve(),
            output_dir=tmp_path / "audit",
        )
    finally:
        os.chdir(previous_cwd)

    rows = _read_rows(outputs["csv"])
    by_image = {row["image_id"]: row for row in rows}
    assert by_image["XY01"]["validation_ready_image"] == "True"
    assert "status_manifest_reference_path_mismatch" not in by_image["XY01"]["blocking_reasons"]
    assert "status_manifest_annotation_panel_path_mismatch" not in by_image["XY01"]["blocking_reasons"]


def _write_annotation_package(
    tmp_path: Path,
    *,
    xy01_status: str = "complete_non_empty",
    xy01_non_empty: bool = True,
    xy03_status: str = "not_started",
    manifest_paths: str = "absolute",
) -> Path:
    package_dir = tmp_path / "manual_validation" / "package"
    reference_dir = package_dir / "reference_masks_to_fill"
    annotation_dir = package_dir / "annotation_panels_raw_only"
    reference_dir.mkdir(parents=True)
    annotation_dir.mkdir(parents=True)
    manifest_rows = []
    status_rows = []
    specs = [
        ("XY01", xy01_status, xy01_non_empty),
        ("XY02", "confirmed_empty", False),
        ("XY03", xy03_status, False),
    ]
    for image_id, status, non_empty in specs:
        reference = np.zeros((12, 14), dtype=np.uint32)
        if non_empty:
            reference[2:6, 3:8] = 5
        reference_path = reference_dir / f"{image_id}_manual_reference_labels.tif"
        tifffile.imwrite(reference_path, reference)
        candidate_path = package_dir / "candidate_masks" / f"{image_id}_candidate_labels.tif"
        candidate_path.parent.mkdir(exist_ok=True)
        tifffile.imwrite(candidate_path, np.zeros_like(reference, dtype=np.uint32))
        annotation_path = annotation_dir / f"{image_id}_manual_annotation_panel.png"
        Image.new("RGB", (80, 48), color=(30, 30, 30)).save(annotation_path)
        manifest_rows.append(
            {
                "image_id": image_id,
                "source_id": image_id,
                "validation_task": "asma_associated_region",
                "ch2_path": "",
                "ch4_path": "",
                "candidate_mask_path": _manifest_path(
                    candidate_path,
                    package_dir=package_dir,
                    mode=manifest_paths,
                ),
                "nuclei_mask_path": "",
                "manual_reference_mask_path": _manifest_path(
                    reference_path,
                    package_dir=package_dir,
                    mode=manifest_paths,
                ),
                "annotation_panel_path": _manifest_path(
                    annotation_path,
                    package_dir=package_dir,
                    mode=manifest_paths,
                ),
                "guide_panel_path": "",
                "method": "",
                "foreground_method": "",
                "dapi_positive_nucleus_count": "",
                "candidate_integrated_raw": "",
                "candidate_intensity_per_DAPI_positive_nucleus": "",
                "qc_status": "",
                "qc_flags": "",
            }
        )
        status_rows.append(
            {
                "image_id": image_id,
                "manual_reference_mask_path": str(reference_path.resolve()),
                "annotation_panel_path": str(annotation_path.resolve()),
                "status": status,
                "labeler": "tester" if status != "not_started" else "",
                "completed_date": "2026-06-29" if status != "not_started" else "",
                "notes": "",
            }
        )
    _write_csv(package_dir / "manual_validation_manifest.csv", manifest_rows)
    _write_csv(package_dir / "manual_labeling_status.csv", status_rows)
    return package_dir


def _manifest_path(path: Path, *, package_dir: Path, mode: str) -> str:
    if mode == "absolute":
        return str(path)
    if mode == "package_relative":
        return str(path.relative_to(package_dir))
    if mode == "repo_relative":
        return str(path.relative_to(package_dir.parent.parent))
    raise ValueError(f"unknown manifest path mode: {mode}")


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))
