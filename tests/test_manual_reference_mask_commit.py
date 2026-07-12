from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
import tifffile

from dapi_norm.manual_reference_commit import commit_manual_reference_mask


def test_commit_manual_reference_mask_writes_labels_and_updates_status(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    edited_labels = np.zeros((8, 9), dtype=np.uint16)
    edited_labels[2:5, 3:6] = 7
    edited_path = tmp_path / "edited_XY01_labels.tif"
    tifffile.imwrite(edited_path, edited_labels)

    result = commit_manual_reference_mask(
        package_dir=package_dir,
        image_id="xy01",
        labels_path=edited_path,
        labeler="tester",
        completed_date="2026-06-29",
        notes="reviewed in napari",
    )

    reference_path = package_dir / "reference_masks_to_fill" / "XY01_manual_reference_labels.tif"
    committed = tifffile.imread(reference_path)
    assert committed.dtype == np.uint32
    assert np.array_equal(committed, edited_labels.astype(np.uint32))
    assert result["image_id"] == "XY01"
    assert result["status"] == "complete_non_empty"
    assert result["positive_label_count"] == 1

    status_row = _read_status_by_image(package_dir)["XY01"]
    assert status_row["status"] == "complete_non_empty"
    assert status_row["labeler"] == "tester"
    assert status_row["completed_date"] == "2026-06-29"
    assert status_row["notes"] == "reviewed in napari"


def test_commit_manual_reference_mask_can_confirm_empty_reference_explicitly(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    edited_path = tmp_path / "empty_XY01_labels.tif"
    tifffile.imwrite(edited_path, np.zeros((8, 9), dtype=np.uint32))

    result = commit_manual_reference_mask(
        package_dir=package_dir,
        image_id="XY01",
        labels_path=edited_path,
        status="confirmed_empty",
        completed_date="2026-06-29",
    )

    assert result["status"] == "confirmed_empty"
    status_row = _read_status_by_image(package_dir)["XY01"]
    assert status_row["status"] == "confirmed_empty"


def test_commit_manual_reference_mask_rejects_empty_auto_commit(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    edited_path = tmp_path / "empty_XY01_labels.tif"
    tifffile.imwrite(edited_path, np.zeros((8, 9), dtype=np.uint32))

    with pytest.raises(ValueError, match="empty labels require status='confirmed_empty'"):
        commit_manual_reference_mask(
            package_dir=package_dir,
            image_id="XY01",
            labels_path=edited_path,
        )


def test_commit_manual_reference_mask_rejects_shape_mismatch(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    edited_path = tmp_path / "wrong_shape.tif"
    tifffile.imwrite(edited_path, np.ones((4, 4), dtype=np.uint32))

    with pytest.raises(ValueError, match="shape mismatch"):
        commit_manual_reference_mask(
            package_dir=package_dir,
            image_id="XY01",
            labels_path=edited_path,
        )


def test_commit_manual_reference_mask_rejects_disconnected_single_label(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    edited_labels = np.zeros((8, 9), dtype=np.uint32)
    edited_labels[1:3, 1:3] = 2
    edited_labels[6:8, 6:8] = 2
    edited_path = tmp_path / "disconnected_single_label.tif"
    tifffile.imwrite(edited_path, edited_labels)

    with pytest.raises(ValueError, match="invalid instance labels"):
        commit_manual_reference_mask(
            package_dir=package_dir,
            image_id="XY01",
            labels_path=edited_path,
        )


def test_commit_manual_reference_mask_rejects_labels_that_exceed_uint32(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    edited_labels = np.zeros((8, 9), dtype=np.uint64)
    edited_labels[1:3, 1:3] = 1
    edited_labels[5:7, 5:7] = np.iinfo(np.uint32).max + 2
    edited_path = tmp_path / "overflow_label.tif"
    tifffile.imwrite(edited_path, edited_labels)

    with pytest.raises(ValueError, match="exceeds uint32"):
        commit_manual_reference_mask(
            package_dir=package_dir,
            image_id="XY01",
            labels_path=edited_path,
        )


def test_commit_manual_reference_mask_rejects_boolean_masks(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    edited_labels = np.zeros((8, 9), dtype=bool)
    edited_labels[2:5, 3:6] = True
    edited_path = tmp_path / "boolean_labels.tif"
    tifffile.imwrite(edited_path, edited_labels)

    with pytest.raises(ValueError, match="must use an integer dtype, not bool"):
        commit_manual_reference_mask(
            package_dir=package_dir,
            image_id="XY01",
            labels_path=edited_path,
        )


def test_commit_manual_reference_mask_accepts_single_connected_label_one(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    edited_labels = np.zeros((8, 9), dtype=np.uint8)
    edited_labels[2:5, 3:6] = 1
    edited_path = tmp_path / "single_object_label_one.tif"
    tifffile.imwrite(edited_path, edited_labels)

    result = commit_manual_reference_mask(
        package_dir=package_dir,
        image_id="XY01",
        labels_path=edited_path,
    )

    assert result["status"] == "complete_non_empty"
    assert result["positive_label_count"] == 1


def test_commit_manual_reference_mask_loads_npz_manual_reference_layer(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    labels = np.zeros((8, 9), dtype=np.uint32)
    labels[2:6, 2:6] = 3
    bundle_path = tmp_path / "XY01_annotation_layers.npz"
    np.savez_compressed(bundle_path, manual_reference_labels=labels, ch2=np.zeros_like(labels))

    result = commit_manual_reference_mask(
        package_dir=package_dir,
        image_id="XY01",
        labels_path=bundle_path,
        labeler="tester",
        completed_date="2026-06-29",
    )

    assert result["status"] == "complete_non_empty"
    reference_path = package_dir / "reference_masks_to_fill" / "XY01_manual_reference_labels.tif"
    assert np.array_equal(tifffile.imread(reference_path), labels)


def _write_package(tmp_path: Path) -> Path:
    package_dir = tmp_path / "manual_validation" / "package"
    reference_dir = package_dir / "reference_masks_to_fill"
    annotation_dir = package_dir / "annotation_panels_raw_only"
    candidate_dir = package_dir / "candidate_masks"
    for path in [reference_dir, annotation_dir, candidate_dir]:
        path.mkdir(parents=True)
    reference_path = reference_dir / "XY01_manual_reference_labels.tif"
    annotation_path = annotation_dir / "XY01_manual_annotation_panel.png"
    candidate_path = candidate_dir / "XY01_candidate_labels.tif"
    tifffile.imwrite(reference_path, np.zeros((8, 9), dtype=np.uint32))
    tifffile.imwrite(candidate_path, np.zeros((8, 9), dtype=np.uint32))
    annotation_path.write_bytes(b"annotation")
    manifest_row = {
        "image_id": "XY01",
        "source_id": "XY01",
        "validation_task": "asma_associated_region",
        "ch2_path": "",
        "ch4_path": "",
        "candidate_mask_path": str(candidate_path),
        "nuclei_mask_path": "",
        "manual_reference_mask_path": str(reference_path),
        "annotation_panel_path": str(annotation_path),
        "guide_panel_path": "",
        "method": "candidate",
        "foreground_method": "synthetic",
        "dapi_positive_nucleus_count": "1",
        "candidate_integrated_raw": "10",
        "candidate_intensity_per_DAPI_positive_nucleus": "10",
        "qc_status": "reviewable_not_validated",
        "qc_flags": "not_validated_whole_cell_mask",
    }
    status_row = {
        "image_id": "XY01",
        "manual_reference_mask_path": str(reference_path.resolve()),
        "annotation_panel_path": str(annotation_path.resolve()),
        "status": "not_started",
        "labeler": "",
        "completed_date": "",
        "notes": "",
    }
    _write_csv(package_dir / "manual_validation_manifest.csv", [manifest_row])
    _write_csv(package_dir / "manual_labeling_status.csv", [status_row])
    return package_dir


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_status_by_image(package_dir: Path) -> dict[str, dict[str, str]]:
    with (package_dir / "manual_labeling_status.csv").open(newline="", encoding="utf-8") as handle:
        return {row["image_id"]: row for row in csv.DictReader(handle)}
