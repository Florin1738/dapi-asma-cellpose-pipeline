from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
import tifffile

from dapi_norm.manual_annotation_handoff import prepare_manual_annotation_handoff


def test_prepare_manual_annotation_handoff_writes_layer_bundles_without_touching_masks(
    tmp_path: Path,
):
    package_dir = _write_package(tmp_path)
    manual_path = package_dir / "reference_masks_to_fill" / "XY01_manual_reference_labels.tif"
    before = tifffile.imread(manual_path).copy()

    outputs = prepare_manual_annotation_handoff(
        package_dir=package_dir,
        output_dir=tmp_path / "handoff",
    )

    assert outputs["manifest"].exists()
    assert outputs["readme"].exists()
    assert np.array_equal(tifffile.imread(manual_path), before)
    rows = _read_rows(outputs["manifest"])
    assert len(rows) == 1
    row = rows[0]
    assert row["image_id"] == "XY01"
    assert Path(row["layer_bundle_path"]).exists()
    assert Path(row["manual_reference_mask_path"]) == manual_path.resolve()
    bundle = np.load(row["layer_bundle_path"])
    assert set(bundle.files) >= {
        "ch2",
        "ch4",
        "candidate_labels",
        "nuclei_labels",
        "manual_reference_labels",
    }
    assert bundle["ch2"].shape == (8, 9)
    assert bundle["candidate_labels"].dtype == np.uint32
    assert np.count_nonzero(bundle["manual_reference_labels"]) == 0
    readme = outputs["readme"].read_text(encoding="utf-8")
    assert "napari" in readme.lower()
    assert "Do not edit candidate_labels" in readme
    assert "<initials>" not in readme
    assert "YOUR_INITIALS" in readme


def test_prepare_manual_annotation_handoff_rejects_stale_status_paths(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    stale_manual_path = tmp_path / "stale_manual_reference_labels.tif"
    tifffile.imwrite(stale_manual_path, np.zeros((8, 9), dtype=np.uint32))
    status_path = package_dir / "manual_labeling_status.csv"
    rows = _read_rows(status_path)
    rows[0]["manual_reference_mask_path"] = str(stale_manual_path)
    _write_csv(status_path, rows)

    with pytest.raises(ValueError, match="status/manifest path mismatch"):
        prepare_manual_annotation_handoff(
            package_dir=package_dir,
            output_dir=tmp_path / "handoff",
        )


def test_prepare_manual_annotation_handoff_rejects_duplicate_status_image_ids(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    status_path = package_dir / "manual_labeling_status.csv"
    rows = _read_rows(status_path)
    stale_manual_path = tmp_path / "stale_manual_reference_labels.tif"
    tifffile.imwrite(stale_manual_path, np.zeros((8, 9), dtype=np.uint32))
    stale_row = dict(rows[0])
    stale_row["manual_reference_mask_path"] = str(stale_manual_path)
    _write_csv(status_path, [stale_row, rows[0]])

    with pytest.raises(ValueError, match="duplicate image_id rows in manual_labeling_status.csv: XY01"):
        prepare_manual_annotation_handoff(
            package_dir=package_dir,
            output_dir=tmp_path / "handoff",
        )


def test_prepare_manual_annotation_handoff_rejects_duplicate_manifest_image_ids(tmp_path: Path):
    package_dir = _write_package(tmp_path)
    manifest_path = package_dir / "manual_validation_manifest.csv"
    rows = _read_rows(manifest_path)
    _write_csv(manifest_path, [rows[0], rows[0]])

    with pytest.raises(ValueError, match="duplicate image_id rows in manual_validation_manifest.csv: XY01"):
        prepare_manual_annotation_handoff(
            package_dir=package_dir,
            output_dir=tmp_path / "handoff",
        )


def _write_package(tmp_path: Path) -> Path:
    package_dir = tmp_path / "manual_validation" / "package"
    reference_dir = package_dir / "reference_masks_to_fill"
    annotation_dir = package_dir / "annotation_panels_raw_only"
    guide_dir = package_dir / "guide_panels"
    data_dir = tmp_path / "data"
    masks_dir = tmp_path / "masks"
    for path in [reference_dir, annotation_dir, guide_dir, data_dir, masks_dir]:
        path.mkdir(parents=True)
    ch2 = np.arange(72, dtype=np.uint16).reshape(8, 9)
    ch4 = np.flipud(ch2)
    candidate = np.zeros((8, 9), dtype=np.uint32)
    candidate[2:5, 3:6] = 1
    nuclei = np.zeros((8, 9), dtype=np.uint32)
    nuclei[3:4, 4:5] = 7
    manual = np.zeros((8, 9), dtype=np.uint32)
    ch2_path = data_dir / "XY01_CH2.tif"
    ch4_path = data_dir / "XY01_CH4.tif"
    candidate_path = masks_dir / "XY01_candidate_labels.tif"
    nuclei_path = masks_dir / "XY01_nuclei_labels.tif"
    manual_path = reference_dir / "XY01_manual_reference_labels.tif"
    for path, arr in [
        (ch2_path, ch2),
        (ch4_path, ch4),
        (candidate_path, candidate),
        (nuclei_path, nuclei),
        (manual_path, manual),
    ]:
        tifffile.imwrite(path, arr)
    annotation_path = annotation_dir / "XY01_manual_annotation_panel.png"
    guide_path = guide_dir / "XY01_manual_validation_guide.png"
    annotation_path.write_bytes(b"annotation")
    guide_path.write_bytes(b"guide")
    manifest_row = {
        "image_id": "XY01",
        "source_id": "XY01",
        "validation_task": "asma_associated_region",
        "ch2_path": str(ch2_path),
        "ch4_path": str(ch4_path),
        "candidate_mask_path": str(candidate_path),
        "nuclei_mask_path": str(nuclei_path),
        "manual_reference_mask_path": str(manual_path),
        "annotation_panel_path": str(annotation_path),
        "guide_panel_path": str(guide_path),
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
        "manual_reference_mask_path": str(manual_path.resolve()),
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


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))
