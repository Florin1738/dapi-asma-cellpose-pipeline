from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
import tifffile

from dapi_norm.manual_reference_bulk_import import import_raw_annotation_labels


def test_import_raw_annotation_labels_commits_non_empty_and_confirmed_empty(
    tmp_path: Path,
):
    package_dir, manifest_path = _write_package_and_raw_export(tmp_path)
    xy01_labels = np.zeros((8, 9), dtype=np.uint32)
    xy01_labels[2:5, 3:6] = 4
    tifffile.imwrite(
        package_dir / "raw_annotation_exports" / "XY01" / "XY01_editable_manual_reference_labels.tif",
        xy01_labels,
    )

    outputs = import_raw_annotation_labels(
        package_dir=package_dir,
        raw_export_manifest_path=manifest_path,
        output_dir=tmp_path / "import",
        labeler="tester",
        completed_date="2026-06-29",
        notes="bulk import test",
        confirm_empty_ids={"XY02"},
    )

    rows = _read_rows(outputs["summary"])
    by_image = {row["image_id"]: row for row in rows}
    assert by_image["XY01"]["action"] == "committed_non_empty"
    assert by_image["XY01"]["status_after"] == "complete_non_empty"
    assert by_image["XY02"]["action"] == "committed_confirmed_empty"
    assert by_image["XY02"]["status_after"] == "confirmed_empty"
    assert np.array_equal(
        tifffile.imread(package_dir / "reference_masks_to_fill" / "XY01_manual_reference_labels.tif"),
        xy01_labels,
    )
    status = _status_by_image(package_dir)
    assert status["XY01"]["labeler"] == "tester"
    assert status["XY02"]["status"] == "confirmed_empty"


def test_import_raw_annotation_labels_require_all_decisions_does_not_partially_commit(
    tmp_path: Path,
):
    package_dir, manifest_path = _write_package_and_raw_export(tmp_path)
    xy01_labels = np.zeros((8, 9), dtype=np.uint32)
    xy01_labels[2:5, 3:6] = 4
    tifffile.imwrite(
        package_dir / "raw_annotation_exports" / "XY01" / "XY01_editable_manual_reference_labels.tif",
        xy01_labels,
    )

    with pytest.raises(ValueError, match="unconfirmed empty labels"):
        import_raw_annotation_labels(
            package_dir=package_dir,
            raw_export_manifest_path=manifest_path,
            output_dir=tmp_path / "import",
            labeler="tester",
            completed_date="2026-06-29",
            require_all_decisions=True,
        )

    status = _status_by_image(package_dir)
    assert status["XY01"]["status"] == "not_started"
    assert status["XY02"]["status"] == "not_started"
    assert not (tmp_path / "import" / "bulk_import_summary.csv").exists()


def test_import_raw_annotation_labels_rejects_duplicate_raw_manifest_image_ids_before_commit(
    tmp_path: Path,
):
    package_dir, manifest_path = _write_package_and_raw_export(tmp_path)
    rows = _read_rows(manifest_path)
    rows.append(dict(rows[0]))
    _write_csv(manifest_path, rows)
    xy01_labels = np.zeros((8, 9), dtype=np.uint32)
    xy01_labels[2:5, 3:6] = 4
    tifffile.imwrite(
        package_dir / "raw_annotation_exports" / "XY01" / "XY01_editable_manual_reference_labels.tif",
        xy01_labels,
    )

    with pytest.raises(ValueError, match="duplicate image_id rows in raw annotation export manifest: XY01"):
        import_raw_annotation_labels(
            package_dir=package_dir,
            raw_export_manifest_path=manifest_path,
            output_dir=tmp_path / "import",
            labeler="tester",
            completed_date="2026-06-29",
        )

    status = _status_by_image(package_dir)
    assert status["XY01"]["status"] == "not_started"
    assert not (tmp_path / "import" / "bulk_import_summary.csv").exists()


def test_import_raw_annotation_labels_rejects_reference_mismatch_before_commit(
    tmp_path: Path,
):
    package_dir, manifest_path = _write_package_and_raw_export(tmp_path)
    rows = _read_rows(manifest_path)
    rows[0]["authoritative_reference_mask_path"] = str(
        (package_dir / "reference_masks_to_fill" / "XY02_manual_reference_labels.tif").resolve()
    )
    _write_csv(manifest_path, rows)
    xy01_labels = np.zeros((8, 9), dtype=np.uint32)
    xy01_labels[2:5, 3:6] = 4
    tifffile.imwrite(
        package_dir / "raw_annotation_exports" / "XY01" / "XY01_editable_manual_reference_labels.tif",
        xy01_labels,
    )

    with pytest.raises(ValueError, match="XY01 raw export authoritative reference path mismatch"):
        import_raw_annotation_labels(
            package_dir=package_dir,
            raw_export_manifest_path=manifest_path,
            output_dir=tmp_path / "import",
            labeler="tester",
            completed_date="2026-06-29",
        )

    status = _status_by_image(package_dir)
    assert status["XY01"]["status"] == "not_started"
    assert not (tmp_path / "import" / "bulk_import_summary.csv").exists()


def test_import_raw_annotation_labels_rejects_candidate_included_manifest_before_commit(
    tmp_path: Path,
):
    package_dir, manifest_path = _write_package_and_raw_export(tmp_path)
    rows = _read_rows(manifest_path)
    rows[0]["candidate_mask_path_included"] = "True"
    _write_csv(manifest_path, rows)
    xy01_labels = np.zeros((8, 9), dtype=np.uint32)
    xy01_labels[2:5, 3:6] = 4
    tifffile.imwrite(
        package_dir / "raw_annotation_exports" / "XY01" / "XY01_editable_manual_reference_labels.tif",
        xy01_labels,
    )

    with pytest.raises(ValueError, match="XY01 raw export manifest is not raw-only"):
        import_raw_annotation_labels(
            package_dir=package_dir,
            raw_export_manifest_path=manifest_path,
            output_dir=tmp_path / "import",
            labeler="tester",
            completed_date="2026-06-29",
        )

    status = _status_by_image(package_dir)
    assert status["XY01"]["status"] == "not_started"
    assert not (tmp_path / "import" / "bulk_import_summary.csv").exists()


def test_import_raw_annotation_labels_validates_all_committable_labels_before_writing(
    tmp_path: Path,
):
    package_dir, manifest_path = _write_package_and_raw_export(tmp_path)
    xy01_labels = np.zeros((8, 9), dtype=np.uint32)
    xy01_labels[2:5, 3:6] = 4
    tifffile.imwrite(
        package_dir / "raw_annotation_exports" / "XY01" / "XY01_editable_manual_reference_labels.tif",
        xy01_labels,
    )
    xy02_bad_shape = np.zeros((7, 9), dtype=np.uint32)
    xy02_bad_shape[2:5, 3:6] = 8
    tifffile.imwrite(
        package_dir / "raw_annotation_exports" / "XY02" / "XY02_editable_manual_reference_labels.tif",
        xy02_bad_shape,
    )

    with pytest.raises(ValueError, match="XY02 shape mismatch"):
        import_raw_annotation_labels(
            package_dir=package_dir,
            raw_export_manifest_path=manifest_path,
            output_dir=tmp_path / "import",
            labeler="tester",
            completed_date="2026-06-29",
        )

    status = _status_by_image(package_dir)
    assert status["XY01"]["status"] == "not_started"
    assert status["XY02"]["status"] == "not_started"
    assert np.count_nonzero(
        tifffile.imread(package_dir / "reference_masks_to_fill" / "XY01_manual_reference_labels.tif")
    ) == 0
    assert not (tmp_path / "import" / "bulk_import_summary.csv").exists()


def _write_package_and_raw_export(tmp_path: Path) -> tuple[Path, Path]:
    package_dir = tmp_path / "manual_validation" / "package"
    reference_dir = package_dir / "reference_masks_to_fill"
    annotation_dir = package_dir / "annotation_panels_raw_only"
    raw_export_dir = package_dir / "raw_annotation_exports"
    for path in [reference_dir, annotation_dir, raw_export_dir]:
        path.mkdir(parents=True)

    manifest_rows = []
    status_rows = []
    raw_rows = []
    for image_id in ["XY01", "XY02"]:
        reference_path = reference_dir / f"{image_id}_manual_reference_labels.tif"
        tifffile.imwrite(reference_path, np.zeros((8, 9), dtype=np.uint32))
        annotation_path = annotation_dir / f"{image_id}_manual_annotation_panel.png"
        annotation_path.write_bytes(b"annotation")
        image_export_dir = raw_export_dir / image_id
        image_export_dir.mkdir()
        editable_path = image_export_dir / f"{image_id}_editable_manual_reference_labels.tif"
        tifffile.imwrite(editable_path, np.zeros((8, 9), dtype=np.uint32))
        ch2_path = image_export_dir / f"{image_id}_CH2_raw.tif"
        ch4_path = image_export_dir / f"{image_id}_CH4_raw.tif"
        tifffile.imwrite(ch2_path, np.zeros((8, 9), dtype=np.uint16))
        tifffile.imwrite(ch4_path, np.zeros((8, 9), dtype=np.uint16))
        manifest_rows.append(
            {
                "image_id": image_id,
                "source_id": image_id,
                "validation_task": "asma_associated_region",
                "ch2_path": str(ch2_path),
                "ch4_path": str(ch4_path),
                "candidate_mask_path": "",
                "nuclei_mask_path": "",
                "manual_reference_mask_path": str(reference_path),
                "annotation_panel_path": str(annotation_path),
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
                "status": "not_started",
                "labeler": "",
                "completed_date": "",
                "notes": "",
            }
        )
        raw_rows.append(
            {
                "image_id": image_id,
                "status": "not_started",
                "ch2_export_path": str(ch2_path),
                "ch4_export_path": str(ch4_path),
                "editable_reference_labels_path": str(editable_path),
                "authoritative_reference_mask_path": str(reference_path.resolve()),
                "commit_command": "",
                "confirm_empty_command": "",
                "candidate_mask_path_included": "False",
            }
        )
    _write_csv(package_dir / "manual_validation_manifest.csv", manifest_rows)
    _write_csv(package_dir / "manual_labeling_status.csv", status_rows)
    raw_manifest_path = raw_export_dir / "raw_annotation_export_manifest.csv"
    _write_csv(raw_manifest_path, raw_rows)
    return package_dir, raw_manifest_path


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _status_by_image(package_dir: Path) -> dict[str, dict[str, str]]:
    return {
        row["image_id"]: row
        for row in _read_rows(package_dir / "manual_labeling_status.csv")
    }
