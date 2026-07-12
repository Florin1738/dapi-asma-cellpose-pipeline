from __future__ import annotations

import csv
import shlex
from pathlib import Path

import numpy as np
import tifffile

from dapi_norm.raw_annotation_export import prepare_raw_annotation_export


def test_prepare_raw_annotation_export_writes_raw_only_editor_files(tmp_path: Path):
    package_dir = _write_package(tmp_path)

    outputs = prepare_raw_annotation_export(
        package_dir=package_dir,
        output_dir=package_dir / "raw_annotation_exports",
    )

    assert outputs["manifest"].exists()
    assert outputs["readme"].exists()
    rows = _read_rows(outputs["manifest"])
    assert len(rows) == 1
    row = rows[0]
    assert row["image_id"] == "XY01"
    assert row["status"] == "not_started"
    assert Path(row["ch2_export_path"]).exists()
    assert Path(row["ch4_export_path"]).exists()
    assert Path(row["editable_reference_labels_path"]).exists()
    assert row["candidate_mask_path_included"] == "False"
    assert "candidate" not in row["ch2_export_path"].lower()
    assert "candidate" not in row["ch4_export_path"].lower()

    assert np.array_equal(tifffile.imread(row["ch2_export_path"]), np.arange(72, dtype=np.uint16).reshape(8, 9))
    assert np.array_equal(tifffile.imread(row["ch4_export_path"]), np.flipud(np.arange(72, dtype=np.uint16).reshape(8, 9)))
    editable = tifffile.imread(row["editable_reference_labels_path"])
    assert editable.dtype == np.uint32
    assert np.count_nonzero(editable) == 0

    readme = outputs["readme"].read_text(encoding="utf-8")
    assert "candidate masks are intentionally not exported" in readme
    assert "scripts/commit_manual_reference_mask.py" in readme
    assert "--image-id XY01" in readme
    assert "Do not manually overwrite reference_masks_to_fill/" in readme


def test_prepare_raw_annotation_export_writes_shell_safe_command_templates(
    tmp_path: Path,
):
    package_dir = _write_package(tmp_path, package_name="package with spaces")
    output_dir = package_dir / "raw annotation exports"

    outputs = prepare_raw_annotation_export(
        package_dir=package_dir,
        output_dir=output_dir,
    )

    row = _read_rows(outputs["manifest"])[0]
    for column in ["commit_command", "confirm_empty_command"]:
        command = row[column]
        assert "<" not in command
        assert ">" not in command
        args = shlex.split(command)
        assert args[:2] == [".venv/bin/python", "scripts/commit_manual_reference_mask.py"]
        assert args[args.index("--package") + 1] == str(package_dir)
        assert args[args.index("--image-id") + 1] == "XY01"
        assert args[args.index("--labels") + 1] == row["editable_reference_labels_path"]
        assert args[args.index("--labeler") + 1] == "YOUR_INITIALS"


def test_prepare_raw_annotation_export_refuses_stale_candidate_files(
    tmp_path: Path,
):
    package_dir = _write_package(tmp_path)
    output_dir = package_dir / "raw_annotation_exports"
    stale_dir = output_dir / "XY01"
    stale_dir.mkdir(parents=True)
    (stale_dir / "XY01_candidate_labels.tif").write_bytes(b"stale candidate")

    try:
        prepare_raw_annotation_export(
            package_dir=package_dir,
            output_dir=output_dir,
        )
    except ValueError as exc:
        assert "candidate-looking files" in str(exc)
    else:
        raise AssertionError("expected stale candidate file failure")


def test_prepare_raw_annotation_export_does_not_partially_write_after_validation_failure(
    tmp_path: Path,
):
    package_dir = _write_package(tmp_path, image_ids=("XY01", "XY02"))
    rows = _read_rows(package_dir / "manual_labeling_status.csv")
    rows[1]["manual_reference_mask_path"] = str(package_dir / "reference_masks_to_fill" / "XY02_stale.tif")
    _write_csv(package_dir / "manual_labeling_status.csv", rows)
    output_dir = package_dir / "raw_annotation_exports"

    try:
        prepare_raw_annotation_export(
            package_dir=package_dir,
            output_dir=output_dir,
        )
    except ValueError as exc:
        assert "status/manifest path mismatch" in str(exc)
    else:
        raise AssertionError("expected validation failure")

    assert not output_dir.exists()


def test_prepare_raw_annotation_export_refuses_missing_manifest_status_pair(
    tmp_path: Path,
):
    package_dir = _write_package(tmp_path)
    (package_dir / "manual_labeling_status.csv").write_text(
        "image_id,manual_reference_mask_path,annotation_panel_path,status,labeler,completed_date,notes\n",
        encoding="utf-8",
    )

    try:
        prepare_raw_annotation_export(
            package_dir=package_dir,
            output_dir=package_dir / "raw_annotation_exports",
        )
    except ValueError as exc:
        assert "missing status rows" in str(exc)
    else:
        raise AssertionError("expected missing status row failure")


def test_prepare_raw_annotation_export_rejects_stale_status_reference_path(
    tmp_path: Path,
):
    package_dir = _write_package(tmp_path)
    stale_path = package_dir / "reference_masks_to_fill" / "XY01_stale_labels.tif"
    tifffile.imwrite(stale_path, np.zeros((8, 9), dtype=np.uint32))
    rows = _read_rows(package_dir / "manual_labeling_status.csv")
    rows[0]["manual_reference_mask_path"] = str(stale_path)
    _write_csv(package_dir / "manual_labeling_status.csv", rows)

    try:
        prepare_raw_annotation_export(
            package_dir=package_dir,
            output_dir=package_dir / "raw_annotation_exports",
        )
    except ValueError as exc:
        assert "status/manifest path mismatch" in str(exc)
    else:
        raise AssertionError("expected stale status reference path failure")


def _write_package(
    tmp_path: Path,
    *,
    package_name: str = "package",
    image_ids: tuple[str, ...] = ("XY01",),
) -> Path:
    package_dir = tmp_path / "manual_validation" / package_name
    reference_dir = package_dir / "reference_masks_to_fill"
    annotation_dir = package_dir / "annotation_panels_raw_only"
    guide_dir = package_dir / "guide_panels"
    data_dir = tmp_path / "data"
    masks_dir = tmp_path / "masks"
    for path in [reference_dir, annotation_dir, guide_dir, data_dir, masks_dir]:
        path.mkdir(parents=True)
    manifest_rows = []
    status_rows = []
    for index, image_id in enumerate(image_ids):
        ch2 = (np.arange(72, dtype=np.uint16).reshape(8, 9) + index).astype(np.uint16)
        ch4 = np.flipud(ch2)
        manual = np.zeros((8, 9), dtype=np.uint32)
        candidate = np.ones((8, 9), dtype=np.uint32)
        ch2_path = data_dir / f"{image_id}_CH2.tif"
        ch4_path = data_dir / f"{image_id}_CH4.tif"
        manual_path = reference_dir / f"{image_id}_manual_reference_labels.tif"
        candidate_path = masks_dir / f"{image_id}_candidate_labels.tif"
        annotation_path = annotation_dir / f"{image_id}_manual_annotation_panel.png"
        guide_path = guide_dir / f"{image_id}_manual_validation_guide.png"
        for path, arr in [
            (ch2_path, ch2),
            (ch4_path, ch4),
            (manual_path, manual),
            (candidate_path, candidate),
        ]:
            tifffile.imwrite(path, arr)
        annotation_path.write_bytes(b"annotation")
        guide_path.write_bytes(b"guide")
        manifest_rows.append(
            {
                "image_id": image_id,
                "source_id": image_id,
                "validation_task": "asma_associated_region",
                "ch2_path": str(ch2_path),
                "ch4_path": str(ch4_path),
                "candidate_mask_path": str(candidate_path),
                "nuclei_mask_path": "",
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
        )
        status_rows.append(
            {
                "image_id": image_id,
                "manual_reference_mask_path": str(manual_path.resolve()),
                "annotation_panel_path": str(annotation_path.resolve()),
                "status": "not_started",
                "labeler": "",
                "completed_date": "",
                "notes": "",
            }
        )
    _write_csv(package_dir / "manual_validation_manifest.csv", manifest_rows)
    _write_csv(package_dir / "manual_labeling_status.csv", status_rows)
    return package_dir


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))
