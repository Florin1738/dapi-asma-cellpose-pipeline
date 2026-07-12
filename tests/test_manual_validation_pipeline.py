from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import tifffile

from dapi_norm.manual_validation_pipeline import run_manual_validation_pipeline


def test_run_manual_validation_pipeline_stops_when_annotation_audit_is_not_ready(
    tmp_path: Path,
):
    package_dir = _write_package(tmp_path, statuses={"XY01": "complete_non_empty", "XY02": "not_started"})
    candidate_dirs = _write_candidate_methods(tmp_path)

    outputs = run_manual_validation_pipeline(
        package_dir=package_dir,
        candidate_dirs={"perfect": candidate_dirs["perfect"]},
        output_dir=tmp_path / "pipeline",
        iou_threshold=0.5,
    )

    assert outputs["gate_report"].exists()
    assert outputs["audit_csv"].exists()
    assert outputs["validation_ready"] is False
    assert "Validation report was not run" in outputs["gate_report"].read_text(encoding="utf-8")
    assert "XY02" in outputs["gate_report"].read_text(encoding="utf-8")
    assert "status_not_complete" in outputs["gate_report"].read_text(encoding="utf-8")
    assert not (tmp_path / "pipeline" / "validation_report").exists()


def test_run_manual_validation_pipeline_flags_stale_report_when_blocked(
    tmp_path: Path,
):
    package_dir = _write_package(tmp_path, statuses={"XY01": "complete_non_empty", "XY02": "not_started"})
    candidate_dirs = _write_candidate_methods(tmp_path)
    stale_report_dir = tmp_path / "pipeline" / "validation_report"
    stale_report_dir.mkdir(parents=True)
    (stale_report_dir / "method_validation_summary.csv").write_text("stale\n", encoding="utf-8")

    outputs = run_manual_validation_pipeline(
        package_dir=package_dir,
        candidate_dirs={"perfect": candidate_dirs["perfect"]},
        output_dir=tmp_path / "pipeline",
        iou_threshold=0.5,
    )

    gate_report = outputs["gate_report"].read_text(encoding="utf-8")
    assert outputs["validation_ready"] is False
    assert "Stale validation report directory exists" in gate_report
    assert "do not use those metrics" in gate_report


def test_run_manual_validation_pipeline_runs_report_after_annotation_audit_passes(
    tmp_path: Path,
):
    package_dir = _write_package(tmp_path, statuses={"XY01": "complete_non_empty", "XY02": "complete_non_empty"})
    candidate_dirs = _write_candidate_methods(tmp_path)

    outputs = run_manual_validation_pipeline(
        package_dir=package_dir,
        candidate_dirs=candidate_dirs,
        output_dir=tmp_path / "pipeline",
        iou_threshold=0.5,
        min_precision=0.9,
        min_recall=0.9,
        min_f1=0.9,
        min_mean_iou=0.9,
    )

    assert outputs["validation_ready"] is True
    assert outputs["gate_report"].exists()
    assert outputs["method_summary"].exists()
    assert outputs["image_summary"].exists()
    assert outputs["validation_report"].exists()
    report = outputs["gate_report"].read_text(encoding="utf-8")
    assert "Validation report was run" in report
    rows = _read_rows(outputs["method_summary"])
    by_method = {row["candidate_method"]: row for row in rows}
    assert by_method["perfect"]["passes_acceptance_criteria"] == "True"
    assert by_method["missed"]["passes_acceptance_criteria"] == "False"
    assert (tmp_path / "pipeline" / "validation_report" / "overlays" / "perfect_contact_sheet.png").exists()


def _write_package(tmp_path: Path, *, statuses: dict[str, str]) -> Path:
    package_dir = tmp_path / "manual_validation" / "package"
    reference_dir = package_dir / "reference_masks_to_fill"
    annotation_dir = package_dir / "annotation_panels_raw_only"
    image_dir = tmp_path / "images"
    candidate_shape_dir = package_dir / "candidate_shape_masks"
    for path in [reference_dir, annotation_dir, image_dir, candidate_shape_dir]:
        path.mkdir(parents=True)

    manifest_rows = []
    status_rows = []
    object_slices = {
        "XY01": (slice(1, 3), slice(1, 3)),
        "XY02": (slice(5, 7), slice(5, 7)),
    }
    for image_id, object_slice in object_slices.items():
        ch2 = np.zeros((10, 10), dtype=np.uint16)
        ch2[object_slice] = 1000
        ch2_path = image_dir / f"{image_id}_CH2.tif"
        tifffile.imwrite(ch2_path, ch2)

        reference = np.zeros((10, 10), dtype=np.uint32)
        if statuses[image_id] == "complete_non_empty":
            reference[object_slice] = 1
        reference_path = reference_dir / f"{image_id}_manual_reference_labels.tif"
        tifffile.imwrite(reference_path, reference)

        shape_mask_path = candidate_shape_dir / f"{image_id}_shape_labels.tif"
        tifffile.imwrite(shape_mask_path, np.zeros((10, 10), dtype=np.uint32))
        annotation_path = annotation_dir / f"{image_id}_manual_annotation_panel.png"
        annotation_path.write_bytes(b"placeholder")
        manifest_rows.append(
            {
                "image_id": image_id,
                "source_id": image_id,
                "validation_task": "asma_associated_region",
                "ch2_path": str(ch2_path),
                "ch4_path": str(ch2_path),
                "candidate_mask_path": str(shape_mask_path),
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
                "status": statuses[image_id],
                "labeler": "tester" if statuses[image_id] != "not_started" else "",
                "completed_date": "2026-06-29" if statuses[image_id] != "not_started" else "",
                "notes": "",
            }
        )
    _write_csv(package_dir / "manual_validation_manifest.csv", manifest_rows)
    _write_csv(package_dir / "manual_labeling_status.csv", status_rows)
    return package_dir


def _write_candidate_methods(tmp_path: Path) -> dict[str, Path]:
    candidate_root = tmp_path / "candidates"
    perfect_dir = candidate_root / "perfect"
    missed_dir = candidate_root / "missed"
    perfect_dir.mkdir(parents=True)
    missed_dir.mkdir(parents=True)
    for image_id, object_slice in {
        "XY01": (slice(1, 3), slice(1, 3)),
        "XY02": (slice(5, 7), slice(5, 7)),
    }.items():
        perfect = np.zeros((10, 10), dtype=np.uint32)
        perfect[object_slice] = 7
        tifffile.imwrite(perfect_dir / f"{image_id}_candidate_labels.tif", perfect)

    missed_xy01 = np.zeros((10, 10), dtype=np.uint32)
    missed_xy01[1:3, 1:3] = 1
    missed_xy02 = np.zeros((10, 10), dtype=np.uint32)
    tifffile.imwrite(missed_dir / "XY01_candidate_labels.tif", missed_xy01)
    tifffile.imwrite(missed_dir / "XY02_candidate_labels.tif", missed_xy02)
    return {"perfect": perfect_dir, "missed": missed_dir}


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))
