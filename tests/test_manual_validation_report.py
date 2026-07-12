from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
import tifffile

from dapi_norm.manual_validation_report import (
    _object_match_status_masks,
    _status_overlay_masks,
    parse_candidate_specs,
    run_manual_validation_report,
)


def test_run_manual_validation_report_ranks_candidates_and_writes_overlays(
    tmp_path: Path,
):
    reference_dir, status_path, manifest_path = _write_reference_package(tmp_path)
    candidate_dirs = _write_candidate_methods(tmp_path)

    outputs = run_manual_validation_report(
        candidate_dirs=candidate_dirs,
        reference_dir=reference_dir,
        completion_status_path=status_path,
        manifest_path=manifest_path,
        output_dir=tmp_path / "validation_report",
        iou_threshold=0.5,
        min_precision=0.9,
        min_recall=0.9,
        min_f1=0.9,
        min_mean_iou=0.9,
    )

    assert outputs["method_summary"].exists()
    assert outputs["image_summary"].exists()
    assert outputs["report"].exists()
    rows = _read_rows(outputs["method_summary"])
    by_method = {row["candidate_method"]: row for row in rows}
    assert by_method["perfect"]["passes_acceptance_criteria"] == "True"
    assert by_method["perfect"]["micro_precision"] == "1.0"
    assert by_method["perfect"]["micro_recall"] == "1.0"
    assert by_method["missed"]["passes_acceptance_criteria"] == "False"
    assert by_method["missed"]["micro_recall"] == "0.5"
    assert (
        tmp_path
        / "validation_report"
        / "overlays"
        / "perfect"
        / "XY01_candidate_vs_reference_overlay.png"
    ).exists()
    assert (
        tmp_path
        / "validation_report"
        / "overlays"
        / "missed"
        / "XY02_candidate_vs_reference_overlay.png"
    ).exists()
    assert (tmp_path / "validation_report" / "overlays" / "missed_contact_sheet.png").exists()
    image_rows = _read_rows(outputs["image_summary"])
    assert image_rows
    assert all(row["overlay_path"] for row in image_rows)
    assert all(Path(row["overlay_path"]).exists() for row in image_rows)
    assert (
        tmp_path
        / "validation_report"
        / "per_candidate"
        / "perfect"
        / "manual_mask_validation_summary.csv"
    ).exists()
    report = outputs["report"].read_text(encoding="utf-8")
    assert "Validation run completed against manual/reference masks" in report
    assert "Status: Validated against manual/reference masks" not in report
    assert "This validates only the supplied manual-reference task" in report


def test_object_match_status_masks_split_true_positive_false_positive_and_false_negative():
    candidate = np.zeros((8, 8), dtype=np.uint32)
    reference = np.zeros((8, 8), dtype=np.uint32)
    candidate[1:3, 1:3] = 2
    reference[1:3, 1:3] = 8
    candidate[4:6, 1:3] = 3
    reference[4:6, 5:7] = 9
    matches = [
        {
            "candidate_label": "2",
            "reference_label": "8",
            "match_status": "true_positive",
            "matched": "True",
        },
        {
            "candidate_label": "3",
            "reference_label": "",
            "match_status": "false_positive",
            "matched": "False",
        },
        {
            "candidate_label": "",
            "reference_label": "9",
            "match_status": "false_negative",
            "matched": "False",
        },
    ]

    masks = _object_match_status_masks(
        candidate_labels=candidate,
        reference_labels=reference,
        match_rows=matches,
    )

    assert np.array_equal(masks["true_positive"], (candidate == 2) | (reference == 8))
    assert np.array_equal(masks["false_positive"], candidate == 3)
    assert np.array_equal(masks["false_negative"], reference == 9)
    assert not np.any(
        masks["true_positive"] & masks["false_positive"]
    )
    assert not np.any(
        masks["true_positive"] & masks["false_negative"]
    )


def test_status_overlay_masks_do_not_draw_fp_or_fn_boundaries_on_matched_tp():
    candidate = np.zeros((8, 8), dtype=np.uint32)
    reference = np.zeros((8, 8), dtype=np.uint32)
    candidate[2:6, 2:6] = 2
    reference[2:6, 3:7] = 8
    matches = [
        {
            "candidate_label": "2",
            "reference_label": "8",
            "match_status": "true_positive",
            "matched": "True",
        }
    ]

    masks = _object_match_status_masks(
        candidate_labels=candidate,
        reference_labels=reference,
        match_rows=matches,
    )
    overlay_masks = _status_overlay_masks(masks)

    assert np.any(overlay_masks["true_positive_boundary"])
    assert not np.any(overlay_masks["false_positive_fill"])
    assert not np.any(overlay_masks["false_negative_fill"])
    assert not np.any(overlay_masks["false_positive_boundary"])
    assert not np.any(overlay_masks["false_negative_boundary"])


def test_status_overlay_boundaries_stay_inside_their_own_status_regions():
    true_positive = np.zeros((8, 8), dtype=bool)
    false_positive = np.zeros((8, 8), dtype=bool)
    false_negative = np.zeros((8, 8), dtype=bool)
    true_positive[2:6, 2:4] = True
    false_positive[2:6, 4:6] = True

    overlay_masks = _status_overlay_masks(
        {
            "true_positive": true_positive,
            "false_positive": false_positive,
            "false_negative": false_negative,
        }
    )

    assert np.any(overlay_masks["true_positive_boundary"])
    assert np.any(overlay_masks["false_positive_boundary"])
    assert not np.any(overlay_masks["true_positive_boundary"] & false_positive)
    assert not np.any(overlay_masks["false_positive_boundary"] & true_positive)


def test_status_overlay_masks_split_fp_fn_overlap_as_below_threshold_mismatch():
    true_positive = np.zeros((8, 8), dtype=bool)
    false_positive = np.zeros((8, 8), dtype=bool)
    false_negative = np.zeros((8, 8), dtype=bool)
    false_positive[2:6, 2:6] = True
    false_negative[4:7, 4:7] = True

    overlay_masks = _status_overlay_masks(
        {
            "true_positive": true_positive,
            "false_positive": false_positive,
            "false_negative": false_negative,
        }
    )

    expected_overlap = false_positive & false_negative
    assert np.array_equal(overlay_masks["below_threshold_overlap_fill"], expected_overlap)
    assert not np.any(overlay_masks["false_positive_fill"] & expected_overlap)
    assert not np.any(overlay_masks["false_negative_fill"] & expected_overlap)
    assert not np.any(
        overlay_masks["false_positive_fill"] & overlay_masks["false_negative_fill"]
    )
    assert not np.any(
        overlay_masks["false_positive_fill"] & overlay_masks["below_threshold_overlap_fill"]
    )
    assert not np.any(
        overlay_masks["false_negative_fill"] & overlay_masks["below_threshold_overlap_fill"]
    )
    assert np.any(overlay_masks["below_threshold_overlap_boundary"])


def test_run_manual_validation_report_refuses_not_started_status(tmp_path: Path):
    reference_dir, status_path, manifest_path = _write_reference_package(
        tmp_path,
        status="not_started",
    )
    candidate_dirs = _write_candidate_methods(tmp_path)
    output_dir = tmp_path / "validation_report"

    with pytest.raises(ValueError, match="not marked complete"):
        run_manual_validation_report(
            candidate_dirs=candidate_dirs,
            reference_dir=reference_dir,
            completion_status_path=status_path,
            manifest_path=manifest_path,
            output_dir=output_dir,
            iou_threshold=0.5,
        )
    assert not output_dir.exists()


def test_parse_candidate_specs_requires_named_unique_directories(tmp_path: Path):
    first = tmp_path / "first"
    second = tmp_path / "second"

    parsed = parse_candidate_specs([f"first={first}", f"second={second}"])

    assert parsed == {"first": first, "second": second}
    with pytest.raises(ValueError, match="name=directory"):
        parse_candidate_specs([str(first)])
    with pytest.raises(ValueError, match="Duplicate candidate method"):
        parse_candidate_specs([f"first={first}", f"first={second}"])


def test_run_manual_validation_report_refuses_extra_not_started_status_row(
    tmp_path: Path,
):
    reference_dir, status_path, manifest_path = _write_reference_package(tmp_path)
    with status_path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "image_id",
                "manual_reference_mask_path",
                "annotation_panel_path",
                "status",
                "labeler",
                "completed_date",
                "notes",
            ],
        )
        writer.writerow(
            {
                "image_id": "XY99",
                "manual_reference_mask_path": str((reference_dir / "XY99_manual_reference_labels.tif").resolve()),
                "annotation_panel_path": "",
                "status": "not_started",
                "labeler": "",
                "completed_date": "",
                "notes": "extra incomplete row",
            }
        )
    candidate_dirs = _write_candidate_methods(tmp_path)

    with pytest.raises(ValueError, match="not marked complete"):
        run_manual_validation_report(
            candidate_dirs=candidate_dirs,
            reference_dir=reference_dir,
            completion_status_path=status_path,
            manifest_path=manifest_path,
            output_dir=tmp_path / "validation_report",
            iou_threshold=0.5,
        )


def test_run_manual_validation_report_requires_manifest_ch2_for_each_image(
    tmp_path: Path,
):
    reference_dir, status_path, manifest_path = _write_reference_package(tmp_path)
    rows = _read_rows(manifest_path)
    rows = [row for row in rows if row["image_id"] != "XY02"]
    _write_csv(manifest_path, rows)
    candidate_dirs = _write_candidate_methods(tmp_path)

    with pytest.raises(ValueError, match="manifest is missing CH2 path"):
        run_manual_validation_report(
            candidate_dirs=candidate_dirs,
            reference_dir=reference_dir,
            completion_status_path=status_path,
            manifest_path=manifest_path,
            output_dir=tmp_path / "validation_report",
            iou_threshold=0.5,
        )


def test_run_manual_validation_report_renders_overlay_for_absent_candidate_tiff(
    tmp_path: Path,
):
    reference_dir, status_path, manifest_path = _write_reference_package(tmp_path)
    candidate_dirs = _write_candidate_methods(tmp_path)
    (candidate_dirs["missed"] / "XY02_candidate_labels.tif").unlink()

    outputs = run_manual_validation_report(
        candidate_dirs={"missed": candidate_dirs["missed"]},
        reference_dir=reference_dir,
        completion_status_path=status_path,
        manifest_path=manifest_path,
        output_dir=tmp_path / "validation_report",
        iou_threshold=0.5,
    )

    image_rows = _read_rows(outputs["image_summary"])
    by_image = {row["image_id"]: row for row in image_rows}
    assert by_image["XY02"]["false_negatives"] == "1"
    assert Path(by_image["XY02"]["overlay_path"]).exists()
    assert (
        tmp_path
        / "validation_report"
        / "debug"
        / "missing_candidate_masks"
        / "missed"
        / "XY02_missing_candidate_as_empty_labels.tif"
    ).exists()


def _write_reference_package(
    tmp_path: Path,
    *,
    status: str = "complete_non_empty",
) -> tuple[Path, Path, Path]:
    package = tmp_path / "manual_validation" / "package"
    reference_dir = package / "reference_masks_to_fill"
    annotation_dir = package / "annotation_panels_raw_only"
    image_dir = tmp_path / "images"
    reference_dir.mkdir(parents=True)
    annotation_dir.mkdir(parents=True)
    image_dir.mkdir(parents=True)

    ch2_paths: dict[str, Path] = {}
    rows = []
    status_rows = []
    for image_id, object_slice in {
        "XY01": (slice(1, 3), slice(1, 3)),
        "XY02": (slice(5, 7), slice(5, 7)),
    }.items():
        ch2 = np.zeros((10, 10), dtype=np.uint16)
        ch2[object_slice] = 1000
        ch2_path = image_dir / f"{image_id}_CH2.tif"
        tifffile.imwrite(ch2_path, ch2)
        ch2_paths[image_id] = ch2_path
        reference = np.zeros((10, 10), dtype=np.uint32)
        reference[object_slice] = 1
        reference_path = reference_dir / f"{image_id}_manual_reference_labels.tif"
        tifffile.imwrite(reference_path, reference)
        annotation_path = annotation_dir / f"{image_id}_manual_annotation_panel.png"
        annotation_path.write_bytes(b"placeholder")
        rows.append(
            {
                "image_id": image_id,
                "source_id": image_id,
                "validation_task": "asma_associated_region",
                "ch2_path": str(ch2_path),
                "ch4_path": str(ch2_path),
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
                "status": status,
                "labeler": "tester",
                "completed_date": "2026-06-29",
                "notes": "",
            }
        )

    manifest_path = package / "manual_validation_manifest.csv"
    _write_csv(manifest_path, rows)
    status_path = package / "manual_labeling_status.csv"
    _write_csv(status_path, status_rows)
    return reference_dir, status_path, manifest_path


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
