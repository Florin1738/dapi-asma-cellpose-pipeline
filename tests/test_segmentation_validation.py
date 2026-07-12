from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
import tifffile
import yaml

from dapi_norm.pi_simple_summary import ImagePair
from dapi_norm.seeded_regions import run_seeded_region_batch
from dapi_norm.segmentation_validation import (
    evaluate_instance_mask_iou,
    run_manual_mask_validation,
    validate_seeded_region_outputs,
    write_instance_mask_validation_tables,
)


def test_validate_seeded_region_outputs_checks_formulas_artifacts_and_qc_status(tmp_path: Path):
    output_dir = _write_seeded_run(tmp_path)

    result = validate_seeded_region_outputs(output_dir)

    assert result["summary_rows"] == 1
    assert result["formulas_match"] is True
    assert result["masks_exist"] is True
    assert result["qc_panels_exist"] is True
    assert result["manual_ground_truth_available"] is False


def test_validate_seeded_region_outputs_rejects_corrupted_background_corrected_sum(
    tmp_path: Path,
):
    output_dir = _write_seeded_run(tmp_path)
    _rewrite_summary_values(
        output_dir / "summaries" / "seeded_region_image_metrics.csv",
        {
            "seeded_region_integrated_background_corrected": "999",
            "seeded_region_intensity_per_DAPI_positive_nucleus": "999",
        },
    )

    with pytest.raises(ValueError, match="seeded_region_integrated_background_corrected"):
        validate_seeded_region_outputs(output_dir)


def test_validate_seeded_region_outputs_rejects_missing_summary_rows(tmp_path: Path):
    output_dir = _write_seeded_run(tmp_path)
    config_path = output_dir / "logs" / "config_resolved.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["images_processed"] = 2
    config["image_records"].append({**config["image_records"][0], "source_id": "XY02", "location": "XY02"})
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    with pytest.raises(ValueError, match="missing rows"):
        validate_seeded_region_outputs(output_dir)


def test_validate_seeded_region_outputs_rejects_missing_scientific_warning(tmp_path: Path):
    output_dir = _write_seeded_run(tmp_path)
    _rewrite_summary_values(
        output_dir / "summaries" / "seeded_region_image_metrics.csv",
        {"warnings": "exploratory_asma_associated_region_not_validated_cell_mask"},
    )

    with pytest.raises(ValueError, match="required warning"):
        validate_seeded_region_outputs(output_dir)


def test_validate_seeded_region_outputs_rejects_seeded_labels_not_present_in_dapi_mask(
    tmp_path: Path,
):
    output_dir = _write_seeded_run(tmp_path)
    rows = _read_csv_rows(output_dir / "summaries" / "seeded_region_image_metrics.csv")
    mask_path = Path(rows[0]["mask_path"])
    mask = tifffile.imread(mask_path)
    mask[8:10, 8:10] = 99
    tifffile.imwrite(mask_path, mask)

    with pytest.raises(ValueError, match="not present in DAPI"):
        validate_seeded_region_outputs(output_dir)


def test_evaluate_instance_mask_iou_reports_manual_mask_validation_metrics():
    candidate = np.zeros((8, 8), dtype=np.uint32)
    candidate[1:3, 1:3] = 1
    candidate[5:7, 5:7] = 2
    reference = np.zeros((8, 8), dtype=np.uint32)
    reference[1:3, 1:3] = 10
    reference[1:3, 5:7] = 20

    summary, matches = evaluate_instance_mask_iou(
        image_id="synthetic",
        candidate_mask=candidate,
        reference_mask=reference,
        iou_threshold=0.5,
    )

    assert summary["image_id"] == "synthetic"
    assert summary["candidate_count"] == 2
    assert summary["reference_count"] == 2
    assert summary["n_predicted"] == 2
    assert summary["n_manual"] == 2
    assert summary["true_positives"] == 1
    assert summary["matched_count"] == 1
    assert summary["false_positives"] == 1
    assert summary["false_positive_count"] == 1
    assert summary["false_negatives"] == 1
    assert summary["false_negative_count"] == 1
    assert summary["precision"] == 0.5
    assert summary["recall"] == 0.5
    assert summary["f1"] == 0.5
    assert summary["mean_iou_matched"] == 1.0
    assert summary["count_error"] == 0
    assert summary["count_error_percent"] == 0.0
    assert matches == [
        {
            "image_id": "synthetic",
            "candidate_label": 1,
            "reference_label": 10,
            "iou": 1.0,
            "match_status": "true_positive",
            "matched": True,
        },
        {
            "image_id": "synthetic",
            "candidate_label": 2,
            "reference_label": "",
            "iou": 0.0,
            "match_status": "false_positive",
            "matched": False,
        },
        {
            "image_id": "synthetic",
            "candidate_label": "",
            "reference_label": 20,
            "iou": 0.0,
            "match_status": "false_negative",
            "matched": False,
        },
    ]


def test_evaluate_instance_mask_iou_rejects_zero_threshold():
    mask = np.zeros((4, 4), dtype=np.uint32)

    with pytest.raises(ValueError, match="greater than 0"):
        evaluate_instance_mask_iou(
            image_id="synthetic",
            candidate_mask=mask,
            reference_mask=mask,
            iou_threshold=0.0,
        )


@pytest.mark.parametrize(
    ("bad_mask", "expected_message"),
    [
        (np.zeros((4, 4), dtype=np.float32), "integer"),
        (np.zeros((4, 4, 3), dtype=np.uint8), "2-D"),
        (-np.ones((4, 4), dtype=np.int16), "negative"),
    ],
)
def test_evaluate_instance_mask_iou_rejects_invalid_label_masks(
    bad_mask: np.ndarray,
    expected_message: str,
):
    valid = np.zeros((4, 4), dtype=np.uint32)

    with pytest.raises(ValueError, match=expected_message):
        evaluate_instance_mask_iou(
            image_id="synthetic",
            candidate_mask=bad_mask,
            reference_mask=valid,
            iou_threshold=0.5,
        )


def test_evaluate_instance_mask_iou_rejects_disconnected_binary_instance_mask():
    binary = np.zeros((8, 8), dtype=np.uint8)
    binary[1:3, 1:3] = 1
    binary[5:7, 5:7] = 1
    valid = np.zeros((8, 8), dtype=np.uint32)

    with pytest.raises(ValueError, match="binary"):
        evaluate_instance_mask_iou(
            image_id="synthetic",
            candidate_mask=binary,
            reference_mask=valid,
            iou_threshold=0.5,
        )


def test_write_instance_mask_validation_tables_exports_summary_and_matches(tmp_path: Path):
    summary = {
        "image_id": "synthetic",
        "roi_id": "full_image",
        "iou_threshold": 0.5,
        "n_manual": 1,
        "n_predicted": 1,
        "candidate_count": 1,
        "reference_count": 1,
        "true_positives": 1,
        "matched_count": 1,
        "false_positives": 0,
        "false_positive_count": 0,
        "false_negatives": 0,
        "false_negative_count": 0,
        "precision": 1.0,
        "recall": 1.0,
        "f1": 1.0,
        "mean_iou_matched": 0.75,
        "count_error": 0,
        "count_error_percent": 0.0,
    }
    matches = [
        {
            "image_id": "synthetic",
            "candidate_label": 1,
            "reference_label": 2,
            "iou": 0.75,
            "match_status": "true_positive",
            "matched": True,
        }
    ]

    write_instance_mask_validation_tables(
        output_dir=tmp_path / "validation",
        summaries=[summary],
        matches=matches,
    )

    assert (tmp_path / "validation" / "manual_mask_validation_summary.csv").exists()
    assert (tmp_path / "validation" / "manual_mask_validation_matches.csv").exists()


def test_run_manual_mask_validation_reads_masks_and_writes_iou_outputs(tmp_path: Path):
    candidate = np.zeros((8, 8), dtype=np.uint32)
    candidate[1:3, 1:3] = 1
    reference = np.zeros((8, 8), dtype=np.uint32)
    reference[1:3, 1:3] = 2
    candidate_path = tmp_path / "candidate" / "XY01_candidate_labels.tif"
    reference_path = tmp_path / "manual" / "XY01_manual_labels.tif"
    candidate_path.parent.mkdir()
    reference_path.parent.mkdir()
    tifffile.imwrite(candidate_path, candidate)
    tifffile.imwrite(reference_path, reference)

    summaries = run_manual_mask_validation(
        candidate_mask_paths={"XY01": candidate_path},
        reference_mask_paths={"XY01": reference_path},
        output_dir=tmp_path / "validation",
        iou_threshold=0.5,
    )

    assert summaries[0]["image_id"] == "XY01"
    assert summaries[0]["matched_count"] == 1
    assert (tmp_path / "validation" / "manual_mask_validation_summary.csv").exists()
    assert (tmp_path / "validation" / "manual_mask_validation_matches.csv").exists()
    assert (tmp_path / "validation" / "logs" / "config_resolved.yaml").exists()
    assert (tmp_path / "validation" / "logs" / "run_log.txt").exists()


def test_run_manual_mask_validation_counts_reference_without_candidate_as_false_negative(
    tmp_path: Path,
):
    candidate = np.zeros((8, 8), dtype=np.uint32)
    candidate[1:3, 1:3] = 1
    reference_xy01 = np.zeros((8, 8), dtype=np.uint32)
    reference_xy01[1:3, 1:3] = 2
    reference_xy02 = np.zeros((8, 8), dtype=np.uint32)
    reference_xy02[4:6, 4:6] = 3
    candidate_path = tmp_path / "candidate" / "XY01_candidate_labels.tif"
    reference_path_xy01 = tmp_path / "manual" / "XY01_manual_labels.tif"
    reference_path_xy02 = tmp_path / "manual" / "XY02_manual_labels.tif"
    candidate_path.parent.mkdir()
    reference_path_xy01.parent.mkdir()
    tifffile.imwrite(candidate_path, candidate)
    tifffile.imwrite(reference_path_xy01, reference_xy01)
    tifffile.imwrite(reference_path_xy02, reference_xy02)

    summaries = run_manual_mask_validation(
        candidate_mask_paths={"XY01": candidate_path},
        reference_mask_paths={"XY01": reference_path_xy01, "XY02": reference_path_xy02},
        output_dir=tmp_path / "validation",
        iou_threshold=0.5,
    )

    by_id = {summary["image_id"]: summary for summary in summaries}
    assert by_id["XY01"]["matched_count"] == 1
    assert by_id["XY02"]["candidate_count"] == 0
    assert by_id["XY02"]["reference_count"] == 1
    assert by_id["XY02"]["false_negative_count"] == 1
    assert by_id["XY02"]["recall"] == 0.0


def test_run_manual_mask_validation_rejects_all_empty_reference_masks(tmp_path: Path):
    candidate = np.zeros((8, 8), dtype=np.uint32)
    candidate[1:3, 1:3] = 1
    reference_xy01 = np.zeros((8, 8), dtype=np.uint32)
    reference_xy02 = np.zeros((8, 8), dtype=np.uint32)
    candidate_path = tmp_path / "candidate" / "XY01_candidate_labels.tif"
    reference_path_xy01 = tmp_path / "manual" / "XY01_manual_labels.tif"
    reference_path_xy02 = tmp_path / "manual" / "XY02_manual_labels.tif"
    candidate_path.parent.mkdir()
    reference_path_xy01.parent.mkdir()
    tifffile.imwrite(candidate_path, candidate)
    tifffile.imwrite(reference_path_xy01, reference_xy01)
    tifffile.imwrite(reference_path_xy02, reference_xy02)

    with pytest.raises(ValueError, match="all reference masks are empty"):
        run_manual_mask_validation(
            candidate_mask_paths={"XY01": candidate_path},
            reference_mask_paths={"XY01": reference_path_xy01, "XY02": reference_path_xy02},
            output_dir=tmp_path / "validation",
            iou_threshold=0.5,
        )


def test_run_manual_mask_validation_rejects_incomplete_reference_status(tmp_path: Path):
    candidate = np.zeros((8, 8), dtype=np.uint32)
    candidate[1:3, 1:3] = 1
    reference = np.zeros((8, 8), dtype=np.uint32)
    reference[1:3, 1:3] = 2
    candidate_path = tmp_path / "candidate" / "XY01_candidate_labels.tif"
    reference_path = tmp_path / "manual" / "XY01_manual_labels.tif"
    candidate_path.parent.mkdir()
    reference_path.parent.mkdir()
    tifffile.imwrite(candidate_path, candidate)
    tifffile.imwrite(reference_path, reference)

    with pytest.raises(ValueError, match="not marked complete"):
        run_manual_mask_validation(
            candidate_mask_paths={"XY01": candidate_path},
            reference_mask_paths={"XY01": reference_path},
            output_dir=tmp_path / "validation",
            iou_threshold=0.5,
            reference_completion_status={"XY01": "not_started"},
        )


def test_run_manual_mask_validation_requires_empty_references_confirmed_empty(
    tmp_path: Path,
):
    candidate_xy01 = np.zeros((8, 8), dtype=np.uint32)
    candidate_xy01[1:3, 1:3] = 1
    reference_xy01 = np.zeros((8, 8), dtype=np.uint32)
    reference_xy01[1:3, 1:3] = 2
    reference_xy02 = np.zeros((8, 8), dtype=np.uint32)
    candidate_path_xy01 = tmp_path / "candidate" / "XY01_candidate_labels.tif"
    reference_path_xy01 = tmp_path / "manual" / "XY01_manual_labels.tif"
    reference_path_xy02 = tmp_path / "manual" / "XY02_manual_labels.tif"
    candidate_path_xy01.parent.mkdir()
    reference_path_xy01.parent.mkdir()
    tifffile.imwrite(candidate_path_xy01, candidate_xy01)
    tifffile.imwrite(reference_path_xy01, reference_xy01)
    tifffile.imwrite(reference_path_xy02, reference_xy02)

    with pytest.raises(ValueError, match="empty reference masks must be marked confirmed_empty"):
        run_manual_mask_validation(
            candidate_mask_paths={"XY01": candidate_path_xy01},
            reference_mask_paths={
                "XY01": reference_path_xy01,
                "XY02": reference_path_xy02,
            },
            output_dir=tmp_path / "validation",
            iou_threshold=0.5,
            reference_completion_status={
                "XY01": "complete_non_empty",
                "XY02": "complete_non_empty",
            },
        )


def test_run_manual_mask_validation_rejects_status_for_different_reference_mask(
    tmp_path: Path,
):
    candidate = np.zeros((8, 8), dtype=np.uint32)
    candidate[1:3, 1:3] = 1
    reference = np.zeros((8, 8), dtype=np.uint32)
    reference[1:3, 1:3] = 2
    candidate_path = tmp_path / "candidate" / "XY01_candidate_labels.tif"
    reference_path = tmp_path / "manual" / "XY01_manual_labels.tif"
    stale_reference_path = tmp_path / "other_package" / "XY01_manual_labels.tif"
    candidate_path.parent.mkdir()
    reference_path.parent.mkdir()
    tifffile.imwrite(candidate_path, candidate)
    tifffile.imwrite(reference_path, reference)

    with pytest.raises(ValueError, match="status mask path does not match"):
        run_manual_mask_validation(
            candidate_mask_paths={"XY01": candidate_path},
            reference_mask_paths={"XY01": reference_path},
            output_dir=tmp_path / "validation",
            iou_threshold=0.5,
            reference_completion_status={"XY01": "complete_non_empty"},
            reference_completion_mask_paths={"XY01": stale_reference_path},
        )


def _write_seeded_run(tmp_path: Path) -> Path:
    dataset = tmp_path / "dataset"
    xy01 = dataset / "XY01"
    xy01.mkdir(parents=True)
    ch2 = np.zeros((28, 28), dtype=np.uint16) + 100
    ch2[8:22, 7:23] = 1000
    ch4 = np.zeros((28, 28), dtype=np.uint16)
    ch4[12:15, 13:16] = 2000
    ch2_path = xy01 / "sample_XY01_CH2.tif"
    ch4_path = xy01 / "sample_XY01_CH4.tif"
    tifffile.imwrite(ch2_path, ch2)
    tifffile.imwrite(ch4_path, ch4)
    nuclei_mask = np.zeros((28, 28), dtype=np.uint32)
    nuclei_mask[12:15, 13:16] = 1
    nucleus_mask_path = tmp_path / "counts" / "masks" / "XY01_CH4_fake_labels.tif"
    nucleus_mask_path.parent.mkdir(parents=True)
    tifffile.imwrite(nucleus_mask_path, nuclei_mask)
    pair = ImagePair(location="XY01", source_id="XY01", ch2_path=ch2_path, ch4_path=ch4_path)
    output_dir = tmp_path / "seeded"
    run_seeded_region_batch(
        image_pairs=[pair],
        mask_lookup={"XY01": nucleus_mask_path},
        output_dir=output_dir,
        foreground_method="otsu",
        background_value=100,
    )
    return output_dir


def _rewrite_summary_values(path: Path, updates: dict[str, str]) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    for row in rows:
        row.update(updates)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))
