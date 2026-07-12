from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from dapi_norm.method_review_panels import (
    MethodReviewRecord,
    summarize_method_overlap,
    write_method_review_package,
)
from scripts.render_method_review_panels import _metrics_lookup


def test_summarize_method_overlap_reports_union_jaccard_and_method_only_area():
    propagation = np.zeros((8, 8), dtype=np.uint32)
    cellpose = np.zeros((8, 8), dtype=np.uint32)
    propagation[1:5, 1:5] = 1
    cellpose[3:7, 3:7] = 2

    summary = summarize_method_overlap(
        image_id="XY01",
        nuclei_mask=np.zeros((8, 8), dtype=np.uint32),
        propagation_labels=propagation,
        cellpose_labels=cellpose,
        propagation_metrics={"qc_status": "reviewable_not_validated"},
        cellpose_metrics={"qc_status": "needs_manual_review"},
        crop_box=(0, 0, 8, 8),
    )

    assert summary["image_id"] == "XY01"
    assert summary["propagation_region_area_px"] == 16
    assert summary["cellpose_region_area_px"] == 16
    assert summary["both_region_area_px"] == 4
    assert summary["propagation_only_area_px"] == 12
    assert summary["cellpose_only_area_px"] == 12
    assert summary["method_region_jaccard"] == pytest.approx(4 / 28)
    assert summary["interpretation"] == "qualitative_qc_only_not_manual_validation"


def test_write_method_review_package_writes_panels_and_summary(tmp_path: Path):
    record = _synthetic_record()

    outputs = write_method_review_package(
        records=[record],
        output_dir=tmp_path / "review",
        crop_size=12,
    )

    assert outputs["full_field_panel"].exists()
    assert outputs["crop_panel"].exists()
    assert outputs["summary_csv"].exists()
    assert outputs["readme"].exists()
    rows = _read_rows(outputs["summary_csv"])
    assert len(rows) == 1
    assert rows[0]["image_id"] == "XY01"
    assert rows[0]["interpretation"] == "qualitative_qc_only_not_manual_validation"
    assert "manual/reference masks" in outputs["readme"].read_text(encoding="utf-8")


def test_write_method_review_package_rejects_shape_mismatch(tmp_path: Path):
    record = _synthetic_record(
        cellpose_labels=np.zeros((9, 8), dtype=np.uint32),
    )

    with pytest.raises(ValueError, match="shape mismatch"):
        write_method_review_package(
            records=[record],
            output_dir=tmp_path / "review",
        )


def test_metrics_lookup_rejects_duplicate_image_keys(tmp_path: Path):
    metrics_path = tmp_path / "metrics.csv"
    _write_csv(
        metrics_path,
        [
            {"image_id": "XY01", "source_id": "XY01", "mask_path": "a.tif"},
            {"image_id": "XY01", "source_id": "XY01", "mask_path": "b.tif"},
        ],
    )

    with pytest.raises(ValueError, match="Duplicate metrics key XY01"):
        _metrics_lookup(metrics_path)


def _synthetic_record(cellpose_labels: np.ndarray | None = None) -> MethodReviewRecord:
    ch2 = np.zeros((24, 24), dtype=np.uint16)
    ch2[3:16, 3:16] = 180
    ch2[12:22, 12:22] = 90
    ch4 = np.zeros((24, 24), dtype=np.uint16)
    ch4[5:7, 5:7] = 255
    ch4[17:19, 17:19] = 255
    nuclei = np.zeros((24, 24), dtype=np.uint32)
    nuclei[5:7, 5:7] = 1
    nuclei[17:19, 17:19] = 2
    propagation = np.zeros((24, 24), dtype=np.uint32)
    propagation[3:16, 3:16] = 1
    if cellpose_labels is None:
        cellpose_labels = np.zeros((24, 24), dtype=np.uint32)
        cellpose_labels[10:22, 10:22] = 1
    return MethodReviewRecord(
        image_id="XY01",
        ch2_image=ch2,
        ch4_image=ch4,
        nuclei_mask=nuclei,
        propagation_labels=propagation,
        cellpose_labels=cellpose_labels,
        propagation_metrics={
            "seeded_region_intensity_per_DAPI_positive_nucleus": 123.0,
            "qc_status": "reviewable_not_validated",
            "qc_flags": "not_validated_whole_cell_mask",
        },
        cellpose_metrics={
            "target_integrated_intensity_per_DAPI_positive_nucleus": 456.0,
            "qc_status": "needs_manual_review",
            "qc_flags": "not_validated_whole_cell_mask",
        },
    )


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
