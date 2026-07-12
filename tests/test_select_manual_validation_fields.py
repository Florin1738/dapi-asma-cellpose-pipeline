from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import tifffile

from scripts.select_manual_validation_fields import build_validation_feature_records, _write_readme


def test_build_validation_feature_records_reads_pairs_masks_and_method_jaccard(tmp_path: Path):
    input_root, counts_root, method_summary = _write_synthetic_feature_inputs(tmp_path)

    records = build_validation_feature_records(
        input_root=input_root,
        nuclei_counts_root=counts_root,
        method_review_summary=method_summary,
    )

    assert len(records) == 1
    record = records[0]
    assert record.image_id == "XY01"
    assert record.source_id == "XY01"
    assert record.dapi_positive_nucleus_count == 2
    assert record.target_integrated_raw == pytest.approx(65635.0)
    assert record.target_integrated_raw_per_DAPI_positive_nucleus == pytest.approx(32817.5)
    assert record.target_saturation_fraction == pytest.approx(1 / 6)
    assert record.dapi_saturation_fraction == pytest.approx(1 / 6)
    assert record.method_region_jaccard == pytest.approx(0.25)


def test_build_validation_feature_records_rejects_duplicate_method_summary_ids(
    tmp_path: Path,
):
    input_root, counts_root, method_summary = _write_synthetic_feature_inputs(tmp_path)
    method_summary.write_text(
        "image_id,method_region_jaccard\nXY01,0.25\nXY01,0.5\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Duplicate image_id"):
        build_validation_feature_records(
            input_root=input_root,
            nuclei_counts_root=counts_root,
            method_review_summary=method_summary,
        )


def test_selection_readme_includes_audit_context_and_reason_table(tmp_path: Path):
    input_root, counts_root, method_summary = _write_synthetic_feature_inputs(tmp_path)
    records = build_validation_feature_records(
        input_root=input_root,
        nuclei_counts_root=counts_root,
        method_review_summary=method_summary,
    )
    records[0].selection_reasons = "must_include;low_raw_target_integrated"

    readme_path = tmp_path / "selection" / "README.md"
    _write_readme(
        readme_path,
        input_root=input_root,
        nuclei_counts_root=counts_root,
        output_dir=tmp_path / "selection",
        method_review_summary=method_summary,
        max_images=16,
        per_bucket=2,
        must_include=["XY01"],
        all_records=records,
        selected=records,
    )

    text = readme_path.read_text(encoding="utf-8")
    assert "Candidate fields evaluated: `1`" in text
    assert "Maximum selected fields requested: `16`" in text
    assert "```bash" in text
    assert "scripts/select_manual_validation_fields.py" in text
    assert "| XY01 | `must_include;low_raw_target_integrated`" in text


def _write_synthetic_feature_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    input_root = tmp_path / "input"
    xy_dir = input_root / "XY01"
    xy_dir.mkdir(parents=True)
    ch2 = np.array([[0, 65535, 10], [20, 30, 40]], dtype=np.uint16)
    ch4 = np.array([[0, 1, 2], [3, 4, 65535]], dtype=np.uint16)
    tifffile.imwrite(xy_dir / "sample_XY01_CH2.tif", ch2, photometric="minisblack")
    tifffile.imwrite(xy_dir / "sample_XY01_CH4.tif", ch4, photometric="minisblack")

    counts_root = tmp_path / "counts"
    masks_dir = counts_root / "masks"
    summaries_dir = counts_root / "summaries"
    masks_dir.mkdir(parents=True)
    summaries_dir.mkdir(parents=True)
    mask = np.array([[0, 1, 1], [2, 2, 0]], dtype=np.uint32)
    mask_path = masks_dir / "XY01_CH4_cpsam_v2_labels.tif"
    tifffile.imwrite(mask_path, mask, photometric="minisblack")
    (summaries_dir / "nucleus_counts.csv").write_text(
        "image_id,mask_path\nXY01,masks/XY01_CH4_cpsam_v2_labels.tif\n",
        encoding="utf-8",
    )

    method_summary = tmp_path / "method_comparison_review_summary.csv"
    method_summary.write_text(
        "image_id,method_region_jaccard\nXY01,0.25\n",
        encoding="utf-8",
    )
    return input_root, counts_root, method_summary
