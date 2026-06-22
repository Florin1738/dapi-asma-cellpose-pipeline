from __future__ import annotations

from pathlib import Path
import csv

import numpy as np
import pytest
import tifffile

from dapi_norm.target_normalization import run_target_normalization
from dapi_norm.target_validation import validate_target_outputs


def test_validate_target_outputs_checks_formula_and_artifacts(tmp_path: Path):
    input_root, counts_dir = _write_dataset_and_counts(tmp_path)
    output_dir = tmp_path / "target"
    run_target_normalization(
        input_root=input_root,
        counts_dir=counts_dir,
        output_dir=output_dir,
        target_channel_id="CH2",
        dapi_channel_id="CH4",
        background_percentile=0,
    )

    result = validate_target_outputs(output_dir)

    assert result["summary_rows"] == 1
    assert result["plots_exist"] is True
    assert result["qc_overlays_exist"] is True
    assert result["formulas_match"] is True


def test_validate_target_outputs_rejects_corrupted_normalized_value(tmp_path: Path):
    input_root, counts_dir = _write_dataset_and_counts(tmp_path)
    output_dir = tmp_path / "target"
    run_target_normalization(
        input_root=input_root,
        counts_dir=counts_dir,
        output_dir=output_dir,
        target_channel_id="CH2",
        dapi_channel_id="CH4",
        background_percentile=0,
    )
    updates = {"target_integrated_intensity_per_DAPI_positive_nucleus": "999"}
    _rewrite_summary_values(output_dir / "summaries" / "image_level_summary.csv", updates)
    _rewrite_summary_values(output_dir / "summaries" / "well_level_summary.csv", updates)

    with pytest.raises(ValueError, match="normalized intensity"):
        validate_target_outputs(output_dir)


def test_validate_target_outputs_rejects_mask_count_mismatch(tmp_path: Path):
    input_root, counts_dir = _write_dataset_and_counts(tmp_path)
    output_dir = tmp_path / "target"
    run_target_normalization(
        input_root=input_root,
        counts_dir=counts_dir,
        output_dir=output_dir,
        target_channel_id="CH2",
        dapi_channel_id="CH4",
        background_percentile=0,
    )
    summary_path = output_dir / "summaries" / "image_level_summary.csv"
    _rewrite_summary_values(
        summary_path,
        {"filtered_nucleus_count": "1", "target_integrated_intensity_per_DAPI_positive_nucleus": "60"},
    )
    _rewrite_summary_values(
        output_dir / "summaries" / "well_level_summary.csv",
        {"filtered_nucleus_count": "1", "target_integrated_intensity_per_DAPI_positive_nucleus": "60"},
    )

    with pytest.raises(ValueError, match="mask label count"):
        validate_target_outputs(output_dir)


def test_validate_target_outputs_recomputes_background_from_method(tmp_path: Path):
    input_root, counts_dir = _write_dataset_and_counts(tmp_path)
    output_dir = tmp_path / "target"
    run_target_normalization(
        input_root=input_root,
        counts_dir=counts_dir,
        output_dir=output_dir,
        target_channel_id="CH2",
        dapi_channel_id="CH4",
        background_percentile=10,
    )
    updates = {
        "background_value_per_px": "0",
        "target_integrated_background_corrected": "100",
        "target_integrated_intensity_per_DAPI_positive_nucleus": "50",
    }
    _rewrite_summary_values(output_dir / "summaries" / "image_level_summary.csv", updates)
    _rewrite_summary_values(output_dir / "summaries" / "well_level_summary.csv", updates)

    with pytest.raises(ValueError, match="background_value_per_px"):
        validate_target_outputs(output_dir)


def test_validate_target_outputs_rejects_well_summary_divergence(tmp_path: Path):
    input_root, counts_dir = _write_dataset_and_counts(tmp_path)
    output_dir = tmp_path / "target"
    run_target_normalization(
        input_root=input_root,
        counts_dir=counts_dir,
        output_dir=output_dir,
        target_channel_id="CH2",
        dapi_channel_id="CH4",
        background_percentile=0,
    )
    _rewrite_summary_values(
        output_dir / "summaries" / "well_level_summary.csv",
        {"target_integrated_intensity_per_DAPI_positive_nucleus": "999"},
    )

    with pytest.raises(ValueError, match="well_level_summary"):
        validate_target_outputs(output_dir)


def test_validate_target_outputs_resolves_repo_relative_paths_from_absolute_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    project_dir = tmp_path / "project"
    input_root, counts_dir = _write_dataset_and_counts(project_dir)
    output_dir = project_dir / "output" / "target_normalization"
    run_target_normalization(
        input_root=input_root,
        counts_dir=counts_dir,
        output_dir=output_dir,
        target_channel_id="CH2",
        dapi_channel_id="CH4",
        background_percentile=0,
    )
    _rewrite_target_paths_relative_to_project(output_dir, project_dir)

    monkeypatch.chdir(tmp_path)

    result = validate_target_outputs(output_dir.resolve())

    assert result["summary_rows"] == 1
    assert result["formulas_match"] is True


def _write_dataset_and_counts(tmp_path: Path) -> tuple[Path, Path]:
    input_root = tmp_path / "dataset"
    xy01 = input_root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    target = np.array([[10, 20], [30, 40]], dtype=np.uint16)
    dapi = np.array([[0, 65535], [10, 20]], dtype=np.uint16)
    tifffile.imwrite(xy01 / "sample_XY01_CH2.tif", target)
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", dapi)

    counts_dir = tmp_path / "counts"
    summaries_dir = counts_dir / "summaries"
    masks_dir = counts_dir / "masks"
    qc_dir = counts_dir / "qc"
    summaries_dir.mkdir(parents=True)
    masks_dir.mkdir()
    qc_dir.mkdir()
    mask_path = masks_dir / "XY01_CH4_fake_labels.tif"
    qc_path = qc_dir / "XY01_CH4_fake_montage.png"
    tifffile.imwrite(mask_path, np.array([[0, 1], [2, 0]], dtype=np.uint32))
    qc_path.write_bytes(b"not-a-real-png")
    (summaries_dir / "nucleus_counts.csv").write_text(
        "image_id,input_path,backend,model_name,channel_id,candidate_stain,"
        "channel_identity_confirmed,nucleus_count,mask_path,qc_montage_path,warnings\n"
        f"XY01,input/sample_XY01_CH4.tif,cellpose,fake,CH4,candidate_DAPI,False,2,"
        f"{mask_path},{qc_path},channel_identity_unconfirmed\n",
        encoding="utf-8",
    )
    return input_root, counts_dir


def _rewrite_target_paths_relative_to_project(output_dir: Path, project_dir: Path) -> None:
    summary_path = output_dir / "summaries" / "image_level_summary.csv"
    with summary_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    for row in rows:
        row["input_path"] = str(Path(row["input_path"]).relative_to(project_dir))
        row["mask_path"] = str(Path(row["mask_path"]).relative_to(project_dir))
        row["qc_overlay_path"] = str(Path(row["qc_overlay_path"]).relative_to(project_dir))
    for csv_path in [
        output_dir / "summaries" / "image_level_summary.csv",
        output_dir / "summaries" / "well_level_summary.csv",
    ]:
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)


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
