from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import tifffile

from dapi_norm.target_normalization import run_target_normalization


def test_run_target_normalization_writes_image_summary_with_background_corrected_endpoint(
    tmp_path: Path,
):
    input_root = tmp_path / "dataset"
    xy01 = input_root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    target = np.array([[10, 20], [30, 40]], dtype=np.uint16)
    dapi = np.array([[0, 65535], [10, 20]], dtype=np.uint16)
    tifffile.imwrite(xy01 / "sample_XY01_CH2.tif", target)
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", dapi)
    counts_dir = _write_count_outputs(tmp_path, nucleus_count=2, mask=np.array([[0, 1], [2, 0]]))

    output_dir = tmp_path / "target"
    rows = run_target_normalization(
        input_root=input_root,
        counts_dir=counts_dir,
        output_dir=output_dir,
        target_channel_id="CH2",
        dapi_channel_id="CH4",
        background_percentile=0,
    )

    assert len(rows) == 1
    row = rows[0]
    assert row["image_id"] == "XY01"
    assert row["well_id"] == "XY01"
    assert row["target_channel"] == "CH2"
    assert row["dapi_channel"] == "CH4"
    assert row["raw_nucleus_count"] == 2
    assert row["filtered_nucleus_count"] == 2
    assert row["target_area_px"] == 4
    assert row["target_integrated_raw"] == 100.0
    assert row["target_mean_raw"] == 25.0
    assert row["background_method"] == "percentile_0"
    assert row["background_value_per_px"] == 10.0
    assert row["target_integrated_background_corrected"] == 60.0
    assert row["target_integrated_intensity_per_DAPI_positive_nucleus"] == 30.0
    assert row["dapi_saturation_fraction"] == 0.25
    assert row["target_saturation_fraction"] == 0.0
    assert "channel_identity_unconfirmed" in row["warnings"]

    summary_csv = output_dir / "summaries" / "image_level_summary.csv"
    assert summary_csv.exists()
    assert "target_integrated_intensity_per_DAPI_positive_nucleus" in summary_csv.read_text()


def test_run_target_normalization_writes_visual_qc_artifacts(tmp_path: Path):
    input_root = tmp_path / "dataset"
    xy01 = input_root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    target = np.zeros((16, 16), dtype=np.uint16)
    target[4:12, 5:13] = 1000
    dapi = np.zeros((16, 16), dtype=np.uint16)
    dapi[5:8, 6:9] = 2000
    tifffile.imwrite(xy01 / "sample_XY01_CH2.tif", target)
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", dapi)
    mask = np.zeros((16, 16), dtype=np.uint32)
    mask[5:8, 6:9] = 1
    counts_dir = _write_count_outputs(tmp_path, nucleus_count=1, mask=mask)

    output_dir = tmp_path / "target"
    run_target_normalization(
        input_root=input_root,
        counts_dir=counts_dir,
        output_dir=output_dir,
        target_channel_id="CH2",
        dapi_channel_id="CH4",
    )

    assert (output_dir / "qc" / "XY01_CH2_target_with_CH4_nucleus_outlines.png").exists()
    assert (output_dir / "plots" / "normalized_intensity_by_well.png").exists()
    assert (output_dir / "plots" / "target_integrated_vs_nucleus_count.png").exists()
    assert (output_dir / "qc_contact_sheet.png").exists()
    assert (output_dir / "logs" / "config_resolved.yaml").exists()


def test_run_target_normalization_rejects_counts_from_wrong_channel(tmp_path: Path):
    input_root = tmp_path / "dataset"
    xy01 = input_root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    tifffile.imwrite(xy01 / "sample_XY01_CH2.tif", np.ones((4, 4), dtype=np.uint16))
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", np.ones((4, 4), dtype=np.uint16))
    counts_dir = _write_count_outputs(
        tmp_path, nucleus_count=1, mask=np.ones((4, 4), dtype=np.uint32), channel_id="CH2"
    )

    with pytest.raises(ValueError, match="Expected CH4 count row"):
        run_target_normalization(
            input_root=input_root,
            counts_dir=counts_dir,
            output_dir=tmp_path / "target",
            target_channel_id="CH2",
            dapi_channel_id="CH4",
        )


def test_run_target_normalization_rejects_count_input_from_wrong_directory(tmp_path: Path):
    input_root = tmp_path / "dataset"
    xy01 = input_root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    tifffile.imwrite(xy01 / "sample_XY01_CH2.tif", np.ones((4, 4), dtype=np.uint16))
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", np.ones((4, 4), dtype=np.uint16))
    wrong_root = tmp_path / "other_dataset"
    wrong_xy01 = wrong_root / "sample_02" / "XY01"
    wrong_xy01.mkdir(parents=True)
    wrong_dapi_path = wrong_xy01 / "sample_XY01_CH4.tif"
    tifffile.imwrite(wrong_dapi_path, np.ones((4, 4), dtype=np.uint16))
    counts_dir = _write_count_outputs(
        tmp_path,
        nucleus_count=1,
        mask=np.ones((4, 4), dtype=np.uint32),
        input_path=wrong_dapi_path,
    )

    with pytest.raises(ValueError, match="Count row for XY01 points to"):
        run_target_normalization(
            input_root=input_root,
            counts_dir=counts_dir,
            output_dir=tmp_path / "target",
            target_channel_id="CH2",
            dapi_channel_id="CH4",
        )


def _write_count_outputs(
    tmp_path: Path,
    *,
    nucleus_count: int,
    mask: np.ndarray,
    channel_id: str = "CH4",
    backend: str = "cellpose",
    input_path: Path | None = None,
) -> Path:
    counts_dir = tmp_path / "counts"
    summaries_dir = counts_dir / "summaries"
    masks_dir = counts_dir / "masks"
    qc_dir = counts_dir / "qc"
    summaries_dir.mkdir(parents=True)
    masks_dir.mkdir()
    qc_dir.mkdir()
    mask_path = masks_dir / "XY01_CH4_fake_labels.tif"
    qc_path = qc_dir / "XY01_CH4_fake_montage.png"
    tifffile.imwrite(mask_path, mask.astype(np.uint32), photometric="minisblack")
    qc_path.write_bytes(b"not-a-real-png")
    (summaries_dir / "nucleus_counts.csv").write_text(
        "image_id,input_path,backend,model_name,channel_id,candidate_stain,"
        "channel_identity_confirmed,nucleus_count,mask_path,qc_montage_path,warnings\n"
        f"XY01,{input_path or f'input/sample_XY01_{channel_id}.tif'},{backend},fake,{channel_id},candidate_DAPI,False,{nucleus_count},"
        f"{mask_path},{qc_path},channel_identity_unconfirmed\n",
        encoding="utf-8",
    )
    return counts_dir
