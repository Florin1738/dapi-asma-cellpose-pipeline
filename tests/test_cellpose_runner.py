from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import tifffile
import yaml

from dapi_norm.cellpose_runner import count_labels, run_nuclei_count_batch


def test_count_labels_counts_objects_not_max_label():
    labels = np.array([[0, 2], [5, 0]], dtype=np.uint16)

    assert count_labels(labels) == 2


def test_run_nuclei_count_batch_writes_masks_qc_and_tables_with_fake_segmenter(tmp_path: Path):
    input_root = tmp_path / "dataset"
    xy01 = input_root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    image = np.zeros((24, 24, 3), dtype=np.uint16)
    image[4:8, 5:9, 2] = 2000
    image[12:17, 14:19, 2] = 2500
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", image, photometric="rgb")

    def fake_segmenter(_image: np.ndarray) -> np.ndarray:
        labels = np.zeros((24, 24), dtype=np.uint16)
        labels[4:8, 5:9] = 1
        labels[12:17, 14:19] = 2
        return labels

    output_dir = tmp_path / "out"
    summaries = run_nuclei_count_batch(
        input_root=input_root,
        output_dir=output_dir,
        channel_id="CH4",
        model_name="fake_model",
        segmenter=fake_segmenter,
    )

    assert summaries[0]["image_id"] == "XY01"
    assert summaries[0]["nucleus_count"] == 2
    assert (output_dir / "masks" / "XY01_CH4_fake_model_labels.tif").exists()
    assert (output_dir / "qc" / "XY01_CH4_fake_model_montage.png").exists()
    assert (output_dir / "qc_contact_sheet.png").exists()
    assert (output_dir / "summaries" / "nucleus_counts.csv").exists()
    assert (output_dir / "summaries" / "per_nucleus_locations.csv").exists()
    assert (output_dir / "logs" / "config_resolved.yaml").exists()
    assert (output_dir / "logs" / "run_log.txt").exists()
    counts_csv = (output_dir / "summaries" / "nucleus_counts.csv").read_text(encoding="utf-8")
    assert "XY01" in counts_csv
    assert "channel_identity_unconfirmed" in counts_csv
    assert "kept_after_filtering" in (
        output_dir / "summaries" / "per_nucleus_locations.csv"
    ).read_text(encoding="utf-8")
    config = yaml.safe_load((output_dir / "logs" / "config_resolved.yaml").read_text())
    assert config["channel_id"] == "CH4"
    assert config["model"]["name"] == "fake_model"
    assert config["segmentation_parameters"]["flow_threshold"] == 0.4
    assert config["image_inputs"][0]["extracted_plane_shape"] == [24, 24]


def test_run_nuclei_count_batch_records_confirmed_channel_identity(tmp_path: Path):
    input_root = tmp_path / "dataset"
    xy01 = input_root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    image = np.zeros((12, 12), dtype=np.uint16)
    image[3:6, 4:7] = 2000
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", image)

    def fake_segmenter(_image: np.ndarray) -> np.ndarray:
        labels = np.zeros((12, 12), dtype=np.uint16)
        labels[3:6, 4:7] = 1
        return labels

    output_dir = tmp_path / "out"
    summaries = run_nuclei_count_batch(
        input_root=input_root,
        output_dir=output_dir,
        channel_id="CH4",
        model_name="fake_model",
        segmenter=fake_segmenter,
        channel_identity_confirmed=True,
    )

    assert summaries[0]["channel_identity_confirmed"] is True
    assert summaries[0]["warnings"] == ""
    counts_csv = (output_dir / "summaries" / "nucleus_counts.csv").read_text(encoding="utf-8")
    assert "candidate_DAPI,True,1" in counts_csv
    assert "channel_identity_unconfirmed" not in counts_csv
    config = yaml.safe_load((output_dir / "logs" / "config_resolved.yaml").read_text())
    assert config["channel_identity_confirmed"] is True
    assert config["image_inputs"][0]["channel_identity_confirmed"] is True
    run_log = (output_dir / "logs" / "run_log.txt").read_text(encoding="utf-8")
    assert "channel_identity_confirmed: true" in run_log


def test_run_nuclei_count_batch_rejects_wrong_shape_segmenter_output(tmp_path: Path):
    input_root = tmp_path / "dataset"
    xy01 = input_root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    image = np.zeros((24, 24, 3), dtype=np.uint16)
    image[4:8, 5:9, 2] = 2000
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", image, photometric="rgb")

    def bad_segmenter(_image: np.ndarray) -> np.ndarray:
        return np.zeros((12, 12), dtype=np.uint16)

    with pytest.raises(ValueError, match="same shape"):
        run_nuclei_count_batch(
            input_root=input_root,
            output_dir=tmp_path / "out",
            channel_id="CH4",
            model_name="fake_model",
            segmenter=bad_segmenter,
        )


def test_run_nuclei_count_batch_preserves_large_label_ids(tmp_path: Path):
    input_root = tmp_path / "dataset"
    xy01 = input_root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    image = np.zeros((24, 24, 3), dtype=np.uint16)
    image[4:8, 5:9, 2] = 2000
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", image, photometric="rgb")

    def large_label_segmenter(_image: np.ndarray) -> np.ndarray:
        labels = np.zeros((24, 24), dtype=np.uint32)
        labels[4:8, 5:9] = 70_000
        return labels

    output_dir = tmp_path / "out"
    summaries = run_nuclei_count_batch(
        input_root=input_root,
        output_dir=output_dir,
        channel_id="CH4",
        model_name="fake_model",
        segmenter=large_label_segmenter,
    )

    written_mask = tifffile.imread(output_dir / "masks" / "XY01_CH4_fake_model_labels.tif")
    assert summaries[0]["nucleus_count"] == 1
    assert written_mask.max() == 70_000
