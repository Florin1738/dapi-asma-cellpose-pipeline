from __future__ import annotations

from pathlib import Path
import csv

import numpy as np
import pytest
import tifffile
import yaml

from dapi_norm.cellpose_runner import run_nuclei_count_batch
from dapi_norm.count_validation import validate_count_outputs


def test_validate_count_outputs_checks_generated_batch_artifacts(tmp_path: Path):
    input_root = tmp_path / "dataset"
    xy01 = input_root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    image = np.zeros((20, 20, 3), dtype=np.uint16)
    image[3:6, 4:7, 2] = 2000
    image[12:16, 13:17, 2] = 2500
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", image, photometric="rgb")

    def fake_segmenter(_image: np.ndarray) -> np.ndarray:
        labels = np.zeros((20, 20), dtype=np.uint32)
        labels[3:6, 4:7] = 1
        labels[12:16, 13:17] = 2
        return labels

    output_dir = tmp_path / "out"
    run_nuclei_count_batch(
        input_root=input_root,
        output_dir=output_dir,
        channel_id="CH4",
        model_name="fake_model",
        segmenter=fake_segmenter,
    )

    result = validate_count_outputs(output_dir)

    assert result["summary_rows"] == 1
    assert result["total_nucleus_count"] == 2
    assert result["per_nucleus_rows"] == 2
    assert result["mask_counts_match_csv"] is True


def test_validate_count_outputs_rejects_mismatched_per_nucleus_rows(tmp_path: Path):
    output_dir = tmp_path / "out"
    run_nuclei_count_batch(
        input_root=_write_one_image_dataset(tmp_path),
        output_dir=output_dir,
        channel_id="CH4",
        model_name="fake_model",
        segmenter=lambda image: np.ones(image.shape, dtype=np.uint32),
    )
    per_nucleus_path = output_dir / "summaries" / "per_nucleus_locations.csv"
    per_nucleus_path.write_text(per_nucleus_path.read_text(encoding="utf-8").splitlines()[0] + "\n")

    with pytest.raises(ValueError, match="per-nucleus rows"):
        validate_count_outputs(output_dir)


def test_validate_count_outputs_resolves_paths_from_absolute_output_dir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    project_dir = tmp_path / "project"
    output_dir = project_dir / "output" / "cellpose_counts"
    run_nuclei_count_batch(
        input_root=_write_one_image_dataset(project_dir),
        output_dir=output_dir,
        channel_id="CH4",
        model_name="fake_model",
        segmenter=lambda image: np.ones(image.shape, dtype=np.uint32),
    )
    _rewrite_summary_paths_relative_to_project(output_dir, project_dir)

    monkeypatch.chdir(tmp_path)

    result = validate_count_outputs(output_dir.resolve())

    assert result["summary_rows"] == 1
    assert result["mask_counts_match_csv"] is True


def test_validate_count_outputs_rejects_wrong_shaped_mask(tmp_path: Path):
    output_dir = tmp_path / "out"
    run_nuclei_count_batch(
        input_root=_write_one_image_dataset(tmp_path),
        output_dir=output_dir,
        channel_id="CH4",
        model_name="fake_model",
        segmenter=lambda image: np.ones(image.shape, dtype=np.uint32),
    )
    summary_path = output_dir / "summaries" / "nucleus_counts.csv"
    with summary_path.open(encoding="utf-8", newline="") as handle:
        summary_row = next(csv.DictReader(handle))
    mask_path = Path(summary_row["mask_path"])
    tifffile.imwrite(mask_path, np.ones((1, 1), dtype=np.uint32), photometric="minisblack")

    with pytest.raises(ValueError, match="shape"):
        validate_count_outputs(output_dir)


def test_validate_count_outputs_rejects_missing_image_shape_metadata(tmp_path: Path):
    output_dir = tmp_path / "out"
    run_nuclei_count_batch(
        input_root=_write_one_image_dataset(tmp_path),
        output_dir=output_dir,
        channel_id="CH4",
        model_name="fake_model",
        segmenter=lambda image: np.ones(image.shape, dtype=np.uint32),
    )
    config_path = output_dir / "logs" / "config_resolved.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["image_inputs"] = []
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    with pytest.raises(ValueError, match="No image shape metadata"):
        validate_count_outputs(output_dir)


def _write_one_image_dataset(tmp_path: Path) -> Path:
    input_root = tmp_path / "dataset"
    xy01 = input_root / "sample_02" / "XY01"
    xy01.mkdir(parents=True, exist_ok=True)
    image = np.zeros((10, 10, 3), dtype=np.uint16)
    image[..., 2] = 100
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", image, photometric="rgb")
    return input_root


def _rewrite_summary_paths_relative_to_project(output_dir: Path, project_dir: Path) -> None:
    summary_path = output_dir / "summaries" / "nucleus_counts.csv"
    with summary_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    for row in rows:
        row["mask_path"] = str(Path(row["mask_path"]).relative_to(project_dir))
        row["qc_montage_path"] = str(Path(row["qc_montage_path"]).relative_to(project_dir))
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
