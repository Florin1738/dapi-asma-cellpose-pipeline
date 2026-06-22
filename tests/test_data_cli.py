from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import tifffile

from dapi_norm.data_cli import run_data_inventory


def test_run_data_inventory_writes_reports_and_selected_previews(tmp_path: Path):
    root = tmp_path / "dataset"
    xy01 = root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    tifffile.imwrite(xy01 / "sample_XY01_CH2.tif", np.full((16, 16), 10, dtype=np.uint16))
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", np.full((16, 16), 100, dtype=np.uint16))

    output = tmp_path / "out"
    run_data_inventory(root=root, output=output, preview_positions=["XY01"])

    assert (output / "image_manifest.csv").exists()
    assert (output / "dataset_summary.md").exists()
    assert (output / "previews" / "XY01_channels_preview.png").exists()


def test_run_data_inventory_fails_for_missing_root(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="Dataset root does not exist"):
        run_data_inventory(root=tmp_path / "missing", output=tmp_path / "out")


def test_run_data_inventory_fails_for_unknown_requested_preview_position(tmp_path: Path):
    root = tmp_path / "dataset"
    xy01 = root / "sample_02" / "XY01"
    xy01.mkdir(parents=True)
    tifffile.imwrite(xy01 / "sample_XY01_CH2.tif", np.full((16, 16), 10, dtype=np.uint16))

    with pytest.raises(ValueError, match="Requested preview positions not found: XY99"):
        run_data_inventory(root=root, output=tmp_path / "out", preview_positions=["XY99"])
