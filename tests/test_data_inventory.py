from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

import pytest

from dapi_norm.data_inventory import inventory_dataset, parse_image_filename, write_inventory_reports


def write_tiff(path: Path, value: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    image = np.full((6, 8), value, dtype=np.uint16)
    tifffile.imwrite(path, image)


def test_parse_image_filename_identifies_sample_position_and_channel():
    parsed = parse_image_filename("ApYYM20AGGSMA_XY07_CH4.tif")

    assert parsed.sample_id == "ApYYM20AGGSMA"
    assert parsed.position_id == "XY07"
    assert parsed.channel_id == "CH4"
    assert parsed.kind == "channel"


def test_parse_image_filename_identifies_overlay():
    parsed = parse_image_filename("ApYYM20AGGSMA_XY07_Overlay.tif")

    assert parsed.sample_id == "ApYYM20AGGSMA"
    assert parsed.position_id == "XY07"
    assert parsed.channel_id is None
    assert parsed.kind == "overlay"


def test_inventory_groups_channels_and_ignores_sidecars(tmp_path: Path):
    data_root = tmp_path / "ApYYM20AGGSMA_02"
    write_tiff(data_root / "ApYYM20AGGSMA_02" / "XY01" / "ApYYM20AGGSMA_XY01_CH2.tif", 20)
    write_tiff(data_root / "ApYYM20AGGSMA_02" / "XY01" / "ApYYM20AGGSMA_XY01_CH4.tif", 40)
    write_tiff(
        data_root / "ApYYM20AGGSMA_02" / "XY01" / "ApYYM20AGGSMA_XY01_Overlay.tif",
        60,
    )
    (data_root / "ApYYM20AGGSMA_02" / "H02.lnk").write_text("shortcut", encoding="utf-8")
    (data_root / "__ApYYM20AGGSMA_02" / "XY01").mkdir(parents=True)
    (data_root / "__ApYYM20AGGSMA_02" / "XY01" / "_H02_Error.txt").write_text(
        "cannot download", encoding="utf-8"
    )

    inventory = inventory_dataset(data_root)

    assert inventory.root == data_root
    assert inventory.file_counts_by_suffix == {".lnk": 1, ".tif": 3, ".txt": 1}
    assert list(inventory.positions) == ["XY01"]
    position = inventory.positions["XY01"]
    assert sorted(position.channels) == ["CH2", "CH4"]
    assert position.overlay is not None
    assert position.channels["CH2"].width == 8
    assert position.channels["CH2"].height == 6
    assert position.channels["CH2"].dtype == "uint16"
    assert position.channels["CH2"].mean_intensity == 20.0
    assert inventory.sidecar_count == 2


def test_inventory_extracts_active_sample_from_pseudocolor_rgb_tiff(tmp_path: Path):
    data_root = tmp_path / "ApYYM20AGGSMA_02"
    image = np.zeros((6, 8, 3), dtype=np.uint16)
    image[..., 2] = 123
    write_path = data_root / "ApYYM20AGGSMA_02" / "XY01" / "ApYYM20AGGSMA_XY01_CH4.tif"
    write_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(write_path, image, photometric="rgb")

    inventory = inventory_dataset(data_root)

    info = inventory.positions["XY01"].channels["CH4"]
    assert info.width == 8
    assert info.height == 6
    assert info.dtype == "uint16"
    assert info.mean_intensity == 123.0


def test_inventory_rejects_multi_active_rgb_tiff_for_intensity_stats(tmp_path: Path):
    data_root = tmp_path / "ApYYM20AGGSMA_02"
    image = np.zeros((6, 8, 3), dtype=np.uint16)
    image[..., 0] = 50
    image[..., 2] = 123
    write_path = data_root / "ApYYM20AGGSMA_02" / "XY01" / "ApYYM20AGGSMA_XY01_CH2.tif"
    write_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(write_path, image, photometric="rgb")

    with pytest.raises(ValueError, match="multi-active RGB"):
        inventory_dataset(data_root)


def test_write_inventory_reports_creates_manifest_and_markdown(tmp_path: Path):
    data_root = tmp_path / "ApYYM20AGGSMA_02"
    write_tiff(data_root / "ApYYM20AGGSMA_02" / "XY01" / "ApYYM20AGGSMA_XY01_CH2.tif", 20)
    write_tiff(data_root / "ApYYM20AGGSMA_02" / "XY01" / "ApYYM20AGGSMA_XY01_CH4.tif", 40)
    inventory = inventory_dataset(data_root)

    output_dir = tmp_path / "reports"
    write_inventory_reports(inventory, output_dir)

    manifest = output_dir / "image_manifest.csv"
    interpretation = output_dir / "channel_interpretation_manifest.csv"
    summary = output_dir / "dataset_summary.md"
    assert manifest.exists()
    assert interpretation.exists()
    assert summary.exists()
    manifest_text = manifest.read_text(encoding="utf-8")
    assert "position_id,kind,channel_id,path,width,height,dtype" in manifest_text
    assert "XY01,channel,CH2" in manifest_text
    interpretation_text = interpretation.read_text(encoding="utf-8")
    assert "position_id,source_file,filename_channel,rgb_component,candidate_stain" in (
        interpretation_text
    )
    assert "requires_confirmation" in interpretation_text
    assert "XY01,ApYYM20AGGSMA_02/XY01/ApYYM20AGGSMA_XY01_CH2.tif,CH2,red,candidate_target,,unconfirmed,unconfirmed,true" in (
        interpretation_text
    )
    assert "XY01,ApYYM20AGGSMA_02/XY01/ApYYM20AGGSMA_XY01_CH4.tif,CH4,blue,candidate_DAPI,,unconfirmed,unconfirmed,true" in (
        interpretation_text
    )
    summary_text = summary.read_text(encoding="utf-8")
    assert "# Dataset Summary" in summary_text
    assert "Positions: 1" in summary_text
    assert "Channel IDs: CH2, CH4" in summary_text
