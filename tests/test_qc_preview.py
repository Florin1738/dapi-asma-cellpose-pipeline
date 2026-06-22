from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from dapi_norm.qc_preview import compose_channel_rgb, generate_position_preview


def test_generate_position_preview_writes_png(tmp_path: Path):
    ch2_path = tmp_path / "sample_XY01_CH2.tif"
    ch4_path = tmp_path / "sample_XY01_CH4.tif"
    out_path = tmp_path / "preview.png"

    y, x = np.ogrid[:32, :32]
    ch2 = np.zeros((32, 32), dtype=np.uint16)
    ch2[(y - 12) ** 2 + (x - 12) ** 2 < 16] = 2000
    ch2[(y - 22) ** 2 + (x - 22) ** 2 < 9] = 1800
    ch4 = np.full((32, 32), 200, dtype=np.uint16)
    ch4[8:24, 8:24] = 1200

    tifffile.imwrite(ch2_path, ch2)
    tifffile.imwrite(ch4_path, ch4)

    generate_position_preview(
        position_id="XY01",
        channel_paths={"CH2": ch2_path, "CH4": ch4_path},
        output_path=out_path,
    )

    assert out_path.exists()
    assert out_path.stat().st_size > 1000


def test_compose_channel_rgb_maps_ch2_to_red_and_ch4_to_blue():
    ch2 = np.ones((2, 2), dtype=np.float32)
    ch4 = np.full((2, 2), 0.5, dtype=np.float32)

    composite = compose_channel_rgb({"CH2": ch2, "CH4": ch4})

    np.testing.assert_allclose(composite[..., 0], 1.0)
    np.testing.assert_allclose(composite[..., 1], 0.0)
    np.testing.assert_allclose(composite[..., 2], 0.5)


def test_generate_position_preview_accepts_pseudocolor_rgb_tiffs(tmp_path: Path):
    ch2_path = tmp_path / "sample_XY01_CH2.tif"
    ch4_path = tmp_path / "sample_XY01_CH4.tif"
    out_path = tmp_path / "preview.png"

    ch2 = np.zeros((32, 32, 3), dtype=np.uint16)
    ch4 = np.zeros((32, 32, 3), dtype=np.uint16)
    ch2[..., 0] = 500
    ch4[..., 2] = 1200
    tifffile.imwrite(ch2_path, ch2, photometric="rgb")
    tifffile.imwrite(ch4_path, ch4, photometric="rgb")

    generate_position_preview(
        position_id="XY01",
        channel_paths={"CH2": ch2_path, "CH4": ch4_path},
        output_path=out_path,
    )

    assert out_path.exists()
    assert out_path.stat().st_size > 1000
