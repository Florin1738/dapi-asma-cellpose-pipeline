from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile

from dapi_norm.centroid_overlays import render_centroid_overlays


def test_render_centroid_overlays_writes_overlay_and_contact_sheet(tmp_path: Path):
    image_dir = tmp_path / "data" / "plate 1" / "RunA" / "XY01"
    image_dir.mkdir(parents=True)
    image_path = image_dir / "sample_XY01_CH4.tif"
    image = np.zeros((32, 32), dtype=np.uint16)
    image[10, 12] = 500
    image[20, 22] = 700
    tifffile.imwrite(image_path, image)

    summary_dir = tmp_path / "counts" / "Plate_1" / "RunA" / "summaries"
    summary_dir.mkdir(parents=True)
    (summary_dir / "nucleus_counts.csv").write_text(
        "image_id,input_path,backend,model_name,channel_id,candidate_stain,"
        "channel_identity_confirmed,nucleus_count,mask_path,qc_montage_path,warnings\n"
        f"XY01,{image_path},cellpose,cpsam_v2,CH4,candidate_DAPI,True,"
        "2,masks/XY01.tif,qc/XY01.png,\n",
        encoding="utf-8",
    )
    (summary_dir / "per_nucleus_locations.csv").write_text(
        "image_id,input_path,nucleus_id,x_centroid,y_centroid,area_px,bbox_min_row,"
        "bbox_min_col,bbox_max_row,bbox_max_col,touches_border,kept_after_filtering\n"
        f"XY01,{image_path},1,12,10,20,8,10,12,14,False,True\n"
        f"XY01,{image_path},2,22,20,22,18,20,22,24,False,True\n",
        encoding="utf-8",
    )

    result = render_centroid_overlays(
        counts_root=tmp_path / "counts",
        output_dir=tmp_path / "overlays",
    )

    assert result == {"overlay_count": 1, "contact_sheet_count": 1}
    assert (tmp_path / "overlays" / "Plate_1" / "RunA" / "XY01_CH4_green_centroids.png").exists()
    assert (tmp_path / "overlays" / "Plate_1" / "RunA" / "contact_sheet_green_centroids.png").exists()


def test_render_centroid_overlays_writes_overlay_for_zero_count_image(tmp_path: Path):
    image_dir = tmp_path / "data" / "plate 1" / "RunA" / "XY02"
    image_dir.mkdir(parents=True)
    image_path = image_dir / "sample_XY02_CH4.tif"
    tifffile.imwrite(image_path, np.zeros((32, 32), dtype=np.uint16))

    summary_dir = tmp_path / "counts" / "Plate_1" / "RunA" / "summaries"
    summary_dir.mkdir(parents=True)
    (summary_dir / "nucleus_counts.csv").write_text(
        "image_id,input_path,backend,model_name,channel_id,candidate_stain,"
        "channel_identity_confirmed,nucleus_count,mask_path,qc_montage_path,warnings\n"
        f"XY02,{image_path},cellpose,cpsam_v2,CH4,candidate_DAPI,True,"
        "0,masks/XY02.tif,qc/XY02.png,\n",
        encoding="utf-8",
    )
    (summary_dir / "per_nucleus_locations.csv").write_text(
        "image_id,input_path,nucleus_id,x_centroid,y_centroid,area_px,bbox_min_row,"
        "bbox_min_col,bbox_max_row,bbox_max_col,touches_border,kept_after_filtering\n",
        encoding="utf-8",
    )

    result = render_centroid_overlays(
        counts_root=tmp_path / "counts",
        output_dir=tmp_path / "overlays",
    )

    assert result == {"overlay_count": 1, "contact_sheet_count": 1}
    assert (tmp_path / "overlays" / "Plate_1" / "RunA" / "XY02_CH4_green_centroids.png").exists()
