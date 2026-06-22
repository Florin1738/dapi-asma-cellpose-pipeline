from __future__ import annotations

from pathlib import Path

import numpy as np

from dapi_norm.nuclei_outputs import summarize_labeled_mask, write_nuclei_count_tables


def test_summarize_labeled_mask_returns_count_and_centroids():
    labels = np.zeros((6, 8), dtype=np.uint16)
    labels[1:3, 2:4] = 1
    labels[3:6, 5:8] = 2

    summary, rows = summarize_labeled_mask(
        image_id="XY01",
        input_path=Path("input/XY01_CH4.tif"),
        mask=labels,
        backend="cellpose",
        model_name="cpsam_v2",
        channel_id="CH4",
        candidate_stain="candidate_DAPI",
    )

    assert summary["image_id"] == "XY01"
    assert summary["nucleus_count"] == 2
    assert summary["backend"] == "cellpose"
    assert summary["model_name"] == "cpsam_v2"
    assert summary["channel_id"] == "CH4"
    assert summary["candidate_stain"] == "candidate_DAPI"
    assert rows[0]["nucleus_id"] == 1
    assert rows[0]["area_px"] == 4
    assert rows[0]["x_centroid"] == 2.5
    assert rows[0]["y_centroid"] == 1.5
    assert rows[0]["bbox_min_row"] == 1
    assert rows[0]["bbox_min_col"] == 2
    assert rows[0]["bbox_max_row"] == 3
    assert rows[0]["bbox_max_col"] == 4
    assert rows[0]["touches_border"] is False
    assert rows[0]["kept_after_filtering"] is True
    assert rows[1]["touches_border"] is True


def test_write_nuclei_count_tables_uses_stable_columns(tmp_path: Path):
    summary_rows = [
        {
            "image_id": "XY01",
            "input_path": "input/XY01_CH4.tif",
            "backend": "cellpose",
            "model_name": "cpsam_v2",
            "channel_id": "CH4",
            "candidate_stain": "candidate_DAPI",
            "channel_identity_confirmed": False,
            "nucleus_count": 2,
            "mask_path": "masks/XY01.tif",
            "qc_montage_path": "qc/XY01.png",
            "warnings": "channel_identity_unconfirmed",
        }
    ]
    nucleus_rows = [
        {
            "image_id": "XY01",
            "input_path": "input/XY01_CH4.tif",
            "nucleus_id": 1,
            "x_centroid": 2.5,
            "y_centroid": 1.5,
            "area_px": 4,
            "bbox_min_row": 1,
            "bbox_min_col": 2,
            "bbox_max_row": 3,
            "bbox_max_col": 4,
            "touches_border": False,
            "kept_after_filtering": True,
        }
    ]

    write_nuclei_count_tables(
        output_dir=tmp_path,
        summary_rows=summary_rows,
        nucleus_rows=nucleus_rows,
    )

    summary_csv = (tmp_path / "nucleus_counts.csv").read_text(encoding="utf-8")
    per_nucleus_csv = (tmp_path / "per_nucleus_locations.csv").read_text(encoding="utf-8")
    assert summary_csv.startswith(
        "image_id,input_path,backend,model_name,channel_id,candidate_stain,"
        "channel_identity_confirmed,nucleus_count,mask_path,qc_montage_path,warnings\n"
    )
    assert "XY01,input/XY01_CH4.tif,cellpose,cpsam_v2,CH4,candidate_DAPI,False,2" in (
        summary_csv
    )
    assert per_nucleus_csv.startswith(
        "image_id,input_path,nucleus_id,x_centroid,y_centroid,area_px,"
        "bbox_min_row,bbox_min_col,bbox_max_row,bbox_max_col,touches_border,"
        "kept_after_filtering\n"
    )
