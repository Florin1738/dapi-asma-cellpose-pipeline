from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import tifffile
import yaml

from dapi_norm.cellpose_plate_summary import (
    CELLPOSE_PLATE_SUMMARY_COLUMNS,
    build_cellpose_plate_summary,
    write_cellpose_plate_summary_csv,
    write_cellpose_plate_summary_markdown,
)


def test_build_cellpose_plate_summary_merges_runs_and_preserves_source_id(tmp_path: Path):
    runs_root = tmp_path / "runs"
    _write_fake_run(
        runs_root / "Plate_1" / "RunA",
        plate="Plate 1",
        run_name="RunA",
        image_id="XY01",
        ch2=np.array([[1, 2], [3, 4]], dtype=np.uint16),
        nuclei_count=2,
        candidate_area=3,
        candidate_raw=9,
        candidate_per_nucleus=4.5,
        anchored_area=2,
        anchored_raw=7,
        anchored_background_corrected=6,
        anchored_per_nucleus=3.5,
        anchored_background_corrected_per_nucleus=3.0,
        no_dapi_objects=1,
    )
    _write_fake_run(
        runs_root / "Plate_1" / "RunB",
        plate="Plate 1",
        run_name="RunB",
        image_id="XY01",
        ch2=np.array([[10, 20], [30, 40]], dtype=np.uint16),
        nuclei_count=4,
        candidate_area=4,
        candidate_raw=100,
        candidate_per_nucleus=25,
        anchored_area=3,
        anchored_raw=70,
        anchored_background_corrected=60,
        anchored_per_nucleus=17.5,
        anchored_background_corrected_per_nucleus=15.0,
        no_dapi_objects=2,
    )
    _write_fake_run(
        runs_root / "Plate_2" / "RunC",
        plate="Plate 2",
        run_name="RunC",
        image_id="XY02",
        ch2=np.array([[5, 5], [5, 5]], dtype=np.uint16),
        nuclei_count=1,
        candidate_area=1,
        candidate_raw=5,
        candidate_per_nucleus=5,
        anchored_area=1,
        anchored_raw=5,
        anchored_background_corrected=4,
        anchored_per_nucleus=5,
        anchored_background_corrected_per_nucleus=4.0,
        no_dapi_objects=0,
    )

    rows = build_cellpose_plate_summary(runs_root)

    assert [row["plate"] for row in rows] == ["Plate 1", "Plate 1", "Plate 2"]
    assert [row["source_id"] for row in rows] == ["RunA/XY01", "RunB/XY01", "RunC/XY02"]
    assert [row["location"] for row in rows] == ["XY01", "XY01", "XY02"]
    assert rows[0]["target_channel_id"] == "CH2"
    assert rows[0]["dapi_channel_id"] == "CH4"
    assert rows[0]["whole_field_target_integrated_raw"] == 10.0
    assert rows[0]["whole_field_target_integrated_raw_per_DAPI_positive_nucleus"] == 5.0
    assert rows[0]["whole_field_ch2_integrated_raw"] == 10.0
    assert rows[0]["whole_field_ch2_integrated_raw_per_DAPI_positive_nucleus"] == 5.0
    assert rows[1]["whole_field_ch2_integrated_raw"] == 100.0
    assert rows[0]["cellpose_masked_ch2_integrated_raw"] == 9.0
    assert rows[0]["cellpose_masked_target_integrated_raw"] == 9.0
    assert rows[0]["cellpose_masked_ch2_integrated_raw_per_DAPI_positive_nucleus"] == 4.5
    assert rows[0]["cellpose_masked_target_integrated_raw_per_DAPI_positive_nucleus"] == 4.5
    assert rows[0]["cellpose_masked_ch2_integrated_background_corrected"] == 8.0
    assert rows[0]["cellpose_masked_target_integrated_background_corrected"] == 8.0
    assert (
        rows[0]["cellpose_masked_ch2_integrated_background_corrected_per_DAPI_positive_nucleus"]
        == 4.0
    )
    assert (
        rows[0]["cellpose_masked_target_integrated_background_corrected_per_DAPI_positive_nucleus"]
        == 4.0
    )
    assert rows[0]["dapi_anchored_cellpose_ch2_integrated_raw"] == 7.0
    assert rows[0]["dapi_anchored_cellpose_target_integrated_raw"] == 7.0
    assert rows[0]["dapi_anchored_cellpose_ch2_integrated_raw_per_DAPI_positive_nucleus"] == 3.5
    assert rows[0]["dapi_anchored_cellpose_target_integrated_raw_per_DAPI_positive_nucleus"] == 3.5
    assert rows[0]["dapi_anchored_cellpose_ch2_integrated_background_corrected"] == 6.0
    assert rows[0]["dapi_anchored_cellpose_target_integrated_background_corrected"] == 6.0
    assert (
        rows[0]["dapi_anchored_cellpose_ch2_integrated_background_corrected_per_DAPI_positive_nucleus"]
        == 3.0
    )
    assert (
        rows[0][
            "dapi_anchored_cellpose_target_integrated_background_corrected_per_DAPI_positive_nucleus"
        ]
        == 3.0
    )
    assert rows[1]["no_dapi_cellpose_object_count_excluded_in_anchored_variant"] == 2
    assert rows[0]["qc_status"] == "needs_manual_review"
    assert rows[0]["qc_flags"] == "synthetic_flag"
    assert rows[0]["source_warnings"] == "synthetic_warning"
    assert rows[0]["source_qc_panel_path"].endswith("qc/XY01.png")
    assert rows[0]["source_excluded_signal_check_path"] == ""
    assert rows[0]["target_path"].endswith("XY01_CH2.tif")
    assert rows[0]["dapi_path"].endswith("XY01_CH4.tif")
    assert rows[0]["ch4_path"].endswith("XY01_CH4.tif")


def test_write_cellpose_plate_summary_outputs_csv_and_markdown(tmp_path: Path):
    rows = [
        {
            "plate": "Plate 1",
            "source_id": "RunA/XY01",
            "location": "XY01",
            "whole_field_ch2_integrated_raw": 10.0,
            "whole_field_ch2_integrated_raw_per_DAPI_positive_nucleus": 5.0,
            "dapi_positive_nucleus_count": 2,
            "cellpose_masked_ch2_integrated_raw": 9.0,
            "cellpose_masked_ch2_integrated_raw_per_DAPI_positive_nucleus": 4.5,
            "cellpose_masked_ch2_integrated_background_corrected": 8.0,
            "cellpose_masked_ch2_integrated_background_corrected_per_DAPI_positive_nucleus": 4.0,
            "dapi_anchored_cellpose_ch2_integrated_raw": 7.0,
            "dapi_anchored_cellpose_ch2_integrated_raw_per_DAPI_positive_nucleus": 3.5,
            "dapi_anchored_cellpose_ch2_integrated_background_corrected": 6.0,
            "dapi_anchored_cellpose_ch2_integrated_background_corrected_per_DAPI_positive_nucleus": 3.0,
            "cellpose_masked_area_px": 3,
            "dapi_anchored_cellpose_masked_area_px": 2,
            "cellpose_masked_area_per_DAPI_positive_nucleus": 1.5,
            "dapi_anchored_cellpose_masked_area_per_DAPI_positive_nucleus": 1.0,
            "no_dapi_cellpose_object_count_excluded_in_anchored_variant": 1,
            "cellpose_object_count": 2,
            "qc_status": "needs_manual_review",
            "qc_flags": "synthetic_flag",
            "source_warnings": "synthetic_warning",
            "cellpose_mask_path": "masks/XY01.tif",
            "dapi_nuclei_mask_path": "counts/XY01.tif",
            "source_qc_panel_path": "qc/XY01.png",
            "source_excluded_signal_check_path": "",
            "ch2_path": "data/XY01_CH2.tif",
            "ch4_path": "data/XY01_CH4.tif",
        }
    ]

    csv_path = write_cellpose_plate_summary_csv(tmp_path / "summary.csv", rows)
    md_path = write_cellpose_plate_summary_markdown(tmp_path / "summary.md", rows)

    with csv_path.open(newline="", encoding="utf-8") as handle:
        csv_rows = list(csv.DictReader(handle))
    assert csv_rows[0]["source_id"] == "RunA/XY01"
    assert list(csv_rows[0].keys()) == CELLPOSE_PLATE_SUMMARY_COLUMNS
    text = md_path.read_text(encoding="utf-8")
    assert "| Plate 1 | RunA/XY01 | XY01 |" in text
    assert "dapi_anchored_cellpose_ch2_integrated_raw_per_DAPI_positive_nucleus" in text


def _write_fake_run(
    run_dir: Path,
    *,
    plate: str,
    run_name: str,
    image_id: str,
    ch2: np.ndarray,
    nuclei_count: int,
    candidate_area: int,
    candidate_raw: float,
    candidate_per_nucleus: float,
    anchored_area: int,
    anchored_raw: float,
    anchored_background_corrected: float,
    anchored_per_nucleus: float,
    anchored_background_corrected_per_nucleus: float,
    no_dapi_objects: int,
) -> None:
    summaries = run_dir / "summaries"
    logs = run_dir / "logs"
    image_dir = run_dir / "images"
    summaries.mkdir(parents=True)
    logs.mkdir()
    image_dir.mkdir()
    ch2_path = image_dir / f"{image_id}_CH2.tif"
    ch4_path = image_dir / f"{image_id}_CH4.tif"
    tifffile.imwrite(ch2_path, ch2)
    tifffile.imwrite(ch4_path, np.where(ch2 > 0, 100, 0).astype(np.uint16))
    metrics_path = summaries / "cellpose_cell_region_image_metrics.csv"
    metrics_path.write_text(
        "image_id,source_id,dapi_positive_nucleus_count,cellpose_object_count,"
        "candidate_region_area_px,target_integrated_raw_in_cellpose_region,"
        "target_integrated_background_corrected_in_cellpose_region,"
        "target_integrated_intensity_per_DAPI_positive_nucleus,"
        "dapi_anchored_candidate_region_area_px,dapi_anchored_positive_area_per_DAPI_positive_nucleus,"
        "dapi_anchored_target_integrated_raw,dapi_anchored_target_integrated_background_corrected,"
        "dapi_anchored_target_integrated_intensity_per_DAPI_positive_nucleus,"
        "dapi_anchored_excluded_no_dapi_object_count,qc_status,qc_flags,warnings,"
        "mask_path,nuclei_mask_path,qc_panel_path,excluded_signal_check_path\n"
        f"{image_id},{image_id},{nuclei_count},{no_dapi_objects + 1},{candidate_area},"
        f"{candidate_raw},{candidate_raw - 1},{candidate_per_nucleus},{anchored_area},{anchored_area / nuclei_count},"
        f"{anchored_raw},{anchored_background_corrected},{anchored_background_corrected_per_nucleus},"
        f"{no_dapi_objects},needs_manual_review,synthetic_flag,synthetic_warning,"
        f"masks/{image_id}.tif,counts/{image_id}.tif,qc/{image_id}.png,\n",
        encoding="utf-8",
    )
    config = {
        "image_inputs": [
            {
                "source_id": image_id,
                "location": image_id,
                "ch2_path": str(ch2_path),
                "ch4_path": str(ch4_path),
                "nuclei_mask_path": f"counts/{image_id}.tif",
                "cellpose_mask_path": f"masks/{image_id}.tif",
            }
        ],
        "plate": plate,
        "run_name": run_name,
    }
    (logs / "config_resolved.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
