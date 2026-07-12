from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import tifffile

from dapi_norm.cellpose_endpoint_figures import (
    render_cellpose_endpoint_figures,
    render_dapi_nuclei_qc_pages,
)


def test_render_cellpose_endpoint_figures_writes_endpoint_plots_and_overlay_pages(tmp_path: Path):
    summary_path = tmp_path / "summary.csv"
    rows = []
    for index in range(1, 5):
        rows.append(_write_record(tmp_path, image_id=f"XY{index:02d}", index=index))
    _write_csv(summary_path, rows)

    outputs = render_cellpose_endpoint_figures(
        summary_csv=summary_path,
        output_dir=tmp_path / "figures",
        panel_page_size=3,
    )

    assert outputs["metric_contrast"].exists()
    assert outputs["masking_effect"].exists()
    assert outputs["plate_summary"].exists()
    assert outputs["representative_cell_segmentation"].exists()
    assert len(outputs["overlay_pages"]) == 2
    assert outputs["overlay_pages"][0].name.endswith("page_01.png")
    assert outputs["overlay_index"].exists()
    with outputs["overlay_index"].open(newline="", encoding="utf-8") as handle:
        index_rows = list(csv.DictReader(handle))
    assert len(index_rows) == 4
    assert index_rows[0]["source_id"].startswith("RunA/XY")
    assert index_rows[0]["page"].endswith("page_01.png")
    assert outputs["captions_markdown"].exists()
    captions = outputs["captions_markdown"].read_text(encoding="utf-8")
    assert "Figure 1" in captions
    assert "DAPI-positive nucleus" in captions


def test_render_cellpose_endpoint_figures_labels_non_default_target_channel(tmp_path: Path):
    summary_path = tmp_path / "summary.csv"
    rows = []
    for index in range(1, 3):
        row = _write_record(tmp_path, image_id=f"XY{index:02d}", index=index)
        row["target_channel_id"] = "CH1"
        row["dapi_channel_id"] = "CH4"
        row["target_path"] = row["ch2_path"]
        row["dapi_path"] = row["ch4_path"]
        row["whole_field_target_integrated_raw"] = row["whole_field_ch2_integrated_raw"]
        row["whole_field_target_integrated_raw_per_DAPI_positive_nucleus"] = row[
            "whole_field_ch2_integrated_raw_per_DAPI_positive_nucleus"
        ]
        row["cellpose_masked_target_integrated_raw"] = row["cellpose_masked_ch2_integrated_raw"]
        row["cellpose_masked_target_integrated_raw_per_DAPI_positive_nucleus"] = row[
            "cellpose_masked_ch2_integrated_raw_per_DAPI_positive_nucleus"
        ]
        row["cellpose_masked_target_integrated_background_corrected"] = row[
            "cellpose_masked_ch2_integrated_background_corrected"
        ]
        row["cellpose_masked_target_integrated_background_corrected_per_DAPI_positive_nucleus"] = row[
            "cellpose_masked_ch2_integrated_background_corrected_per_DAPI_positive_nucleus"
        ]
        row["dapi_anchored_cellpose_target_integrated_raw"] = row[
            "dapi_anchored_cellpose_ch2_integrated_raw"
        ]
        row["dapi_anchored_cellpose_target_integrated_raw_per_DAPI_positive_nucleus"] = row[
            "dapi_anchored_cellpose_ch2_integrated_raw_per_DAPI_positive_nucleus"
        ]
        row["dapi_anchored_cellpose_target_integrated_background_corrected"] = row[
            "dapi_anchored_cellpose_ch2_integrated_background_corrected"
        ]
        row[
            "dapi_anchored_cellpose_target_integrated_background_corrected_per_DAPI_positive_nucleus"
        ] = row["dapi_anchored_cellpose_ch2_integrated_background_corrected_per_DAPI_positive_nucleus"]
        rows.append(row)
    _write_csv(summary_path, rows)

    outputs = render_cellpose_endpoint_figures(
        summary_csv=summary_path,
        output_dir=tmp_path / "figures",
        panel_page_size=2,
    )

    captions = outputs["captions_markdown"].read_text(encoding="utf-8")
    assert "CH1/aSMA" in captions
    assert "raw CH2/aSMA" not in captions


def test_render_dapi_nuclei_qc_pages_writes_plate_split_pages_and_index(tmp_path: Path):
    summary_path = tmp_path / "summary.csv"
    rows = []
    for index in range(1, 5):
        rows.append(_write_record(tmp_path, image_id=f"XY{index:02d}", index=index))
    _write_csv(summary_path, rows)

    outputs = render_dapi_nuclei_qc_pages(
        summary_csv=summary_path,
        output_dir=tmp_path / "dapi_pages",
        page_size=2,
    )

    assert outputs["field_count"] == 4
    assert outputs["plate_count"] == 2
    assert len(outputs["pages"]) == 2
    assert (tmp_path / "dapi_pages" / "Plate_1" / "dapi_nuclei_overlay_page_01.png").exists()
    assert (tmp_path / "dapi_pages" / "Plate_2" / "dapi_nuclei_overlay_page_01.png").exists()
    assert outputs["index"].exists()
    with outputs["index"].open(newline="", encoding="utf-8") as handle:
        index_rows = list(csv.DictReader(handle))
    assert len(index_rows) == 4
    assert index_rows[0]["dapi_positive_nucleus_count"] == "2"


def _write_record(tmp_path: Path, *, image_id: str, index: int) -> dict[str, str]:
    image_dir = tmp_path / "images"
    image_dir.mkdir(exist_ok=True)
    ch2 = np.zeros((30, 40), dtype=np.uint16)
    ch2[4:20, 5:25] = 1000 * index
    ch4 = np.zeros((30, 40), dtype=np.uint16)
    ch4[8:10, 8:10] = 5000
    ch4[15:17, 20:22] = 5000
    nuclei = np.zeros((30, 40), dtype=np.uint16)
    nuclei[8:10, 8:10] = 1
    nuclei[15:17, 20:22] = 2
    labels = np.zeros((30, 40), dtype=np.uint16)
    labels[5:18, 6:24] = 1
    labels[21:26, 28:35] = 2
    ch2_path = image_dir / f"{image_id}_CH2.tif"
    ch4_path = image_dir / f"{image_id}_CH4.tif"
    nuclei_path = image_dir / f"{image_id}_nuclei.tif"
    labels_path = image_dir / f"{image_id}_cellpose.tif"
    tifffile.imwrite(ch2_path, ch2)
    tifffile.imwrite(ch4_path, ch4)
    tifffile.imwrite(nuclei_path, nuclei)
    tifffile.imwrite(labels_path, labels)
    nuclei_count = 2
    whole_raw = float(np.sum(ch2))
    cellpose_raw = float(np.sum(ch2[labels > 0]))
    anchored_raw = float(np.sum(ch2[labels == 1]))
    return {
        "plate": "Plate 1" if index <= 2 else "Plate 2",
        "source_id": f"RunA/{image_id}",
        "location": image_id,
        "whole_field_ch2_integrated_raw": str(whole_raw),
        "whole_field_ch2_integrated_raw_per_DAPI_positive_nucleus": str(whole_raw / nuclei_count),
        "dapi_positive_nucleus_count": str(nuclei_count),
        "cellpose_masked_ch2_integrated_raw": str(cellpose_raw),
        "cellpose_masked_ch2_integrated_raw_per_DAPI_positive_nucleus": str(
            cellpose_raw / nuclei_count
        ),
        "cellpose_masked_ch2_integrated_background_corrected": str(cellpose_raw),
        "cellpose_masked_ch2_integrated_background_corrected_per_DAPI_positive_nucleus": str(
            cellpose_raw / nuclei_count
        ),
        "dapi_anchored_cellpose_ch2_integrated_raw": str(anchored_raw),
        "dapi_anchored_cellpose_ch2_integrated_raw_per_DAPI_positive_nucleus": str(
            anchored_raw / nuclei_count
        ),
        "dapi_anchored_cellpose_ch2_integrated_background_corrected": str(anchored_raw),
        "dapi_anchored_cellpose_ch2_integrated_background_corrected_per_DAPI_positive_nucleus": str(
            anchored_raw / nuclei_count
        ),
        "cellpose_masked_area_px": str(int(np.count_nonzero(labels))),
        "dapi_anchored_cellpose_masked_area_px": str(int(np.count_nonzero(labels == 1))),
        "cellpose_masked_area_per_DAPI_positive_nucleus": "134.5",
        "dapi_anchored_cellpose_masked_area_per_DAPI_positive_nucleus": "117.0",
        "no_dapi_cellpose_object_count_excluded_in_anchored_variant": "1",
        "cellpose_object_count": "2",
        "cellpose_mask_path": str(labels_path),
        "dapi_nuclei_mask_path": str(nuclei_path),
        "ch2_path": str(ch2_path),
        "ch4_path": str(ch4_path),
    }


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
