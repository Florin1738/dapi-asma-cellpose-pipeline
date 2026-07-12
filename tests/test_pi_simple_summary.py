from __future__ import annotations

from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pytest
import tifffile

from dapi_norm.pi_simple_summary import (
    HEADERS,
    PiSummaryRow,
    build_pi_summary,
    find_image_pairs,
    load_count_csv,
    load_counts_by_plate,
    write_pi_plots,
    write_pi_workbook,
)


def test_build_pi_summary_uses_raw_ch2_intensity_and_count_ratio(tmp_path: Path):
    xy01 = tmp_path / "Plate 1" / "RunA" / "XY01"
    xy01.mkdir(parents=True)
    tifffile.imwrite(xy01 / "sample_XY01_CH2.tif", np.array([[1, 2], [3, 4]], dtype=np.uint16))
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", np.array([[5, 0], [0, 5]], dtype=np.uint16))

    counts_csv = tmp_path / "nucleus_counts.csv"
    counts_csv.write_text(
        "image_id,input_path,backend,model_name,channel_id,candidate_stain,"
        "channel_identity_confirmed,nucleus_count,mask_path,qc_montage_path,warnings\n"
        "XY01,input/sample_XY01_CH4.tif,cellpose,cpsam_v2,CH4,candidate_DAPI,True,"
        "2,masks/XY01.tif,qc/XY01.png,\n",
        encoding="utf-8",
    )

    summary = build_pi_summary(input_root=tmp_path, counts_by_plate={"Plate 1": load_count_csv(counts_csv)})

    assert list(summary) == ["Plate 1", "Plate 2"]
    assert summary["Plate 1"][0].location == "XY01"
    assert summary["Plate 1"][0].asma_intensity == 10.0
    assert summary["Plate 1"][0].nuclei_count == 2
    assert summary["Plate 1"][0].ratio == 5.0
    assert summary["Plate 2"] == []


def test_find_image_pairs_accepts_explicit_target_and_dapi_channels(tmp_path: Path):
    xy01 = tmp_path / "Experiment" / "XY01"
    xy01.mkdir(parents=True)
    ch1_path = xy01 / "sample_XY01_CH1.tif"
    ch2_path = xy01 / "sample_XY01_CH2.tif"
    ch4_path = xy01 / "sample_XY01_CH4.tif"
    ch1_path.touch()
    ch2_path.touch()
    ch4_path.touch()

    pairs = find_image_pairs(tmp_path, target_channel="CH1", dapi_channel="CH4")

    assert len(pairs) == 1
    assert pairs[0].ch2_path == ch1_path
    assert pairs[0].ch4_path == ch4_path
    assert pairs[0].target_channel_id == "CH1"
    assert pairs[0].dapi_channel_id == "CH4"


def test_build_pi_summary_falls_back_to_plate_1_when_no_plate_folders(tmp_path: Path):
    xy02 = tmp_path / "Experiment" / "Nested" / "XY02"
    xy02.mkdir(parents=True)
    tifffile.imwrite(xy02 / "sample_XY02_CH2.tif", np.array([[10, 20], [0, 0]], dtype=np.uint16))
    tifffile.imwrite(xy02 / "sample_XY02_CH4.tif", np.array([[1, 1], [0, 0]], dtype=np.uint16))

    summary = build_pi_summary(input_root=tmp_path, counts_by_plate={"Plate 1": {"XY02": 5}})

    assert [row.location for row in summary["Plate 1"]] == ["XY02"]
    assert summary["Plate 1"][0].asma_intensity == 30.0
    assert summary["Plate 1"][0].nuclei_count == 5
    assert summary["Plate 2"] == []


def test_build_pi_summary_matches_duplicate_xy_labels_by_run_source(tmp_path: Path):
    for run_name, ch2_value in [("RunA", 5), ("RunB", 7)]:
        xy01 = tmp_path / "Plate 1" / run_name / "XY01"
        xy01.mkdir(parents=True)
        tifffile.imwrite(
            xy01 / f"{run_name}_XY01_CH2.tif",
            np.full((2, 2), ch2_value, dtype=np.uint16),
        )
        tifffile.imwrite(
            xy01 / f"{run_name}_XY01_CH4.tif",
            np.full((2, 2), 1, dtype=np.uint16),
        )

    summary = build_pi_summary(
        input_root=tmp_path,
        counts_by_plate={"Plate 1": {"RUNA/XY01": 2, "RUNB/XY01": 4}},
    )

    rows = summary["Plate 1"]
    assert [row.location for row in rows] == ["XY01", "XY01"]
    assert [row.source_id for row in rows] == ["RunA/XY01", "RunB/XY01"]
    assert [row.asma_intensity for row in rows] == [20.0, 28.0]
    assert [row.nuclei_count for row in rows] == [2, 4]
    assert [row.ratio for row in rows] == [10.0, 7.0]


def test_build_pi_summary_rejects_bare_counts_for_duplicate_xy_labels(tmp_path: Path):
    for run_name in ["RunA", "RunB"]:
        xy01 = tmp_path / "Plate 1" / run_name / "XY01"
        xy01.mkdir(parents=True)
        tifffile.imwrite(xy01 / f"{run_name}_XY01_CH2.tif", np.ones((2, 2), dtype=np.uint16))
        tifffile.imwrite(xy01 / f"{run_name}_XY01_CH4.tif", np.ones((2, 2), dtype=np.uint16))

    with pytest.raises(ValueError, match="Unprefixed count key XY01 is ambiguous"):
        build_pi_summary(input_root=tmp_path, counts_by_plate={"Plate 1": {"XY01": 99}})


def test_load_counts_by_plate_loads_per_run_count_outputs(tmp_path: Path):
    for run_name, count in [("RunA", 2), ("RunB", 4)]:
        summary_dir = tmp_path / "Plate_1" / run_name / "summaries"
        summary_dir.mkdir(parents=True)
        (summary_dir / "nucleus_counts.csv").write_text(
            "image_id,input_path,backend,model_name,channel_id,candidate_stain,"
            "channel_identity_confirmed,nucleus_count,mask_path,qc_montage_path,warnings\n"
            f"XY01,input/{run_name}_XY01_CH4.tif,cellpose,cpsam_v2,CH4,candidate_DAPI,"
            f"True,{count},masks/XY01.tif,qc/XY01.png,\n",
            encoding="utf-8",
        )

    counts = load_counts_by_plate(tmp_path)

    assert counts["Plate 1"] == {"RUNA/XY01": 2, "RUNB/XY01": 4}


def test_write_pi_workbook_creates_plate_sheets_with_requested_headers_and_ratio_formula(
    tmp_path: Path,
):
    workbook_path = tmp_path / "summary.xlsx"
    summary = {
        "Plate 1": [PiSummaryRow(location="XY01", asma_intensity=100.0, nuclei_count=4)],
        "Plate 2": [],
    }
    write_pi_workbook(workbook_path, summary)

    with ZipFile(workbook_path) as zf:
        workbook_xml = zf.read("xl/workbook.xml").decode("utf-8")
        sheet1_xml = zf.read("xl/worksheets/sheet1.xml").decode("utf-8")
        sheet2_xml = zf.read("xl/worksheets/sheet2.xml").decode("utf-8")

    assert 'name="Plate 1"' in workbook_xml
    assert 'name="Plate 2"' in workbook_xml
    for header in HEADERS:
        assert header in sheet1_xml
        assert header in sheet2_xml
    assert '<f>IF(C2&gt;0,B2/C2,"")</f>' in sheet1_xml
    assert "<v>25</v>" in sheet1_xml


def test_write_pi_plots_creates_summary_plots_for_non_empty_plate(tmp_path: Path):
    summary = {
        "Plate 1": [
            PiSummaryRow(location="XY01", asma_intensity=100.0, nuclei_count=4),
            PiSummaryRow(location="XY02", asma_intensity=60.0, nuclei_count=3),
        ],
        "Plate 2": [],
    }

    created = write_pi_plots(tmp_path, summary)

    assert tmp_path / "Plate_1_aSMA_intensity_by_location.png" in created
    assert tmp_path / "Plate_1_ratio_by_location.png" in created
    assert (tmp_path / "Plate_1_aSMA_intensity_by_location.png").exists()
    assert (tmp_path / "Plate_1_ratio_by_location.png").exists()
    assert not (tmp_path / "Plate_2_aSMA_intensity_by_location.png").exists()
