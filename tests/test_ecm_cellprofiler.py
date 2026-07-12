from __future__ import annotations

from pathlib import Path
from zipfile import ZipFile

import numpy as np
import tifffile

from dapi_norm.ecm_cellprofiler import (
    ECM_ENDPOINT,
    discover_ecm_records,
    measure_ecm_from_mask,
    robust_background_threshold_preview,
    write_ecm_workbook,
)
from scripts.run_cellprofiler_ecm_quantification import resolve_selected_deviations


def test_discover_ecm_records_finds_ch1_ch4_acquisition_pairs(tmp_path: Path):
    xy01 = tmp_path / "Drug 1-3" / "Run A" / "XY01"
    xy01.mkdir(parents=True)
    ch1 = xy01 / "RunA_XY01_CH1.tif"
    ch2 = xy01 / "RunA_XY01_CH2.tif"
    ch4 = xy01 / "RunA_XY01_CH4.tif"
    ch1.touch()
    ch2.touch()
    ch4.touch()

    records = discover_ecm_records(tmp_path)

    assert len(records) == 1
    assert records[0].acquisition_name == "Drug 1-3/Run A"
    assert records[0].location == "XY01"
    assert records[0].ecm_path == ch1.resolve()
    assert records[0].dapi_path == ch4.resolve()
    assert records[0].staged_name.endswith("__XY01_ECM.tif")


def test_measure_ecm_from_mask_uses_mask_negative_background_without_dapi_normalization(
    tmp_path: Path,
):
    ecm_path = tmp_path / "XY01_CH1.tif"
    dapi_path = tmp_path / "XY01_CH4.tif"
    tifffile.imwrite(ecm_path, np.zeros((20, 20), dtype=np.uint16))
    tifffile.imwrite(dapi_path, np.zeros((20, 20), dtype=np.uint16))
    image = np.full((20, 20), 10, dtype=np.uint16)
    image[0:2, 0:2] = 110
    mask = np.zeros((20, 20), dtype=np.uint16)
    mask[0:2, 0:2] = 1
    record = discover_ecm_records(_make_single_record_root(tmp_path))[0]

    measurement = measure_ecm_from_mask(
        record=record,
        image=image,
        mask=mask,
        dapi_count=2,
        threshold_deviations=3,
        ecm_channel_id="CH1",
        dapi_channel_id="CH4",
        mask_path=tmp_path / "mask.tif",
        overlay_path=tmp_path / "overlay.png",
        qc_panel_path=tmp_path / "qc",
        root_for_paths=tmp_path,
    )

    assert measurement.ecm_background_value_per_px == 10.0
    assert measurement.ecm_positive_area_px == 4
    assert measurement.ecm_positive_integrated_raw == 440.0
    assert measurement.ecm_positive_integrated_background_corrected == 400.0
    assert measurement.as_row()[ECM_ENDPOINT] == 400.0


def test_robust_background_threshold_preview_responds_to_deviation_parameter():
    image = np.concatenate([np.full(100, 10), np.full(10, 50)]).astype(np.uint16)

    low = robust_background_threshold_preview(image, deviations=1)
    high = robust_background_threshold_preview(image, deviations=5)

    assert high > low


def test_write_ecm_workbook_contains_method_sheet_and_endpoint(tmp_path: Path):
    workbook = tmp_path / "ecm.xlsx"
    rows = [
        {
            "acquisition": "Run A",
            "location": "XY01",
            "ecm_channel_id": "CH7",
            "dapi_channel_id": "CH9",
            "ecm_positive_integrated_background_corrected": 400.0,
            "dapi_positive_nucleus_count": 2,
            "ecm_positive_area_fraction": 0.01,
            "ecm_background_value_per_px": 10.0,
            "qc_status": "pass",
            "qc_flags": "",
        }
    ]

    write_ecm_workbook(workbook, rows)

    with ZipFile(workbook) as zf:
        workbook_xml = zf.read("xl/workbook.xml").decode("utf-8")
        sheet1 = zf.read("xl/worksheets/sheet1.xml").decode("utf-8")
        sheet2 = zf.read("xl/worksheets/sheet2.xml").decode("utf-8")

    assert 'name="Run A"' in workbook_xml
    assert 'name="Method"' in workbook_xml
    assert "ECM integrated intensity" in sheet1
    assert "DAPI count (context only)" in sheet1
    assert "ECM intensity / DAPI-positive nucleus" not in sheet1
    assert "CH7" in sheet2
    assert "CH9" in sheet2
    assert ECM_ENDPOINT in sheet2
    assert "context/QC only; not used for ECM normalization" in sheet2


def test_resolve_selected_deviations_honors_explicit_override():
    rows = [{"candidate_deviations": 2.0}, {"candidate_deviations": 3.0}]

    selected = resolve_selected_deviations(rows, selected_deviations=2.0)

    assert selected == 2.0


def test_resolve_selected_deviations_rejects_unswept_override():
    rows = [{"candidate_deviations": 2.0}, {"candidate_deviations": 3.0}]

    try:
        resolve_selected_deviations(rows, selected_deviations=4.0)
    except ValueError as exc:
        assert "--selected-deviations 4" in str(exc)
        assert "2, 3" in str(exc)
    else:
        raise AssertionError("Expected invalid selected deviation override to fail")


def _make_single_record_root(tmp_path: Path) -> Path:
    xy01 = tmp_path / "Drug 1-3" / "Run A" / "XY01"
    xy01.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(xy01 / "RunA_XY01_CH1.tif", np.zeros((2, 2), dtype=np.uint16))
    tifffile.imwrite(xy01 / "RunA_XY01_CH4.tif", np.zeros((2, 2), dtype=np.uint16))
    return tmp_path
