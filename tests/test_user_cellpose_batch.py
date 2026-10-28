from __future__ import annotations

import csv
from pathlib import Path
from zipfile import ZipFile

import pytest
import yaml

from dapi_norm.user_cellpose_batch import (
    PRIMARY_ENDPOINT,
    _acquisition_log_record,
    _pairs_for_acquisition,
    _prepare_output_root,
    discover_acquisitions,
    user_facing_workbook_filename,
    workbook_filename_for_assay_folder,
    write_pi_style_workbook_from_rows,
    write_user_summary_csv,
)


def test_discover_acquisitions_finds_leaf_acquisition_folders(tmp_path: Path):
    root = tmp_path / "plates"
    _touch_pair(root / "plate 1" / "Run A" / "XY01", "RunA_XY01")
    _touch_pair(root / "plate 1" / "Run A" / "reports" / "XY99", "RunA_XY99")
    _touch_pair(root / "plate 1" / "Run B" / "XY01", "RunB_XY01")
    _touch_pair(root / "plate 2" / "Run C" / "XY05", "RunC_XY05")
    _touch_pair(root / "output" / "Ignored" / "XY99", "Ignored_XY99")

    acquisitions = discover_acquisitions(root)

    assert [acquisition.output_key for acquisition in acquisitions] == [
        "Plate_1/Run_A",
        "Plate_1/Run_B",
        "Plate_2/Run_C",
    ]
    assert [acquisition.image_count for acquisition in acquisitions] == [1, 1, 1]
    assert acquisitions[0].display_name == "plate 1/Run A"
    assert [pair.source_id for pair in _pairs_for_acquisition(root / "plate 1" / "Run A", max_images=None)] == [
        "XY01"
    ]


def test_discover_acquisitions_accepts_single_acquisition_root(tmp_path: Path):
    root = tmp_path / "Acquisition One"
    _touch_pair(root / "XY02", "Sample_XY02")

    acquisitions = discover_acquisitions(root)

    assert len(acquisitions) == 1
    assert acquisitions[0].plate_name == "Plate 1"
    assert acquisitions[0].run_name == "Acquisition_One"
    assert acquisitions[0].display_name == "Acquisition One"
    assert acquisitions[0].image_count == 1


def test_discover_acquisitions_accepts_explicit_channel_mapping(tmp_path: Path):
    root = tmp_path / "plates"
    xy01 = root / "drug group" / "Run A" / "XY01"
    xy01.mkdir(parents=True)
    (xy01 / "RunA_XY01_CH1.tif").touch()
    (xy01 / "RunA_XY01_CH2.tif").touch()
    (xy01 / "RunA_XY01_CH4.tif").touch()

    acquisitions = discover_acquisitions(root, target_channel_id="CH2", dapi_channel_id="CH4")
    pairs = _pairs_for_acquisition(
        root / "drug group" / "Run A",
        max_images=None,
        target_channel_id="CH2",
        dapi_channel_id="CH4",
    )

    assert [acquisition.output_key for acquisition in acquisitions] == ["Plate_1/Run_A"]
    assert acquisitions[0].image_count == 1
    assert pairs[0].target_channel_id == "CH2"
    assert pairs[0].dapi_channel_id == "CH4"
    assert pairs[0].ch2_path.name.endswith("_CH2.tif")
    assert pairs[0].ch4_path.name.endswith("_CH4.tif")


def test_user_facing_workbook_filename_uses_assay_folder_name(tmp_path: Path):
    assay_root = tmp_path / "003-7b TW aSMA"
    _touch_pair(assay_root / "XY01", "RunA_XY01")

    acquisition = discover_acquisitions(assay_root)[0]

    assert user_facing_workbook_filename([acquisition]) == (
        "aSMA Quantification Workbook - 003-7b TW aSMA.xlsx"
    )


def test_user_facing_workbook_filename_marks_combined_runs_and_sanitizes_invalid_characters(
    tmp_path: Path,
):
    _touch_pair(tmp_path / "Run A" / "XY01", "RunA_XY01")
    _touch_pair(tmp_path / "Run B" / "XY01", "RunB_XY01")

    acquisitions = discover_acquisitions(tmp_path)

    assert user_facing_workbook_filename(acquisitions) == (
        "aSMA Quantification Workbook - Combined 2 Assays.xlsx"
    )
    assert workbook_filename_for_assay_folder('assay:01/CH2') == (
        "aSMA Quantification Workbook - assay_01_CH2.xlsx"
    )
    with pytest.raises(ValueError, match="valid workbook filename"):
        workbook_filename_for_assay_folder(". ")


def test_acquisition_log_record_is_yaml_serializable(tmp_path: Path):
    root = tmp_path / "Acquisition One"
    _touch_pair(root / "XY02", "Sample_XY02")
    acquisition = discover_acquisitions(root)[0]

    dumped = yaml.safe_dump({"acquisitions": [_acquisition_log_record(acquisition)]})

    assert "Acquisition One" in dumped


def test_user_summary_csv_exposes_scientific_endpoint(tmp_path: Path):
    rows = [_fake_summary_row()]

    csv_path = write_user_summary_csv(tmp_path / "summary.csv", rows)

    with csv_path.open(newline="", encoding="utf-8") as handle:
        exported = list(csv.DictReader(handle))
    assert exported[0]["target_integrated_intensity"] == "250.0"
    assert exported[0]["normalization_denominator_count"] == "5"
    assert exported[0][PRIMARY_ENDPOINT] == "50.0"
    assert exported[0]["dapi_anchored_cellpose_object_count"] == "3"
    assert exported[0]["target_channel_id"] == "CH2"
    assert exported[0]["dapi_channel_id"] == "CH4"
    assert exported[0]["target_path"] == "XY01_CH2.tif"
    assert exported[0]["dapi_path"] == "XY01_CH4.tif"
    assert "cellpose_object_count" not in exported[0]
    assert "DAPI-normalized" not in csv_path.read_text(encoding="utf-8")


def test_prepare_output_root_overwrite_clears_managed_outputs_only(tmp_path: Path):
    output_root = tmp_path / "run"
    stale_regions = output_root / "cellpose_cell_regions" / "Plate_9" / "Stale"
    stale_regions.mkdir(parents=True)
    (stale_regions / "old.csv").write_text("stale", encoding="utf-8")
    keep_file = output_root / "operator_note.txt"
    keep_file.write_text("keep", encoding="utf-8")

    _prepare_output_root(output_root, overwrite=True)

    assert not (output_root / "cellpose_cell_regions").exists()
    assert keep_file.read_text(encoding="utf-8") == "keep"


def test_pi_style_workbook_from_rows_is_standalone_xlsx(tmp_path: Path):
    plate_3_row = _fake_summary_row()
    plate_3_row["plate"] = "Plate 3"
    plate_3_row["source_id"] = "RunC/XY03"
    plate_3_row["field_id"] = "XY03"
    plate_3_row["well_marker"] = "XY03/_A01"
    workbook_path = write_pi_style_workbook_from_rows(
        tmp_path / "summary.xlsx",
        [_fake_summary_row(), plate_3_row],
    )

    with ZipFile(workbook_path) as zf:
        names = set(zf.namelist())
        assert "xl/workbook.xml" in names
        assert "xl/worksheets/sheet1.xml" in names
        assert "xl/worksheets/sheet2.xml" in names
        workbook_xml = zf.read("xl/workbook.xml").decode("utf-8")
        sheet_xml = zf.read("xl/worksheets/sheet1.xml").decode("utf-8")
        sheet2_xml = zf.read("xl/worksheets/sheet2.xml").decode("utf-8")
        metadata_xml = zf.read("xl/worksheets/sheet3.xml").decode("utf-8")
    assert "Plate 3" in workbook_xml
    assert "Channel Mapping" in workbook_xml
    assert "RunA/XY01" in sheet_xml
    assert "RunC/XY03" in sheet2_xml
    assert "target_channel_id" in metadata_xml
    assert "CH2" in metadata_xml
    assert "dapi_channel_id" in metadata_xml
    assert "CH4" in metadata_xml
    assert 'IF(C2&gt;0,B2/C2,"")' in sheet_xml


def test_pi_style_workbook_expands_location_column_for_long_source_ids(tmp_path: Path):
    long_location_row = _fake_summary_row()
    long_location_row["source_id"] = "0087cTW_01/XY82"
    long_location_row["field_id"] = "XY82"
    long_location_row["well_marker"] = "XY82/_A01"
    workbook_path = write_pi_style_workbook_from_rows(
        tmp_path / "summary.xlsx",
        [long_location_row],
    )

    with ZipFile(workbook_path) as zf:
        sheet_xml = zf.read("xl/worksheets/sheet1.xml").decode("utf-8")

    assert '<col min="1" max="1" width="17" customWidth="1"/>' in sheet_xml


def _touch_pair(xy_dir: Path, stem: str) -> None:
    xy_dir.mkdir(parents=True)
    (xy_dir / f"{stem}_CH2.tif").touch()
    (xy_dir / f"{stem}_CH4.tif").touch()


def _fake_summary_row() -> dict[str, object]:
    return {
        "well_id": "A01",
        "field_in_well": 1,
        "fields_in_well": 1,
        "treatment_mapping_status": "not_provided",
        "assay_folder": "RunA",
        "field_id": "XY01",
        "well_marker": "XY01/_A01",
        "well_marker_sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        "plate_map_source": "", "plate_map_sha256": "",
        "compound": "", "batch": "", "dose": "", "condition": "",
        "tgf_beta": "", "serum": "", "cells_seeded": "",
        "plate": "Plate 1",
        "location": "XY01",
        "source_id": "RunA/XY01",
        "dapi_anchored_cellpose_ch2_integrated_background_corrected": 250.0,
        "dapi_positive_nucleus_count": 5,
        "dapi_anchored_cellpose_ch2_integrated_background_corrected_per_DAPI_positive_nucleus": 50.0,
        "qc_status": "reviewable_not_validated",
        "qc_flags": "",
        "cellpose_object_count": 4,
        "no_dapi_cellpose_object_count_excluded_in_anchored_variant": 1,
        "dapi_anchored_cellpose_masked_area_px": 120,
        "cellpose_mask_path": "masks/XY01.tif",
        "dapi_nuclei_mask_path": "counts/XY01.tif",
        "ch2_path": "XY01_CH2.tif",
        "ch4_path": "XY01_CH4.tif",
    }
