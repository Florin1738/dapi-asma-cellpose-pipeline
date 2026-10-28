from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock
import csv
import pytest

from dapi_norm.well_mapping import (
    TREATMENT_COLUMNS, load_acquisition_mapping, mapping_for_acquisitions,
    attach_mapping, write_mapping_artifacts, validate_mapping_rows,
)
from dapi_norm.user_cellpose_batch import run_user_cellpose_batch, write_pi_style_workbook_from_rows
from dapi_norm.parallel_user_cellpose_batch import build_per_acquisition_runs


def marker(root, field, well):
    folder = root / field
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"_{well}").touch()
    for channel in ("CH2", "CH4"):
        (folder / f"sample_{field}_{channel}.tif").touch()


def plate_map(root, wells, **changes):
    path = root / "plate_map.csv"
    rows = [dict(assay_folder=root.name, well_id=well,
                 **{col: "n/a" for col in TREATMENT_COLUMNS}) for well in wells]
    for row in rows:
        row.update(changes)
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=("assay_folder", "well_id", *TREATMENT_COLUMNS))
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_irregular_scan_and_repeats_are_preserved(tmp_path):
    marker(tmp_path, "XY01", "G11")
    marker(tmp_path, "XY02", "B02")
    marker(tmp_path, "XY99", "B02")
    rows = load_acquisition_mapping(tmp_path, ["XY99", "XY02", "XY01"])
    assert rows["XY01"]["well_id"] == "G11"
    assert rows["XY02"]["fields_in_well"] == 2
    assert rows["XY02"]["field_in_well"] == 1
    assert rows["XY99"]["field_in_well"] == 2
    assert rows["XY01"]["treatment_mapping_status"] == "not_provided"
    write_mapping_artifacts(tmp_path / "export", {str(tmp_path): rows})
    assert (tmp_path / "export/image_to_well.csv").is_file()
    assert '"H12"' in (tmp_path / "export/mapping_validation.json").read_text()


@pytest.mark.parametrize("mode", ["missing", "conflict", "duplicate_alias", "invalid"])
def test_bad_marker_rejected(tmp_path, mode):
    marker(tmp_path, "XY01", "A01")
    if mode == "missing":
        (tmp_path / "XY01/_A01").unlink()
    elif mode == "conflict":
        (tmp_path / "XY01/_B02").touch()
    elif mode == "duplicate_alias":
        (tmp_path / "XY01/_A1").touch()
    else:
        (tmp_path / "XY01/_A01").rename(tmp_path / "XY01/_I01")
    with pytest.raises(ValueError, match="exactly one well marker"):
        load_acquisition_mapping(tmp_path, ["XY01"])


def test_duplicate_field_rejected(tmp_path):
    marker(tmp_path, "XY01", "A01")
    with pytest.raises(ValueError, match="Duplicate image field"):
        load_acquisition_mapping(tmp_path, ["XY01", "XY01"])


def test_valid_plate_map_is_joined_by_well_not_row(tmp_path):
    marker(tmp_path, "XY01", "H12")
    marker(tmp_path, "XY02", "B02")
    path = plate_map(tmp_path, ["B02", "H12"], compound="test drug", dose="5 uM")
    records = load_acquisition_mapping(tmp_path, ["XY01", "XY02"], path)
    assert records["XY01"]["compound"] == "test drug"
    assert records["XY01"]["treatment_mapping_status"] == "provided"
    assert len(records["XY01"]["plate_map_sha256"]) == 64
    rows = [dict(target_path=str(tmp_path / "XY01/sample_XY01_CH2.tif"))]
    attach_mapping(rows, {str(tmp_path.resolve()): records})
    assert rows[0]["well_id"] == "H12"


@pytest.mark.parametrize("wells,changes,message", [
    (["A01", "A01"], {}, "Duplicate plate-map"),
    (["A02"], {}, "Plate map missing"),
    (["A01"], {"dose": ""}, "Incomplete treatment"),
    (["A01"], {"assay_folder": "wrong"}, "Plate map missing"),
    (["A13"], {}, "Invalid plate-map"),
])
def test_bad_treatment_maps_rejected(tmp_path, wells, changes, message):
    marker(tmp_path, "XY01", "A01")
    path = plate_map(tmp_path, wells, **changes)
    with pytest.raises(ValueError, match=message):
        load_acquisition_mapping(tmp_path, ["XY01"], path)


def test_mapping_rejected_before_model_or_output_creation(tmp_path, monkeypatch):
    raw = tmp_path / "inputs"
    marker(raw, "XY01", "A01")
    (raw / "XY01/_A01").unlink()
    model = Mock(side_effect=AssertionError("Model must not run"))
    monkeypatch.setattr("dapi_norm.user_cellpose_batch.run_nuclei_count_batch", model)
    with pytest.raises(ValueError, match="well marker"):
        run_user_cellpose_batch(input_root=raw, output_root=tmp_path / "outputs")
    model.assert_not_called()
    assert not (tmp_path / "outputs").exists()
    with pytest.raises(ValueError, match="well marker"):
        build_per_acquisition_runs(input_root=raw, output_root=tmp_path / "outputs",
            python_executable=Path("python"), runner_script=Path("runner.py"))


def test_workbook_export_rejects_unmapped_result(tmp_path):
    with pytest.raises(ValueError, match="Verified well mapping"):
        write_pi_style_workbook_from_rows(tmp_path / "bad.xlsx", [{"source_id": "XY01"}])
    assert not (tmp_path / "bad.xlsx").exists()


def test_shared_treatment_map_rejects_ambiguous_assays(tmp_path):
    acquisitions = [SimpleNamespace(input_root=tmp_path / p / "Assay") for p in ("P1", "P2")]
    with pytest.raises(ValueError, match="ambiguous"):
        mapping_for_acquisitions(acquisitions, tmp_path / "map.csv")


@pytest.mark.parametrize('change', [
    {'well_id':'B02'}, {'well_marker_sha256':''}, {'fields_in_well':0},
    {'compound':'guessed drug'}, {'treatment_mapping_status':'provided'},
    {'field_id':'XY02'}, {'plate_map_source':'unknown.csv'},
])
def test_export_rejects_inconsistent_provenance(tmp_path, change):
    marker(tmp_path,'XY01','A01')
    row = dict(load_acquisition_mapping(tmp_path,['XY01'])['XY01'],source_id='Assay/XY01')
    row.update(change)
    with pytest.raises(ValueError):
        validate_mapping_rows([row])


def test_export_accepts_complete_marker_and_treatment_provenance(tmp_path):
    marker(tmp_path,'XY01','A01')
    plate_map(tmp_path,['A01'])
    row = dict(load_acquisition_mapping(tmp_path,['XY01'])['XY01'],source_id='Assay/XY01')
    validate_mapping_rows([row])


def test_duplicate_channel_is_rejected_before_planning(tmp_path):
    marker(tmp_path,'XY01','A01')
    (tmp_path/'XY01/second_XY01_CH2.tif').touch()
    with pytest.raises(ValueError,match='Duplicate channel CH2'):
        build_per_acquisition_runs(input_root=tmp_path,output_root=tmp_path/'output',
            python_executable=Path('python'),runner_script=Path('runner.py'))


def test_shared_plate_map_is_forwarded_as_absolute_path(tmp_path, monkeypatch):
    raw = tmp_path/'Assay'
    marker(raw,'XY01','A01')
    path = plate_map(raw,['A01'])
    monkeypatch.chdir(tmp_path)
    run = build_per_acquisition_runs(input_root=raw,output_root=tmp_path/'output',
        python_executable=Path('python'),runner_script=Path('runner.py'),
        plate_map=Path('Assay/plate_map.csv'))[0]
    assert run.command[run.command.index('--plate-map')+1] == str(path.resolve())
