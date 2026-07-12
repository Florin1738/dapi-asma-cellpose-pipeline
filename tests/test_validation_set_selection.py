from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from dapi_norm.validation_set_selection import (
    ValidationFeatureRecord,
    render_validation_selection_panel,
    select_validation_set,
    write_selection_csvs,
)


def test_select_validation_set_preserves_must_include_and_covers_extremes():
    records = [
        _record("XY01", raw=10, nuclei=100, saturation=0.00),
        _record("XY02", raw=20, nuclei=20, saturation=0.01),
        _record("XY03", raw=30, nuclei=60, saturation=0.02),
        _record("XY04", raw=40, nuclei=5, saturation=0.10),
        _record("XY05", raw=50, nuclei=200, saturation=0.00),
        _record("XY06", raw=60, nuclei=30, saturation=0.03, method_jaccard=0.2),
        _record("XY07", raw=70, nuclei=70, saturation=0.04),
        _record("XY08", raw=80, nuclei=80, saturation=0.00),
    ]

    selected = select_validation_set(
        records,
        max_images=6,
        per_bucket=1,
        must_include=["XY03"],
    )

    selected_ids = [row.image_id for row in selected]
    assert selected_ids[0] == "XY03"
    assert "XY01" in selected_ids
    assert "XY08" in selected_ids
    assert "XY04" in selected_ids
    assert "XY06" in selected_ids
    reasons_by_id = {row.image_id: row.selection_reasons for row in selected}
    assert "must_include" in reasons_by_id["XY03"]
    assert "low_raw_target_integrated" in reasons_by_id["XY01"]
    assert "high_raw_target_integrated" in reasons_by_id["XY08"]
    assert "low_dapi_positive_nucleus_count" in reasons_by_id["XY04"]
    assert "high_method_disagreement" in reasons_by_id["XY06"]


def test_select_validation_set_deduplicates_reasons_for_same_image():
    records = [
        _record("XY01", raw=10, nuclei=10, saturation=0.0),
        _record("XY02", raw=100, nuclei=100, saturation=0.5),
    ]

    selected = select_validation_set(records, max_images=2, per_bucket=1, must_include=["XY02"])

    xy02 = next(row for row in selected if row.image_id == "XY02")
    reasons = xy02.selection_reasons.split(";")
    assert reasons.count("must_include") == 1
    assert len(reasons) == len(set(reasons))


def test_write_selection_csvs_and_panel(tmp_path: Path):
    selected = [
        _record("XY01", raw=10, nuclei=10, saturation=0.0),
        _record("XY02", raw=100, nuclei=100, saturation=0.5),
    ]
    selected[0].selection_reasons = "low_raw_target_integrated"
    selected[1].selection_reasons = "high_raw_target_integrated"

    outputs = write_selection_csvs(
        all_records=selected,
        selected_records=selected,
        output_dir=tmp_path,
    )
    render_validation_selection_panel(
        selected_records=selected,
        output_path=tmp_path / "panel.png",
        image_loader=_synthetic_image_loader,
        nuclei_loader=_synthetic_nuclei_loader,
    )

    assert outputs["all_features"].exists()
    assert outputs["selected"].exists()
    assert (tmp_path / "panel.png").exists()
    rows = _read_rows(outputs["selected"])
    assert rows[0]["selection_reasons"] == "low_raw_target_integrated"


def _record(
    image_id: str,
    *,
    raw: float,
    nuclei: int,
    saturation: float,
    method_jaccard: float | None = None,
) -> ValidationFeatureRecord:
    return ValidationFeatureRecord(
        image_id=image_id,
        source_id=image_id,
        ch2_path=Path(f"{image_id}_CH2.tif"),
        ch4_path=Path(f"{image_id}_CH4.tif"),
        nuclei_mask_path=Path(f"{image_id}_nuclei.tif"),
        dapi_positive_nucleus_count=nuclei,
        target_integrated_raw=raw,
        target_integrated_raw_per_DAPI_positive_nucleus=raw / nuclei if nuclei else float("nan"),
        target_saturation_fraction=saturation,
        dapi_saturation_fraction=0.0,
        method_region_jaccard=method_jaccard,
    )


def _synthetic_image_loader(_path: Path) -> np.ndarray:
    image = np.zeros((12, 12), dtype=np.uint16)
    image[2:8, 2:8] = 100
    return image


def _synthetic_nuclei_loader(_path: Path) -> np.ndarray:
    mask = np.zeros((12, 12), dtype=np.uint32)
    mask[3:5, 3:5] = 1
    mask[8:10, 8:10] = 2
    return mask


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))
