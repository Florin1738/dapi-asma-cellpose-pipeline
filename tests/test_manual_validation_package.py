from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
import tifffile
import yaml

from dapi_norm.manual_validation_package import prepare_manual_validation_package


def test_prepare_manual_validation_package_writes_manifest_guides_and_blank_masks(
    tmp_path: Path,
    monkeypatch,
):
    project_root = tmp_path / "project"
    input_root = project_root / "data" / "run"
    xy01 = input_root / "XY01"
    seeded_run = project_root / "output" / "seeded_asma_regions" / "run1"
    summary_dir = seeded_run / "summaries"
    logs_dir = seeded_run / "logs"
    seeded_mask_path = seeded_run / "masks" / "XY01_candidate_labels.tif"
    nuclei_mask_path = project_root / "output" / "counts" / "masks" / "XY01_nuclei.tif"
    xy01.mkdir(parents=True)
    summary_dir.mkdir(parents=True)
    logs_dir.mkdir(parents=True)
    seeded_mask_path.parent.mkdir(parents=True)
    nuclei_mask_path.parent.mkdir(parents=True)
    ch2 = np.zeros((32, 32), dtype=np.uint16)
    ch2[8:24, 9:25] = 2000
    ch4 = np.zeros((32, 32), dtype=np.uint16)
    ch4[14:17, 15:18] = 3000
    candidate = np.zeros((32, 32), dtype=np.uint32)
    candidate[8:24, 9:25] = 1
    nuclei = np.zeros((32, 32), dtype=np.uint32)
    nuclei[14:17, 15:18] = 1
    tifffile.imwrite(xy01 / "sample_XY01_CH2.tif", ch2)
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", ch4)
    tifffile.imwrite(seeded_mask_path, candidate)
    tifffile.imwrite(nuclei_mask_path, nuclei)
    (summary_dir / "seeded_region_image_metrics.csv").write_text(
        "\n".join(
            [
                "image_id,source_id,method,foreground_method,dapi_positive_nucleus_count,"
                "seeded_region_integrated_raw,seeded_region_intensity_per_DAPI_positive_nucleus,"
                "mask_path,qc_status,qc_flags",
                "XY01,XY01,seeded_intensity_propagation,otsu,1,512000.0,512000.0,"
                "output/seeded_asma_regions/run1/masks/XY01_candidate_labels.tif,"
                "reviewable_not_validated,not_validated_whole_cell_mask",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (logs_dir / "config_resolved.yaml").write_text(
        yaml.safe_dump(
            {
                "output_dir": "output/seeded_asma_regions/run1",
                "image_records": [
                    {
                        "source_id": "XY01",
                        "location": "XY01",
                        "nuclei_mask_path": "output/counts/masks/XY01_nuclei.tif",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    output_dir = project_root / "manual_validation" / "package1"
    monkeypatch.chdir(tmp_path)

    outputs = prepare_manual_validation_package(
        input_root=input_root,
        seeded_run_dir=seeded_run,
        output_dir=output_dir,
        positions=["XY01"],
        iou_threshold=0.5,
        task="asma_associated_region",
    )

    assert outputs["manifest"].exists()
    assert outputs["status"].exists()
    assert outputs["readme"].exists()
    assert outputs["annotation_dir"].joinpath("XY01_manual_annotation_panel.png").exists()
    assert outputs["guide_dir"].joinpath("XY01_manual_validation_guide.png").exists()
    blank_mask_path = outputs["reference_mask_dir"] / "XY01_manual_reference_labels.tif"
    blank = tifffile.imread(blank_mask_path)
    assert blank.shape == (32, 32)
    assert blank.dtype == np.uint32
    assert np.count_nonzero(blank) == 0

    with outputs["manifest"].open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["image_id"] == "XY01"
    assert rows[0]["manual_reference_mask_path"] == str(blank_mask_path)
    assert rows[0]["annotation_panel_path"].endswith("XY01_manual_annotation_panel.png")
    assert rows[0]["candidate_mask_path"].endswith("XY01_candidate_labels.tif")
    assert rows[0]["validation_task"] == "asma_associated_region"

    readme = outputs["readme"].read_text(encoding="utf-8")
    assert "scripts/validate_manual_instance_masks.py" in readme
    assert "--iou-threshold 0.5" in readme
    assert "Do not report precision, recall, F1, or IoU until manual masks are filled" in readme
    assert "Primary instance unit" in readme
    assert "Do not trace the candidate overlay" in readme
    assert "manual_labeling_status.csv" in readme
    assert "--completion-status" in readme

    with outputs["status"].open(newline="", encoding="utf-8") as handle:
        status_rows = list(csv.DictReader(handle))
    assert status_rows[0]["image_id"] == "XY01"
    assert status_rows[0]["status"] == "not_started"
    assert Path(status_rows[0]["manual_reference_mask_path"]).is_absolute()
    assert Path(status_rows[0]["manual_reference_mask_path"]).resolve() == blank_mask_path.resolve()


def test_prepare_manual_validation_package_refuses_to_overwrite_existing_manual_labels(
    tmp_path: Path,
    monkeypatch,
):
    project_root, input_root, seeded_run = _write_tiny_seeded_package_source(tmp_path)
    output_dir = project_root / "manual_validation" / "package1"
    monkeypatch.chdir(tmp_path)
    prepare_manual_validation_package(
        input_root=input_root,
        seeded_run_dir=seeded_run,
        output_dir=output_dir,
        positions=["XY01"],
    )
    manual_mask_path = output_dir / "reference_masks_to_fill" / "XY01_manual_reference_labels.tif"
    filled = tifffile.imread(manual_mask_path)
    filled[2:4, 2:4] = 7
    tifffile.imwrite(manual_mask_path, filled)

    with pytest.raises(FileExistsError, match="manual reference mask already exists"):
        prepare_manual_validation_package(
            input_root=input_root,
            seeded_run_dir=seeded_run,
            output_dir=output_dir,
            positions=["XY01"],
        )

    assert np.max(tifffile.imread(manual_mask_path)) == 7


def test_prepare_manual_validation_package_force_overwrites_placeholders(
    tmp_path: Path,
    monkeypatch,
):
    project_root, input_root, seeded_run = _write_tiny_seeded_package_source(tmp_path)
    output_dir = project_root / "manual_validation" / "package1"
    monkeypatch.chdir(tmp_path)
    outputs = prepare_manual_validation_package(
        input_root=input_root,
        seeded_run_dir=seeded_run,
        output_dir=output_dir,
        positions=["XY01"],
    )
    manual_mask_path = output_dir / "reference_masks_to_fill" / "XY01_manual_reference_labels.tif"
    filled = tifffile.imread(manual_mask_path)
    filled[2:4, 2:4] = 7
    tifffile.imwrite(manual_mask_path, filled)
    outputs["status"].write_text(
        "\n".join(
            [
                "image_id,manual_reference_mask_path,annotation_panel_path,status,labeler,completed_date,notes",
                f"XY01,{manual_mask_path},example.png,complete_non_empty,tester,2026-06-29,temporary",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    prepare_manual_validation_package(
        input_root=input_root,
        seeded_run_dir=seeded_run,
        output_dir=output_dir,
        positions=["XY01"],
        force_overwrite_reference_masks=True,
    )

    assert np.count_nonzero(tifffile.imread(manual_mask_path)) == 0
    with (output_dir / "manual_labeling_status.csv").open(newline="", encoding="utf-8") as handle:
        status_rows = list(csv.DictReader(handle))
    assert status_rows[0]["status"] == "not_started"


def _write_tiny_seeded_package_source(tmp_path: Path) -> tuple[Path, Path, Path]:
    project_root = tmp_path / "project"
    input_root = project_root / "data" / "run"
    xy01 = input_root / "XY01"
    seeded_run = project_root / "output" / "seeded_asma_regions" / "run1"
    summary_dir = seeded_run / "summaries"
    logs_dir = seeded_run / "logs"
    seeded_mask_path = seeded_run / "masks" / "XY01_candidate_labels.tif"
    nuclei_mask_path = project_root / "output" / "counts" / "masks" / "XY01_nuclei.tif"
    xy01.mkdir(parents=True)
    summary_dir.mkdir(parents=True)
    logs_dir.mkdir(parents=True)
    seeded_mask_path.parent.mkdir(parents=True)
    nuclei_mask_path.parent.mkdir(parents=True)
    ch2 = np.zeros((32, 32), dtype=np.uint16)
    ch2[8:24, 9:25] = 2000
    ch4 = np.zeros((32, 32), dtype=np.uint16)
    ch4[14:17, 15:18] = 3000
    candidate = np.zeros((32, 32), dtype=np.uint32)
    candidate[8:24, 9:25] = 1
    nuclei = np.zeros((32, 32), dtype=np.uint32)
    nuclei[14:17, 15:18] = 1
    tifffile.imwrite(xy01 / "sample_XY01_CH2.tif", ch2)
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", ch4)
    tifffile.imwrite(seeded_mask_path, candidate)
    tifffile.imwrite(nuclei_mask_path, nuclei)
    (summary_dir / "seeded_region_image_metrics.csv").write_text(
        "\n".join(
            [
                "image_id,source_id,method,foreground_method,dapi_positive_nucleus_count,"
                "seeded_region_integrated_raw,seeded_region_intensity_per_DAPI_positive_nucleus,"
                "mask_path,qc_status,qc_flags",
                "XY01,XY01,seeded_intensity_propagation,otsu,1,512000.0,512000.0,"
                "output/seeded_asma_regions/run1/masks/XY01_candidate_labels.tif,"
                "reviewable_not_validated,not_validated_whole_cell_mask",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (logs_dir / "config_resolved.yaml").write_text(
        yaml.safe_dump(
            {
                "output_dir": "output/seeded_asma_regions/run1",
                "image_records": [
                    {
                        "source_id": "XY01",
                        "location": "XY01",
                        "nuclei_mask_path": "output/counts/masks/XY01_nuclei.tif",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    return project_root, input_root, seeded_run
