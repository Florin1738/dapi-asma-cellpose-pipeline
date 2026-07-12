from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from scripts.render_seeded_region_crops import (
    _largest_label_crop_box,
    _load_nuclei_mask_for_pair,
    main,
)


def test_largest_label_crop_box_falls_back_to_nuclei_when_no_seeded_region():
    labels = np.zeros((100, 100), dtype=np.uint32)
    nuclei = np.zeros((100, 100), dtype=np.uint32)
    nuclei[8:12, 78:82] = 1

    y0, x0, y1, x1 = _largest_label_crop_box(
        labels,
        nuclei_mask=nuclei,
        crop_size=20,
    )

    assert y0 <= 10 < y1
    assert x0 <= 80 < x1


def test_load_nuclei_mask_resolves_relative_paths_from_recorded_output_root(
    tmp_path: Path,
    monkeypatch,
):
    project_root = tmp_path / "project"
    seeded_run_dir = project_root / "output" / "seeded_asma_regions" / "run1"
    logs_dir = seeded_run_dir / "logs"
    mask_path = project_root / "output" / "counts" / "masks" / "XY01_labels.tif"
    logs_dir.mkdir(parents=True)
    mask_path.parent.mkdir(parents=True)
    expected = np.zeros((8, 8), dtype=np.uint32)
    expected[2:4, 3:5] = 1
    tifffile.imwrite(mask_path, expected)
    (logs_dir / "config_resolved.yaml").write_text(
        yaml.safe_dump(
            {
                "output_dir": "output/seeded_asma_regions/run1",
                "image_records": [
                    {
                        "source_id": "XY01",
                        "location": "XY01",
                        "nuclei_mask_path": "output/counts/masks/XY01_labels.tif",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.chdir(tmp_path)

    loaded = _load_nuclei_mask_for_pair(seeded_run_dir, "XY01")

    assert np.array_equal(loaded, expected)


def test_render_seeded_region_crops_resolves_relative_mask_path_from_outside_project(
    tmp_path: Path,
    monkeypatch,
):
    project_root = tmp_path / "project"
    input_root = project_root / "data" / "run"
    xy01 = input_root / "XY01"
    seeded_run_dir = project_root / "output" / "seeded_asma_regions" / "run1"
    summary_dir = seeded_run_dir / "summaries"
    logs_dir = seeded_run_dir / "logs"
    seeded_mask_path = seeded_run_dir / "masks" / "XY01_labels.tif"
    nuclei_mask_path = project_root / "output" / "counts" / "masks" / "XY01_nuclei.tif"
    xy01.mkdir(parents=True)
    summary_dir.mkdir(parents=True)
    logs_dir.mkdir(parents=True)
    seeded_mask_path.parent.mkdir(parents=True)
    nuclei_mask_path.parent.mkdir(parents=True)
    ch2 = np.zeros((24, 24), dtype=np.uint16)
    ch2[6:18, 6:18] = 2000
    ch4 = np.zeros((24, 24), dtype=np.uint16)
    ch4[10:13, 10:13] = 3000
    labels = np.zeros((24, 24), dtype=np.uint32)
    labels[6:18, 6:18] = 1
    nuclei = np.zeros((24, 24), dtype=np.uint32)
    nuclei[10:13, 10:13] = 1
    tifffile.imwrite(xy01 / "sample_XY01_CH2.tif", ch2)
    tifffile.imwrite(xy01 / "sample_XY01_CH4.tif", ch4)
    tifffile.imwrite(seeded_mask_path, labels)
    tifffile.imwrite(nuclei_mask_path, nuclei)
    (summary_dir / "seeded_region_image_metrics.csv").write_text(
        "\n".join(
            [
                "image_id,source_id,mask_path,qc_status,qc_flags",
                "XY01,XY01,output/seeded_asma_regions/run1/masks/XY01_labels.tif,"
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
    output_path = tmp_path / "review.png"
    monkeypatch.chdir(tmp_path)

    main(
        input_root=input_root,
        seeded_run_dir=seeded_run_dir,
        output_path=output_path,
        positions="XY01",
        crop_size=20,
    )

    assert output_path.exists()
    assert output_path.stat().st_size > 0
