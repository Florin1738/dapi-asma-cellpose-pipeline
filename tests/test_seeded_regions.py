from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image
import tifffile
import yaml

from dapi_norm.pi_simple_summary import ImagePair
from dapi_norm.seeded_regions import (
    _format_percent,
    _non_seeded_area_rgb,
    _retained_region_overlay_rgb,
    build_ch2_foreground_mask,
    load_mask_lookup_from_counts_root,
    measure_seeded_region_image,
    segment_seeded_regions_propagation,
    run_seeded_region_batch,
    score_seeded_region_qc,
    segment_seeded_regions_random_walker,
    segment_seeded_regions,
    write_seeded_region_crop_panel,
    write_seeded_region_qc_panel,
)


def test_build_ch2_foreground_mask_uses_data_driven_threshold_and_cleanup():
    image = np.zeros((20, 20), dtype=np.uint16) + 100
    image[5:15, 6:16] = 1000

    mask, stats = build_ch2_foreground_mask(image, method="otsu", min_size=16)

    assert mask[10, 10]
    assert not mask[0, 0]
    assert stats.method == "otsu"
    assert 100 < stats.threshold < 1000
    assert stats.foreground_area_px == 100
    assert stats.image_area_px == 400


def test_seeded_intensity_watershed_partitions_foreground_from_nuclei_without_radius():
    ch2 = np.zeros((30, 40), dtype=np.uint16) + 50
    ch2[8:22, 5:35] = 500
    nuclei = np.zeros((30, 40), dtype=np.uint32)
    nuclei[12:15, 10:13] = 1
    nuclei[12:15, 27:30] = 2
    foreground = ch2 > 200

    labels, stats = segment_seeded_regions(ch2, nuclei, foreground)

    assert set(np.unique(labels)) == {0, 1, 2}
    assert labels[13, 11] == 1
    assert labels[13, 28] == 2
    assert labels[0, 0] == 0
    assert stats.method == "seeded_intensity_watershed"
    assert stats.nucleus_labels == 2
    assert stats.foreground_components_with_seed == 1


def test_seeded_intensity_random_walker_partitions_foreground_from_nuclei_without_radius():
    ch2 = np.zeros((30, 40), dtype=np.uint16) + 50
    ch2[8:22, 5:35] = 500
    ch2[:, 19:21] = 150
    nuclei = np.zeros((30, 40), dtype=np.uint32)
    nuclei[12:15, 10:13] = 7
    nuclei[12:15, 27:30] = 11
    foreground = ch2 > 100

    labels, stats = segment_seeded_regions_random_walker(ch2, nuclei, foreground, beta=30.0)

    assert set(np.unique(labels)) == {0, 7, 11}
    assert labels[13, 11] == 7
    assert labels[13, 28] == 11
    assert labels[0, 0] == 0
    assert stats.method == "seeded_intensity_random_walker"
    assert stats.nucleus_labels == 2
    assert stats.foreground_components_with_seed == 1


def test_seeded_intensity_random_walker_handles_uniform_components_without_solver_warning(recwarn):
    ch2 = np.zeros((24, 32), dtype=np.uint16) + 50
    ch2[6:18, 5:27] = 500
    nuclei = np.zeros((24, 32), dtype=np.uint32)
    nuclei[10:13, 9:12] = 1
    nuclei[10:13, 21:24] = 2
    foreground = ch2 > 100

    labels, stats = segment_seeded_regions_random_walker(ch2, nuclei, foreground, beta=90.0)

    assert set(np.unique(labels)) == {0, 1, 2}
    assert labels[11, 10] == 1
    assert labels[11, 22] == 2
    assert stats.method == "seeded_intensity_random_walker"
    assert len(recwarn) == 0


def test_seeded_intensity_propagation_uses_cellprofiler_style_seeded_growth():
    ch2 = np.zeros((36, 48), dtype=np.uint16) + 50
    ch2[8:28, 5:43] = 800
    ch2[14:22, 22:26] = 100
    ch2[2:6, 2:8] = 900
    nuclei = np.zeros((36, 48), dtype=np.uint32)
    nuclei[17:20, 11:14] = 3
    nuclei[17:20, 35:38] = 9
    foreground = ch2 > 200

    labels, stats = segment_seeded_regions_propagation(
        ch2,
        nuclei,
        foreground,
        regularization_factor=0.05,
    )

    assert stats.method == "seeded_intensity_propagation"
    assert stats.nucleus_labels == 2
    assert stats.foreground_components == 2
    assert stats.foreground_components_with_seed == 1
    assert labels[18, 12] == 3
    assert labels[18, 36] == 9
    assert labels[0, 0] == 0
    assert labels[3, 4] == 0
    assert not np.any((labels > 0) & ~foreground)
    assert np.count_nonzero(labels > 0) > np.count_nonzero(nuclei > 0)


def test_seeded_intensity_propagation_does_not_measure_nuclei_outside_ch2_foreground():
    ch2 = np.zeros((24, 24), dtype=np.uint16) + 50
    nuclei = np.zeros((24, 24), dtype=np.uint32)
    nuclei[9:12, 9:12] = 1
    foreground = np.zeros((24, 24), dtype=bool)

    labels, stats = segment_seeded_regions_propagation(
        ch2,
        nuclei,
        foreground,
        regularization_factor=0.05,
    )

    assert np.count_nonzero(labels) == 0
    assert stats.foreground_area_px == 0
    assert stats.seeded_region_area_px == 0
    assert stats.unseeded_foreground_area_px == 0


def test_seeded_intensity_propagation_does_not_assign_foreground_without_overlapping_seed():
    ch2 = np.zeros((24, 24), dtype=np.uint16) + 50
    ch2[10:13, 10:13] = 1000
    nuclei = np.zeros((24, 24), dtype=np.uint32)
    nuclei[9, 10:13] = 1
    foreground = ch2 > 200

    labels, stats = segment_seeded_regions_propagation(
        ch2,
        nuclei,
        foreground,
        regularization_factor=0.05,
    )

    assert np.count_nonzero(labels) == 0
    assert stats.foreground_area_px == 9
    assert stats.seeded_region_area_px == 0
    assert stats.unseeded_foreground_area_px == 9
    assert stats.foreground_components == 1
    assert stats.foreground_components_with_seed == 0


def test_measure_seeded_regions_reports_non_seeded_area_and_per_nucleus_endpoint():
    ch2 = np.array([[10, 10, 10], [10, 100, 100], [10, 100, 100]], dtype=np.uint16)
    labels = np.array([[0, 0, 0], [0, 1, 1], [0, 1, 1]], dtype=np.uint32)
    nuclei = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 0]], dtype=np.uint32)

    row = measure_seeded_region_image("XY01", ch2, nuclei, labels, background_value=10)

    assert row["image_id"] == "XY01"
    assert row["method"] == "seeded_intensity_watershed"
    assert row["dapi_positive_nucleus_count"] == 1
    assert row["seeded_region_area_px"] == 4
    assert row["non_seeded_area_px"] == 5
    assert row["seeded_region_integrated_raw"] == 400.0
    assert row["seeded_region_integrated_background_corrected"] == 360.0
    assert row["seeded_region_intensity_per_DAPI_positive_nucleus"] == 360.0
    assert "low_nucleus_count_qc_required" in row["warnings"]


def test_write_seeded_region_qc_panel_creates_visual_artifact(tmp_path: Path):
    ch2 = np.zeros((24, 24), dtype=np.uint16)
    ch2[6:18, 5:19] = 2000
    ch4 = np.zeros((24, 24), dtype=np.uint16)
    ch4[9:12, 10:13] = 3000
    nuclei = np.zeros((24, 24), dtype=np.uint32)
    nuclei[9:12, 10:13] = 1
    foreground = ch2 > 100
    labels = np.zeros((24, 24), dtype=np.uint32)
    labels[6:18, 5:19] = 1

    output_path = tmp_path / "qc" / "XY01_seeded_region_qc.png"
    write_seeded_region_qc_panel(
        image_id="XY01",
        ch2_image=ch2,
        ch4_image=ch4,
        nuclei_mask=nuclei,
        foreground_mask=foreground,
        seeded_labels=labels,
        output_path=output_path,
        metrics={
            "seeded_region_integrated_raw": 336000.0,
            "seeded_region_intensity_per_DAPI_positive_nucleus": 336000.0,
        },
    )

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_non_seeded_area_overlay_marks_only_unseeded_foreground_magenta():
    ch2_scaled = np.ones((3, 3), dtype=np.float32) * 0.5
    region_mask = np.zeros((3, 3), dtype=bool)
    foreground_mask = np.zeros((3, 3), dtype=bool)
    region_mask[1, 1] = True
    foreground_mask[1, 1] = True
    foreground_mask[0, 1] = True

    rgb = _non_seeded_area_rgb(ch2_scaled, region_mask, foreground_mask)

    assert np.allclose(rgb[1, 1], [0.0, 0.75, 0.18])
    assert np.allclose(rgb[0, 1], [0.9, 0.0, 0.9])
    assert np.all(rgb[2, 2] < 0.5)


def test_retained_region_overlay_fills_segmentation_light_green_without_hiding_ch2():
    ch2_scaled = np.ones((3, 3), dtype=np.float32) * 0.5
    region_mask = np.zeros((3, 3), dtype=bool)
    region_mask[1, 1] = True

    rgb = _retained_region_overlay_rgb(ch2_scaled, region_mask)

    assert rgb[1, 1, 1] > rgb[1, 1, 0]
    assert rgb[1, 1, 1] > rgb[1, 1, 2]
    assert 0.2 < rgb[1, 1, 0] < 0.8
    assert np.allclose(rgb[0, 0], [0.5, 0.5, 0.5])


def test_format_percent_marks_tiny_nonzero_values():
    assert _format_percent(0.0) == "0.0%"
    assert _format_percent(0.0002) == "<0.1%"
    assert _format_percent(0.0123) == "1.2%"


def test_write_seeded_region_crop_panel_creates_boundary_crop_artifact(tmp_path: Path):
    ch2 = np.zeros((32, 32), dtype=np.uint16)
    ch2[8:24, 9:25] = 2000
    ch4 = np.zeros((32, 32), dtype=np.uint16)
    ch4[14:17, 15:18] = 3000
    nuclei = np.zeros((32, 32), dtype=np.uint32)
    nuclei[14:17, 15:18] = 1
    labels = np.zeros((32, 32), dtype=np.uint32)
    labels[8:24, 9:25] = 1

    output_path = tmp_path / "crops" / "crop_panel.png"
    write_seeded_region_crop_panel(
        crops=[
            {
                "image_id": "XY01",
                "ch2_image": ch2,
                "ch4_image": ch4,
                "nuclei_mask": nuclei,
                "seeded_labels": labels,
                "box": (4, 4, 28, 28),
                "caption": "reviewable_not_validated | synthetic crop",
            }
        ],
        output_path=output_path,
    )

    assert output_path.exists()
    assert output_path.stat().st_size > 0
    with Image.open(output_path) as image:
        arr = np.asarray(image.convert("RGB"))
        width, height = image.size
    assert width / height > 4.0
    right_caption_panel = arr[:, int(arr.shape[1] * 0.78) :, :]
    assert np.count_nonzero(np.any(right_caption_panel < 245, axis=2)) > 100


def test_run_seeded_region_batch_writes_summary_masks_qc_and_log(tmp_path: Path):
    dataset = tmp_path / "dataset"
    xy01 = dataset / "XY01"
    xy01.mkdir(parents=True)
    ch2 = np.zeros((28, 28), dtype=np.uint16) + 100
    ch2[8:22, 7:23] = 1000
    ch4 = np.zeros((28, 28), dtype=np.uint16)
    ch4[12:15, 13:16] = 2000
    ch2_path = xy01 / "sample_XY01_CH2.tif"
    ch4_path = xy01 / "sample_XY01_CH4.tif"
    tifffile.imwrite(ch2_path, ch2)
    tifffile.imwrite(ch4_path, ch4)
    nuclei_mask = np.zeros((28, 28), dtype=np.uint32)
    nuclei_mask[12:15, 13:16] = 1
    nucleus_mask_path = tmp_path / "counts" / "masks" / "XY01_CH4_fake_labels.tif"
    nucleus_mask_path.parent.mkdir(parents=True)
    tifffile.imwrite(nucleus_mask_path, nuclei_mask)
    pair = ImagePair(location="XY01", source_id="XY01", ch2_path=ch2_path, ch4_path=ch4_path)

    rows = run_seeded_region_batch(
        image_pairs=[pair],
        mask_lookup={"XY01": nucleus_mask_path},
        output_dir=tmp_path / "out",
        foreground_method="otsu",
        background_value=100,
    )

    assert len(rows) == 1
    assert rows[0]["method"] == "seeded_intensity_watershed"
    assert rows[0]["foreground_method"] == "otsu"
    assert rows[0]["dapi_positive_nucleus_count"] == 1
    assert (tmp_path / "out" / "masks" / "XY01_seeded_intensity_watershed_labels.tif").exists()
    assert (tmp_path / "out" / "summaries" / "seeded_region_image_metrics.csv").exists()
    assert (tmp_path / "out" / "qc" / "XY01_seeded_region_qc.png").exists()
    assert (tmp_path / "out" / "qc_contact_sheet.png").exists()
    config_path = tmp_path / "out" / "logs" / "config_resolved.yaml"
    assert config_path.exists()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    assert config["method"]["name"] == "seeded_intensity_watershed"
    assert config["method"]["whole_cell_claim"] is False
    summary_text = (tmp_path / "out" / "summaries" / "seeded_region_image_metrics.csv").read_text(
        encoding="utf-8"
    )
    assert "non_seeded_area_px" in summary_text
    assert "excluded_background_area_px" not in summary_text


def test_run_seeded_region_batch_can_use_random_walker_comparator(tmp_path: Path):
    dataset = tmp_path / "dataset"
    xy01 = dataset / "XY01"
    xy01.mkdir(parents=True)
    ch2 = np.zeros((28, 28), dtype=np.uint16) + 100
    ch2[8:22, 7:23] = 1000
    ch4 = np.zeros((28, 28), dtype=np.uint16)
    ch4[12:15, 10:13] = 2000
    ch4[12:15, 18:21] = 2000
    ch2_path = xy01 / "sample_XY01_CH2.tif"
    ch4_path = xy01 / "sample_XY01_CH4.tif"
    tifffile.imwrite(ch2_path, ch2)
    tifffile.imwrite(ch4_path, ch4)
    nuclei_mask = np.zeros((28, 28), dtype=np.uint32)
    nuclei_mask[12:15, 10:13] = 1
    nuclei_mask[12:15, 18:21] = 2
    nucleus_mask_path = tmp_path / "counts" / "masks" / "XY01_CH4_fake_labels.tif"
    nucleus_mask_path.parent.mkdir(parents=True)
    tifffile.imwrite(nucleus_mask_path, nuclei_mask)
    pair = ImagePair(location="XY01", source_id="XY01", ch2_path=ch2_path, ch4_path=ch4_path)

    rows = run_seeded_region_batch(
        image_pairs=[pair],
        mask_lookup={"XY01": nucleus_mask_path},
        output_dir=tmp_path / "out",
        foreground_method="otsu",
        background_value=100,
        segmentation_method="random_walker",
        random_walker_beta=30.0,
    )

    assert rows[0]["method"] == "seeded_intensity_random_walker"
    assert (tmp_path / "out" / "masks" / "XY01_seeded_intensity_random_walker_labels.tif").exists()
    config_path = tmp_path / "out" / "logs" / "config_resolved.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    assert config["method"]["name"] == "seeded_intensity_random_walker"
    assert config["method"]["random_walker_beta"] == 30.0


def test_run_seeded_region_batch_can_use_cellprofiler_style_propagation(tmp_path: Path):
    dataset = tmp_path / "dataset"
    xy01 = dataset / "XY01"
    xy01.mkdir(parents=True)
    ch2 = np.zeros((28, 28), dtype=np.uint16) + 100
    ch2[8:22, 7:23] = 1000
    ch4 = np.zeros((28, 28), dtype=np.uint16)
    ch4[12:15, 10:13] = 2000
    ch4[12:15, 18:21] = 2000
    ch2_path = xy01 / "sample_XY01_CH2.tif"
    ch4_path = xy01 / "sample_XY01_CH4.tif"
    tifffile.imwrite(ch2_path, ch2)
    tifffile.imwrite(ch4_path, ch4)
    nuclei_mask = np.zeros((28, 28), dtype=np.uint32)
    nuclei_mask[12:15, 10:13] = 1
    nuclei_mask[12:15, 18:21] = 2
    nucleus_mask_path = tmp_path / "counts" / "masks" / "XY01_CH4_fake_labels.tif"
    nucleus_mask_path.parent.mkdir(parents=True)
    tifffile.imwrite(nucleus_mask_path, nuclei_mask)
    pair = ImagePair(location="XY01", source_id="XY01", ch2_path=ch2_path, ch4_path=ch4_path)

    rows = run_seeded_region_batch(
        image_pairs=[pair],
        mask_lookup={"XY01": nucleus_mask_path},
        output_dir=tmp_path / "out",
        foreground_method="otsu",
        background_value=100,
        segmentation_method="propagation",
        propagation_regularization=0.05,
    )

    assert rows[0]["method"] == "seeded_intensity_propagation"
    assert (tmp_path / "out" / "masks" / "XY01_seeded_intensity_propagation_labels.tif").exists()
    config_path = tmp_path / "out" / "logs" / "config_resolved.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    assert config["method"]["name"] == "seeded_intensity_propagation"
    assert config["method"]["propagation_regularization"] == 0.05


def test_seeded_region_metrics_warn_when_foreground_component_lacks_dapi_seed():
    ch2 = np.zeros((10, 10), dtype=np.uint16)
    labels = np.zeros((10, 10), dtype=np.uint32)
    labels[1:3, 1:3] = 1
    nuclei = np.zeros((10, 10), dtype=np.uint32)
    nuclei[1:2, 1:2] = 1
    from dapi_norm.seeded_regions import ForegroundStats, SeededSegmentationStats

    row = measure_seeded_region_image(
        "XY01",
        ch2,
        nuclei,
        labels,
        foreground_stats=ForegroundStats(
            method="otsu",
            threshold=1.0,
            image_area_px=100,
            foreground_area_px=10,
            foreground_fraction=0.1,
            min_size=1,
            fill_holes=True,
        ),
        segmentation_stats=SeededSegmentationStats(
            method="seeded_intensity_watershed",
            nucleus_labels=1,
            foreground_area_px=10,
            seeded_region_area_px=4,
            unseeded_foreground_area_px=6,
            foreground_components=2,
            foreground_components_with_seed=1,
        ),
    )

    assert row["unseeded_foreground_area_px"] == 6
    assert "foreground_components_without_dapi_seed_excluded" in row["warnings"]


def test_score_seeded_region_qc_flags_artifact_like_low_nucleus_high_area_case():
    row = {
        "dapi_positive_nucleus_count": 3,
        "seeded_region_fraction": 0.25,
        "unseeded_foreground_fraction": 0.17,
        "foreground_components": 50,
        "foreground_components_with_seed": 2,
    }

    status, flags = score_seeded_region_qc(row)

    assert status == "reject_qc_failure"
    assert "low_nucleus_count" in flags
    assert "large_seeded_area_with_low_nucleus_count" in flags
    assert "high_unseeded_foreground_fraction" in flags
    assert "low_fraction_of_foreground_components_with_dapi_seed" in flags


def test_score_seeded_region_qc_marks_cleaner_field_reviewable_not_validated():
    row = {
        "dapi_positive_nucleus_count": 120,
        "seeded_region_fraction": 0.28,
        "unseeded_foreground_fraction": 0.02,
        "foreground_components": 50,
        "foreground_components_with_seed": 48,
    }

    status, flags = score_seeded_region_qc(row)

    assert status == "reviewable_not_validated"
    assert flags == ["not_validated_whole_cell_mask"]


def test_score_seeded_region_qc_marks_sizeable_unseeded_foreground_for_manual_review():
    row = {
        "dapi_positive_nucleus_count": 57,
        "seeded_region_fraction": 0.04,
        "unseeded_foreground_fraction": 0.13,
        "foreground_components": 43,
        "foreground_components_with_seed": 32,
    }

    status, flags = score_seeded_region_qc(row)

    assert status == "needs_manual_review"
    assert "sizeable_unseeded_target_foreground" in flags


def test_score_seeded_region_qc_rejects_near_full_field_region_even_with_many_nuclei():
    row = {
        "dapi_positive_nucleus_count": 250,
        "seeded_region_fraction": 0.94,
        "unseeded_foreground_fraction": 0.01,
        "foreground_components": 20,
        "foreground_components_with_seed": 20,
    }

    status, flags = score_seeded_region_qc(row)

    assert status == "reject_qc_failure"
    assert "near_full_field_seeded_region_fraction" in flags


def test_score_seeded_region_qc_rejects_low_nucleus_fields_for_method_comparison():
    row = {
        "dapi_positive_nucleus_count": 3,
        "seeded_region_fraction": 0.0,
        "unseeded_foreground_fraction": 0.0,
        "foreground_components": 0,
        "foreground_components_with_seed": 0,
    }

    status, flags = score_seeded_region_qc(row)

    assert status == "reject_qc_failure"
    assert "low_nucleus_count" in flags


def test_seeded_region_metrics_handle_zero_nuclei_without_finite_per_nucleus_value():
    ch2 = np.ones((5, 5), dtype=np.uint16) * 100
    labels = np.zeros((5, 5), dtype=np.uint32)
    labels[1:4, 1:4] = 1
    nuclei = np.zeros((5, 5), dtype=np.uint32)

    row = measure_seeded_region_image("XY00", ch2, nuclei, labels, background_value=10)

    assert row["dapi_positive_nucleus_count"] == 0
    assert np.isnan(row["seeded_region_intensity_per_DAPI_positive_nucleus"])
    assert "zero_nucleus_count" in row["warnings"]


def test_seeded_batch_records_channel_extraction_and_validation_status(tmp_path: Path):
    dataset = tmp_path / "dataset"
    xy01 = dataset / "XY01"
    xy01.mkdir(parents=True)
    ch2 = np.zeros((16, 16), dtype=np.uint16)
    ch2[4:12, 4:12] = 1000
    ch4 = np.zeros((16, 16), dtype=np.uint16)
    ch4[7:9, 7:9] = 2000
    ch2_path = xy01 / "sample_XY01_CH2.tif"
    ch4_path = xy01 / "sample_XY01_CH4.tif"
    tifffile.imwrite(ch2_path, ch2)
    tifffile.imwrite(ch4_path, ch4)
    nuclei_mask_path = tmp_path / "counts" / "masks" / "XY01_CH4_fake_labels.tif"
    nuclei_mask_path.parent.mkdir(parents=True)
    nuclei = np.zeros((16, 16), dtype=np.uint32)
    nuclei[7:9, 7:9] = 1
    tifffile.imwrite(nuclei_mask_path, nuclei)

    run_seeded_region_batch(
        image_pairs=[
            ImagePair(location="XY01", source_id="XY01", ch2_path=ch2_path, ch4_path=ch4_path)
        ],
        mask_lookup={"XY01": nuclei_mask_path},
        output_dir=tmp_path / "out",
        foreground_method="otsu",
    )

    config = yaml.safe_load((tmp_path / "out" / "logs" / "config_resolved.yaml").read_text())
    assert config["channel_extraction"]["ch2_channel"] == "CH2/aSMA"
    assert config["channel_extraction"]["ch4_channel"] == "CH4/DAPI"
    assert config["channel_extraction"]["z_projection"] == "none"
    assert config["validation_status"]["manual_ground_truth_available"] is False
    assert config["validation_status"]["precision_recall_f1_allowed"] is False


def test_load_mask_lookup_from_counts_root_reads_count_summary_and_resolves_paths(tmp_path: Path):
    counts_root = tmp_path / "counts"
    summaries = counts_root / "summaries"
    masks = counts_root / "masks"
    summaries.mkdir(parents=True)
    masks.mkdir()
    mask_path = masks / "XY01_CH4_fake_labels.tif"
    tifffile.imwrite(mask_path, np.ones((4, 4), dtype=np.uint32))
    (summaries / "nucleus_counts.csv").write_text(
        "image_id,mask_path\n"
        "XY01,masks/XY01_CH4_fake_labels.tif\n",
        encoding="utf-8",
    )

    lookup = load_mask_lookup_from_counts_root(counts_root)

    assert lookup["XY01"] == mask_path
    assert lookup["xy01"] == mask_path


def test_load_mask_lookup_from_counts_root_rejects_duplicate_discovered_xy_masks(tmp_path: Path):
    first = tmp_path / "run_a" / "masks" / "XY01_labels.tif"
    second = tmp_path / "run_b" / "masks" / "XY01_labels.tif"
    first.parent.mkdir(parents=True)
    second.parent.mkdir(parents=True)
    tifffile.imwrite(first, np.zeros((4, 4), dtype=np.uint32))
    tifffile.imwrite(second, np.zeros((4, 4), dtype=np.uint32))

    import pytest

    with pytest.raises(ValueError, match="Duplicate nuclei masks"):
        load_mask_lookup_from_counts_root(tmp_path)
