from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import yaml

from dapi_norm.cellpose_cell_regions import (
    build_cellpose_cell_input,
    measure_cellpose_cell_region_image,
    run_cellpose_cell_region_batch,
    score_cellpose_cell_region_qc,
)
from dapi_norm.pi_simple_summary import ImagePair


def test_build_cellpose_cell_input_stacks_ch2_and_ch4_as_channel_first():
    ch2 = np.ones((5, 6), dtype=np.uint16) * 10
    ch4 = np.ones((5, 6), dtype=np.uint16) * 20

    stacked = build_cellpose_cell_input(ch2, ch4)

    assert stacked.shape == (2, 5, 6)
    assert stacked.dtype == np.float32
    assert np.all(stacked[0] == 10)
    assert np.all(stacked[1] == 20)


def test_measure_cellpose_cell_region_image_reports_exploratory_per_nucleus_endpoint():
    ch2 = np.array([[10, 10, 10], [10, 100, 100], [10, 100, 100]], dtype=np.uint16)
    nuclei = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 2]], dtype=np.uint32)
    cell_labels = np.array([[0, 0, 0], [0, 7, 7], [0, 7, 7]], dtype=np.uint32)

    row = measure_cellpose_cell_region_image(
        image_id="XY01",
        ch2_image=ch2,
        nuclei_mask=nuclei,
        cell_labels=cell_labels,
        background_value=10,
        source_id="XY01",
        mask_path=Path("masks/XY01.tif"),
        qc_panel_path=Path("qc/XY01.png"),
    )

    assert row["method"] == "cellpose_ch2_ch4_candidate_asma_associated_region"
    assert row["cellpose_object_count"] == 1
    assert row["dapi_positive_nucleus_count"] == 2
    assert row["normalization_denominator_count"] == 2
    assert row["nuclei_filtering_applied"] is False
    assert row["nuclei_filtering_policy"] == "none_count_nonzero_labels_in_supplied_mask"
    assert row["candidate_region_area_px"] == 4
    assert row["outside_candidate_region_area_px"] == 5
    assert row["target_integrated_raw_in_cellpose_region"] == 400.0
    assert row["target_integrated_background_corrected_in_cellpose_region"] == 360.0
    assert row["target_integrated_intensity_per_DAPI_positive_nucleus"] == 180.0
    assert row["target_integrated_intensity_per_cellpose_object"] == 360.0
    assert row["dapi_nuclei_with_centroid_inside_cellpose_region"] == 2
    assert row["cellpose_objects_with_dapi_centroid"] == 1
    assert row["cellpose_objects_without_dapi_centroid"] == 0
    assert row["cellpose_objects_with_multiple_dapi_centroids"] == 1
    assert row["background_method"] == "constant_value_10"
    assert row["excluded_region_integrated_background_corrected"] == 0.0
    assert row["excluded_region_background_corrected_fraction"] == 0.0
    assert "ch2_asma_used_as_candidate_cytoplasm_channel" in row["warnings"]


def test_measure_cellpose_cell_region_image_reports_dapi_anchored_variant():
    ch2 = np.array(
        [
            [10, 10, 10, 10, 10],
            [10, 100, 100, 10, 80],
            [10, 100, 100, 10, 80],
            [10, 10, 10, 200, 200],
            [10, 10, 10, 200, 200],
        ],
        dtype=np.uint16,
    )
    nuclei = np.array(
        [
            [0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 2, 0, 0],
            [0, 0, 0, 3, 0],
            [0, 0, 0, 0, 0],
        ],
        dtype=np.uint32,
    )
    cell_labels = np.array(
        [
            [0, 0, 0, 0, 0],
            [0, 7, 7, 0, 8],
            [0, 7, 7, 0, 8],
            [0, 0, 0, 9, 9],
            [0, 0, 0, 9, 9],
        ],
        dtype=np.uint32,
    )

    row = measure_cellpose_cell_region_image(
        image_id="XY01",
        ch2_image=ch2,
        nuclei_mask=nuclei,
        cell_labels=cell_labels,
        background_value=10,
        source_id="XY01",
    )

    assert row["cellpose_object_count"] == 3
    assert row["cellpose_objects_without_dapi_centroid"] == 1
    assert row["cellpose_objects_with_multiple_dapi_centroids"] == 1
    assert row["candidate_region_area_px"] == 10
    assert row["target_integrated_raw_in_cellpose_region"] == 1360.0

    assert row["dapi_anchored_cellpose_object_count"] == 2
    assert row["dapi_anchored_excluded_no_dapi_object_count"] == 1
    assert row["dapi_anchored_candidate_region_area_px"] == 8
    assert row["dapi_anchored_candidate_region_fraction"] == 8 / 25
    assert row["dapi_anchored_positive_area_per_DAPI_positive_nucleus"] == 8 / 3
    assert row["dapi_anchored_target_integrated_raw"] == 1200.0
    assert row["dapi_anchored_target_integrated_background_corrected"] == 1120.0
    assert row["dapi_anchored_target_integrated_intensity_per_DAPI_positive_nucleus"] == 1120.0 / 3


def test_score_cellpose_cell_region_qc_rejects_empty_or_near_full_field_masks():
    empty_status, empty_flags = score_cellpose_cell_region_qc(
        {
            "cellpose_object_count": 0,
            "candidate_region_fraction": 0.0,
            "dapi_positive_nucleus_count": 100,
            "dapi_nuclei_centroid_coverage_fraction": 0.0,
        }
    )
    full_status, full_flags = score_cellpose_cell_region_qc(
        {
            "cellpose_object_count": 10,
            "candidate_region_fraction": 0.92,
            "dapi_positive_nucleus_count": 100,
            "dapi_nuclei_centroid_coverage_fraction": 0.8,
        }
    )

    assert empty_status == "reject_qc_failure"
    assert "zero_cellpose_object_count" in empty_flags
    assert full_status == "reject_qc_failure"
    assert "near_full_field_cellpose_region_fraction" in full_flags


def test_score_cellpose_cell_region_qc_marks_object_mismatch_for_manual_review():
    status, flags = score_cellpose_cell_region_qc(
        {
            "cellpose_object_count": 10,
            "candidate_region_fraction": 0.25,
            "dapi_positive_nucleus_count": 10,
            "dapi_nuclei_centroid_coverage_fraction": 0.9,
            "cellpose_objects_without_dapi_centroid_fraction": 0.3,
            "cellpose_objects_with_multiple_dapi_centroids_fraction": 0.2,
            "excluded_region_background_corrected_fraction": 0.35,
        }
    )

    assert status == "needs_manual_review"
    assert "candidate_objects_without_dapi_centroid" in flags
    assert "candidate_objects_with_multiple_dapi_centroids" in flags
    assert "sizeable_background_corrected_ch2_outside_cellpose_region" in flags


def test_score_cellpose_cell_region_qc_rejects_majority_excluded_ch2():
    status, flags = score_cellpose_cell_region_qc(
        {
            "cellpose_object_count": 10,
            "candidate_region_fraction": 0.25,
            "dapi_positive_nucleus_count": 10,
            "dapi_nuclei_centroid_coverage_fraction": 0.9,
            "excluded_region_background_corrected_fraction": 0.7,
        }
    )

    assert status == "reject_qc_failure"
    assert "majority_background_corrected_ch2_outside_cellpose_region" in flags


def test_run_cellpose_cell_region_batch_writes_masks_qc_summary_and_logs(tmp_path: Path):
    dataset = tmp_path / "dataset"
    xy01 = dataset / "XY01"
    xy01.mkdir(parents=True)
    ch2 = np.zeros((24, 24), dtype=np.uint16) + 50
    ch2[6:18, 5:19] = 1000
    ch4 = np.zeros((24, 24), dtype=np.uint16)
    ch4[10:13, 9:12] = 2000
    ch4[10:13, 14:17] = 2000
    ch2_path = xy01 / "sample_XY01_CH2.tif"
    ch4_path = xy01 / "sample_XY01_CH4.tif"
    tifffile.imwrite(ch2_path, ch2)
    tifffile.imwrite(ch4_path, ch4)
    nuclei_mask = np.zeros((24, 24), dtype=np.uint32)
    nuclei_mask[10:13, 9:12] = 1
    nuclei_mask[10:13, 14:17] = 2
    nuclei_mask_path = tmp_path / "counts" / "masks" / "XY01_CH4_labels.tif"
    nuclei_mask_path.parent.mkdir(parents=True)
    tifffile.imwrite(nuclei_mask_path, nuclei_mask)

    def fake_segmenter(image: np.ndarray) -> np.ndarray:
        assert image.shape == (2, 24, 24)
        labels = np.zeros((24, 24), dtype=np.uint32)
        labels[6:18, 5:19] = 1
        return labels

    rows = run_cellpose_cell_region_batch(
        image_pairs=[ImagePair(location="XY01", source_id="XY01", ch2_path=ch2_path, ch4_path=ch4_path)],
        mask_lookup={"XY01": nuclei_mask_path},
        output_dir=tmp_path / "out",
        model_name="fake_model",
        segmenter=fake_segmenter,
        background_value=50,
    )

    assert rows[0]["cellpose_object_count"] == 1
    assert rows[0]["dapi_positive_nucleus_count"] == 2
    assert rows[0]["normalization_denominator_count"] == 2
    assert rows[0]["cellpose_objects_with_multiple_dapi_centroids"] == 1
    assert rows[0]["excluded_signal_check_path"]
    assert rows[0]["nuclei_mask_path"] == str(nuclei_mask_path)
    assert (tmp_path / "out" / "masks" / "XY01_cellpose_ch2_ch4_fake_model_labels.tif").exists()
    assert (tmp_path / "out" / "qc" / "XY01_cellpose_cell_region_qc.png").exists()
    assert (
        tmp_path / "out" / "qc" / "XY01_cellpose_cell_region_qc_excluded_signal_check.png"
    ).exists()
    assert (tmp_path / "out" / "qc_contact_sheet.png").exists()
    assert (tmp_path / "out" / "summaries" / "cellpose_cell_region_image_metrics.csv").exists()
    summary_text = (
        tmp_path / "out" / "summaries" / "cellpose_cell_region_image_metrics.csv"
    ).read_text(encoding="utf-8")
    assert "cellpose_object_count" in summary_text
    assert "cellpose_cell_count" not in summary_text
    assert "candidate_region_area_px" in summary_text
    assert "cell_region_area_px" not in summary_text
    assert "outside_candidate_region_area_px" in summary_text
    assert "non_cell_region_area_px" not in summary_text
    assert "nuclei_mask_path" in summary_text
    assert "normalization_denominator_count" in summary_text
    assert "excluded_signal_check_path" in summary_text
    config = yaml.safe_load((tmp_path / "out" / "logs" / "config_resolved.yaml").read_text())
    assert config["method"]["name"] == "cellpose_ch2_ch4_candidate_asma_associated_region"
    assert config["method"]["whole_cell_claim"] is False
    assert config["channel_extraction"]["ch2_channel"] == "CH2/aSMA candidate cytoplasm"
    assert config["validation_status"]["whole_cell_segmentation_validated"] is False
    assert config["measurement"]["background_method"] == "constant_value_50"
    assert (
        config["measurement"]["dapi_anchored_region_definition"]
        == "Cellpose objects retained only when at least one DAPI nucleus centroid falls inside the object"
    )


def test_run_cellpose_cell_region_batch_logs_explicit_channel_mapping(tmp_path: Path):
    dataset = tmp_path / "dataset"
    xy01 = dataset / "XY01"
    xy01.mkdir(parents=True)
    target = np.zeros((16, 16), dtype=np.uint16)
    target[4:12, 4:12] = 1000
    dapi = np.zeros((16, 16), dtype=np.uint16)
    dapi[7:9, 7:9] = 2000
    target_path = xy01 / "sample_XY01_CH1.tif"
    dapi_path = xy01 / "sample_XY01_CH4.tif"
    tifffile.imwrite(target_path, target)
    tifffile.imwrite(dapi_path, dapi)
    nuclei_mask = np.zeros((16, 16), dtype=np.uint32)
    nuclei_mask[7:9, 7:9] = 1
    nuclei_mask_path = tmp_path / "counts" / "masks" / "XY01_CH4_labels.tif"
    nuclei_mask_path.parent.mkdir(parents=True)
    tifffile.imwrite(nuclei_mask_path, nuclei_mask)

    def fake_segmenter(image: np.ndarray) -> np.ndarray:
        assert image.shape == (2, 16, 16)
        labels = np.zeros((16, 16), dtype=np.uint32)
        labels[4:12, 4:12] = 1
        return labels

    rows = run_cellpose_cell_region_batch(
        image_pairs=[
            ImagePair(
                location="XY01",
                source_id="XY01",
                ch2_path=target_path,
                ch4_path=dapi_path,
                target_channel_id="CH1",
                dapi_channel_id="CH4",
            )
        ],
        mask_lookup={"XY01": nuclei_mask_path},
        output_dir=tmp_path / "out",
        model_name="fake_model",
        segmenter=fake_segmenter,
        target_channel_id="CH1",
        dapi_channel_id="CH4",
        background_value=0,
        write_internal_qc=False,
    )

    assert rows[0]["target_channel_id"] == "CH1"
    assert rows[0]["dapi_channel_id"] == "CH4"
    assert rows[0]["target_path"] == str(target_path)
    assert rows[0]["dapi_path"] == str(dapi_path)
    assert (tmp_path / "out" / "masks" / "XY01_cellpose_ch1_ch4_fake_model_labels.tif").exists()
    config = yaml.safe_load((tmp_path / "out" / "logs" / "config_resolved.yaml").read_text())
    assert config["method"]["name"] == "cellpose_ch1_ch4_candidate_asma_associated_region"
    assert config["channel_extraction"]["target_channel_id"] == "CH1"
    assert config["channel_extraction"]["dapi_channel_id"] == "CH4"
    assert config["image_inputs"][0]["target_path"] == str(target_path)


def test_run_cellpose_cell_region_batch_can_skip_internal_qc_rendering(tmp_path: Path):
    dataset = tmp_path / "dataset"
    xy01 = dataset / "XY01"
    xy01.mkdir(parents=True)
    ch2 = np.zeros((12, 12), dtype=np.uint16)
    ch2[3:9, 3:9] = 1000
    ch4 = np.zeros((12, 12), dtype=np.uint16)
    ch4[5:7, 5:7] = 2000
    ch2_path = xy01 / "sample_XY01_CH2.tif"
    ch4_path = xy01 / "sample_XY01_CH4.tif"
    tifffile.imwrite(ch2_path, ch2)
    tifffile.imwrite(ch4_path, ch4)
    nuclei_mask = np.zeros((12, 12), dtype=np.uint32)
    nuclei_mask[5:7, 5:7] = 1
    nuclei_mask_path = tmp_path / "counts" / "masks" / "XY01_CH4_labels.tif"
    nuclei_mask_path.parent.mkdir(parents=True)
    tifffile.imwrite(nuclei_mask_path, nuclei_mask)

    def fake_segmenter(_image: np.ndarray) -> np.ndarray:
        labels = np.zeros((12, 12), dtype=np.uint32)
        labels[3:9, 3:9] = 1
        return labels

    rows = run_cellpose_cell_region_batch(
        image_pairs=[ImagePair(location="XY01", source_id="XY01", ch2_path=ch2_path, ch4_path=ch4_path)],
        mask_lookup={"XY01": nuclei_mask_path},
        output_dir=tmp_path / "out",
        model_name="fake_model",
        segmenter=fake_segmenter,
        write_internal_qc=False,
    )

    assert rows[0]["qc_panel_path"] == ""
    assert rows[0]["excluded_signal_check_path"] == ""
    assert (tmp_path / "out" / "masks" / "XY01_cellpose_ch2_ch4_fake_model_labels.tif").exists()
    assert not (tmp_path / "out" / "qc_contact_sheet.png").exists()
    config = yaml.safe_load((tmp_path / "out" / "logs" / "config_resolved.yaml").read_text())
    assert config["outputs"]["internal_qc_generated"] is False
