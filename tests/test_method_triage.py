from __future__ import annotations

import csv
from pathlib import Path

from dapi_norm.method_triage import (
    _triage_plot_style,
    build_method_triage_rows,
    triage_region_restricted_status,
    write_method_triage_outputs,
)


def test_triage_region_restricted_status_never_accepts_without_manual_validation():
    assert (
        triage_region_restricted_status(["reviewable_not_validated"], manual_validation_available=False)
        == "not_validated_manual_validation_required"
    )
    assert (
        triage_region_restricted_status(
            ["reviewable_not_validated", "reject_qc_failure"],
            manual_validation_available=False,
        )
        == "mixed_region_restricted_qc_rejection_not_validated"
    )
    assert (
        triage_region_restricted_status(
            ["reviewable_not_validated", "needs_manual_review"],
            manual_validation_available=False,
        )
        == "manual_review_flagged_not_validated"
    )
    assert (
        triage_region_restricted_status(["reject_qc_failure"], manual_validation_available=False)
        == "all_region_restricted_methods_rejected"
    )
    assert (
        triage_region_restricted_status([], manual_validation_available=False)
        == "region_restricted_sources_missing"
    )
    assert (
        triage_region_restricted_status(["reviewable_not_validated"], manual_validation_available=True)
        == "manual_validation_available_review_required"
    )


def test_triage_plot_style_labels_qc_statuses_for_standalone_png():
    assert (
        _triage_plot_style("not_validated_manual_validation_required")["annotation"]
        == "not\nvalidated"
    )
    assert (
        _triage_plot_style("mixed_region_restricted_qc_rejection_not_validated")[
            "annotation"
        ]
        == "mixed\nQC reject"
    )
    assert (
        _triage_plot_style("manual_review_flagged_not_validated")["annotation"]
        == "manual\nreview"
    )
    assert (
        _triage_plot_style("all_region_restricted_methods_rejected")["annotation"]
        == "all\nrejected"
    )


def test_build_method_triage_rows_combines_whole_field_seeded_and_cellpose_metrics(
    tmp_path: Path,
):
    pi = tmp_path / "pi.csv"
    pi.write_text(
        "\n".join(
            [
                "loc,raw,nuclei,raw_per_nuc",
                "XY22,1000,10,100",
                "XY23,600,10,60",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    seeded = tmp_path / "seeded.csv"
    seeded.write_text(
        "\n".join(
            [
                "image_id,seeded_region_intensity_per_DAPI_positive_nucleus,"
                "seeded_region_fraction,unseeded_foreground_fraction,qc_status,qc_flags",
                "XY22,80,0.5,0.02,reviewable_not_validated,not_validated_whole_cell_mask",
                "XY23,10,0.1,0.2,needs_manual_review,sizeable_unseeded_target_foreground",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    cellpose = tmp_path / "cellpose.csv"
    cellpose.write_text(
        "\n".join(
            [
                "image_id,target_integrated_intensity_per_DAPI_positive_nucleus,"
                "candidate_region_fraction,excluded_region_background_corrected_fraction,"
                "qc_status,qc_flags",
                "XY22,90,0.6,0.1,reviewable_not_validated,not_validated_whole_cell_mask",
                "XY23,20,0.2,0.7,reject_qc_failure,majority_background_corrected_ch2_outside_cellpose_region",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    propagation = tmp_path / "propagation.csv"
    propagation.write_text(
        "\n".join(
            [
                "image_id,seeded_region_intensity_per_DAPI_positive_nucleus,"
                "seeded_region_fraction,unseeded_foreground_fraction,qc_status,qc_flags",
                "XY22,85,0.55,0.03,reviewable_not_validated,not_validated_whole_cell_mask",
                "XY23,15,0.12,0.18,needs_manual_review,sizeable_unseeded_target_foreground",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = build_method_triage_rows(
        pi_metrics_path=pi,
        seeded_watershed_path=seeded,
        seeded_random_walker_path=seeded,
        seeded_propagation_path=propagation,
        cellpose_paths=[cellpose],
    )

    assert [row["image_id"] for row in rows] == ["XY22", "XY23"]
    xy22 = rows[0]
    assert xy22["whole_field_raw_per_nucleus"] == 100.0
    assert xy22["seeded_watershed_per_nucleus"] == 80.0
    assert xy22["seeded_propagation_per_nucleus"] == 85.0
    assert xy22["cellpose_candidate_per_nucleus"] == 90.0
    assert xy22["accepted_region_restricted_method"] is False
    assert xy22["triage_status"] == "not_validated_manual_validation_required"
    xy23 = rows[1]
    assert xy23["cellpose_qc_status"] == "reject_qc_failure"
    assert xy23["triage_status"] == "mixed_region_restricted_qc_rejection_not_validated"
    assert xy23["any_region_restricted_qc_reject"] is True


def test_build_method_triage_rows_can_filter_positions_and_report_missing_sources(
    tmp_path: Path,
):
    pi = tmp_path / "pi.csv"
    pi.write_text(
        "\n".join(
            [
                "loc,raw,nuclei,raw_per_nuc",
                "XY22,1000,10,100",
                "XY23,600,10,60",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    seeded = tmp_path / "seeded.csv"
    seeded.write_text(
        "\n".join(
            [
                "image_id,seeded_region_intensity_per_DAPI_positive_nucleus,"
                "seeded_region_fraction,unseeded_foreground_fraction,qc_status,qc_flags",
                "XY22,80,0.5,0.02,reviewable_not_validated,not_validated_whole_cell_mask",
                "XY23,10,0.1,0.2,needs_manual_review,sizeable_unseeded_target_foreground",
                "XY95,5,0.05,0.3,reject_qc_failure,extra_position",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    cellpose = tmp_path / "cellpose.csv"
    cellpose.write_text(
        "\n".join(
            [
                "image_id,target_integrated_intensity_per_DAPI_positive_nucleus,"
                "candidate_region_fraction,excluded_region_background_corrected_fraction,"
                "qc_status,qc_flags",
                "XY22,90,0.6,0.1,reviewable_not_validated,not_validated_whole_cell_mask",
                "XY40,2,0.001,0.95,reject_qc_failure,majority_background_corrected_ch2_outside_cellpose_region",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = build_method_triage_rows(
        pi_metrics_path=pi,
        seeded_watershed_path=seeded,
        seeded_random_walker_path=seeded,
        seeded_propagation_path=None,
        cellpose_paths=[cellpose],
        positions=["XY22", "XY23", "XY40"],
    )

    assert [row["image_id"] for row in rows] == ["XY22", "XY23", "XY40"]
    assert rows[0]["missing_sources"] == ""
    assert rows[1]["missing_sources"] == "cellpose"
    assert rows[2]["missing_sources"] == "whole_field;seeded_watershed;seeded_random_walker"


def test_build_method_triage_rows_treats_missing_method_files_as_missing_sources(
    tmp_path: Path,
):
    pi = tmp_path / "pi.csv"
    pi.write_text(
        "\n".join(
            [
                "loc,raw,nuclei,raw_per_nuc",
                "XY22,1000,10,100",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    seeded = tmp_path / "seeded.csv"
    seeded.write_text(
        "\n".join(
            [
                "image_id,seeded_region_intensity_per_DAPI_positive_nucleus,"
                "seeded_region_fraction,unseeded_foreground_fraction,qc_status,qc_flags",
                "XY22,80,0.5,0.02,reviewable_not_validated,not_validated_whole_cell_mask",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    rows = build_method_triage_rows(
        pi_metrics_path=pi,
        seeded_watershed_path=seeded,
        seeded_random_walker_path=tmp_path / "missing_random_walker.csv",
        seeded_propagation_path=tmp_path / "missing_propagation.csv",
        cellpose_paths=[tmp_path / "missing_cellpose.csv"],
        positions=["XY22"],
    )

    assert rows[0]["missing_sources"] == "seeded_random_walker;seeded_propagation;cellpose"
    assert rows[0]["triage_status"] == "not_validated_manual_validation_required"
    assert rows[0]["any_region_restricted_qc_reject"] is False


def test_write_method_triage_outputs_writes_csv_markdown_and_plot(tmp_path: Path):
    rows = [
        {
            "image_id": "XY22",
            "whole_field_raw_per_nucleus": 100.0,
            "seeded_watershed_per_nucleus": 80.0,
            "seeded_random_walker_per_nucleus": 80.0,
            "seeded_propagation_per_nucleus": 85.0,
            "cellpose_candidate_per_nucleus": 90.0,
            "triage_status": "not_validated_manual_validation_required",
            "accepted_region_restricted_method": False,
            "any_region_restricted_qc_reject": False,
            "cellpose_qc_status": "reviewable_not_validated",
            "seeded_watershed_qc_status": "reviewable_not_validated",
            "seeded_random_walker_qc_status": "reviewable_not_validated",
            "seeded_propagation_qc_status": "reviewable_not_validated",
            "missing_sources": "",
        },
        {
            "image_id": "XY40",
            "whole_field_raw_per_nucleus": float("nan"),
            "seeded_watershed_per_nucleus": float("nan"),
            "seeded_random_walker_per_nucleus": float("nan"),
            "seeded_propagation_per_nucleus": float("nan"),
            "cellpose_candidate_per_nucleus": 2.0,
            "triage_status": "all_region_restricted_methods_rejected",
            "accepted_region_restricted_method": False,
            "any_region_restricted_qc_reject": True,
            "cellpose_qc_status": "reject_qc_failure",
            "seeded_watershed_qc_status": "",
            "seeded_random_walker_qc_status": "",
            "seeded_propagation_qc_status": "",
            "missing_sources": "whole_field;seeded_watershed;seeded_random_walker;seeded_propagation",
        }
    ]

    outputs = write_method_triage_outputs(rows, tmp_path / "out")

    assert outputs["csv"].exists()
    assert outputs["report"].exists()
    assert outputs["plot"].exists()
    with outputs["csv"].open(newline="", encoding="utf-8") as handle:
        csv_rows = list(csv.DictReader(handle))
    assert csv_rows[0]["triage_status"] == "not_validated_manual_validation_required"
    xy40 = csv_rows[1]
    assert xy40["whole_field_raw_per_nucleus"] == ""
    assert xy40["seeded_watershed_per_nucleus"] == ""
    assert xy40["seeded_random_walker_per_nucleus"] == ""
    assert xy40["seeded_propagation_per_nucleus"] == ""
    assert xy40["cellpose_candidate_per_nucleus"] == "2.0"
    report = outputs["report"].read_text()
    assert "No region-restricted method is accepted" in report
    assert "| XY40 | NA | NA | NA | NA | 2.000e+00 | yes |" in report
