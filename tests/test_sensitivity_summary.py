from __future__ import annotations

from pathlib import Path

from dapi_norm.sensitivity_summary import (
    build_sensitivity_rows,
    summarize_image_stability,
    summarize_runs,
    write_sensitivity_outputs,
)


def _write_seeded_csv(path: Path, rows: list[tuple[str, float, float, str]]) -> None:
    path.write_text(
        "\n".join(
            [
                "image_id,method,foreground_method,seeded_region_intensity_per_DAPI_positive_nucleus,"
                "seeded_region_fraction,unseeded_foreground_fraction,qc_status,qc_flags",
                *[
                    f"{image_id},seeded_intensity_propagation,{run_method},{value},0.1,0.0,{qc},flags"
                    for image_id, value, _region_fraction, qc in rows
                    for run_method in ["otsu"]
                ],
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_build_sensitivity_rows_reads_labeled_seeded_metric_tables(tmp_path: Path):
    table = tmp_path / "otsu.csv"
    _write_seeded_csv(
        table,
        [
            ("XY22", 30.0, 0.2, "reviewable_not_validated"),
            ("XY23", 10.0, 0.1, "needs_manual_review"),
            ("XY40", 0.0, 0.0, "reject_qc_failure"),
        ],
    )

    rows = build_sensitivity_rows(
        table_specs=[("prop_otsu", table)],
        positions=["XY22", "XY23", "XY40"],
    )

    assert [row["run_id"] for row in rows] == ["prop_otsu", "prop_otsu", "prop_otsu"]
    assert rows[0]["image_id"] == "XY22"
    assert rows[0]["per_nucleus"] == 30.0
    assert rows[0]["qc_status"] == "reviewable_not_validated"


def test_summarize_runs_flags_ordering_challenge_rejection_and_manual_review(tmp_path: Path):
    table_a = tmp_path / "a.csv"
    table_b = tmp_path / "b.csv"
    _write_seeded_csv(
        table_a,
        [
            ("XY22", 30.0, 0.2, "reviewable_not_validated"),
            ("XY23", 10.0, 0.1, "needs_manual_review"),
            ("XY24", 5.0, 0.1, "reviewable_not_validated"),
            ("XY40", 0.0, 0.0, "reject_qc_failure"),
        ],
    )
    _write_seeded_csv(
        table_b,
        [
            ("XY22", 5.0, 0.2, "reviewable_not_validated"),
            ("XY23", 10.0, 0.1, "reviewable_not_validated"),
            ("XY24", 30.0, 0.1, "reviewable_not_validated"),
            ("XY40", 2.0, 0.0, "reviewable_not_validated"),
        ],
    )
    rows = build_sensitivity_rows(table_specs=[("good", table_a), ("bad", table_b)])

    summary = summarize_runs(
        rows,
        ordered_positions=["XY22", "XY23", "XY24"],
        challenge_positions=["XY40"],
    )

    good = next(row for row in summary if row["run_id"] == "good")
    bad = next(row for row in summary if row["run_id"] == "bad")
    assert good["expected_order_preserved"] is True
    assert good["challenge_all_rejected"] is True
    assert good["challenge_all_zero_per_nucleus"] is True
    assert good["manual_review_present"] is True
    assert bad["expected_order_preserved"] is False
    assert bad["challenge_all_rejected"] is False
    assert bad["challenge_all_zero_per_nucleus"] is False


def test_summarize_image_stability_reports_range_and_cv():
    rows = [
        {"run_id": "a", "image_id": "XY22", "per_nucleus": 10.0, "qc_status": "reviewable"},
        {"run_id": "b", "image_id": "XY22", "per_nucleus": 20.0, "qc_status": "reviewable"},
        {"run_id": "a", "image_id": "XY40", "per_nucleus": 0.0, "qc_status": "reject_qc_failure"},
    ]

    summary = summarize_image_stability(rows)

    xy22 = next(row for row in summary if row["image_id"] == "XY22")
    assert xy22["n_runs"] == 2
    assert xy22["min_per_nucleus"] == 10.0
    assert xy22["max_per_nucleus"] == 20.0
    assert xy22["coefficient_of_variation"] > 0
    xy40 = next(row for row in summary if row["image_id"] == "XY40")
    assert xy40["all_runs_rejected"] is True


def test_write_sensitivity_outputs_writes_csvs_report_and_plot(tmp_path: Path):
    rows = [
        {
            "run_id": "prop_otsu",
            "image_id": "XY22",
            "per_nucleus": 30.0,
            "region_fraction": 0.2,
            "unseeded_foreground_fraction": 0.0,
            "qc_status": "reviewable_not_validated",
            "qc_flags": "not_validated_whole_cell_mask",
        },
        {
            "run_id": "prop_otsu",
            "image_id": "XY40",
            "per_nucleus": 0.0,
            "region_fraction": 0.0,
            "unseeded_foreground_fraction": 0.0,
            "qc_status": "reject_qc_failure",
            "qc_flags": "low_nucleus_count",
        },
    ]

    outputs = write_sensitivity_outputs(
        rows,
        output_dir=tmp_path / "out",
        ordered_positions=["XY22", "XY23", "XY24"],
        challenge_positions=["XY40"],
    )

    assert outputs["long_csv"].exists()
    assert outputs["run_summary_csv"].exists()
    assert outputs["image_summary_csv"].exists()
    assert outputs["report"].exists()
    assert outputs["plot"].exists()
    assert "exploratory robustness diagnostic" in outputs["report"].read_text(encoding="utf-8")
