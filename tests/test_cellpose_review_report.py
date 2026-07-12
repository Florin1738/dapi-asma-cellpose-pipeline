from __future__ import annotations

import csv
from pathlib import Path

import pytest
from typer.testing import CliRunner

from dapi_norm.cellpose_review_report import write_cellpose_review_report
from scripts.render_cellpose_review_report import app


def test_write_cellpose_review_report_summarizes_statuses_and_links_panels(tmp_path: Path):
    run_dir = tmp_path / "cellpose_run"
    metrics_path = run_dir / "summaries" / "cellpose_cell_region_image_metrics.csv"
    qc_dir = run_dir / "qc"
    qc_dir.mkdir(parents=True)
    for image_id in ["XY22", "XY40"]:
        (qc_dir / f"{image_id}_cellpose_cell_region_qc.png").write_bytes(b"png")
        (qc_dir / f"{image_id}_cellpose_cell_region_qc_excluded_signal_check.png").write_bytes(b"png")
    _write_csv(
        metrics_path,
        [
            _metrics_row(
                image_id="XY22",
                qc_status="reviewable_not_validated",
                qc_flags="not_validated_whole_cell_mask",
                per_nucleus="1.5e7",
                region_fraction="0.61",
                dapi_coverage="0.95",
                excluded_fraction="0.17",
                panel_path=str(qc_dir / "XY22_cellpose_cell_region_qc.png"),
                excluded_path=str(qc_dir / "XY22_cellpose_cell_region_qc_excluded_signal_check.png"),
            ),
            _metrics_row(
                image_id="XY40",
                qc_status="reject_qc_failure",
                qc_flags="not_validated_whole_cell_mask;very_low_cellpose_region_fraction_with_dapi_nuclei",
                per_nucleus="3.5e6",
                region_fraction="0.0008",
                dapi_coverage="0.67",
                excluded_fraction="0.99",
                panel_path=str(qc_dir / "XY40_cellpose_cell_region_qc.png"),
                excluded_path=str(qc_dir / "XY40_cellpose_cell_region_qc_excluded_signal_check.png"),
            ),
        ],
    )

    outputs = write_cellpose_review_report(
        metrics_csv=metrics_path,
        output_dir=run_dir / "review_report",
        title="Plate 1 selected-15 Cellpose candidate review",
    )

    html = outputs["html"].read_text(encoding="utf-8")
    summary = outputs["summary"].read_text(encoding="utf-8")
    assert "Plate 1 selected-15 Cellpose candidate review" in html
    assert "Qualitative QC only" in html
    assert "reject_qc_failure" in html
    assert "reviewable_not_validated" in html
    assert "XY40_cellpose_cell_region_qc.png" in html
    assert "XY40_cellpose_cell_region_qc_excluded_signal_check.png" in html
    assert html.index("XY40") < html.index("XY22")
    assert "reject_qc_failure: 1" in summary
    assert "reviewable_not_validated: 1" in summary
    assert "No precision, recall, F1, or IoU" in summary


def test_write_cellpose_review_report_rejects_missing_qc_panel(tmp_path: Path):
    metrics_path = tmp_path / "metrics.csv"
    _write_csv(
        metrics_path,
        [
            _metrics_row(
                image_id="XY22",
                qc_status="reviewable_not_validated",
                qc_flags="not_validated_whole_cell_mask",
                per_nucleus="1.5e7",
                region_fraction="0.61",
                dapi_coverage="0.95",
                excluded_fraction="0.17",
                panel_path=str(tmp_path / "missing.png"),
                excluded_path=str(tmp_path / "missing_excluded.png"),
            )
        ],
    )

    with pytest.raises(FileNotFoundError, match="QC panel missing for XY22"):
        write_cellpose_review_report(
            metrics_csv=metrics_path,
            output_dir=tmp_path / "review_report",
        )


def test_render_cellpose_review_report_cli(tmp_path: Path):
    runner = CliRunner()
    run_dir = tmp_path / "cellpose_run"
    metrics_path = run_dir / "summaries" / "cellpose_cell_region_image_metrics.csv"
    qc_dir = run_dir / "qc"
    qc_dir.mkdir(parents=True)
    (qc_dir / "XY22_cellpose_cell_region_qc.png").write_bytes(b"png")
    (qc_dir / "XY22_cellpose_cell_region_qc_excluded_signal_check.png").write_bytes(b"png")
    _write_csv(
        metrics_path,
        [
            _metrics_row(
                image_id="XY22",
                qc_status="reviewable_not_validated",
                qc_flags="not_validated_whole_cell_mask",
                per_nucleus="1.5e7",
                region_fraction="0.61",
                dapi_coverage="0.95",
                excluded_fraction="0.17",
                panel_path=str(qc_dir / "XY22_cellpose_cell_region_qc.png"),
                excluded_path=str(qc_dir / "XY22_cellpose_cell_region_qc_excluded_signal_check.png"),
            )
        ],
    )

    result = runner.invoke(
        app,
        [
            "--metrics",
            str(metrics_path),
            "--output",
            str(run_dir / "review_report"),
            "--title",
            "CLI review",
        ],
    )

    assert result.exit_code == 0
    assert "wrote" in result.stdout
    assert (run_dir / "review_report" / "index.html").exists()
    assert (run_dir / "review_report" / "README.md").exists()


def _metrics_row(
    *,
    image_id: str,
    qc_status: str,
    qc_flags: str,
    per_nucleus: str,
    region_fraction: str,
    dapi_coverage: str,
    excluded_fraction: str,
    panel_path: str,
    excluded_path: str,
) -> dict[str, str]:
    return {
        "image_id": image_id,
        "source_id": image_id,
        "qc_status": qc_status,
        "qc_flags": qc_flags,
        "dapi_positive_nucleus_count": "10",
        "cellpose_object_count": "9",
        "candidate_region_fraction": region_fraction,
        "dapi_nuclei_centroid_coverage_fraction": dapi_coverage,
        "excluded_region_background_corrected_fraction": excluded_fraction,
        "target_integrated_intensity_per_DAPI_positive_nucleus": per_nucleus,
        "qc_panel_path": panel_path,
        "excluded_signal_check_path": excluded_path,
    }


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
