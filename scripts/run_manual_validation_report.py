#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.manual_validation_report import (
    parse_candidate_specs,
    run_manual_validation_report,
)

app = typer.Typer(
    help=(
        "Run manual/reference mask validation for one or more candidate methods, "
        "aggregate method-level metrics, and render candidate-vs-reference overlays."
    )
)


@app.command()
def main(
    candidate: list[str] = typer.Option(
        ...,
        "--candidate",
        help="Candidate method as name=mask_dir. Repeat for multiple methods.",
    ),
    reference_dir: Path = typer.Option(
        ...,
        "--reference-dir",
        help="Directory of completed manual/reference label masks.",
    ),
    completion_status: Path = typer.Option(
        ...,
        "--completion-status",
        help="manual_labeling_status.csv for the reference package.",
    ),
    manifest: Path = typer.Option(
        ...,
        "--manifest",
        help="manual_validation_manifest.csv with CH2 paths used in required overlays.",
    ),
    output_dir: Path = typer.Option(..., "--output", help="Validation report output directory."),
    iou_threshold: float = typer.Option(
        0.5,
        "--iou-threshold",
        min=1e-12,
        max=1.0,
        help="Object match IoU threshold.",
    ),
    min_precision: float = typer.Option(0.8, "--min-precision", min=0.0, max=1.0),
    min_recall: float = typer.Option(0.8, "--min-recall", min=0.0, max=1.0),
    min_f1: float = typer.Option(0.8, "--min-f1", min=0.0, max=1.0),
    min_mean_iou: float = typer.Option(0.5, "--min-mean-iou", min=0.0, max=1.0),
) -> None:
    try:
        outputs = run_manual_validation_report(
            candidate_dirs=parse_candidate_specs(candidate),
            reference_dir=reference_dir,
            completion_status_path=completion_status,
            manifest_path=manifest,
            output_dir=output_dir,
            iou_threshold=iou_threshold,
            min_precision=min_precision,
            min_recall=min_recall,
            min_f1=min_f1,
            min_mean_iou=min_mean_iou,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"manual validation report failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"wrote {outputs['method_summary']}")
    typer.echo(f"wrote {outputs['image_summary']}")
    typer.echo(f"wrote {outputs['report']}")
    typer.echo(f"wrote overlays under {outputs['overlay_dir']}")


if __name__ == "__main__":
    app()
