#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.manual_validation_pipeline import run_manual_validation_pipeline
from dapi_norm.manual_validation_report import parse_candidate_specs

app = typer.Typer(
    help=(
        "Run the manual-validation gate: audit annotation readiness first, "
        "then run quantitative candidate-vs-reference validation only if ready."
    )
)


@app.command()
def main(
    package_dir: Path = typer.Option(
        ...,
        "--package",
        help="Manual validation package directory.",
    ),
    candidate: list[str] = typer.Option(
        ...,
        "--candidate",
        help="Candidate method as name=mask_dir. Repeat for multiple methods.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output",
        help="Pipeline output directory.",
    ),
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
        outputs = run_manual_validation_pipeline(
            package_dir=package_dir,
            candidate_dirs=parse_candidate_specs(candidate),
            output_dir=output_dir,
            iou_threshold=iou_threshold,
            min_precision=min_precision,
            min_recall=min_recall,
            min_f1=min_f1,
            min_mean_iou=min_mean_iou,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"manual validation pipeline failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"validation_ready={outputs['validation_ready']}")
    typer.echo(f"wrote {outputs['gate_report']}")
    typer.echo(f"wrote {outputs['audit_report']}")
    if outputs["validation_ready"]:
        typer.echo(f"wrote {outputs['validation_report']}")
        typer.echo(f"wrote overlays under {outputs['overlay_dir']}")


if __name__ == "__main__":
    app()
