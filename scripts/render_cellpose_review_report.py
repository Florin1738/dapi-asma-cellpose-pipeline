#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.cellpose_review_report import write_cellpose_review_report

app = typer.Typer(
    help="Render an HTML visual QC report for Cellpose CH2+CH4 candidate-region outputs."
)


@app.command()
def main(
    metrics_csv: Path = typer.Option(
        ...,
        "--metrics",
        help="Cellpose candidate-region image metrics CSV.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output",
        help="Output directory for index.html and README.md.",
    ),
    title: str = typer.Option(
        "Cellpose CH2+CH4 Candidate Region Visual Review",
        "--title",
        help="Report title.",
    ),
) -> None:
    try:
        outputs = write_cellpose_review_report(
            metrics_csv=metrics_csv,
            output_dir=output_dir,
            title=title,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"cellpose review report failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"wrote {outputs['html']}")
    typer.echo(f"wrote {outputs['summary']}")


if __name__ == "__main__":
    app()
