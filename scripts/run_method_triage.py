#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.method_triage import build_method_triage_rows, write_method_triage_outputs

app = typer.Typer(
    help="Build an exploratory cross-method triage report for aSMA region-restricted methods."
)


@app.command()
def main(
    pi_metrics: Path | None = typer.Option(
        None,
        "--pi-metrics",
        help="PI/whole-field metrics CSV with loc, raw, nuclei, raw_per_nuc columns.",
    ),
    seeded_watershed: Path | None = typer.Option(
        None,
        "--seeded-watershed",
        help="Seeded watershed image metrics CSV.",
    ),
    seeded_random_walker: Path | None = typer.Option(
        None,
        "--seeded-random-walker",
        help="Seeded random-walker image metrics CSV.",
    ),
    seeded_propagation: Path | None = typer.Option(
        None,
        "--seeded-propagation",
        help="CellProfiler-style seeded propagation image metrics CSV.",
    ),
    cellpose: list[Path] | None = typer.Option(
        None,
        "--cellpose",
        help="Cellpose candidate-region metrics CSV. Pass multiple times to merge runs.",
    ),
    output_dir: Path = typer.Option(..., "--output", help="Output directory."),
    positions: str | None = typer.Option(
        None,
        "--positions",
        help="Optional comma-separated image positions to include, e.g. XY22,XY23,XY24.",
    ),
    manual_validation_available: bool = typer.Option(
        False,
        "--manual-validation-available/--no-manual-validation-available",
        help=(
            "Set only when manual ground truth is available for review; "
            "this never auto-accepts a method."
        ),
    ),
) -> None:
    rows = build_method_triage_rows(
        pi_metrics_path=pi_metrics,
        seeded_watershed_path=seeded_watershed,
        seeded_random_walker_path=seeded_random_walker,
        seeded_propagation_path=seeded_propagation,
        cellpose_paths=cellpose,
        manual_validation_available=manual_validation_available,
        positions=_parse_positions(positions),
    )
    outputs = write_method_triage_outputs(rows, output_dir)
    typer.echo(f"wrote {outputs['csv']}")
    typer.echo(f"wrote {outputs['report']}")
    typer.echo(f"wrote {outputs['plot']}")


def _parse_positions(positions: str | None) -> list[str] | None:
    if positions is None:
        return None
    parsed = [position.strip() for position in positions.split(",") if position.strip()]
    return parsed or None


if __name__ == "__main__":
    app()
