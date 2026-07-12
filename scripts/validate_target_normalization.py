#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.target_validation import validate_target_outputs

app = typer.Typer(help="Validate target per-DAPI-positive-nucleus output artifacts.")


@app.command()
def main(
    output_dir: Path = typer.Option(
        ...,
        "--output",
        help="Target per-DAPI-positive-nucleus output directory.",
    ),
) -> None:
    result = validate_target_outputs(output_dir)
    typer.echo(f"summary_rows={result['summary_rows']}")
    typer.echo(f"plots_exist={result['plots_exist']}")
    typer.echo(f"qc_overlays_exist={result['qc_overlays_exist']}")
    typer.echo(f"formulas_match={result['formulas_match']}")


if __name__ == "__main__":
    app()
