#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.segmentation_validation import validate_seeded_region_outputs

app = typer.Typer(help="Validate seeded aSMA-region output artifacts and formulas.")


@app.command()
def main(
    output_dir: Path = typer.Option(..., "--output", help="Seeded-region output directory."),
) -> None:
    result = validate_seeded_region_outputs(output_dir)
    for key, value in result.items():
        typer.echo(f"{key}: {value}")


if __name__ == "__main__":
    app()
