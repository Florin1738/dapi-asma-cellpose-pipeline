#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.cellpose_endpoint_figures import render_representative_cell_segmentation_panel

app = typer.Typer(
    help="Render a representative Cellpose retained-region panel without rebuilding all QC pages."
)


@app.command()
def main(
    summary_csv: Path = typer.Option(
        ...,
        "--summary",
        help="Merged full-plate Cellpose endpoint summary CSV.",
    ),
    output_path: Path = typer.Option(
        ...,
        "--output",
        help="Output PNG path for the representative panel.",
    ),
    fields_per_plate: int = typer.Option(
        4,
        "--fields-per-plate",
        help="Number of low-to-high endpoint examples to select per plate.",
    ),
) -> None:
    panel = render_representative_cell_segmentation_panel(
        summary_csv=summary_csv,
        output_path=output_path,
        per_plate=fields_per_plate,
    )
    typer.echo(f"representative_cell_segmentation={panel}")


if __name__ == "__main__":
    app()
