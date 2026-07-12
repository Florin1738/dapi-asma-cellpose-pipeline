#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.cellpose_endpoint_figures import render_cellpose_endpoint_figures

app = typer.Typer(help="Render Cellpose CH2+CH4 endpoint plots and visual overlay panels.")


@app.command()
def main(
    summary_csv: Path = typer.Option(
        ...,
        "--summary",
        help="Merged full-plate Cellpose endpoint summary CSV.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output",
        help="Directory for plots, captions, and overlay pages.",
    ),
    panel_page_size: int = typer.Option(
        12,
        "--panel-page-size",
        help="Number of image overlays per visual panel page.",
    ),
    selected_fields_per_plate: int = typer.Option(
        12,
        "--selected-fields-per-plate",
        help="Stratified fields per plate for the aligned endpoint bar plots.",
    ),
    max_overlay_images: int | None = typer.Option(
        None,
        "--max-overlay-images",
        help="Optional cap for overlay pages. Omit to render every field.",
    ),
) -> None:
    outputs = render_cellpose_endpoint_figures(
        summary_csv=summary_csv,
        output_dir=output_dir,
        panel_page_size=panel_page_size,
        selected_fields_per_plate=selected_fields_per_plate,
        max_overlay_images=max_overlay_images,
    )
    typer.echo(f"metric_contrast={outputs['metric_contrast']}")
    typer.echo(f"masking_effect={outputs['masking_effect']}")
    typer.echo(f"plate_summary={outputs['plate_summary']}")
    typer.echo(f"representative_cell_segmentation={outputs['representative_cell_segmentation']}")
    typer.echo(f"overlay_pages={len(outputs['overlay_pages'])}")
    typer.echo(f"captions={outputs['captions_markdown']}")


if __name__ == "__main__":
    app()
