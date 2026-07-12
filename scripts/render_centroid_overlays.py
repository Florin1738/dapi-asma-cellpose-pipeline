#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.centroid_overlays import render_centroid_overlays

app = typer.Typer(help="Render DAPI-only centroid overlays from Cellpose per-nucleus CSVs.")


@app.command()
def main(
    counts_root: Path = typer.Option(
        Path("output/pi_simple_summary/cellpose_counts"),
        "--counts",
        help="Root containing per-run Cellpose outputs.",
    ),
    output_dir: Path = typer.Option(
        Path("output/pi_simple_summary/centroid_overlays"),
        "--output",
        help="Directory for green-centroid DAPI overlays.",
    ),
    contact_sheet_limit: int = typer.Option(
        12,
        "--contact-sheet-limit",
        min=1,
        help="Maximum overlays per run contact sheet.",
    ),
) -> None:
    outputs = render_centroid_overlays(
        counts_root=counts_root,
        output_dir=output_dir,
        contact_sheet_limit=contact_sheet_limit,
    )
    typer.echo(f"overlays={outputs['overlay_count']}")
    typer.echo(f"contact_sheets={outputs['contact_sheet_count']}")


if __name__ == "__main__":
    app()
