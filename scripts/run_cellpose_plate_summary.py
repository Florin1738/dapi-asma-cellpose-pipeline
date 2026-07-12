#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.cellpose_plate_summary import (
    build_cellpose_plate_summary,
    write_cellpose_plate_summary_csv,
    write_cellpose_plate_summary_markdown,
)

app = typer.Typer(help="Merge per-acquisition Cellpose CH2+CH4 outputs into plate-level tables.")


@app.command()
def main(
    runs_root: Path = typer.Option(..., "--runs-root", help="Root containing Plate_*/Run*/ Cellpose outputs."),
    output_dir: Path = typer.Option(..., "--output", help="Directory for merged CSV/Markdown tables."),
) -> None:
    rows = build_cellpose_plate_summary(runs_root)
    csv_path = write_cellpose_plate_summary_csv(
        output_dir / "cellpose_full_plate_endpoint_summary.csv",
        rows,
    )
    md_path = write_cellpose_plate_summary_markdown(
        output_dir / "cellpose_full_plate_endpoint_summary.md",
        rows,
    )
    typer.echo(f"rows={len(rows)}")
    typer.echo(f"csv={csv_path}")
    typer.echo(f"markdown={md_path}")


if __name__ == "__main__":
    app()
