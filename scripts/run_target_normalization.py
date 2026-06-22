#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.target_normalization import run_target_normalization

app = typer.Typer(help="Measure candidate target intensity normalized by Cellpose nuclei counts.")


@app.command()
def main(
    input_root: Path = typer.Option(..., "--input", help="Dataset root to process."),
    counts_dir: Path = typer.Option(..., "--counts", help="Cellpose count output directory."),
    output_dir: Path = typer.Option(..., "--output", help="Target-normalization output directory."),
    target_channel_id: str = typer.Option("CH2", "--target-channel", help="Candidate target channel."),
    dapi_channel_id: str = typer.Option("CH4", "--dapi-channel", help="Candidate DAPI channel."),
    background_percentile: float = typer.Option(
        10, "--background-percentile", help="Percentile background estimate."
    ),
) -> None:
    rows = run_target_normalization(
        input_root=input_root,
        counts_dir=counts_dir,
        output_dir=output_dir,
        target_channel_id=target_channel_id,
        dapi_channel_id=dapi_channel_id,
        background_percentile=background_percentile,
    )
    for row in rows:
        typer.echo(
            f"{row['well_id']}: "
            f"normalized={row['target_integrated_intensity_per_DAPI_positive_nucleus']:.3f} "
            f"nuclei={row['filtered_nucleus_count']}"
        )


if __name__ == "__main__":
    app()
