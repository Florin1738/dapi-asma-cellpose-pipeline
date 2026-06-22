#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.cellpose_runner import run_nuclei_count_batch

app = typer.Typer(help="Run Cellpose nuclei counting on candidate DAPI channel images.")


@app.command()
def main(
    input_root: Path = typer.Option(..., "--input", help="Dataset root to process."),
    output_dir: Path = typer.Option(..., "--output", help="Output directory."),
    channel_id: str = typer.Option("CH4", "--channel", help="Candidate DAPI channel ID."),
    model_name: str = typer.Option("cpsam_v2", "--model", help="Cellpose model name."),
    max_images: int | None = typer.Option(None, "--max-images", help="Limit images for pilot runs."),
    gpu: bool = typer.Option(True, "--gpu/--cpu", help="Request GPU/MPS acceleration."),
) -> None:
    summaries = run_nuclei_count_batch(
        input_root=input_root,
        output_dir=output_dir,
        channel_id=channel_id,
        model_name=model_name,
        gpu=gpu,
        max_images=max_images,
    )
    for summary in summaries:
        typer.echo(f"{summary['image_id']}: nucleus_count={summary['nucleus_count']}")


if __name__ == "__main__":
    app()

