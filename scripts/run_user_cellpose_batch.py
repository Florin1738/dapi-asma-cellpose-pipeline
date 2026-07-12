#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.user_cellpose_batch import discover_acquisitions, run_user_cellpose_batch

app = typer.Typer(
    help=(
        "Run the low-friction Cellpose DAPI nuclei plus target/aSMA Cellpose region pipeline "
        "over one folder of plate/acquisition data."
    )
)


@app.command()
def main(
    input_root: Path = typer.Option(
        ...,
        "--input",
        help="Folder containing plate/acquisition data with XY## folders and channel TIFFs.",
    ),
    output_root: Path = typer.Option(
        ...,
        "--output",
        help="New output folder for this complete run.",
    ),
    model_name: str = typer.Option(
        "cpsam_v2",
        "--model",
        help="Cellpose model name or local model path.",
    ),
    gpu: bool = typer.Option(
        True,
        "--gpu/--cpu",
        help="Request GPU/MPS/CUDA acceleration when available.",
    ),
    background_value: float = typer.Option(
        0.0,
        "--background-value",
        help="Scalar CH2 background value subtracted inside retained Cellpose regions.",
    ),
    flow_threshold: float = typer.Option(
        0.4,
        "--flow-threshold",
        help="Cellpose flow threshold for CH2+CH4 region segmentation.",
    ),
    cellprob_threshold: float = typer.Option(
        0.0,
        "--cellprob-threshold",
        help="Cellpose cell probability threshold for CH2+CH4 region segmentation.",
    ),
    diameter: float | None = typer.Option(
        None,
        "--diameter",
        help="Optional Cellpose region diameter. Leave unset for model/default behavior.",
    ),
    target_channel: str = typer.Option(
        "CH2",
        "--target-channel",
        help="Channel to measure as target/aSMA intensity, e.g. CH2.",
    ),
    dapi_channel: str = typer.Option(
        "CH4",
        "--dapi-channel",
        help="Channel to segment/count as DAPI-positive nuclei, e.g. CH4.",
    ),
    max_images_per_acquisition: int | None = typer.Option(
        None,
        "--max-images-per-acquisition",
        help="Limit images per acquisition for smoke tests.",
    ),
    render_figures: bool = typer.Option(
        True,
        "--render-figures/--skip-figures",
        help="Render final QC overlay pages and summary figures.",
    ),
    max_overlay_images: int | None = typer.Option(
        None,
        "--max-overlay-images",
        help="Optional cap for final overlay pages.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Only list discovered acquisition folders; do not run Cellpose.",
    ),
    overwrite: bool = typer.Option(
        False,
        "--overwrite",
        help="Allow writing into a non-empty output folder.",
    ),
) -> None:
    acquisitions = discover_acquisitions(
        input_root,
        target_channel_id=target_channel,
        dapi_channel_id=dapi_channel,
    )
    typer.echo(f"Discovered {len(acquisitions)} acquisition folder(s):")
    for acquisition in acquisitions:
        typer.echo(
            f"  - {acquisition.plate_name} / {acquisition.display_name}: "
            f"{acquisition.image_count} image pairs"
        )
    if dry_run:
        return

    result = run_user_cellpose_batch(
        input_root=input_root,
        output_root=output_root,
        model_name=model_name,
        gpu=gpu,
        background_value=background_value,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        diameter=diameter,
        target_channel_id=target_channel,
        dapi_channel_id=dapi_channel,
        max_images_per_acquisition=max_images_per_acquisition,
        render_figures=render_figures,
        max_overlay_images=max_overlay_images,
        overwrite=overwrite,
    )
    typer.echo(f"Processed image rows: {result.rows_processed}")
    typer.echo(f"User summary CSV: {result.outputs.user_summary_csv}")
    typer.echo(f"Workbook: {result.outputs.workbook_path}")
    typer.echo(f"Run summary: {result.outputs.run_summary_html}")


if __name__ == "__main__":
    app()
