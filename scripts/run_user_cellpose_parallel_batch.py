#!/usr/bin/env python
"""Run the existing Cellpose workflow once per acquisition with bounded concurrency."""

from __future__ import annotations

import sys
from pathlib import Path

import typer

from dapi_norm.parallel_user_cellpose_batch import (
    build_per_acquisition_runs,
    execute_parallel_runs,
    validate_new_output_root,
)


app = typer.Typer(
    help=(
        "Schedule isolated Cellpose runs per acquisition folder. Each successful "
        "folder receives its own standard workbook, QC figures, masks, and run log."
    )
)


@app.command()
def main(
    input_root: Path = typer.Option(
        ...,
        "--input",
        help="Parent folder containing one or more valid acquisition folders.",
    ),
    output_root: Path = typer.Option(
        ...,
        "--output",
        help="New or empty root for isolated per-acquisition output bundles.",
    ),
    workers: int = typer.Option(
        1,
        "--workers",
        min=1,
        help="Maximum number of simultaneous Cellpose processes. Start benchmark testing at 2.",
    ),
    model_name: str = typer.Option("cpsam_v2", "--model"),
    gpu: bool = typer.Option(True, "--gpu/--cpu"),
    background_value: float = typer.Option(0.0, "--background-value"),
    flow_threshold: float = typer.Option(0.4, "--flow-threshold"),
    cellprob_threshold: float = typer.Option(0.0, "--cellprob-threshold"),
    diameter: float | None = typer.Option(None, "--diameter"),
    target_channel: str = typer.Option("CH2", "--target-channel"),
    dapi_channel: str = typer.Option("CH4", "--dapi-channel"),
    max_images_per_acquisition: int | None = typer.Option(
        None,
        "--max-images-per-acquisition",
        help="Optional per-folder cap for a controlled benchmark/smoke test.",
    ),
    render_figures: bool = typer.Option(True, "--render-figures/--skip-figures"),
    max_overlay_images: int | None = typer.Option(None, "--max-overlay-images"),
    continue_on_error: bool = typer.Option(
        False,
        "--continue-on-error",
        help="Keep scheduling later acquisitions after a failed child run.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        help="Validate discovery and print per-acquisition destinations without writing files.",
    ),
) -> None:
    project_root = Path(__file__).resolve().parents[1]
    runs = build_per_acquisition_runs(
        input_root=input_root,
        output_root=output_root,
        python_executable=Path(sys.executable),
        runner_script=project_root / "scripts" / "run_user_cellpose_batch.py",
        target_channel_id=target_channel,
        dapi_channel_id=dapi_channel,
        model_name=model_name,
        gpu=gpu,
        background_value=background_value,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        diameter=diameter,
        max_images_per_acquisition=max_images_per_acquisition,
        render_figures=render_figures,
        max_overlay_images=max_overlay_images,
    )
    typer.echo(f"Planned {len(runs)} isolated acquisition run(s); workers={workers}.")
    for run in runs:
        typer.echo(f"  - {run.display_name}: {run.image_count} pairs -> {run.output_root}")
    if dry_run:
        return

    validate_new_output_root(output_root, runs)
    result = execute_parallel_runs(
        runs=runs,
        output_root=output_root,
        project_root=project_root,
        workers=workers,
        continue_on_error=continue_on_error,
    )
    typer.echo(f"Parallel run summary: {result.summary_path}")
    typer.echo(f"Orchestration manifest: {result.manifest_path}")
    if not result.succeeded:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
