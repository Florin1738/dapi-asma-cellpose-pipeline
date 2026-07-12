#!/usr/bin/env python
from __future__ import annotations

from collections import Counter
from pathlib import Path

import typer

from dapi_norm.cellpose_cell_regions import run_cellpose_cell_region_batch
from dapi_norm.pi_simple_summary import find_image_pairs
from dapi_norm.seeded_regions import load_mask_lookup_from_counts_root

app = typer.Typer(
    help=(
        "Run exploratory Cellpose target+DAPI candidate aSMA-associated object segmentation. "
        "This does not produce validated whole-cell masks."
    )
)


@app.command()
def main(
    input_root: Path = typer.Option(
        ...,
        "--input",
        help="Dataset root, preferably a single microscope output/sample folder with XY directories.",
    ),
    counts_root: Path = typer.Option(
        ...,
        "--counts-root",
        help="Existing Cellpose nuclei-count output root with masks and summaries.",
    ),
    output_dir: Path = typer.Option(..., "--output", help="Output directory for masks/QC/metrics."),
    positions: str = typer.Option(
        "",
        "--positions",
        help="Optional comma-separated XY positions to run, e.g. XY22,XY23,XY24.",
    ),
    model_name: str = typer.Option(
        "cpsam_v2",
        "--model-name",
        help="Cellpose pretrained model name or local model path.",
    ),
    gpu: bool = typer.Option(
        True,
        "--gpu/--no-gpu",
        help="Request Cellpose GPU acceleration when available.",
    ),
    background_value: float = typer.Option(
        0.0,
        "--background-value",
        help="Scalar target-channel background value for sensitivity metrics inside Cellpose candidate regions.",
    ),
    target_channel: str = typer.Option(
        "CH2",
        "--target-channel",
        help="Channel to measure as target/aSMA intensity, e.g. CH2.",
    ),
    dapi_channel: str = typer.Option(
        "CH4",
        "--dapi-channel",
        help="Channel used for DAPI nuclei context/mask lookup, e.g. CH4.",
    ),
    flow_threshold: float = typer.Option(
        0.4,
        "--flow-threshold",
        help="Cellpose flow threshold.",
    ),
    cellprob_threshold: float = typer.Option(
        0.0,
        "--cellprob-threshold",
        help="Cellpose cell probability threshold.",
    ),
    diameter: float | None = typer.Option(
        None,
        "--diameter",
        help="Optional Cellpose object diameter. Leave unset for model/default estimation.",
    ),
    max_images: int | None = typer.Option(
        None,
        "--max-images",
        help="Optional cap after position filtering.",
    ),
    internal_qc: bool = typer.Option(
        True,
        "--internal-qc/--skip-internal-qc",
        help=(
            "Render the older per-image Cellpose QC PNGs during segmentation. "
            "Disable when a separate overlay-panel rendering step will be run afterward."
        ),
    ),
) -> None:
    pairs = find_image_pairs(input_root, target_channel=target_channel, dapi_channel=dapi_channel)
    if positions:
        requested = {token.strip().upper() for token in positions.split(",") if token.strip()}
        pairs = [
            pair
            for pair in pairs
            if pair.location.upper() in requested or pair.source_id.upper() in requested
        ]
    if max_images is not None:
        pairs = pairs[:max_images]
    if not pairs:
        raise typer.BadParameter(f"No {target_channel}/{dapi_channel} image pairs matched under {input_root}")
    _reject_ambiguous_location_masks(pairs)

    rows = run_cellpose_cell_region_batch(
        image_pairs=pairs,
        mask_lookup=load_mask_lookup_from_counts_root(counts_root),
        output_dir=output_dir,
        model_name=model_name,
        gpu=gpu,
        background_value=background_value,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        diameter=diameter,
        target_channel_id=target_channel,
        dapi_channel_id=dapi_channel,
        write_internal_qc=internal_qc,
    )
    for row in rows:
        typer.echo(
            f"{row['source_id']}: "
            f"objects={row['cellpose_object_count']} "
            f"area={row['candidate_region_fraction']:.1%} "
            "per_nucleus="
            f"{row['target_integrated_intensity_per_DAPI_positive_nucleus']:.3e} "
            f"qc={row['qc_status']}"
        )


def _reject_ambiguous_location_masks(pairs) -> None:
    counts = Counter(pair.location.upper() for pair in pairs)
    duplicates = sorted(location for location, count in counts.items() if count > 1)
    if duplicates:
        joined = ", ".join(duplicates[:10])
        raise typer.BadParameter(
            "Selected pairs contain duplicate XY locations across microscope output folders "
            f"({joined}). Run this exploratory segmentation per sample folder so XY mask lookup "
            "cannot cross-wire nuclei masks."
        )


if __name__ == "__main__":
    app()
