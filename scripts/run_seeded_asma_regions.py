#!/usr/bin/env python
from __future__ import annotations

from collections import Counter
from pathlib import Path

import typer

from dapi_norm.pi_simple_summary import find_image_pairs
from dapi_norm.seeded_regions import (
    load_mask_lookup_from_counts_root,
    run_seeded_region_batch,
)

app = typer.Typer(
    help=(
        "Run exploratory DAPI-seeded CH2/aSMA-associated region segmentation. "
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
        help="Optional comma-separated XY positions to run, e.g. XY12,XY22,XY23.",
    ),
    foreground_method: str = typer.Option(
        "li",
        "--foreground-method",
        help="CH2 foreground threshold method: li, otsu, triangle, or value:<n>.",
    ),
    background_value: float = typer.Option(
        0.0,
        "--background-value",
        help="Scalar CH2 background value subtracted inside seeded regions for sensitivity metrics.",
    ),
    segmentation_method: str = typer.Option(
        "watershed",
        "--segmentation-method",
        help="Seeded segmentation method: watershed, random_walker, or propagation.",
    ),
    random_walker_beta: float = typer.Option(
        90.0,
        "--random-walker-beta",
        min=1e-12,
        help="Random-walker edge sensitivity parameter; used only with random_walker.",
    ),
    propagation_regularization: float = typer.Option(
        0.05,
        "--propagation-regularization",
        min=1e-12,
        help=(
            "CellProfiler-style propagation regularization/weight; used only with "
            "segmentation_method=propagation."
        ),
    ),
    min_size: int = typer.Option(
        128,
        "--min-size",
        help="Minimum connected CH2 foreground object size in pixels.",
    ),
    max_images: int | None = typer.Option(
        None,
        "--max-images",
        help="Optional cap after position filtering.",
    ),
) -> None:
    pairs = find_image_pairs(input_root)
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
        raise typer.BadParameter(f"No CH2/CH4 image pairs matched under {input_root}")
    _reject_ambiguous_location_masks(pairs)

    mask_lookup = load_mask_lookup_from_counts_root(counts_root)
    rows = run_seeded_region_batch(
        image_pairs=pairs,
        mask_lookup=mask_lookup,
        output_dir=output_dir,
        foreground_method=foreground_method,
        background_value=background_value,
        min_size=min_size,
        segmentation_method=segmentation_method,
        random_walker_beta=random_walker_beta,
        propagation_regularization=propagation_regularization,
    )
    for row in rows:
        typer.echo(
            f"{row['source_id']}: "
            f"area={row['seeded_region_fraction']:.1%} "
            "bgcorr_per_nucleus="
            f"{row['seeded_region_intensity_per_DAPI_positive_nucleus']:.3e} "
            f"warnings={row['warnings']}"
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
