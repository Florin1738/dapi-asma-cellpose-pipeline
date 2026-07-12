#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.manual_validation_package import prepare_manual_validation_package

app = typer.Typer(
    help=(
        "Prepare a manual/reference mask validation package with raw-only annotation "
        "panels, candidate guide panels, blank label masks, manifest, and validation "
        "instructions."
    )
)


@app.command()
def main(
    input_root: Path = typer.Option(..., "--input", help="Image-pair input root."),
    seeded_run_dir: Path = typer.Option(
        ...,
        "--seeded-run",
        help="Seeded-region candidate output directory.",
    ),
    output_dir: Path = typer.Option(..., "--output", help="Manual validation package output."),
    positions: str | None = typer.Option(
        None,
        "--positions",
        help="Optional comma-separated XY positions to include.",
    ),
    iou_threshold: float = typer.Option(
        0.5,
        "--iou-threshold",
        min=1e-12,
        max=1.0,
        help="IoU threshold to record in validation instructions.",
    ),
    task: str = typer.Option(
        "asma_associated_region",
        "--task",
        help="Manual validation task label recorded in manifest.",
    ),
    force_overwrite_reference_masks: bool = typer.Option(
        False,
        "--force-overwrite-reference-masks/--no-force-overwrite-reference-masks",
        help=(
            "Overwrite existing manual/reference TIFFs with blank masks. "
            "Use only before manual labels have been created."
        ),
    ),
) -> None:
    outputs = prepare_manual_validation_package(
        input_root=input_root,
        seeded_run_dir=seeded_run_dir,
        output_dir=output_dir,
        positions=_parse_positions(positions),
        iou_threshold=iou_threshold,
        task=task,
        force_overwrite_reference_masks=force_overwrite_reference_masks,
    )
    typer.echo(f"wrote {outputs['manifest']}")
    typer.echo(f"wrote manual labeling status under {outputs['status']}")
    typer.echo(f"wrote {outputs['readme']}")
    typer.echo(f"wrote raw-only annotation panels under {outputs['annotation_dir']}")
    typer.echo(f"wrote candidate guide panels under {outputs['guide_dir']}")
    typer.echo(f"wrote blank masks under {outputs['reference_mask_dir']}")


def _parse_positions(positions: str | None) -> list[str] | None:
    if positions is None:
        return None
    parsed = [position.strip() for position in positions.split(",") if position.strip()]
    return parsed or None


if __name__ == "__main__":
    app()
