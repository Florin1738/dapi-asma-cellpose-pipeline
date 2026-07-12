#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.manual_annotation_handoff import prepare_manual_annotation_handoff

app = typer.Typer(
    help=(
        "Prepare per-image NPZ layer bundles for manual/reference mask annotation."
    )
)


@app.command()
def main(
    package_dir: Path = typer.Option(
        ...,
        "--package",
        help="Manual validation package directory.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output",
        help="Annotation handoff output directory.",
    ),
) -> None:
    try:
        outputs = prepare_manual_annotation_handoff(
            package_dir=package_dir,
            output_dir=output_dir,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"manual annotation handoff failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"wrote {outputs['manifest']}")
    typer.echo(f"wrote {outputs['readme']}")
    typer.echo(f"wrote layer bundles under {outputs['layer_dir']}")


if __name__ == "__main__":
    app()
