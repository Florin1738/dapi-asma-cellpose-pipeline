#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.raw_annotation_export import prepare_raw_annotation_export

app = typer.Typer(
    help=(
        "Prepare raw-only CH2/CH4 TIFFs and editable manual label scratch TIFFs "
        "for blinded manual/reference annotation."
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
        help="Raw-only annotation export output directory.",
    ),
) -> None:
    try:
        outputs = prepare_raw_annotation_export(
            package_dir=package_dir,
            output_dir=output_dir,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"raw annotation export failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"wrote {outputs['manifest']}")
    typer.echo(f"wrote {outputs['readme']}")


if __name__ == "__main__":
    app()
