#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.manual_annotation_gallery import render_manual_annotation_gallery

app = typer.Typer(help="Render a static HTML worklist for manual/reference annotation.")


@app.command()
def main(
    package_dir: Path = typer.Option(..., "--package", help="Manual validation package directory."),
    output_dir: Path = typer.Option(
        ...,
        "--output",
        help="Output directory for annotation_review_gallery/index.html.",
    ),
) -> None:
    try:
        outputs = render_manual_annotation_gallery(
            package_dir=package_dir,
            output_dir=output_dir,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"manual annotation gallery failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"wrote {outputs['index']}")


if __name__ == "__main__":
    app()
