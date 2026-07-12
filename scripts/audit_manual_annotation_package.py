#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.manual_annotation_audit import run_manual_annotation_audit

app = typer.Typer(
    help=(
        "Audit a manual-validation package for annotation completeness before "
        "running quantitative mask validation."
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
        help="Audit output directory.",
    ),
) -> None:
    try:
        outputs = run_manual_annotation_audit(
            package_dir=package_dir,
            output_dir=output_dir,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"manual annotation audit failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"wrote {outputs['csv']}")
    typer.echo(f"wrote {outputs['report']}")
    typer.echo(f"wrote {outputs['contact_sheet']}")


if __name__ == "__main__":
    app()
