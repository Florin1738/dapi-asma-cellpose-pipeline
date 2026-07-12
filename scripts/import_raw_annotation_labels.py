#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.manual_reference_bulk_import import import_raw_annotation_labels

app = typer.Typer(
    help=(
        "Bulk-import edited manual/reference label TIFFs from a raw annotation export. "
        "Non-empty labels are committed; empty fields require explicit confirmation."
    )
)


@app.command()
def main(
    package_dir: Path = typer.Option(
        ...,
        "--package",
        help="Manual validation package directory.",
    ),
    raw_export_manifest: Path = typer.Option(
        ...,
        "--raw-export-manifest",
        help="raw_annotation_export_manifest.csv produced by prepare_raw_annotation_export.py.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output",
        help="Bulk import summary output directory.",
    ),
    labeler: str = typer.Option(..., "--labeler", help="Name or initials of the labeler."),
    completed_date: str | None = typer.Option(
        None,
        "--completed-date",
        help="Completion date as YYYY-MM-DD. Defaults to commit helper's UTC date.",
    ),
    notes: str = typer.Option("", "--notes", help="Free-text note to write into status CSV."),
    confirm_empty: list[str] = typer.Option(
        [],
        "--confirm-empty",
        help="Image ID to explicitly import as confirmed_empty. Repeat for multiple fields.",
    ),
    require_all_decisions: bool = typer.Option(
        False,
        "--require-all-decisions",
        help="Fail before committing anything if any exported label TIFF is empty and not confirmed.",
    ),
) -> None:
    try:
        outputs = import_raw_annotation_labels(
            package_dir=package_dir,
            raw_export_manifest_path=raw_export_manifest,
            output_dir=output_dir,
            labeler=labeler,
            completed_date=completed_date,
            notes=notes,
            confirm_empty_ids=set(confirm_empty),
            require_all_decisions=require_all_decisions,
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"raw annotation label import failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"wrote {outputs['summary']}")


if __name__ == "__main__":
    app()
