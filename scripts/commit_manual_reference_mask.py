#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.manual_reference_commit import commit_manual_reference_mask

app = typer.Typer(
    help="Commit an edited manual/reference label mask into a manual-validation package."
)


@app.command()
def main(
    package_dir: Path = typer.Option(..., "--package", help="Manual validation package directory."),
    image_id: str = typer.Option(..., "--image-id", help="Image field ID, for example XY22."),
    labels_path: Path = typer.Option(
        ...,
        "--labels",
        help="Edited reference labels as TIFF, or NPZ containing manual_reference_labels.",
    ),
    status: str = typer.Option(
        "auto",
        "--status",
        help="auto, complete_non_empty, or confirmed_empty. Empty masks require confirmed_empty.",
    ),
    labeler: str = typer.Option("", "--labeler", help="Name or initials of the labeler."),
    completed_date: str | None = typer.Option(
        None,
        "--completed-date",
        help="Completion date as YYYY-MM-DD. Defaults to today's UTC date.",
    ),
    notes: str = typer.Option("", "--notes", help="Free-text annotation notes for the status CSV."),
    npz_key: str = typer.Option(
        "manual_reference_labels",
        "--npz-key",
        help="Layer key to read when --labels points to an NPZ bundle.",
    ),
) -> None:
    try:
        result = commit_manual_reference_mask(
            package_dir=package_dir,
            image_id=image_id,
            labels_path=labels_path,
            status=status,
            labeler=labeler,
            completed_date=completed_date,
            notes=notes,
            npz_key=npz_key,
        )
    except (FileNotFoundError, KeyError, ValueError) as exc:
        typer.echo(f"manual reference mask commit failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"image_id={result['image_id']}")
    typer.echo(f"status={result['status']}")
    typer.echo(f"positive_label_count={result['positive_label_count']}")
    typer.echo(f"foreground_area_px={result['foreground_area_px']}")
    typer.echo(f"manual_reference_mask_path={result['manual_reference_mask_path']}")
    typer.echo(f"status_path={result['status_path']}")


if __name__ == "__main__":
    app()
