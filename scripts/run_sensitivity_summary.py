#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.sensitivity_summary import build_sensitivity_rows, write_sensitivity_outputs

app = typer.Typer(
    help=(
        "Build an exploratory robustness summary across region-restricted aSMA methods. "
        "This is not manual validation."
    )
)


@app.command()
def main(
    table: list[str] = typer.Option(
        ...,
        "--table",
        help="Run-labeled seeded metrics CSV as RUN_ID=PATH. Pass multiple times.",
    ),
    output_dir: Path = typer.Option(..., "--output", help="Output directory."),
    positions: str | None = typer.Option(
        None,
        "--positions",
        help="Optional comma-separated image positions to include, e.g. XY22,XY23,XY24.",
    ),
    ordered_positions: str = typer.Option(
        "XY22,XY23,XY24",
        "--ordered-positions",
        help="Comma-separated positions expected to be strictly descending if robust.",
    ),
    challenge_positions: str = typer.Option(
        "XY40,XY41",
        "--challenge-positions",
        help="Comma-separated positions expected to be QC-rejected.",
    ),
) -> None:
    rows = build_sensitivity_rows(
        table_specs=_parse_tables(table),
        positions=_parse_positions(positions),
    )
    outputs = write_sensitivity_outputs(
        rows,
        output_dir=output_dir,
        ordered_positions=_parse_positions(ordered_positions) or [],
        challenge_positions=_parse_positions(challenge_positions) or [],
    )
    typer.echo(f"wrote {outputs['long_csv']}")
    typer.echo(f"wrote {outputs['run_summary_csv']}")
    typer.echo(f"wrote {outputs['image_summary_csv']}")
    typer.echo(f"wrote {outputs['report']}")
    typer.echo(f"wrote {outputs['plot']}")


def _parse_tables(values: list[str]) -> list[tuple[str, Path]]:
    specs: list[tuple[str, Path]] = []
    for value in values:
        if "=" not in value:
            raise typer.BadParameter("--table values must be formatted as RUN_ID=PATH")
        run_id, path = value.split("=", maxsplit=1)
        run_id = run_id.strip()
        if not run_id:
            raise typer.BadParameter("--table RUN_ID cannot be empty")
        specs.append((run_id, Path(path).expanduser()))
    return specs


def _parse_positions(positions: str | None) -> list[str] | None:
    if positions is None:
        return None
    parsed = [position.strip() for position in positions.split(",") if position.strip()]
    return parsed or None


if __name__ == "__main__":
    app()
