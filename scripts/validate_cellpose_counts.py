#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.count_validation import validate_count_outputs

app = typer.Typer(help="Validate Cellpose count-output artifacts.")


@app.command()
def main(
    output_dir: Path = typer.Option(..., "--output", help="Cellpose count output directory."),
) -> None:
    result = validate_count_outputs(output_dir)
    typer.echo(f"summary_rows={result['summary_rows']}")
    typer.echo(f"total_nucleus_count={result['total_nucleus_count']}")
    typer.echo(f"per_nucleus_rows={result['per_nucleus_rows']}")
    typer.echo(f"mask_counts_match_csv={result['mask_counts_match_csv']}")


if __name__ == "__main__":
    app()
