from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.pi_simple_summary import (
    build_pi_summary,
    load_counts_by_plate,
    write_pi_plots,
    write_pi_workbook,
)

app = typer.Typer(
    help=(
        "Create the PI-requested simple aSMA/DAPI Excel workbook: LOCATION, "
        "aSMA intensity, Nuclei Count, Ratio."
    )
)


@app.command()
def main(
    input_root: Path = typer.Option(
        Path("."),
        "--input",
        help="Folder containing Plate 1/Plate 2 folders, or the current single-plate image folder.",
    ),
    counts_root: Path | None = typer.Option(
        Path("output/cellpose_counts"),
        "--counts",
        help=(
            "Cellpose counts folder or nucleus_counts.csv. If omitted or missing, nuclei counts and "
            "ratios are left blank."
        ),
    ),
    output_path: Path = typer.Option(
        Path("output/pi_simple_summary/aSMA_DAPI_nuclei_ratio_summary.xlsx"),
        "--output",
        help="Excel workbook path to write.",
    ),
    plots_dir: Path = typer.Option(
        Path("output/pi_simple_summary/plots"),
        "--plots",
        help="Directory for side plots summarizing aSMA intensity and ratio.",
    ),
) -> None:
    counts = load_counts_by_plate(counts_root if counts_root and counts_root.exists() else None)
    summary = build_pi_summary(input_root=input_root, counts_by_plate=counts)
    write_pi_workbook(output_path, summary)
    plot_paths = write_pi_plots(plots_dir, summary)

    for plate_name, rows in summary.items():
        typer.echo(f"{plate_name}: rows={len(rows)}")
    typer.echo(f"wrote={output_path}")
    typer.echo(f"plots={len(plot_paths)}")


if __name__ == "__main__":
    app()
