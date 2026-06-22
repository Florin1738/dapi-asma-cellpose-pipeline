from __future__ import annotations

from pathlib import Path

import typer

from dapi_norm.data_inventory import inventory_dataset, write_inventory_reports
from dapi_norm.qc_preview import generate_position_preview

app = typer.Typer(help="Inspect local microscopy data organization and generate QC previews.")


def run_data_inventory(
    *,
    root: Path | str,
    output: Path | str,
    preview_positions: list[str] | None = None,
) -> None:
    root_path = Path(root)
    output_path = Path(output)
    inventory = inventory_dataset(root_path)
    write_inventory_reports(inventory, output_path)

    selected_positions = preview_positions
    if selected_positions is None:
        selected_positions = list(inventory.positions)[:4]

    preview_dir = output_path / "previews"
    missing_positions = [
        position_id.upper()
        for position_id in selected_positions
        if position_id.upper() not in inventory.positions
    ]
    if missing_positions:
        raise ValueError(
            "Requested preview positions not found: " + ", ".join(missing_positions)
        )

    for position_id in selected_positions:
        position = inventory.positions.get(position_id.upper())
        channel_paths = {channel: info.path for channel, info in position.channels.items()}
        if channel_paths:
            generate_position_preview(
                position_id=position.position_id,
                channel_paths=channel_paths,
                output_path=preview_dir / f"{position.position_id}_channels_preview.png",
            )


@app.command()
def main(
    root: Path = typer.Option(..., "--root", help="Dataset root to inspect."),
    output: Path = typer.Option(..., "--output", help="Directory for inventory outputs."),
    preview_position: list[str] | None = typer.Option(
        None,
        "--preview-position",
        help="Specific XY position to preview. Repeat for multiple positions.",
    ),
) -> None:
    run_data_inventory(root=root, output=output, preview_positions=preview_position)


if __name__ == "__main__":
    app()
