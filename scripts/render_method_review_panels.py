#!/usr/bin/env python
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import tifffile
import typer

from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.method_review_panels import MethodReviewRecord, write_method_review_package
from dapi_norm.pi_simple_summary import ImagePair, find_image_pairs
from dapi_norm.seeded_regions import load_mask_lookup_from_counts_root

app = typer.Typer(
    help="Render side-by-side visual QC panels comparing candidate aSMA-region methods."
)


@app.command()
def main(
    input_root: Path = typer.Option(..., "--input", help="Image-pair input root."),
    nuclei_counts_root: Path = typer.Option(
        ...,
        "--nuclei-counts-root",
        help="Cellpose DAPI count output root containing nuclei label masks.",
    ),
    propagation_run: Path = typer.Option(
        ...,
        "--propagation-run",
        help="Seeded CellProfiler-style propagation run directory.",
    ),
    cellpose_runs: str = typer.Option(
        ...,
        "--cellpose-runs",
        help="Comma-separated Cellpose candidate-region run directories.",
    ),
    output_dir: Path = typer.Option(..., "--output", help="Review-package output directory."),
    positions: str = typer.Option(
        "XY22,XY23,XY24,XY40,XY41",
        "--positions",
        help="Comma-separated locations/source IDs to include.",
    ),
    crop_size: int = typer.Option(280, "--crop-size", help="Square review crop width/height."),
) -> None:
    requested = [token.strip().upper() for token in positions.split(",") if token.strip()]
    records = _load_review_records(
        input_root=input_root,
        nuclei_counts_root=nuclei_counts_root,
        propagation_run=propagation_run,
        cellpose_runs=[Path(part.strip()) for part in cellpose_runs.split(",") if part.strip()],
        requested_positions=requested,
    )
    outputs = write_method_review_package(records=records, output_dir=output_dir, crop_size=crop_size)
    for key, path in outputs.items():
        typer.echo(f"{key}={path}")


def _load_review_records(
    *,
    input_root: Path,
    nuclei_counts_root: Path,
    propagation_run: Path,
    cellpose_runs: list[Path],
    requested_positions: list[str],
) -> list[MethodReviewRecord]:
    pair_lookup = _pair_lookup(input_root)
    nucleus_lookup = load_mask_lookup_from_counts_root(nuclei_counts_root)
    propagation_metrics = _metrics_lookup(
        propagation_run / "summaries" / "seeded_region_image_metrics.csv"
    )
    cellpose_metrics = _merge_cellpose_metrics(cellpose_runs)

    records: list[MethodReviewRecord] = []
    for position in requested_positions:
        pair = _lookup_pair(pair_lookup, position)
        propagation_row = _lookup_metrics(propagation_metrics, position)
        cellpose_row = _lookup_metrics(cellpose_metrics, position)
        ch2_image, _ = read_primary_intensity_plane(pair.ch2_path)
        ch4_image, _ = read_primary_intensity_plane(pair.ch4_path)
        nuclei_mask = np.asarray(tifffile.imread(_lookup_nucleus_mask(nucleus_lookup, pair)))
        propagation_labels = np.asarray(
            tifffile.imread(_resolve_existing_path(propagation_row["mask_path"], base_dir=propagation_run)),
            dtype=np.uint32,
        )
        cellpose_labels = np.asarray(
            tifffile.imread(_resolve_existing_path(cellpose_row["mask_path"], base_dir=Path.cwd())),
            dtype=np.uint32,
        )
        records.append(
            MethodReviewRecord(
                image_id=position,
                ch2_image=ch2_image,
                ch4_image=ch4_image,
                nuclei_mask=nuclei_mask,
                propagation_labels=propagation_labels,
                cellpose_labels=cellpose_labels,
                propagation_metrics=propagation_row,
                cellpose_metrics=cellpose_row,
            )
        )
    return records


def _pair_lookup(input_root: Path) -> dict[str, ImagePair]:
    lookup: dict[str, ImagePair] = {}
    for pair in find_image_pairs(input_root):
        lookup[pair.location.upper()] = pair
        lookup[pair.source_id.upper()] = pair
    return lookup


def _lookup_pair(lookup: dict[str, ImagePair], position: str) -> ImagePair:
    key = position.upper()
    if key not in lookup:
        raise KeyError(f"No CH2/CH4 image pair found for {position}")
    return lookup[key]


def _lookup_nucleus_mask(mask_lookup: dict[str, Path], pair: ImagePair) -> Path:
    for key in [pair.source_id, pair.location]:
        normalized = key.strip().replace("\\", "/").upper()
        if normalized in mask_lookup:
            return mask_lookup[normalized]
    raise KeyError(f"No DAPI nuclei mask found for {pair.source_id} / {pair.location}")


def _metrics_lookup(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Metrics CSV missing: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        row_keys: set[str] = set()
        for key in [row.get("image_id", ""), row.get("source_id", "")]:
            normalized = key.strip().upper()
            if not normalized or normalized in row_keys:
                continue
            row_keys.add(normalized)
            if normalized in lookup:
                raise ValueError(f"Duplicate metrics key {normalized} in {path}")
            lookup[normalized] = row
    return lookup


def _merge_cellpose_metrics(run_dirs: list[Path]) -> dict[str, dict[str, str]]:
    merged: dict[str, dict[str, str]] = {}
    for run_dir in run_dirs:
        lookup = _metrics_lookup(run_dir / "summaries" / "cellpose_cell_region_image_metrics.csv")
        for key, row in lookup.items():
            if key in merged:
                raise ValueError(f"Duplicate Cellpose metrics for {key} across supplied runs")
            merged[key] = row
    return merged


def _lookup_metrics(lookup: dict[str, dict[str, str]], position: str) -> dict[str, str]:
    key = position.upper()
    if key not in lookup:
        raise KeyError(f"No metrics row found for {position}")
    return lookup[key]


def _resolve_existing_path(value: str, *, base_dir: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    for candidate in [path, base_dir / path, base_dir.parent / path, Path.cwd() / path]:
        if candidate.exists():
            return candidate
    return Path.cwd() / path


if __name__ == "__main__":
    app()
