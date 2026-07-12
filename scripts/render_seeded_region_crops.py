#!/usr/bin/env python
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import tifffile
import typer
import yaml

from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.pi_simple_summary import find_image_pairs
from dapi_norm.seeded_regions import write_seeded_region_crop_panel

app = typer.Typer(help="Render close-up crop panels from seeded aSMA-region outputs.")


@app.command()
def main(
    input_root: Path = typer.Option(..., "--input", help="Image-pair input root."),
    seeded_run_dir: Path = typer.Option(..., "--seeded-run", help="Seeded-region output dir."),
    output_path: Path = typer.Option(..., "--output", help="PNG crop panel path."),
    positions: str = typer.Option(
        "XY22,XY23,XY24,XY40,XY41",
        "--positions",
        help="Comma-separated XY locations/source IDs to include.",
    ),
    crop_size: int = typer.Option(280, "--crop-size", help="Square crop width/height in pixels."),
) -> None:
    requested = [token.strip().upper() for token in positions.split(",") if token.strip()]
    pair_lookup = _pair_lookup(input_root)
    metric_lookup = _metric_lookup(seeded_run_dir / "summaries" / "seeded_region_image_metrics.csv")
    config = _load_seeded_run_config(seeded_run_dir)
    config_path = seeded_run_dir / "logs" / "config_resolved.yaml"
    crops = []
    for position in requested:
        pair = pair_lookup[position]
        metrics = metric_lookup[position]
        ch2_image, _ = read_primary_intensity_plane(pair.ch2_path)
        ch4_image, _ = read_primary_intensity_plane(pair.ch4_path)
        labels_path = _resolve_recorded_path(
            Path(metrics["mask_path"]),
            seeded_run_dir=seeded_run_dir,
            config=config,
            config_path=config_path,
        )
        labels = np.asarray(tifffile.imread(labels_path))
        nuclei_mask = _load_nuclei_mask_for_pair(
            seeded_run_dir,
            position,
            config=config,
            config_path=config_path,
        )
        box = _largest_label_crop_box(labels, nuclei_mask=nuclei_mask, crop_size=crop_size)
        crops.append(
            {
                "image_id": position,
                "ch2_image": ch2_image,
                "ch4_image": ch4_image,
                "nuclei_mask": nuclei_mask,
                "seeded_labels": labels,
                "box": box,
                "caption": f"{metrics['qc_status']} | {metrics['qc_flags']}",
            }
        )
    write_seeded_region_crop_panel(crops=crops, output_path=output_path)
    typer.echo(f"wrote={output_path}")


def _pair_lookup(input_root: Path) -> dict[str, object]:
    pairs = find_image_pairs(input_root)
    lookup = {}
    for pair in pairs:
        lookup[pair.location.upper()] = pair
        lookup[pair.source_id.upper()] = pair
    return lookup


def _metric_lookup(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    lookup = {}
    for row in rows:
        lookup[row["image_id"].upper()] = row
        lookup[row["source_id"].upper()] = row
    return lookup


def _load_seeded_run_config(seeded_run_dir: Path) -> dict:
    config_path = seeded_run_dir / "logs" / "config_resolved.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Seeded run config does not exist: {config_path}")
    return yaml.safe_load(config_path.read_text(encoding="utf-8"))


def _load_nuclei_mask_for_pair(
    seeded_run_dir: Path,
    position: str,
    *,
    config: dict | None = None,
    config_path: Path | None = None,
) -> np.ndarray:
    if config_path is None:
        config_path = seeded_run_dir / "logs" / "config_resolved.yaml"
    if config is None:
        config = _load_seeded_run_config(seeded_run_dir)
    for record in config.get("image_records", []):
        keys = {str(record.get("source_id", "")).upper(), str(record.get("location", "")).upper()}
        if position.upper() not in keys:
            continue
        mask_path = Path(record["nuclei_mask_path"])
        mask_path = _resolve_recorded_path(
            mask_path,
            seeded_run_dir=seeded_run_dir,
            config=config,
            config_path=config_path,
        )
        if not mask_path.exists():
            raise FileNotFoundError(f"Recorded nuclei mask path does not exist: {mask_path}")
        return np.asarray(tifffile.imread(mask_path))
    raise KeyError(f"No image record for {position} in {config_path}")


def _resolve_recorded_path(
    path: Path,
    *,
    seeded_run_dir: Path,
    config: dict,
    config_path: Path,
) -> Path:
    if path.is_absolute():
        return path
    for base in _relative_path_bases(
        seeded_run_dir=seeded_run_dir,
        config=config,
        config_path=config_path,
    ):
        candidate = base / path
        if candidate.exists():
            return candidate
    return Path.cwd() / path


def _relative_path_bases(*, seeded_run_dir: Path, config: dict, config_path: Path) -> list[Path]:
    bases: list[Path] = []
    project_root = _infer_project_root_from_output_dir(seeded_run_dir, config)
    if project_root is not None:
        bases.append(project_root)
    bases.extend([seeded_run_dir, config_path.parent, config_path.parent.parent, Path.cwd()])
    deduped: list[Path] = []
    seen: set[Path] = set()
    for base in bases:
        resolved = base.resolve()
        if resolved not in seen:
            deduped.append(resolved)
            seen.add(resolved)
    return deduped


def _infer_project_root_from_output_dir(seeded_run_dir: Path, config: dict) -> Path | None:
    recorded_output_dir = Path(str(config.get("output_dir", "")))
    if recorded_output_dir.is_absolute() or not recorded_output_dir.parts:
        return None
    run_dir = seeded_run_dir.resolve()
    recorded_parts = recorded_output_dir.parts
    if len(recorded_parts) > len(run_dir.parts):
        return None
    if run_dir.parts[-len(recorded_parts) :] != recorded_parts:
        return None
    return run_dir.parents[len(recorded_parts) - 1]


def _largest_label_crop_box(
    labels: np.ndarray,
    *,
    crop_size: int,
    nuclei_mask: np.ndarray | None = None,
) -> tuple[int, int, int, int]:
    positive = labels[labels > 0]
    if positive.size == 0:
        fallback = None if nuclei_mask is None else np.asarray(nuclei_mask) > 0
        if fallback is not None and np.any(fallback):
            ys, xs = np.nonzero(fallback)
            center_y = int(round(float(np.mean(ys))))
            center_x = int(round(float(np.mean(xs))))
        else:
            center_y = labels.shape[0] // 2
            center_x = labels.shape[1] // 2
    else:
        label_ids, counts = np.unique(positive, return_counts=True)
        label_id = int(label_ids[int(np.argmax(counts))])
        ys, xs = np.nonzero(labels == label_id)
        center_y = int(round(float(np.mean(ys))))
        center_x = int(round(float(np.mean(xs))))
    half = max(8, int(crop_size) // 2)
    y0 = max(0, center_y - half)
    x0 = max(0, center_x - half)
    y1 = min(labels.shape[0], y0 + int(crop_size))
    x1 = min(labels.shape[1], x0 + int(crop_size))
    y0 = max(0, y1 - int(crop_size))
    x0 = max(0, x1 - int(crop_size))
    return y0, x0, y1, x1


if __name__ == "__main__":
    app()
