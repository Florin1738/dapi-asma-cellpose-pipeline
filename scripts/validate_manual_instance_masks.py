#!/usr/bin/env python
from __future__ import annotations

import csv
from pathlib import Path
import re
import sys

import typer

from dapi_norm.segmentation_validation import run_manual_mask_validation

app = typer.Typer(help="Compare candidate instance masks against manual/reference masks by IoU.")


@app.command()
def main(
    candidate_dir: Path = typer.Option(..., "--candidate-dir", help="Directory of candidate label masks."),
    reference_dir: Path = typer.Option(..., "--reference-dir", help="Directory of manual/reference label masks."),
    completion_status: Path | None = typer.Option(
        None,
        "--completion-status",
        help=(
            "Optional manual_labeling_status.csv. If omitted, the validator "
            "uses reference_dir/../manual_labeling_status.csv when present."
        ),
    ),
    output_dir: Path = typer.Option(..., "--output", help="Output directory for validation CSVs."),
    iou_threshold: float = typer.Option(
        0.5,
        "--iou-threshold",
        min=1e-12,
        max=1.0,
        help="Object match IoU threshold; must be greater than 0.",
    ),
) -> None:
    try:
        candidate_lookup = _mask_lookup(candidate_dir)
        reference_lookup = _mask_lookup(reference_dir)
        resolved_completion_status = _resolve_completion_status_path(
            completion_status=completion_status,
            reference_dir=reference_dir,
        )
        completion_lookup, completion_mask_lookup = _completion_status_lookup(
            resolved_completion_status
        )
        summaries = run_manual_mask_validation(
            candidate_mask_paths=candidate_lookup,
            reference_mask_paths=reference_lookup,
            output_dir=output_dir,
            iou_threshold=iou_threshold,
            reference_completion_status=completion_lookup,
            reference_completion_mask_paths=completion_mask_lookup,
            run_metadata={
                "candidate_dir": str(candidate_dir),
                "reference_dir": str(reference_dir),
                "completion_status": str(resolved_completion_status)
                if resolved_completion_status
                else "",
                "argv": sys.argv,
            },
        )
    except (FileNotFoundError, ValueError) as exc:
        typer.echo(f"validation failed: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    for summary in summaries:
        typer.echo(
            f"{summary['image_id']}: "
            f"precision={summary['precision']:.3f} "
            f"recall={summary['recall']:.3f} "
            f"f1={summary['f1']:.3f} "
            f"mean_iou={summary['mean_iou_matched']:.3f}"
        )


def _mask_lookup(root: Path) -> dict[str, Path]:
    if not root.exists():
        raise FileNotFoundError(f"Mask directory does not exist: {root}")
    lookup: dict[str, Path] = {}
    for path in sorted(root.rglob("*.tif")) + sorted(root.rglob("*.tiff")):
        match = re.search(r"(XY\d+)", path.name, flags=re.IGNORECASE)
        if match is None:
            continue
        image_id = match.group(1).upper()
        if image_id in lookup:
            raise ValueError(f"Multiple masks found for {image_id} under {root}")
        lookup[image_id] = path
    if not lookup:
        raise ValueError(f"No XY## TIFF masks found under {root}")
    return lookup


def _resolve_completion_status_path(
    *,
    completion_status: Path | None,
    reference_dir: Path,
) -> Path | None:
    if completion_status is not None:
        return completion_status
    default_path = reference_dir.parent / "manual_labeling_status.csv"
    if default_path.exists():
        return default_path
    return None


def _completion_status_lookup(
    path: Path | None,
) -> tuple[dict[str, str] | None, dict[str, Path] | None]:
    if path is None:
        return None, None
    if not path.exists():
        raise FileNotFoundError(f"Completion status file does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"image_id", "manual_reference_mask_path", "status"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(
                f"{path} is missing required columns: " + ", ".join(missing)
            )
        lookup: dict[str, str] = {}
        mask_paths: dict[str, Path] = {}
        for row in reader:
            image_id = row["image_id"].strip().upper().replace(" ", "")
            if not image_id:
                continue
            if image_id in lookup:
                raise ValueError(f"Duplicate completion status for {image_id} in {path}")
            lookup[image_id] = row["status"].strip().lower()
            mask_paths[image_id] = _resolve_status_mask_path(
                row["manual_reference_mask_path"],
                status_file=path,
            )
    if not lookup:
        raise ValueError(f"No image completion statuses found in {path}")
    return lookup, mask_paths


def _resolve_status_mask_path(value: str, *, status_file: Path) -> Path:
    path = Path(value.strip()).expanduser()
    if path.is_absolute():
        return path
    if path.exists():
        return path
    return status_file.parent / path


if __name__ == "__main__":
    app()
