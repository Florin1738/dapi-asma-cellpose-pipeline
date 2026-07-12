#!/usr/bin/env python
from __future__ import annotations

import csv
from pathlib import Path
import shlex

import numpy as np
import tifffile
import typer

from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.pi_simple_summary import ImagePair, find_image_pairs
from dapi_norm.seeded_regions import load_mask_lookup_from_counts_root
from dapi_norm.validation_set_selection import (
    ValidationFeatureRecord,
    render_validation_selection_panel,
    select_validation_set,
    write_selection_csvs,
)

app = typer.Typer(
    help="Select a stratified set of fields for manual/reference mask validation."
)


@app.command()
def main(
    input_root: Path = typer.Option(..., "--input", help="Single sample folder with XY image pairs."),
    nuclei_counts_root: Path = typer.Option(
        ...,
        "--nuclei-counts-root",
        help="Cellpose DAPI count output root containing nuclei masks.",
    ),
    output_dir: Path = typer.Option(..., "--output", help="Selection output directory."),
    max_images: int = typer.Option(16, "--max-images", min=1, help="Maximum fields to select."),
    per_bucket: int = typer.Option(
        2,
        "--per-bucket",
        min=1,
        help="Fields to draw from each feature bucket before deduplication.",
    ),
    must_include: str = typer.Option(
        "XY22,XY23,XY24,XY40,XY41",
        "--must-include",
        help="Comma-separated fields to include if present.",
    ),
    method_review_summary: Path | None = typer.Option(
        None,
        "--method-review-summary",
        help="Optional method_comparison_review_summary.csv with method Jaccard values.",
    ),
) -> None:
    must_include_ids = _parse_tokens(must_include)
    records = build_validation_feature_records(
        input_root=input_root,
        nuclei_counts_root=nuclei_counts_root,
        method_review_summary=method_review_summary,
    )
    selected = select_validation_set(
        records,
        max_images=max_images,
        per_bucket=per_bucket,
        must_include=must_include_ids,
    )
    outputs = write_selection_csvs(
        all_records=records,
        selected_records=selected,
        output_dir=output_dir,
    )
    panel_path = output_dir / "selected_manual_validation_fields_panel.png"
    render_validation_selection_panel(
        selected_records=selected,
        output_path=panel_path,
        image_loader=_load_primary_plane,
        nuclei_loader=_load_label_mask,
    )
    readme_path = output_dir / "README.md"
    _write_readme(
        readme_path,
        input_root=input_root,
        nuclei_counts_root=nuclei_counts_root,
        output_dir=output_dir,
        method_review_summary=method_review_summary,
        max_images=max_images,
        per_bucket=per_bucket,
        must_include=must_include_ids,
        all_records=records,
        selected=selected,
    )
    typer.echo(f"all_features={outputs['all_features']}")
    typer.echo(f"selected={outputs['selected']}")
    typer.echo(f"panel={panel_path}")
    typer.echo(f"readme={readme_path}")
    typer.echo("selected_positions=" + ",".join(row.image_id for row in selected))


def build_validation_feature_records(
    *,
    input_root: Path,
    nuclei_counts_root: Path,
    method_review_summary: Path | None = None,
) -> list[ValidationFeatureRecord]:
    pairs = find_image_pairs(input_root)
    if not pairs:
        raise ValueError(f"No CH2/CH4 pairs found under {input_root}")
    mask_lookup = load_mask_lookup_from_counts_root(nuclei_counts_root)
    jaccard_lookup = _method_jaccard_lookup(method_review_summary)
    records: list[ValidationFeatureRecord] = []
    for pair in pairs:
        ch2, _ = read_primary_intensity_plane(pair.ch2_path)
        ch4, _ = read_primary_intensity_plane(pair.ch4_path)
        nuclei_path = _lookup_nucleus_mask(mask_lookup, pair)
        nuclei = np.asarray(tifffile.imread(nuclei_path))
        if ch2.shape != ch4.shape or ch2.shape != nuclei.shape:
            raise ValueError(
                f"{pair.source_id} shape mismatch: "
                f"ch2={ch2.shape}, ch4={ch4.shape}, nuclei={nuclei.shape}"
            )
        nucleus_count = _count_nonzero_labels(nuclei)
        raw_integrated = float(np.sum(np.asarray(ch2, dtype=np.float64)))
        records.append(
            ValidationFeatureRecord(
                image_id=pair.location,
                source_id=pair.source_id,
                ch2_path=pair.ch2_path,
                ch4_path=pair.ch4_path,
                nuclei_mask_path=nuclei_path,
                dapi_positive_nucleus_count=nucleus_count,
                target_integrated_raw=raw_integrated,
                target_integrated_raw_per_DAPI_positive_nucleus=(
                    raw_integrated / nucleus_count if nucleus_count else float("nan")
                ),
                target_saturation_fraction=_saturation_fraction(ch2),
                dapi_saturation_fraction=_saturation_fraction(ch4),
                method_region_jaccard=jaccard_lookup.get(pair.location.upper()),
            )
        )
    return records


def _lookup_nucleus_mask(mask_lookup: dict[str, Path], pair: ImagePair) -> Path:
    for key in [pair.source_id, pair.location]:
        normalized = key.strip().replace("\\", "/").upper()
        if normalized in mask_lookup:
            return mask_lookup[normalized]
    raise KeyError(f"No DAPI nuclei mask found for {pair.source_id} / {pair.location}")


def _method_jaccard_lookup(path: Path | None) -> dict[str, float]:
    if path is None:
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    lookup: dict[str, float] = {}
    for row in rows:
        image_id = row["image_id"].strip().upper()
        if image_id in lookup:
            raise ValueError(f"Duplicate image_id in method review summary: {image_id}")
        lookup[image_id] = float(row["method_region_jaccard"])
    return lookup


def _saturation_fraction(image: np.ndarray) -> float:
    arr = np.asarray(image)
    if arr.size == 0:
        return 0.0
    if np.issubdtype(arr.dtype, np.integer):
        saturation_value = np.iinfo(arr.dtype).max
    else:
        saturation_value = float(np.max(arr))
    return float(np.count_nonzero(arr >= saturation_value) / arr.size)


def _count_nonzero_labels(mask: np.ndarray) -> int:
    labels = np.unique(mask)
    return int(np.count_nonzero(labels))


def _load_primary_plane(path: Path) -> np.ndarray:
    image, _ = read_primary_intensity_plane(path)
    return image


def _load_label_mask(path: Path) -> np.ndarray:
    return np.asarray(tifffile.imread(path))


def _parse_tokens(value: str) -> list[str]:
    return [token.strip().upper() for token in value.split(",") if token.strip()]


def _write_readme(
    path: Path,
    *,
    input_root: Path,
    nuclei_counts_root: Path,
    output_dir: Path,
    method_review_summary: Path | None,
    max_images: int,
    per_bucket: int,
    must_include: list[str],
    all_records: list[ValidationFeatureRecord],
    selected: list[ValidationFeatureRecord],
) -> None:
    command_parts = [
        ".venv/bin/python",
        "scripts/select_manual_validation_fields.py",
        "--input",
        str(input_root),
        "--nuclei-counts-root",
        str(nuclei_counts_root),
        "--output",
        str(output_dir),
        "--max-images",
        str(max_images),
        "--per-bucket",
        str(per_bucket),
        "--must-include",
        ",".join(must_include),
    ]
    if method_review_summary is not None:
        command_parts.extend(["--method-review-summary", str(method_review_summary)])
    command = " \\\n  ".join(shlex.quote(part) for part in command_parts)
    method_review_text = str(method_review_summary) if method_review_summary is not None else "not used"
    selected_table = _selected_reason_table(selected)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "# Manual Validation Field Selection",
                "",
                f"Input root: `{input_root}`",
                f"Nuclei counts/masks root: `{nuclei_counts_root}`",
                f"Method-review summary: `{method_review_text}`",
                f"Candidate fields evaluated: `{len(all_records)}`",
                f"Maximum selected fields requested: `{max_images}`",
                f"Fields per feature bucket requested: `{per_bucket}`",
                f"Must-include fields requested: `{','.join(must_include)}`",
                f"Selected fields: `{','.join(row.image_id for row in selected)}`",
                "",
                "Purpose: choose a representative manual/reference labeling set before claiming any segmentation method is validated.",
                "",
                "Selection covers raw CH2/aSMA intensity, CH2 intensity per DAPI-positive nucleus, DAPI-positive nucleus count, CH2 saturation, user-required challenge fields, and optional method-disagreement fields.",
                "",
                "Exact command:",
                "",
                "```bash",
                command,
                "```",
                "",
                "Per-field selection reasons:",
                "",
                selected_table,
                "",
                "This selection is not a validation result. It defines which fields should be manually annotated so candidate masks can later be scored quantitatively.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _selected_reason_table(selected: list[ValidationFeatureRecord]) -> str:
    rows = [
        "| Field | Reasons | Raw CH2 integrated intensity | Raw CH2 per DAPI-positive nucleus | DAPI-positive nuclei | CH2 saturation | Method Jaccard |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for record in selected:
        rows.append(
            "| "
            + " | ".join(
                [
                    record.image_id,
                    f"`{record.selection_reasons}`",
                    _format_scientific(record.target_integrated_raw),
                    _format_scientific(
                        record.target_integrated_raw_per_DAPI_positive_nucleus
                    ),
                    str(record.dapi_positive_nucleus_count),
                    _format_percent(record.target_saturation_fraction),
                    _format_float(record.method_region_jaccard),
                ]
            )
            + " |"
        )
    return "\n".join(rows)


def _format_scientific(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.3e}"


def _format_percent(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{100.0 * value:.3f}%"


def _format_float(value: float | None) -> str:
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.3f}"


if __name__ == "__main__":
    app()
