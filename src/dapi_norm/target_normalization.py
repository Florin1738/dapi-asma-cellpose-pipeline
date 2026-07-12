from __future__ import annotations

import csv
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from skimage.segmentation import find_boundaries
import tifffile
import yaml

from dapi_norm.cellpose_runner import write_qc_contact_sheet
from dapi_norm.data_inventory import inventory_dataset
from dapi_norm.image_arrays import read_primary_intensity_plane

IMAGE_LEVEL_COLUMNS = [
    "image_id",
    "well_id",
    "input_path",
    "backend",
    "dapi_channel",
    "target_channel",
    "raw_nucleus_count",
    "filtered_nucleus_count",
    "target_area_px",
    "target_integrated_raw",
    "target_mean_raw",
    "background_method",
    "background_value_per_px",
    "target_integrated_background_corrected",
    "target_integrated_intensity_per_DAPI_positive_nucleus",
    "dapi_saturation_fraction",
    "target_saturation_fraction",
    "mask_path",
    "qc_overlay_path",
    "warnings",
]


def run_target_normalization(
    *,
    input_root: Path | str,
    counts_dir: Path | str,
    output_dir: Path | str,
    target_channel_id: str = "CH2",
    dapi_channel_id: str = "CH4",
    background_percentile: float = 10,
) -> list[dict[str, Any]]:
    input_path = Path(input_root)
    counts_path = Path(counts_dir)
    output_path = Path(output_dir)
    target_channel = target_channel_id.upper()
    dapi_channel = dapi_channel_id.upper()

    inventory = inventory_dataset(input_path)
    count_rows = _read_count_rows(counts_path / "summaries" / "nucleus_counts.csv")
    rows: list[dict[str, Any]] = []
    qc_paths: list[Path] = []

    qc_dir = output_path / "qc"
    summaries_dir = output_path / "summaries"
    plots_dir = output_path / "plots"
    logs_dir = output_path / "logs"
    qc_dir.mkdir(parents=True, exist_ok=True)
    summaries_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    for position_id, count_row in count_rows.items():
        position = inventory.positions.get(position_id)
        if position is None:
            raise ValueError(f"{position_id} is present in counts but not input inventory")
        if target_channel not in position.channels:
            raise ValueError(f"{position_id} has no {target_channel} target image")
        if dapi_channel not in position.channels:
            raise ValueError(f"{position_id} has no {dapi_channel} DAPI image")

        _validate_count_row(
            count_row=count_row,
            expected_channel=dapi_channel,
            expected_dapi_path=position.channels[dapi_channel].path,
            counts_dir=counts_path,
        )
        target_image, target_pages = read_primary_intensity_plane(
            position.channels[target_channel].path
        )
        dapi_image, dapi_pages = read_primary_intensity_plane(position.channels[dapi_channel].path)
        mask_path = _resolve_reference(count_row["mask_path"], counts_path)
        mask = np.asarray(tifffile.imread(mask_path))
        if mask.shape != target_image.shape:
            raise ValueError(
                f"{position_id} mask shape {mask.shape} does not match target image shape "
                f"{target_image.shape}"
            )

        row = _measure_one_image(
            image_id=position_id,
            target_path=position.channels[target_channel].path,
            target_image=target_image,
            dapi_image=dapi_image,
            count_row=count_row,
            mask_path=mask_path,
            target_channel=target_channel,
            dapi_channel=dapi_channel,
            background_percentile=background_percentile,
        )
        qc_path = qc_dir / f"{position_id}_{target_channel}_target_with_{dapi_channel}_nucleus_outlines.png"
        write_target_overlay(
            target_image=target_image,
            mask=mask,
            output_path=qc_path,
            title=f"{position_id} {target_channel} target with {dapi_channel} nuclei",
        )
        row["qc_overlay_path"] = str(qc_path)
        row["mask_label_count"] = _mask_label_count(mask)
        row["target_page_count"] = int(target_pages)
        row["dapi_page_count"] = int(dapi_pages)
        row["target_extracted_plane_shape"] = [int(target_image.shape[0]), int(target_image.shape[1])]
        row["dapi_extracted_plane_shape"] = [int(dapi_image.shape[0]), int(dapi_image.shape[1])]
        row["count_backend"] = count_row.get("backend", "")
        row["count_model_name"] = count_row.get("model_name", "")
        rows.append(row)
        qc_paths.append(qc_path)

    _write_summary_csv(summaries_dir / "image_level_summary.csv", rows)
    _write_summary_csv(summaries_dir / "well_level_summary.csv", rows)
    write_normalization_plots(rows, plots_dir)
    write_qc_contact_sheet(qc_paths, output_path / "qc_contact_sheet.png")
    _write_run_metadata(
        logs_dir=logs_dir,
        input_root=input_path,
        counts_dir=counts_path,
        output_dir=output_path,
        target_channel=target_channel,
        dapi_channel=dapi_channel,
        background_percentile=background_percentile,
        rows=rows,
    )
    return rows


def _measure_one_image(
    *,
    image_id: str,
    target_path: Path,
    target_image: np.ndarray,
    dapi_image: np.ndarray,
    count_row: dict[str, str],
    mask_path: Path,
    target_channel: str,
    dapi_channel: str,
    background_percentile: float,
) -> dict[str, Any]:
    target_float = target_image.astype(np.float64)
    background_value = float(np.percentile(target_float, background_percentile))
    corrected = np.clip(target_float - background_value, 0, None)
    nucleus_count = int(count_row["nucleus_count"])
    normalized = (
        float(np.sum(corrected)) / nucleus_count if nucleus_count > 0 else float("nan")
    )
    warnings = [
        "channel_identity_unconfirmed",
        "rendered_rgb_source",
        "full_field_target_area",
        "no_nucleus_filtering_applied",
    ]
    if nucleus_count == 0:
        warnings.append("zero_nucleus_count")

    return {
        "image_id": image_id,
        "well_id": image_id,
        "input_path": str(target_path),
        "backend": "cellpose",
        "dapi_channel": dapi_channel,
        "target_channel": target_channel,
        "raw_nucleus_count": nucleus_count,
        "filtered_nucleus_count": nucleus_count,
        "target_area_px": int(target_image.size),
        "target_integrated_raw": float(np.sum(target_float)),
        "target_mean_raw": float(np.mean(target_float)),
        "background_method": f"percentile_{_format_percentile(background_percentile)}",
        "background_value_per_px": background_value,
        "target_integrated_background_corrected": float(np.sum(corrected)),
        "target_integrated_intensity_per_DAPI_positive_nucleus": normalized,
        "dapi_saturation_fraction": _saturation_fraction(dapi_image),
        "target_saturation_fraction": _saturation_fraction(target_image),
        "mask_path": str(mask_path),
        "qc_overlay_path": "",
        "warnings": ";".join(warnings),
    }


def write_target_overlay(
    *,
    target_image: np.ndarray,
    mask: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    scaled = _scale_for_display(target_image)
    overlay = np.dstack([scaled, scaled, scaled])
    boundaries = find_boundaries(mask, mode="outer")
    overlay[boundaries] = [0.0, 1.0, 1.0]

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    axes[0].imshow(scaled, cmap="gray")
    axes[0].set_title(f"{title}\ntarget intensity")
    axes[1].imshow(overlay)
    axes[1].set_title("nucleus outlines")
    for ax in axes:
        ax.axis("off")
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def write_normalization_plots(rows: list[dict[str, Any]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    sorted_rows = sorted(
        rows,
        key=lambda row: float(row["target_integrated_intensity_per_DAPI_positive_nucleus"]),
        reverse=True,
    )
    labels = [row["well_id"] for row in sorted_rows]
    normalized = [
        float(row["target_integrated_intensity_per_DAPI_positive_nucleus"]) / 1_000_000
        for row in sorted_rows
    ]
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    bars = ax.bar(labels, normalized, color="#3b82f6")
    ax.set_title("Candidate target intensity per DAPI-positive nucleus")
    ax.set_ylabel("background-corrected CH2 intensity per CH4 nucleus (millions)")
    ax.set_xlabel("XY position / candidate well")
    ax.tick_params(axis="x", rotation=45)
    for bar, value in zip(bars, normalized, strict=True):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value,
            f"{value:.1f}M",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    fig.savefig(output_dir / "normalized_intensity_by_well.png", dpi=150)
    plt.close(fig)

    nucleus_counts = [int(row["filtered_nucleus_count"]) for row in rows]
    corrected = [float(row["target_integrated_background_corrected"]) / 1_000_000_000 for row in rows]
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ax.scatter(nucleus_counts, corrected, color="#059669", s=48)
    for index, (row, x_value, y_value) in enumerate(
        zip(rows, nucleus_counts, corrected, strict=True)
    ):
        y_offset = 7 if index % 2 == 0 else -12
        ax.annotate(row["well_id"], (x_value, y_value), xytext=(5, y_offset), textcoords="offset points")
    ax.set_title("Target signal versus nucleus count")
    ax.set_xlabel("CH4 Cellpose nucleus count")
    ax.set_ylabel("background-corrected integrated CH2 intensity (billions)")
    ax.margins(x=0.12, y=0.10)
    fig.savefig(output_dir / "target_integrated_vs_nucleus_count.png", dpi=150)
    plt.close(fig)


def _write_run_metadata(
    *,
    logs_dir: Path,
    input_root: Path,
    counts_dir: Path,
    output_dir: Path,
    target_channel: str,
    dapi_channel: str,
    background_percentile: float,
    rows: list[dict[str, Any]],
) -> None:
    logs_dir.mkdir(parents=True, exist_ok=True)
    config = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "input_root": str(input_root),
        "counts_dir": str(counts_dir),
        "output_dir": str(output_dir),
        "target_channel": target_channel,
        "dapi_channel": dapi_channel,
        "measurement": {
            "target_area": "full_image",
            "background_method": "percentile",
            "background_percentile": background_percentile,
            "background_correction": "clip(target - background, lower=0)",
            "primary_endpoint": "target_integrated_intensity_per_DAPI_positive_nucleus",
            "denominator": "filtered_DAPI_positive_nucleus_count",
            "endpoint_formula": (
                "target_integrated_background_corrected / filtered_DAPI_positive_nucleus_count"
            ),
            "axes_policy": (
                "singleton axes are squeezed; 2-D images are used directly; RGB/YXS exports "
                "with exactly one active sample use that active sample"
            ),
            "z_projection": "none",
        },
        "denominator_counts": {
            "source_counts_dir": str(counts_dir),
            "required_backend": "cellpose",
            "required_channel": dapi_channel,
            "segmentation_parameters_source": str(counts_dir / "logs" / "config_resolved.yaml"),
        },
        "channel_identity_confirmed": False,
        "warnings": [
            "channel_identity_unconfirmed",
            "rendered_rgb_source",
            "no_nucleus_filtering_applied",
        ],
        "outputs": {
            "image_level_summary_csv": str(output_dir / "summaries" / "image_level_summary.csv"),
            "well_level_summary_csv": str(output_dir / "summaries" / "well_level_summary.csv"),
            "qc_dir": str(output_dir / "qc"),
            "plots_dir": str(output_dir / "plots"),
            "qc_contact_sheet": str(output_dir / "qc_contact_sheet.png"),
        },
        "image_rows": rows,
    }
    (logs_dir / "config_resolved.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
    )
    (logs_dir / "run_log.txt").write_text(
        "\n".join(
            [
                f"generated_at_utc: {config['generated_at_utc']}",
                f"input_root: {input_root}",
                f"counts_dir: {counts_dir}",
                f"output_dir: {output_dir}",
                f"target_channel: {target_channel}",
                f"dapi_channel: {dapi_channel}",
                f"background_percentile: {background_percentile}",
                f"images_processed: {len(rows)}",
                "channel_identity_confirmed: false",
                "warnings: channel_identity_unconfirmed;rendered_rgb_source",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_summary_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=IMAGE_LEVEL_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _read_count_rows(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Nucleus count summary does not exist: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No rows found in nucleus count summary: {path}")
    return {row["image_id"]: row for row in rows}


def _validate_count_row(
    *,
    count_row: dict[str, str],
    expected_channel: str,
    expected_dapi_path: Path,
    counts_dir: Path,
) -> None:
    image_id = count_row["image_id"]
    if count_row.get("backend") != "cellpose":
        raise ValueError(
            f"Expected cellpose count row for {image_id}, got backend={count_row.get('backend')}"
        )
    if count_row.get("channel_id", "").upper() != expected_channel:
        raise ValueError(
            f"Expected {expected_channel} count row for {image_id}, "
            f"got channel_id={count_row.get('channel_id')}"
        )
    count_input_path = _resolve_reference(count_row.get("input_path", ""), counts_dir)
    expected_resolved = expected_dapi_path.resolve()
    if count_input_path.exists() and count_input_path.resolve() != expected_resolved:
        raise ValueError(
            f"Count row for {image_id} points to {count_input_path}, expected "
            f"{expected_resolved}"
        )
    if not count_input_path.exists() and count_input_path.name != expected_dapi_path.name:
        raise ValueError(
            f"Count row for {image_id} points to unresolved {count_input_path.name}, expected "
            f"{expected_dapi_path.name}"
        )


def _mask_label_count(mask: np.ndarray) -> int:
    return int(np.count_nonzero(np.unique(mask)))


def _resolve_reference(value: str, base_dir: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    candidates = [
        path,
        base_dir / path,
        base_dir.parent / path,
        base_dir.parent.parent / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[-1]


def _scale_for_display(image: np.ndarray) -> np.ndarray:
    low, high = np.percentile(image, [1, 99.8])
    return np.clip((image.astype(np.float32) - low) / max(high - low, 1), 0, 1)


def _saturation_fraction(image: np.ndarray) -> float:
    if np.issubdtype(image.dtype, np.integer):
        saturation_value = np.iinfo(image.dtype).max
        return float(np.count_nonzero(image == saturation_value) / image.size)
    return float(np.count_nonzero(image >= np.nanmax(image)) / image.size)


def _format_percentile(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else str(value).replace(".", "p")
