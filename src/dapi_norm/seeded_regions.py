from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
import re
import textwrap
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage as ndi
from skimage import filters, measure, segmentation
import tifffile
import yaml

from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.pi_simple_summary import ImagePair


SEEDED_REGION_COLUMNS = [
    "image_id",
    "source_id",
    "method",
    "foreground_method",
    "foreground_threshold",
    "background_value_per_px",
    "image_area_px",
    "foreground_area_px",
    "foreground_fraction",
    "seeded_region_area_px",
    "seeded_region_fraction",
    "non_seeded_area_px",
    "non_seeded_area_fraction",
    "dapi_positive_nucleus_count",
    "seeded_region_integrated_raw",
    "seeded_region_mean_raw",
    "seeded_region_integrated_background_corrected",
    "seeded_region_intensity_per_DAPI_positive_nucleus",
    "unseeded_foreground_area_px",
    "unseeded_foreground_fraction",
    "foreground_components",
    "foreground_components_with_seed",
    "qc_status",
    "qc_flags",
    "mask_path",
    "qc_panel_path",
    "warnings",
]

LOW_NUCLEUS_COUNT_WARNING_THRESHOLD = 10
SIZEABLE_UNSEEDED_FOREGROUND_FRACTION_THRESHOLD = 0.10
HIGH_UNSEEDED_FOREGROUND_FRACTION_THRESHOLD = 0.15
LOW_SEEDED_COMPONENT_COVERAGE_THRESHOLD = 0.5
LOW_NUCLEUS_LARGE_AREA_FRACTION_THRESHOLD = 0.2
NEAR_FULL_FIELD_SEEDED_REGION_FRACTION_THRESHOLD = 0.8


@dataclass(frozen=True)
class ForegroundStats:
    method: str
    threshold: float
    image_area_px: int
    foreground_area_px: int
    foreground_fraction: float
    min_size: int
    fill_holes: bool


@dataclass(frozen=True)
class SeededSegmentationStats:
    method: str
    nucleus_labels: int
    foreground_area_px: int
    seeded_region_area_px: int
    unseeded_foreground_area_px: int
    foreground_components: int
    foreground_components_with_seed: int


def build_ch2_foreground_mask(
    image: np.ndarray,
    *,
    method: str = "li",
    min_size: int = 128,
    fill_holes: bool = True,
) -> tuple[np.ndarray, ForegroundStats]:
    arr = np.asarray(image, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"Expected a 2-D CH2 image, got shape {arr.shape}")

    method_key = method.lower().strip()
    threshold = _threshold_for_method(arr, method_key)
    mask = arr > threshold
    if min_size > 1:
        mask = _remove_small_binary_objects(mask, min_size)
    if fill_holes:
        mask = ndi.binary_fill_holes(mask)
    mask = np.asarray(mask, dtype=bool)
    foreground_area = int(np.count_nonzero(mask))
    image_area = int(mask.size)
    return mask, ForegroundStats(
        method=method_key,
        threshold=float(threshold),
        image_area_px=image_area,
        foreground_area_px=foreground_area,
        foreground_fraction=float(foreground_area / image_area) if image_area else 0.0,
        min_size=int(min_size),
        fill_holes=bool(fill_holes),
    )


def segment_seeded_regions(
    ch2_image: np.ndarray,
    nuclei_mask: np.ndarray,
    foreground_mask: np.ndarray,
) -> tuple[np.ndarray, SeededSegmentationStats]:
    ch2 = np.asarray(ch2_image, dtype=np.float64)
    nuclei = np.asarray(nuclei_mask)
    foreground = np.asarray(foreground_mask, dtype=bool)
    if ch2.ndim != 2 or nuclei.ndim != 2 or foreground.ndim != 2:
        raise ValueError("CH2 image, nuclei mask, and foreground mask must all be 2-D")
    if ch2.shape != nuclei.shape or ch2.shape != foreground.shape:
        raise ValueError(
            "CH2 image, nuclei mask, and foreground mask must have matching shapes: "
            f"ch2={ch2.shape}, nuclei={nuclei.shape}, foreground={foreground.shape}"
        )
    if np.any(nuclei < 0):
        raise ValueError("Nuclei mask must not contain negative labels")

    foreground_components = measure.label(foreground, connectivity=1)
    component_ids = np.unique(foreground_components[foreground_components > 0])
    components_with_seed = 0
    seeded_component_mask = np.zeros_like(foreground, dtype=bool)
    for component_id in component_ids:
        component = foreground_components == component_id
        if np.any(nuclei[component] > 0):
            components_with_seed += 1
            seeded_component_mask |= component

    markers = nuclei.astype(np.int64, copy=False)
    elevation = _watershed_elevation(ch2)
    labels = segmentation.watershed(elevation, markers=markers, mask=seeded_component_mask)
    labels = labels.astype(np.uint32, copy=False)
    foreground_area = int(np.count_nonzero(foreground))
    seeded_area = int(np.count_nonzero(labels))
    unseeded_area = int(foreground_area - np.count_nonzero(seeded_component_mask))
    return labels, SeededSegmentationStats(
        method="seeded_intensity_watershed",
        nucleus_labels=_count_nonzero_labels(nuclei),
        foreground_area_px=foreground_area,
        seeded_region_area_px=seeded_area,
        unseeded_foreground_area_px=unseeded_area,
        foreground_components=int(len(component_ids)),
        foreground_components_with_seed=int(components_with_seed),
    )


def segment_seeded_regions_random_walker(
    ch2_image: np.ndarray,
    nuclei_mask: np.ndarray,
    foreground_mask: np.ndarray,
    *,
    beta: float = 90.0,
) -> tuple[np.ndarray, SeededSegmentationStats]:
    ch2 = np.asarray(ch2_image, dtype=np.float64)
    nuclei = np.asarray(nuclei_mask)
    foreground = np.asarray(foreground_mask, dtype=bool)
    if ch2.ndim != 2 or nuclei.ndim != 2 or foreground.ndim != 2:
        raise ValueError("CH2 image, nuclei mask, and foreground mask must all be 2-D")
    if ch2.shape != nuclei.shape or ch2.shape != foreground.shape:
        raise ValueError(
            "CH2 image, nuclei mask, and foreground mask must have matching shapes: "
            f"ch2={ch2.shape}, nuclei={nuclei.shape}, foreground={foreground.shape}"
        )
    if np.any(nuclei < 0):
        raise ValueError("Nuclei mask must not contain negative labels")
    if beta <= 0:
        raise ValueError(f"random walker beta must be > 0, got {beta}")

    foreground_components = measure.label(foreground, connectivity=1)
    component_ids = np.unique(foreground_components[foreground_components > 0])
    labels = np.zeros_like(nuclei, dtype=np.uint32)
    components_with_seed = 0
    seeded_component_area = 0
    for component_id in component_ids:
        component = foreground_components == component_id
        seed_labels = np.unique(nuclei[component])
        seed_labels = seed_labels[seed_labels > 0]
        if len(seed_labels) == 0:
            continue
        components_with_seed += 1
        seeded_component_area += int(np.count_nonzero(component))
        if len(seed_labels) == 1:
            labels[component] = np.uint32(seed_labels[0])
            continue
        _assign_component_by_random_walker(
            output_labels=labels,
            ch2=ch2,
            nuclei=nuclei,
            component=component,
            seed_labels=seed_labels,
            beta=beta,
        )

    foreground_area = int(np.count_nonzero(foreground))
    seeded_area = int(np.count_nonzero(labels))
    return labels, SeededSegmentationStats(
        method="seeded_intensity_random_walker",
        nucleus_labels=_count_nonzero_labels(nuclei),
        foreground_area_px=foreground_area,
        seeded_region_area_px=seeded_area,
        unseeded_foreground_area_px=int(foreground_area - seeded_component_area),
        foreground_components=int(len(component_ids)),
        foreground_components_with_seed=int(components_with_seed),
    )


def segment_seeded_regions_propagation(
    ch2_image: np.ndarray,
    nuclei_mask: np.ndarray,
    foreground_mask: np.ndarray,
    *,
    regularization_factor: float = 0.05,
) -> tuple[np.ndarray, SeededSegmentationStats]:
    """Run CellProfiler-style propagation from DAPI labels inside CH2 foreground.

    The underlying implementation is `centrosome.propagate`, the propagation
    routine used by CellProfiler. The returned labels remain exploratory because
    CH2/aSMA is the biological endpoint and also contributes to the region
    definition.
    """
    ch2 = np.asarray(ch2_image, dtype=np.float64)
    nuclei = np.asarray(nuclei_mask)
    foreground = np.asarray(foreground_mask, dtype=bool)
    if ch2.ndim != 2 or nuclei.ndim != 2 or foreground.ndim != 2:
        raise ValueError("CH2 image, nuclei mask, and foreground mask must all be 2-D")
    if ch2.shape != nuclei.shape or ch2.shape != foreground.shape:
        raise ValueError(
            "CH2 image, nuclei mask, and foreground mask must have matching shapes: "
            f"ch2={ch2.shape}, nuclei={nuclei.shape}, foreground={foreground.shape}"
        )
    if np.any(nuclei < 0):
        raise ValueError("Nuclei mask must not contain negative labels")
    if regularization_factor <= 0:
        raise ValueError(
            "propagation regularization factor must be > 0, "
            f"got {regularization_factor}"
        )

    foreground_components = measure.label(foreground, connectivity=1)
    component_ids = np.unique(foreground_components[foreground_components > 0])
    components_with_seed = 0
    seeded_component_mask = np.zeros_like(foreground, dtype=bool)
    for component_id in component_ids:
        component = foreground_components == component_id
        if np.any(nuclei[component] > 0):
            components_with_seed += 1
            seeded_component_mask |= component

    propagation_mask = seeded_component_mask | (nuclei > 0)
    propagation_image = _propagation_data(ch2)
    labels, _distances = _cellprofiler_propagate(
        propagation_image,
        nuclei.astype(np.int32, copy=False),
        propagation_mask,
        weight=float(regularization_factor),
    )
    labels = labels.astype(np.uint32, copy=False)
    labels[~seeded_component_mask] = 0
    foreground_area = int(np.count_nonzero(foreground))
    seeded_area = int(np.count_nonzero(labels))
    unseeded_area = int(np.count_nonzero(foreground & (labels == 0)))
    return labels, SeededSegmentationStats(
        method="seeded_intensity_propagation",
        nucleus_labels=_count_nonzero_labels(nuclei),
        foreground_area_px=foreground_area,
        seeded_region_area_px=seeded_area,
        unseeded_foreground_area_px=unseeded_area,
        foreground_components=int(len(component_ids)),
        foreground_components_with_seed=int(components_with_seed),
    )


def measure_seeded_region_image(
    image_id: str,
    ch2_image: np.ndarray,
    nuclei_mask: np.ndarray,
    seeded_labels: np.ndarray,
    *,
    background_value: float = 0.0,
    source_id: str | None = None,
    foreground_stats: ForegroundStats | None = None,
    segmentation_stats: SeededSegmentationStats | None = None,
    mask_path: Path | None = None,
    qc_panel_path: Path | None = None,
) -> dict[str, Any]:
    ch2 = np.asarray(ch2_image, dtype=np.float64)
    nuclei = np.asarray(nuclei_mask)
    labels = np.asarray(seeded_labels)
    if ch2.shape != nuclei.shape or ch2.shape != labels.shape:
        raise ValueError(
            "CH2 image, nuclei mask, and seeded labels must have matching shapes: "
            f"ch2={ch2.shape}, nuclei={nuclei.shape}, labels={labels.shape}"
        )

    region_mask = labels > 0
    image_area = int(ch2.size)
    seeded_area = int(np.count_nonzero(region_mask))
    non_seeded_area = int(image_area - seeded_area)
    raw_integrated = float(np.sum(ch2[region_mask]))
    mean_raw = float(np.mean(ch2[region_mask])) if seeded_area else float("nan")
    corrected = np.clip(ch2[region_mask] - float(background_value), 0, None)
    corrected_integrated = float(np.sum(corrected))
    nucleus_count = _count_nonzero_labels(nuclei)
    per_nucleus = corrected_integrated / nucleus_count if nucleus_count else float("nan")
    foreground_area = (
        int(foreground_stats.foreground_area_px)
        if foreground_stats is not None
        else int(segmentation_stats.foreground_area_px) if segmentation_stats is not None else seeded_area
    )
    unseeded_area = (
        int(segmentation_stats.unseeded_foreground_area_px)
        if segmentation_stats is not None
        else int(max(foreground_area - seeded_area, 0))
    )
    warnings = [
        "exploratory_asma_associated_region_not_validated_cell_mask",
        "ch2_endpoint_signal_used_to_define_region",
    ]
    if nucleus_count == 0:
        warnings.append("zero_nucleus_count")
    elif nucleus_count < LOW_NUCLEUS_COUNT_WARNING_THRESHOLD:
        warnings.append("low_nucleus_count_qc_required")
    if unseeded_area > 0:
        warnings.append("foreground_components_without_dapi_seed_excluded")
        warnings.append("unresolved_target_foreground_not_measured")

    row = {
        "image_id": image_id,
        "source_id": source_id or image_id,
        "method": (
            segmentation_stats.method if segmentation_stats is not None else "seeded_intensity_watershed"
        ),
        "foreground_method": foreground_stats.method if foreground_stats is not None else "",
        "foreground_threshold": (
            foreground_stats.threshold if foreground_stats is not None else float("nan")
        ),
        "background_value_per_px": float(background_value),
        "image_area_px": image_area,
        "foreground_area_px": foreground_area,
        "foreground_fraction": float(foreground_area / image_area) if image_area else 0.0,
        "seeded_region_area_px": seeded_area,
        "seeded_region_fraction": float(seeded_area / image_area) if image_area else 0.0,
        "non_seeded_area_px": non_seeded_area,
        "non_seeded_area_fraction": float(non_seeded_area / image_area) if image_area else 0.0,
        "dapi_positive_nucleus_count": nucleus_count,
        "seeded_region_integrated_raw": raw_integrated,
        "seeded_region_mean_raw": mean_raw,
        "seeded_region_integrated_background_corrected": corrected_integrated,
        "seeded_region_intensity_per_DAPI_positive_nucleus": per_nucleus,
        "unseeded_foreground_area_px": unseeded_area,
        "unseeded_foreground_fraction": (
            float(unseeded_area / foreground_area) if foreground_area else 0.0
        ),
        "foreground_components": (
            int(segmentation_stats.foreground_components) if segmentation_stats is not None else 0
        ),
        "foreground_components_with_seed": (
            int(segmentation_stats.foreground_components_with_seed)
            if segmentation_stats is not None
            else 0
        ),
        "mask_path": "" if mask_path is None else str(mask_path),
        "qc_panel_path": "" if qc_panel_path is None else str(qc_panel_path),
        "warnings": ";".join(warnings),
    }
    qc_status, qc_flags = score_seeded_region_qc(row)
    row["qc_status"] = qc_status
    row["qc_flags"] = ";".join(qc_flags)
    return row


def score_seeded_region_qc(row: dict[str, Any]) -> tuple[str, list[str]]:
    flags: list[str] = ["not_validated_whole_cell_mask"]
    nucleus_count = int(row.get("dapi_positive_nucleus_count", 0))
    seeded_fraction = float(row.get("seeded_region_fraction", 0.0))
    unseeded_foreground_fraction = float(row.get("unseeded_foreground_fraction", 0.0))
    foreground_components = int(row.get("foreground_components", 0))
    foreground_components_with_seed = int(row.get("foreground_components_with_seed", 0))

    if nucleus_count == 0:
        flags.append("zero_nucleus_count")
    elif nucleus_count < LOW_NUCLEUS_COUNT_WARNING_THRESHOLD:
        flags.append("low_nucleus_count")
    if (
        nucleus_count < LOW_NUCLEUS_COUNT_WARNING_THRESHOLD
        and seeded_fraction >= LOW_NUCLEUS_LARGE_AREA_FRACTION_THRESHOLD
    ):
        flags.append("large_seeded_area_with_low_nucleus_count")
    if seeded_fraction >= NEAR_FULL_FIELD_SEEDED_REGION_FRACTION_THRESHOLD:
        flags.append("near_full_field_seeded_region_fraction")
    if unseeded_foreground_fraction >= SIZEABLE_UNSEEDED_FOREGROUND_FRACTION_THRESHOLD:
        flags.append("sizeable_unseeded_target_foreground")
    if unseeded_foreground_fraction >= HIGH_UNSEEDED_FOREGROUND_FRACTION_THRESHOLD:
        flags.append("high_unseeded_foreground_fraction")
    if foreground_components > 0:
        seeded_component_fraction = foreground_components_with_seed / foreground_components
        if seeded_component_fraction < LOW_SEEDED_COMPONENT_COVERAGE_THRESHOLD:
            flags.append("low_fraction_of_foreground_components_with_dapi_seed")

    reject_flags = {
        "zero_nucleus_count",
        "low_nucleus_count",
        "large_seeded_area_with_low_nucleus_count",
        "near_full_field_seeded_region_fraction",
        "low_fraction_of_foreground_components_with_dapi_seed",
    }
    if any(flag in flags for flag in reject_flags):
        status = "reject_qc_failure"
    elif len(flags) > 1:
        status = "needs_manual_review"
    else:
        status = "reviewable_not_validated"
    return status, flags


def write_seeded_region_qc_panel(
    *,
    image_id: str,
    ch2_image: np.ndarray,
    ch4_image: np.ndarray,
    nuclei_mask: np.ndarray,
    foreground_mask: np.ndarray,
    seeded_labels: np.ndarray,
    output_path: Path,
    metrics: dict[str, Any],
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ch2_scaled = _scale_for_display(ch2_image)
    ch4_scaled = _scale_for_display(ch4_image)
    labels = np.asarray(seeded_labels)
    foreground = np.asarray(foreground_mask, dtype=bool)
    excluded = ~foreground
    region_boundaries = segmentation.find_boundaries(labels, mode="outer")
    nucleus_boundaries = segmentation.find_boundaries(nuclei_mask, mode="outer")

    seeded_overlay = np.dstack([ch2_scaled, ch2_scaled, ch2_scaled])
    seeded_overlay[excluded] = [0.05, 0.05, 0.05]
    seeded_overlay[region_boundaries] = [0.0, 1.0, 0.2]
    seeded_overlay[nucleus_boundaries] = [0.0, 0.75, 1.0]

    foreground_rgb = np.dstack([ch2_scaled, ch2_scaled, ch2_scaled])
    foreground_rgb[foreground] = np.maximum(foreground_rgb[foreground], [0.85, 0.12, 0.12])
    foreground_rgb[~foreground] = [0.0, 0.0, 0.0]

    ch4_rgb = np.dstack([ch4_scaled, ch4_scaled, ch4_scaled])
    ch4_rgb[nucleus_boundaries] = [0.0, 1.0, 0.2]
    _draw_centroid_crosses(ch4_rgb, nuclei_mask, color=np.array([0.0, 1.0, 0.2]))

    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
    flat_axes = axes.ravel()
    flat_axes[0].imshow(ch2_scaled, cmap="gray")
    flat_axes[0].set_title(f"{image_id} CH2/aSMA raw display")
    flat_axes[1].imshow(ch4_rgb)
    flat_axes[1].set_title(
        f"CH4/DAPI nuclei\ncount={metrics.get('dapi_positive_nucleus_count', 'NA')}"
    )
    flat_axes[2].imshow(foreground_rgb)
    flat_axes[2].set_title(
        "CH2 foreground mask\n"
        f"area={_format_percent(metrics.get('foreground_fraction', float('nan')))}"
    )
    flat_axes[3].imshow(labels, cmap="nipy_spectral")
    flat_axes[3].set_title("DAPI-seeded aSMA regions")
    flat_axes[4].imshow(seeded_overlay)
    flat_axes[4].set_title("Regions over CH2\ncyan nuclei, green region edges")
    flat_axes[5].imshow(_non_seeded_area_rgb(ch2_scaled, labels > 0, foreground))
    flat_axes[5].set_title(
        "Unseeded CH2 foreground in magenta\n"
        "unseeded foreground="
        f"{_format_percent(metrics.get('unseeded_foreground_fraction', float('nan')))}"
    )
    for ax in flat_axes:
        ax.axis("off")
    fig.suptitle(
        "Exploratory aSMA-associated region QC; not a validated whole-cell mask\n"
        f"raw={_format_scientific(metrics.get('seeded_region_integrated_raw'))}, "
        "bg-corrected/nucleus="
        f"{_format_scientific(metrics.get('seeded_region_intensity_per_DAPI_positive_nucleus'))}",
        fontsize=12,
    )
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def write_seeded_region_crop_panel(*, crops: list[dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not crops:
        raise ValueError("Cannot write a crop panel without crops")

    fig, axes = plt.subplots(
        len(crops),
        4,
        figsize=(15.5, max(3.5, 3.4 * len(crops))),
        squeeze=False,
        constrained_layout=True,
        width_ratios=[1.0, 1.0, 1.0, 0.72],
    )
    for row_index, crop in enumerate(crops):
        y0, x0, y1, x1 = crop["box"]
        ch2 = _crop_array(np.asarray(crop["ch2_image"]), y0, x0, y1, x1)
        ch4 = _crop_array(np.asarray(crop["ch4_image"]), y0, x0, y1, x1)
        nuclei = _crop_array(np.asarray(crop["nuclei_mask"]), y0, x0, y1, x1)
        labels = _crop_array(np.asarray(crop["seeded_labels"]), y0, x0, y1, x1)
        ch2_scaled = _scale_for_display(ch2)
        ch4_scaled = _scale_for_display(ch4)

        ch4_rgb = np.dstack([ch4_scaled, ch4_scaled, ch4_scaled])
        nucleus_boundaries = segmentation.find_boundaries(nuclei, mode="outer")
        ch4_rgb[nucleus_boundaries] = [0.0, 1.0, 0.2]
        _draw_centroid_crosses(ch4_rgb, nuclei, color=np.array([0.0, 1.0, 0.2]))

        overlay = _retained_region_overlay_rgb(ch2_scaled, labels > 0)
        region_boundaries = segmentation.find_boundaries(labels, mode="outer")
        overlay[region_boundaries] = [0.0, 1.0, 0.2]
        overlay[nucleus_boundaries] = [0.0, 0.75, 1.0]

        title_prefix = str(crop.get("image_id", f"crop {row_index + 1}"))
        caption = str(crop.get("caption", ""))
        axes[row_index, 0].imshow(ch2_scaled, cmap="gray")
        axes[row_index, 0].set_title(f"{title_prefix} CH2/aSMA")
        axes[row_index, 1].imshow(ch4_rgb)
        axes[row_index, 1].set_title("CH4/DAPI nuclei")
        axes[row_index, 2].imshow(overlay)
        axes[row_index, 2].set_title("Light-green retained region over CH2")
        if caption:
            axes[row_index, 3].text(
                0.0,
                0.98,
                _format_crop_caption(caption),
                ha="left",
                va="top",
                fontsize=8,
                transform=axes[row_index, 3].transAxes,
                color="#212529",
            )
        axes[row_index, 3].text(
            0.0,
            0.08,
            "Exploratory QC only",
            ha="left",
            va="bottom",
            fontsize=8,
            transform=axes[row_index, 3].transAxes,
            color="#495057",
        )
        for col_index in range(4):
            axes[row_index, col_index].axis("off")
    fig.suptitle(
        "Close-up QC crops: exploratory DAPI-seeded aSMA-associated regions",
        fontsize=13,
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def run_seeded_region_batch(
    *,
    image_pairs: Iterable[ImagePair],
    mask_lookup: dict[str, Path],
    output_dir: Path | str,
    foreground_method: str = "li",
    background_value: float = 0.0,
    min_size: int = 128,
    segmentation_method: str = "watershed",
    random_walker_beta: float = 90.0,
    propagation_regularization: float = 0.05,
) -> list[dict[str, Any]]:
    output_path = Path(output_dir)
    masks_dir = output_path / "masks"
    qc_dir = output_path / "qc"
    summaries_dir = output_path / "summaries"
    logs_dir = output_path / "logs"
    masks_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    summaries_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    qc_paths: list[Path] = []
    image_records: list[dict[str, Any]] = []
    for pair in image_pairs:
        mask_path = _mask_path_for_pair(pair, mask_lookup)
        ch2_image, ch2_pages = read_primary_intensity_plane(pair.ch2_path)
        ch4_image, ch4_pages = read_primary_intensity_plane(pair.ch4_path)
        nuclei_mask = np.asarray(tifffile.imread(mask_path))
        if nuclei_mask.shape != ch2_image.shape:
            raise ValueError(
                f"{pair.source_id} nuclei mask shape {nuclei_mask.shape} does not match "
                f"CH2 shape {ch2_image.shape}"
            )
        foreground_mask, foreground_stats = build_ch2_foreground_mask(
            ch2_image,
            method=foreground_method,
            min_size=min_size,
        )
        seeded_labels, segmentation_stats = _segment_seeded_regions_by_method(
            ch2_image,
            nuclei_mask,
            foreground_mask,
            method=segmentation_method,
            random_walker_beta=random_walker_beta,
            propagation_regularization=propagation_regularization,
        )
        safe_id = _safe_output_stem(pair.source_id)
        seeded_mask_path = masks_dir / f"{safe_id}_{segmentation_stats.method}_labels.tif"
        qc_path = qc_dir / f"{safe_id}_seeded_region_qc.png"
        tifffile.imwrite(seeded_mask_path, seeded_labels, photometric="minisblack")
        row = measure_seeded_region_image(
            pair.location,
            ch2_image,
            nuclei_mask,
            seeded_labels,
            background_value=background_value,
            source_id=pair.source_id,
            foreground_stats=foreground_stats,
            segmentation_stats=segmentation_stats,
            mask_path=seeded_mask_path,
            qc_panel_path=qc_path,
        )
        write_seeded_region_qc_panel(
            image_id=pair.source_id,
            ch2_image=ch2_image,
            ch4_image=ch4_image,
            nuclei_mask=nuclei_mask,
            foreground_mask=foreground_mask,
            seeded_labels=seeded_labels,
            output_path=qc_path,
            metrics=row,
        )
        rows.append(row)
        qc_paths.append(qc_path)
        image_records.append(
            {
                "source_id": pair.source_id,
                "location": pair.location,
                "ch2_path": str(pair.ch2_path),
                "ch4_path": str(pair.ch4_path),
                "nuclei_mask_path": str(mask_path),
                "ch2_page_count": int(ch2_pages),
                "ch4_page_count": int(ch4_pages),
                "image_shape": [int(ch2_image.shape[0]), int(ch2_image.shape[1])],
                "foreground": asdict(foreground_stats),
                "segmentation": asdict(segmentation_stats),
            }
        )

    _write_summary_csv(summaries_dir / "seeded_region_image_metrics.csv", rows)
    _write_summary_plots(rows, output_path / "plots")
    write_seeded_region_contact_sheet(qc_paths, output_path / "qc_contact_sheet.png")
    _write_seeded_region_metadata(
        logs_dir=logs_dir,
        output_dir=output_path,
        foreground_method=foreground_method,
        background_value=background_value,
        min_size=min_size,
        segmentation_method=str(rows[0]["method"]) if rows else segmentation_method,
        random_walker_beta=random_walker_beta,
        propagation_regularization=propagation_regularization,
        rows=rows,
        image_records=image_records,
    )
    return rows


def load_mask_lookup_from_counts_root(counts_root: Path | str) -> dict[str, Path]:
    root = Path(counts_root)
    if not root.exists():
        raise FileNotFoundError(f"Counts root does not exist: {root}")
    lookup: dict[str, Path] = {}
    summary_path = root / "summaries" / "nucleus_counts.csv"
    if summary_path.exists():
        with summary_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if "image_id" not in (reader.fieldnames or []) or "mask_path" not in (
                reader.fieldnames or []
            ):
                raise ValueError(f"Expected image_id and mask_path columns in {summary_path}")
            for row in reader:
                image_id = row["image_id"].strip()
                mask_path = _resolve_reference(row["mask_path"], root)
                _add_lookup_aliases(lookup, image_id, mask_path)

    discovered_by_xy: dict[str, list[Path]] = {}
    for mask_path in sorted(root.rglob("*_labels.tif")):
        match = re.search(r"(XY\d+)", mask_path.name, flags=re.IGNORECASE)
        if match is not None:
            discovered_by_xy.setdefault(match.group(1).upper(), []).append(mask_path)
    duplicates = {
        image_id: paths
        for image_id, paths in discovered_by_xy.items()
        if image_id.upper() not in lookup and len(paths) > 1
    }
    if duplicates:
        details = "; ".join(
            f"{image_id}: " + ", ".join(str(path) for path in paths)
            for image_id, paths in sorted(duplicates.items())
        )
        raise ValueError(
            "Duplicate nuclei masks found for unqualified XY image IDs under counts root. "
            "Use a sample-specific counts root or a nucleus_counts.csv summary with exact paths. "
            + details
        )
    for image_id, paths in sorted(discovered_by_xy.items()):
        if image_id.upper() not in lookup:
            _add_lookup_aliases(lookup, image_id, paths[0], overwrite=False)

    if not lookup:
        raise ValueError(f"No nuclei label masks found under {root}")
    return lookup


def _segment_seeded_regions_by_method(
    ch2_image: np.ndarray,
    nuclei_mask: np.ndarray,
    foreground_mask: np.ndarray,
    *,
    method: str,
    random_walker_beta: float,
    propagation_regularization: float,
) -> tuple[np.ndarray, SeededSegmentationStats]:
    method_key = method.lower().strip().replace("-", "_")
    if method_key in {"watershed", "seeded_intensity_watershed"}:
        return segment_seeded_regions(ch2_image, nuclei_mask, foreground_mask)
    if method_key in {"random_walker", "seeded_intensity_random_walker"}:
        return segment_seeded_regions_random_walker(
            ch2_image,
            nuclei_mask,
            foreground_mask,
            beta=random_walker_beta,
        )
    if method_key in {"propagation", "cellprofiler_propagation", "seeded_intensity_propagation"}:
        return segment_seeded_regions_propagation(
            ch2_image,
            nuclei_mask,
            foreground_mask,
            regularization_factor=propagation_regularization,
        )
    raise ValueError(
        "Unsupported seeded-region segmentation method "
        f"{method!r}. Use 'watershed', 'random_walker', or 'propagation'."
    )


def write_seeded_region_contact_sheet(qc_paths: list[Path], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not qc_paths:
        raise ValueError("Cannot write a contact sheet without QC panel paths")
    thumb_width = 760
    padding = 24
    label_height = 34
    columns = 2
    thumbs: list[tuple[str, Image.Image]] = []
    for path in qc_paths:
        with Image.open(path) as image:
            rgb_image = image.convert("RGB")
            scale = thumb_width / rgb_image.width
            thumb = rgb_image.resize((thumb_width, int(rgb_image.height * scale)))
        thumbs.append((path.stem.replace("_seeded_region_qc", ""), thumb))

    cell_height = max(thumb.height for _label, thumb in thumbs) + label_height
    rows = (len(thumbs) + columns - 1) // columns
    sheet = Image.new(
        "RGB",
        (columns * thumb_width + (columns + 1) * padding, rows * cell_height + (rows + 1) * padding),
        "white",
    )
    draw = ImageDraw.Draw(sheet)
    for index, (label, thumb) in enumerate(thumbs):
        row, column = divmod(index, columns)
        x_pos = padding + column * (thumb_width + padding)
        y_pos = padding + row * (cell_height + padding)
        draw.text((x_pos, y_pos), label, fill=(0, 0, 0))
        sheet.paste(thumb, (x_pos, y_pos + label_height))
    sheet.save(output_path)


def _threshold_for_method(arr: np.ndarray, method: str) -> float:
    if np.nanmax(arr) <= np.nanmin(arr):
        return float(np.nanmax(arr))
    if method == "otsu":
        return float(filters.threshold_otsu(arr))
    if method == "li":
        return float(filters.threshold_li(arr))
    if method == "triangle":
        return float(filters.threshold_triangle(arr))
    if method.startswith("value:"):
        return float(method.split(":", 1)[1])
    raise ValueError(
        f"Unsupported foreground method {method!r}. Use 'li', 'otsu', 'triangle', or 'value:<n>'."
    )


def _remove_small_binary_objects(mask: np.ndarray, min_size: int) -> np.ndarray:
    labels, object_count = ndi.label(mask)
    if object_count == 0:
        return np.zeros_like(mask, dtype=bool)
    counts = np.bincount(labels.ravel())
    keep = counts >= min_size
    keep[0] = False
    return keep[labels]


def _watershed_elevation(ch2: np.ndarray) -> np.ndarray:
    arr = np.asarray(ch2, dtype=np.float64)
    low, high = np.percentile(arr, [1, 99.8])
    if high <= low:
        high = float(np.max(arr))
        low = float(np.min(arr))
    if high <= low:
        return np.zeros_like(arr, dtype=np.float64)
    scaled = np.clip((arr - low) / (high - low), 0, 1)
    return 1.0 - filters.gaussian(scaled, sigma=1.0, preserve_range=True)


def _assign_component_by_random_walker(
    *,
    output_labels: np.ndarray,
    ch2: np.ndarray,
    nuclei: np.ndarray,
    component: np.ndarray,
    seed_labels: np.ndarray,
    beta: float,
) -> None:
    ys, xs = np.where(component)
    y0 = int(np.min(ys))
    y1 = int(np.max(ys)) + 1
    x0 = int(np.min(xs))
    x1 = int(np.max(xs)) + 1
    local_component = component[y0:y1, x0:x1]
    local_nuclei = nuclei[y0:y1, x0:x1]
    local_ch2 = ch2[y0:y1, x0:x1]
    marker_lookup = {int(label): index + 1 for index, label in enumerate(seed_labels)}
    markers = np.full(local_component.shape, -1, dtype=np.int32)
    markers[local_component] = 0
    for source_label, marker_label in marker_lookup.items():
        markers[(local_nuclei == source_label) & local_component] = marker_label

    data = _random_walker_data(local_ch2)
    if float(np.nanstd(data[local_component])) == 0.0:
        random_walker_labels = segmentation.watershed(
            np.zeros_like(data),
            markers=np.where(markers > 0, markers, 0),
            mask=local_component,
        )
    else:
        random_walker_labels = segmentation.random_walker(
            data,
            markers,
            beta=float(beta),
            mode="cg_j",
            tol=1e-5,
            prob_tol=0.01,
            channel_axis=None,
        )
    local_output = output_labels[y0:y1, x0:x1]
    for source_label, marker_label in marker_lookup.items():
        local_output[(random_walker_labels == marker_label) & local_component] = np.uint32(
            source_label
        )


def _random_walker_data(ch2: np.ndarray) -> np.ndarray:
    arr = np.asarray(ch2, dtype=np.float64)
    low, high = np.percentile(arr, [1, 99.8])
    if high <= low:
        high = float(np.max(arr))
        low = float(np.min(arr))
    if high <= low:
        return np.zeros_like(arr, dtype=np.float64)
    scaled = np.clip((arr - low) / (high - low), 0, 1)
    return filters.gaussian(scaled, sigma=1.0, preserve_range=True)


def _propagation_data(ch2: np.ndarray) -> np.ndarray:
    return _random_walker_data(ch2)


def _cellprofiler_propagate(
    image: np.ndarray,
    labels: np.ndarray,
    mask: np.ndarray,
    *,
    weight: float,
) -> tuple[np.ndarray, np.ndarray]:
    try:
        from centrosome.propagate import propagate
    except ImportError as exc:  # pragma: no cover - dependency is part of project env
        raise RuntimeError(
            "The CellProfiler-style propagation method requires the 'centrosome' package. "
            "Install project dependencies before running segmentation_method='propagation'."
        ) from exc
    return propagate(image, labels, mask, weight)


def _count_nonzero_labels(mask: np.ndarray) -> int:
    labels = np.unique(mask)
    return int(np.count_nonzero(labels))


def _mask_path_for_pair(pair: ImagePair, mask_lookup: dict[str, Path]) -> Path:
    candidates = [
        pair.source_id,
        pair.source_id.upper(),
        pair.location,
        pair.location.upper(),
        _safe_output_stem(pair.source_id),
    ]
    for candidate in candidates:
        path = mask_lookup.get(candidate)
        if path is not None:
            return Path(path)
    raise KeyError(f"No nuclei mask path found for {pair.source_id}")


def _resolve_reference(value: str, root: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    candidates = [
        root / path,
        root.parent / path,
        root.parent.parent / path,
        path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def _add_lookup_aliases(
    lookup: dict[str, Path], key: str, path: Path, *, overwrite: bool = True
) -> None:
    aliases = {key, key.upper(), key.lower()}
    for alias in aliases:
        if overwrite or alias not in lookup:
            lookup[alias] = path


def _safe_output_stem(value: str) -> str:
    return value.replace("\\", "/").replace("/", "__").replace(" ", "_")


def _scale_for_display(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image, dtype=np.float32)
    low, high = np.percentile(arr, [1, 99.8])
    if high <= low:
        high = float(np.max(arr))
        low = float(np.min(arr))
    if high <= low:
        return np.zeros_like(arr, dtype=np.float32)
    return np.clip((arr - low) / (high - low), 0, 1)


def _non_seeded_area_rgb(
    ch2_scaled: np.ndarray,
    region_mask: np.ndarray,
    foreground_mask: np.ndarray,
) -> np.ndarray:
    rgb = np.dstack([ch2_scaled, ch2_scaled, ch2_scaled]) * 0.35
    foreground = np.asarray(foreground_mask, dtype=bool)
    region = np.asarray(region_mask, dtype=bool)
    unseeded_foreground = foreground & ~region
    rgb[region] = [0.0, 0.75, 0.18]
    rgb[unseeded_foreground] = [0.9, 0.0, 0.9]
    return np.clip(rgb, 0, 1)


def _retained_region_overlay_rgb(ch2_scaled: np.ndarray, region_mask: np.ndarray) -> np.ndarray:
    rgb = np.dstack([ch2_scaled, ch2_scaled, ch2_scaled]).astype(np.float32, copy=True)
    region = np.asarray(region_mask, dtype=bool)
    if np.any(region):
        overlay_color = np.array([0.25, 0.95, 0.35], dtype=np.float32)
        rgb[region] = (0.58 * rgb[region]) + (0.42 * overlay_color)
    return np.clip(rgb, 0, 1)


def _format_crop_caption(caption: str) -> str:
    parts = [part.strip() for part in str(caption).split("|", maxsplit=1)]
    if len(parts) == 1:
        return textwrap.fill(parts[0], width=34)
    status, flags = parts
    return f"QC status:\n{textwrap.fill(status, width=34)}\n\nQC flags:\n{textwrap.fill(flags, width=34)}"


def _crop_array(arr: np.ndarray, y0: int, x0: int, y1: int, x1: int) -> np.ndarray:
    y0 = max(0, int(y0))
    x0 = max(0, int(x0))
    y1 = min(arr.shape[0], int(y1))
    x1 = min(arr.shape[1], int(x1))
    if y1 <= y0 or x1 <= x0:
        raise ValueError(f"Invalid crop box after clipping: {(y0, x0, y1, x1)}")
    return arr[y0:y1, x0:x1]


def _draw_centroid_crosses(rgb: np.ndarray, labels: np.ndarray, *, color: np.ndarray) -> None:
    for region in measure.regionprops(np.asarray(labels, dtype=np.int64)):
        y, x = region.centroid
        yy = int(round(y))
        xx = int(round(x))
        y0 = max(0, yy - 3)
        y1 = min(rgb.shape[0], yy + 4)
        x0 = max(0, xx - 3)
        x1 = min(rgb.shape[1], xx + 4)
        rgb[yy, x0:x1] = color
        rgb[y0:y1, xx] = color


def _format_scientific(value: Any) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(numeric):
        return "NA"
    return f"{numeric:.2e}"


def _format_percent(value: Any) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(numeric):
        return "NA"
    if 0 < abs(numeric) < 0.001:
        return "<0.1%"
    return f"{numeric:.1%}"


def _write_summary_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SEEDED_REGION_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_summary_plots(rows: list[dict[str, Any]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    sorted_rows = sorted(rows, key=lambda row: str(row["source_id"]))
    labels = [str(row["source_id"]) for row in sorted_rows]
    fractions = [float(row["seeded_region_fraction"]) for row in sorted_rows]
    per_nucleus = [
        float(row["seeded_region_intensity_per_DAPI_positive_nucleus"]) / 1_000_000
        for row in sorted_rows
    ]

    fig, ax = plt.subplots(figsize=(max(7, len(rows) * 0.7), 4.5), constrained_layout=True)
    ax.bar(labels, fractions, color="#0f766e")
    ax.set_ylabel("Fraction of image assigned to seeded aSMA region")
    ax.set_title("Exploratory seeded-region retained area")
    ax.tick_params(axis="x", rotation=45)
    fig.savefig(output_dir / "seeded_region_area_fraction.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(max(7, len(rows) * 0.7), 4.5), constrained_layout=True)
    ax.bar(labels, per_nucleus, color="#9333ea")
    ax.set_ylabel("Background-corrected seeded CH2 per nucleus (millions)")
    ax.set_title("QC only: seeded-region aSMA burden per DAPI-positive nucleus")
    ax.tick_params(axis="x", rotation=45)
    for index, row in enumerate(sorted_rows):
        nucleus_count = int(row["dapi_positive_nucleus_count"])
        if nucleus_count < LOW_NUCLEUS_COUNT_WARNING_THRESHOLD:
            ax.text(
                index,
                per_nucleus[index],
                f"n={nucleus_count}\nQC",
                ha="center",
                va="bottom",
                fontsize=8,
                color="#b91c1c",
            )
    fig.savefig(output_dir / "seeded_region_intensity_per_nucleus.png", dpi=150)
    plt.close(fig)


def _write_seeded_region_metadata(
    *,
    logs_dir: Path,
    output_dir: Path,
    foreground_method: str,
    background_value: float,
    min_size: int,
    segmentation_method: str,
    random_walker_beta: float,
    propagation_regularization: float,
    rows: list[dict[str, Any]],
    image_records: list[dict[str, Any]],
) -> None:
    method_config = {
        "name": segmentation_method,
        "whole_cell_claim": False,
        "reported_object": "DAPI-seeded aSMA-associated region",
        "fixed_radius_expansion": False,
        "scientific_warning": (
            "CH2/aSMA is the endpoint signal and is used here to define the region. "
            "The mask is therefore exploratory and must not be interpreted as an "
            "unbiased whole-cell mask without manual validation."
        ),
    }
    if segmentation_method == "seeded_intensity_random_walker":
        method_config["random_walker_beta"] = float(random_walker_beta)
    if segmentation_method == "seeded_intensity_propagation":
        method_config["propagation_regularization"] = float(propagation_regularization)
        method_config["propagation_backend"] = "centrosome.propagate"
    config = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "output_dir": str(output_dir),
        "method": method_config,
        "channel_extraction": {
            "ch2_channel": "CH2/aSMA",
            "ch4_channel": "CH4/DAPI",
            "axes_policy": (
                "singleton axes are squeezed; 2-D images are used directly; RGB/YXS exports "
                "with exactly one active sample use that active sample"
            ),
            "z_projection": "none",
        },
        "foreground": {
            "method": foreground_method,
            "min_size": min_size,
            "fill_holes": True,
        },
        "background": {
            "background_value_per_px": float(background_value),
            "correction": "sum(max(CH2 - background_value_per_px, 0)) inside seeded region",
        },
        "validation_status": {
            "visual_qc_generated": True,
            "manual_ground_truth_available": False,
            "precision_recall_f1_allowed": False,
            "whole_cell_segmentation_validated": False,
        },
        "outputs": {
            "summary_csv": str(output_dir / "summaries" / "seeded_region_image_metrics.csv"),
            "masks_dir": str(output_dir / "masks"),
            "qc_dir": str(output_dir / "qc"),
            "plots_dir": str(output_dir / "plots"),
            "qc_contact_sheet": str(output_dir / "qc_contact_sheet.png"),
        },
        "images_processed": len(rows),
        "image_records": image_records,
    }
    logs_dir.mkdir(parents=True, exist_ok=True)
    (logs_dir / "config_resolved.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
    )
    run_log_lines = [
        f"generated_at_utc: {config['generated_at_utc']}",
        f"method: {segmentation_method}",
        "reported_object: DAPI-seeded aSMA-associated region",
        "whole_cell_claim: false",
        f"foreground_method: {foreground_method}",
    ]
    if segmentation_method == "seeded_intensity_random_walker":
        run_log_lines.append(f"random_walker_beta: {random_walker_beta}")
    if segmentation_method == "seeded_intensity_propagation":
        run_log_lines.append(f"propagation_regularization: {propagation_regularization}")
        run_log_lines.append("propagation_backend: centrosome.propagate")
        run_log_lines.append("reported_labels_clipped_to_ch2_foreground: true")
    run_log_lines.extend(
        [
            f"background_value_per_px: {background_value}",
            f"images_processed: {len(rows)}",
            (
                "warning: CH2/aSMA endpoint signal is used to define regions; "
                "outputs are exploratory, not validated whole-cell masks."
            ),
        ]
    )
    (logs_dir / "run_log.txt").write_text(
        "\n".join(run_log_lines) + "\n",
        encoding="utf-8",
    )
