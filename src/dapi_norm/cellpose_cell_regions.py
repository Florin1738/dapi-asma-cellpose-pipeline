from __future__ import annotations

import csv
from collections.abc import Callable, Iterable
from datetime import UTC, datetime
import importlib.metadata as importlib_metadata
import os
from pathlib import Path
import platform
import re
import sys
import textwrap
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw
from skimage import segmentation
import tifffile
import yaml

from dapi_norm.cellpose_runner import validate_label_mask
from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.pi_simple_summary import ImagePair, _normalize_channel_id

CellRegionSegmenter = Callable[[np.ndarray], np.ndarray]

CELLPOSE_CELL_REGION_METHOD = "cellpose_ch2_ch4_candidate_asma_associated_region"
CELLPOSE_CELL_REGION_COLUMNS = [
    "image_id",
    "source_id",
    "method",
    "model_name",
    "target_channel_id",
    "dapi_channel_id",
    "target_path",
    "dapi_path",
    "background_value_per_px",
    "background_method",
    "image_area_px",
    "cellpose_object_count",
    "candidate_region_area_px",
    "candidate_region_fraction",
    "outside_candidate_region_area_px",
    "outside_candidate_region_fraction",
    "nuclei_mask_path",
    "dapi_positive_nucleus_count",
    "normalization_denominator_count",
    "nuclei_filtering_applied",
    "nuclei_filtering_policy",
    "dapi_nuclei_with_centroid_inside_cellpose_region",
    "dapi_nuclei_centroid_coverage_fraction",
    "cellpose_objects_with_dapi_centroid",
    "cellpose_objects_without_dapi_centroid",
    "cellpose_objects_without_dapi_centroid_fraction",
    "cellpose_objects_with_multiple_dapi_centroids",
    "cellpose_objects_with_multiple_dapi_centroids_fraction",
    "dapi_anchored_cellpose_object_count",
    "dapi_anchored_excluded_no_dapi_object_count",
    "dapi_anchored_candidate_region_area_px",
    "dapi_anchored_candidate_region_fraction",
    "dapi_anchored_positive_area_per_DAPI_positive_nucleus",
    "dapi_anchored_target_integrated_raw",
    "dapi_anchored_target_integrated_background_corrected",
    "dapi_anchored_target_integrated_intensity_per_DAPI_positive_nucleus",
    "target_integrated_raw_in_cellpose_region",
    "target_mean_raw_in_cellpose_region",
    "target_integrated_background_corrected_in_cellpose_region",
    "target_integrated_intensity_per_DAPI_positive_nucleus",
    "target_integrated_intensity_per_cellpose_object",
    "total_image_integrated_background_corrected",
    "excluded_region_integrated_background_corrected",
    "excluded_region_background_corrected_fraction",
    "excluded_signal_display_threshold",
    "qc_status",
    "qc_flags",
    "mask_path",
    "qc_panel_path",
    "excluded_signal_check_path",
    "warnings",
]

LOW_DAPI_CENTROID_COVERAGE_THRESHOLD = 0.5
NEAR_FULL_FIELD_CELL_REGION_FRACTION_THRESHOLD = 0.9
LOW_CELL_REGION_FRACTION_WITH_NUCLEI_THRESHOLD = 0.01
OBJECTS_WITHOUT_DAPI_CENTROID_WARNING_THRESHOLD = 0.10
OBJECTS_WITH_MULTIPLE_DAPI_CENTROIDS_WARNING_THRESHOLD = 0.05
EXCLUDED_BACKGROUND_CORRECTED_FRACTION_WARNING_THRESHOLD = 0.20
EXCLUDED_BACKGROUND_CORRECTED_FRACTION_REJECT_THRESHOLD = 0.50


def build_cellpose_cell_input(ch2_image: np.ndarray, ch4_image: np.ndarray) -> np.ndarray:
    """Stack target and DAPI-context channels as a two-channel image for Cellpose."""
    ch2 = np.asarray(ch2_image)
    ch4 = np.asarray(ch4_image)
    if ch2.ndim != 2 or ch4.ndim != 2:
        raise ValueError(f"Expected 2-D target and DAPI images, got {ch2.shape} and {ch4.shape}")
    if ch2.shape != ch4.shape:
        raise ValueError(
            f"Target and DAPI images must have matching shapes, got {ch2.shape} and {ch4.shape}"
        )
    return np.stack([ch2.astype(np.float32), ch4.astype(np.float32)], axis=0)


def measure_cellpose_cell_region_image(
    *,
    image_id: str,
    ch2_image: np.ndarray,
    nuclei_mask: np.ndarray,
    cell_labels: np.ndarray,
    background_value: float = 0.0,
    source_id: str | None = None,
    model_name: str = "",
    target_channel_id: str = "CH2",
    dapi_channel_id: str = "CH4",
    target_channel_name: str = "aSMA",
    dapi_channel_name: str = "DAPI",
    method_name: str | None = None,
    nuclei_mask_path: Path | None = None,
    mask_path: Path | None = None,
    qc_panel_path: Path | None = None,
    excluded_signal_check_path: Path | None = None,
) -> dict[str, Any]:
    target_channel = _normalize_channel_id(target_channel_id)
    dapi_channel = _normalize_channel_id(dapi_channel_id)
    ch2 = np.asarray(ch2_image, dtype=np.float64)
    nuclei = np.asarray(nuclei_mask)
    labels = np.asarray(cell_labels)
    if ch2.ndim != 2 or nuclei.ndim != 2 or labels.ndim != 2:
        raise ValueError("CH2 image, nuclei mask, and Cellpose labels must all be 2-D")
    if ch2.shape != nuclei.shape or ch2.shape != labels.shape:
        raise ValueError(
            "CH2 image, nuclei mask, and Cellpose labels must have matching shapes: "
            f"ch2={ch2.shape}, nuclei={nuclei.shape}, labels={labels.shape}"
        )
    if np.any(nuclei < 0) or np.any(labels < 0):
        raise ValueError("Nuclei and Cellpose label masks must not contain negative labels")

    region_mask = labels > 0
    image_area = int(ch2.size)
    region_area = int(np.count_nonzero(region_mask))
    raw_integrated = float(np.sum(ch2[region_mask]))
    mean_raw = float(np.mean(ch2[region_mask])) if region_area else float("nan")
    corrected = np.clip(ch2[region_mask] - float(background_value), 0, None)
    corrected_integrated = float(np.sum(corrected))
    corrected_full = np.clip(ch2 - float(background_value), 0, None)
    total_corrected_integrated = float(np.sum(corrected_full))
    excluded_corrected_integrated = float(np.sum(corrected_full[~region_mask]))
    excluded_corrected_fraction = (
        excluded_corrected_integrated / total_corrected_integrated
        if total_corrected_integrated > 0
        else 0.0
    )
    nucleus_count = _count_nonzero_labels(nuclei)
    object_count = _count_nonzero_labels(labels)
    covered_nuclei = _count_nuclei_with_centroid_inside_region(nuclei, region_mask)
    object_occupancy = _cellpose_object_dapi_centroid_occupancy(nuclei, labels)
    anchored_object_ids = _cellpose_object_ids_with_dapi_centroid(nuclei, labels)
    anchored_region_mask = _mask_for_object_ids(labels, anchored_object_ids)
    anchored_region_area = int(np.count_nonzero(anchored_region_mask))
    anchored_raw_integrated = float(np.sum(ch2[anchored_region_mask]))
    anchored_corrected = np.clip(ch2[anchored_region_mask] - float(background_value), 0, None)
    anchored_corrected_integrated = float(np.sum(anchored_corrected))

    warnings = [
        "exploratory_cellpose_output_not_validated_whole_cell_mask",
        "do_not_interpret_as_true_cell_segmentation_without_manual_validation",
    ]
    if target_channel == "CH2":
        warnings.append("ch2_asma_used_as_candidate_cytoplasm_channel")
    else:
        warnings.append(
            f"{target_channel.lower()}_{target_channel_name.lower()}_used_as_candidate_cytoplasm_channel"
        )
    if nucleus_count == 0:
        warnings.append("zero_dapi_positive_nucleus_count")
    if object_count == 0:
        warnings.append("zero_cellpose_object_count")
    if covered_nuclei < nucleus_count:
        warnings.append("some_dapi_nuclei_centroids_not_inside_cellpose_region")
    if object_occupancy["objects_without_dapi_centroid"] > 0:
        warnings.append("some_cellpose_objects_have_no_dapi_centroid")
    if object_occupancy["objects_with_multiple_dapi_centroids"] > 0:
        warnings.append("some_cellpose_objects_have_multiple_dapi_centroids")

    row = {
        "image_id": image_id,
        "source_id": source_id or image_id,
        "method": method_name or _cell_region_method_name(target_channel, dapi_channel),
        "model_name": model_name,
        "target_channel_id": target_channel,
        "dapi_channel_id": dapi_channel,
        "target_path": "",
        "dapi_path": "",
        "background_value_per_px": float(background_value),
        "background_method": _background_method(background_value),
        "image_area_px": image_area,
        "cellpose_object_count": object_count,
        "candidate_region_area_px": region_area,
        "candidate_region_fraction": float(region_area / image_area) if image_area else 0.0,
        "outside_candidate_region_area_px": int(image_area - region_area),
        "outside_candidate_region_fraction": float((image_area - region_area) / image_area)
        if image_area
        else 0.0,
        "nuclei_mask_path": "" if nuclei_mask_path is None else str(nuclei_mask_path),
        "dapi_positive_nucleus_count": nucleus_count,
        "normalization_denominator_count": nucleus_count,
        "nuclei_filtering_applied": False,
        "nuclei_filtering_policy": "none_count_nonzero_labels_in_supplied_mask",
        "dapi_nuclei_with_centroid_inside_cellpose_region": covered_nuclei,
        "dapi_nuclei_centroid_coverage_fraction": float(covered_nuclei / nucleus_count)
        if nucleus_count
        else 0.0,
        "cellpose_objects_with_dapi_centroid": object_occupancy["objects_with_dapi_centroid"],
        "cellpose_objects_without_dapi_centroid": object_occupancy[
            "objects_without_dapi_centroid"
        ],
        "cellpose_objects_without_dapi_centroid_fraction": (
            object_occupancy["objects_without_dapi_centroid"] / object_count
            if object_count
            else 0.0
        ),
        "cellpose_objects_with_multiple_dapi_centroids": object_occupancy[
            "objects_with_multiple_dapi_centroids"
        ],
        "cellpose_objects_with_multiple_dapi_centroids_fraction": (
            object_occupancy["objects_with_multiple_dapi_centroids"] / object_count
            if object_count
            else 0.0
        ),
        "dapi_anchored_cellpose_object_count": len(anchored_object_ids),
        "dapi_anchored_excluded_no_dapi_object_count": int(object_count - len(anchored_object_ids)),
        "dapi_anchored_candidate_region_area_px": anchored_region_area,
        "dapi_anchored_candidate_region_fraction": (
            float(anchored_region_area / image_area) if image_area else 0.0
        ),
        "dapi_anchored_positive_area_per_DAPI_positive_nucleus": (
            anchored_region_area / nucleus_count if nucleus_count else float("nan")
        ),
        "dapi_anchored_target_integrated_raw": anchored_raw_integrated,
        "dapi_anchored_target_integrated_background_corrected": anchored_corrected_integrated,
        "dapi_anchored_target_integrated_intensity_per_DAPI_positive_nucleus": (
            anchored_corrected_integrated / nucleus_count if nucleus_count else float("nan")
        ),
        "target_integrated_raw_in_cellpose_region": raw_integrated,
        "target_mean_raw_in_cellpose_region": mean_raw,
        "target_integrated_background_corrected_in_cellpose_region": corrected_integrated,
        "target_integrated_intensity_per_DAPI_positive_nucleus": (
            corrected_integrated / nucleus_count if nucleus_count else float("nan")
        ),
        "target_integrated_intensity_per_cellpose_object": (
            corrected_integrated / object_count if object_count else float("nan")
        ),
        "total_image_integrated_background_corrected": total_corrected_integrated,
        "excluded_region_integrated_background_corrected": excluded_corrected_integrated,
        "excluded_region_background_corrected_fraction": excluded_corrected_fraction,
        "excluded_signal_display_threshold": _excluded_signal_display_threshold(ch2),
        "mask_path": "" if mask_path is None else str(mask_path),
        "qc_panel_path": "" if qc_panel_path is None else str(qc_panel_path),
        "excluded_signal_check_path": ""
        if excluded_signal_check_path is None
        else str(excluded_signal_check_path),
        "warnings": ";".join(warnings),
    }
    qc_status, qc_flags = score_cellpose_cell_region_qc(row)
    row["qc_status"] = qc_status
    row["qc_flags"] = ";".join(qc_flags)
    return row


def score_cellpose_cell_region_qc(row: dict[str, Any]) -> tuple[str, list[str]]:
    flags: list[str] = ["not_validated_whole_cell_mask"]
    object_count = int(row.get("cellpose_object_count", row.get("cellpose_cell_count", 0)))
    region_fraction = float(
        row.get("candidate_region_fraction", row.get("cell_region_fraction", 0.0))
    )
    nucleus_count = int(row.get("dapi_positive_nucleus_count", 0))
    coverage_fraction = float(row.get("dapi_nuclei_centroid_coverage_fraction", 0.0))
    objects_without_dapi_fraction = float(
        row.get("cellpose_objects_without_dapi_centroid_fraction", 0.0)
    )
    objects_with_multiple_dapi_fraction = float(
        row.get("cellpose_objects_with_multiple_dapi_centroids_fraction", 0.0)
    )
    excluded_corrected_fraction = float(
        row.get("excluded_region_background_corrected_fraction", 0.0)
    )

    if object_count == 0:
        flags.append("zero_cellpose_object_count")
    if nucleus_count == 0:
        flags.append("zero_dapi_positive_nucleus_count")
    if region_fraction >= NEAR_FULL_FIELD_CELL_REGION_FRACTION_THRESHOLD:
        flags.append("near_full_field_cellpose_region_fraction")
    if nucleus_count > 0 and region_fraction <= LOW_CELL_REGION_FRACTION_WITH_NUCLEI_THRESHOLD:
        flags.append("very_low_cellpose_region_fraction_with_dapi_nuclei")
    if nucleus_count > 0 and coverage_fraction < LOW_DAPI_CENTROID_COVERAGE_THRESHOLD:
        flags.append("low_dapi_nuclei_centroid_coverage")
    if objects_without_dapi_fraction > OBJECTS_WITHOUT_DAPI_CENTROID_WARNING_THRESHOLD:
        flags.append("candidate_objects_without_dapi_centroid")
    if objects_with_multiple_dapi_fraction > OBJECTS_WITH_MULTIPLE_DAPI_CENTROIDS_WARNING_THRESHOLD:
        flags.append("candidate_objects_with_multiple_dapi_centroids")
    if excluded_corrected_fraction > EXCLUDED_BACKGROUND_CORRECTED_FRACTION_WARNING_THRESHOLD:
        flags.append("sizeable_background_corrected_ch2_outside_cellpose_region")
    if excluded_corrected_fraction > EXCLUDED_BACKGROUND_CORRECTED_FRACTION_REJECT_THRESHOLD:
        flags.append("majority_background_corrected_ch2_outside_cellpose_region")

    reject_flags = {
        "zero_cellpose_object_count",
        "zero_dapi_positive_nucleus_count",
        "near_full_field_cellpose_region_fraction",
        "very_low_cellpose_region_fraction_with_dapi_nuclei",
        "low_dapi_nuclei_centroid_coverage",
        "majority_background_corrected_ch2_outside_cellpose_region",
    }
    if any(flag in flags for flag in reject_flags):
        return "reject_qc_failure", flags
    if len(flags) > 1:
        return "needs_manual_review", flags
    return "reviewable_not_validated", flags


def run_cellpose_cell_region_batch(
    *,
    image_pairs: Iterable[ImagePair],
    mask_lookup: dict[str, Path],
    output_dir: Path | str,
    model_name: str = "cpsam_v2",
    gpu: bool = True,
    segmenter: CellRegionSegmenter | None = None,
    background_value: float = 0.0,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    diameter: float | None = None,
    target_channel_id: str = "CH2",
    dapi_channel_id: str = "CH4",
    target_channel_name: str = "aSMA",
    dapi_channel_name: str = "DAPI",
    write_internal_qc: bool = True,
) -> list[dict[str, Any]]:
    target_channel = _normalize_channel_id(target_channel_id)
    dapi_channel = _normalize_channel_id(dapi_channel_id)
    if target_channel == dapi_channel:
        raise ValueError("Target and DAPI channels must be different.")
    method_name = _cell_region_method_name(target_channel, dapi_channel)
    output_path = Path(output_dir)
    masks_dir = output_path / "masks"
    qc_dir = output_path / "qc"
    summaries_dir = output_path / "summaries"
    logs_dir = output_path / "logs"
    masks_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    summaries_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    active_segmenter = segmenter
    if active_segmenter is None:
        active_segmenter = _cellpose_cell_segmenter(
            model_name=model_name,
            gpu=gpu,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
            diameter=diameter,
        )

    rows: list[dict[str, Any]] = []
    qc_paths: list[Path] = []
    image_records: list[dict[str, Any]] = []
    for pair in image_pairs:
        pair_target_channel = _normalize_channel_id(pair.target_channel_id)
        pair_dapi_channel = _normalize_channel_id(pair.dapi_channel_id)
        if pair_target_channel != target_channel or pair_dapi_channel != dapi_channel:
            raise ValueError(
                f"{pair.source_id} was discovered as {pair_target_channel}/{pair_dapi_channel}, "
                f"but this run is configured for {target_channel}/{dapi_channel}."
            )
        nuclei_mask_path = _mask_path_for_pair(pair, mask_lookup)
        ch2_image, ch2_pages = read_primary_intensity_plane(pair.ch2_path)
        ch4_image, ch4_pages = read_primary_intensity_plane(pair.ch4_path)
        nuclei_mask = np.asarray(tifffile.imread(nuclei_mask_path))
        if ch2_image.shape != ch4_image.shape:
            raise ValueError(
                f"{pair.source_id} {target_channel} shape {ch2_image.shape} does not match "
                f"{dapi_channel} shape {ch4_image.shape}"
            )
        if nuclei_mask.shape != ch2_image.shape:
            raise ValueError(
                f"{pair.source_id} nuclei mask shape {nuclei_mask.shape} does not match "
                f"{target_channel} shape {ch2_image.shape}"
            )

        cellpose_input = build_cellpose_cell_input(ch2_image, ch4_image)
        cell_labels = validate_label_mask(
            active_segmenter(cellpose_input),
            image_shape=ch2_image.shape,
        )
        safe_id = _safe_output_stem(pair.source_id)
        safe_model = _safe_output_stem(model_name)
        channel_token = f"{target_channel.lower()}_{dapi_channel.lower()}"
        cell_mask_path = masks_dir / f"{safe_id}_cellpose_{channel_token}_{safe_model}_labels.tif"
        qc_path = qc_dir / f"{safe_id}_cellpose_cell_region_qc.png" if write_internal_qc else None
        excluded_signal_check_path = (
            qc_path.with_name(qc_path.stem + "_excluded_signal_check.png")
            if qc_path is not None
            else None
        )
        tifffile.imwrite(cell_mask_path, cell_labels, photometric="minisblack")

        row = measure_cellpose_cell_region_image(
            image_id=pair.location,
            source_id=pair.source_id,
            ch2_image=ch2_image,
            nuclei_mask=nuclei_mask,
            cell_labels=cell_labels,
            background_value=background_value,
            model_name=model_name,
            target_channel_id=target_channel,
            dapi_channel_id=dapi_channel,
            target_channel_name=target_channel_name,
            dapi_channel_name=dapi_channel_name,
            method_name=method_name,
            nuclei_mask_path=nuclei_mask_path,
            mask_path=cell_mask_path,
            qc_panel_path=qc_path,
            excluded_signal_check_path=excluded_signal_check_path,
        )
        row["target_path"] = str(pair.ch2_path)
        row["dapi_path"] = str(pair.ch4_path)
        if write_internal_qc and qc_path is not None:
            write_cellpose_cell_region_qc_panel(
                image_id=pair.source_id,
                ch2_image=ch2_image,
                ch4_image=ch4_image,
                nuclei_mask=nuclei_mask,
                cell_labels=cell_labels,
                output_path=qc_path,
                metrics=row,
                target_channel_id=target_channel,
                dapi_channel_id=dapi_channel,
                target_channel_name=target_channel_name,
                dapi_channel_name=dapi_channel_name,
            )
        rows.append(row)
        if qc_path is not None:
            qc_paths.append(qc_path)
        image_records.append(
            {
                "source_id": pair.source_id,
                "location": pair.location,
                "target_channel_id": target_channel,
                "dapi_channel_id": dapi_channel,
                "target_path": str(pair.ch2_path),
                "dapi_path": str(pair.ch4_path),
                "ch2_path": str(pair.ch2_path),
                "ch4_path": str(pair.ch4_path),
                "nuclei_mask_path": str(nuclei_mask_path),
                "cellpose_mask_path": str(cell_mask_path),
                "qc_panel_path": str(qc_path),
                "excluded_signal_check_path": str(excluded_signal_check_path),
                "ch2_page_count": int(ch2_pages),
                "ch4_page_count": int(ch4_pages),
                "target_page_count": int(ch2_pages),
                "dapi_page_count": int(ch4_pages),
                "image_shape": [int(ch2_image.shape[0]), int(ch2_image.shape[1])],
            }
        )

    _write_summary_csv(summaries_dir / "cellpose_cell_region_image_metrics.csv", rows)
    _write_summary_plots(rows, output_path / "plots")
    if write_internal_qc:
        write_cellpose_cell_region_contact_sheet(qc_paths, output_path / "qc_contact_sheet.png")
    _write_cellpose_cell_region_metadata(
        logs_dir=logs_dir,
        output_dir=output_path,
        model_name=model_name,
        gpu=gpu,
        background_value=background_value,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        diameter=diameter,
        target_channel_id=target_channel,
        dapi_channel_id=dapi_channel,
        target_channel_name=target_channel_name,
        dapi_channel_name=dapi_channel_name,
        method_name=method_name,
        rows=rows,
        image_records=image_records,
        internal_qc_generated=write_internal_qc,
    )
    return rows


def write_cellpose_cell_region_qc_panel(
    *,
    image_id: str,
    ch2_image: np.ndarray,
    ch4_image: np.ndarray,
    nuclei_mask: np.ndarray,
    cell_labels: np.ndarray,
    output_path: Path,
    metrics: dict[str, Any],
    target_channel_id: str = "CH2",
    dapi_channel_id: str = "CH4",
    target_channel_name: str = "aSMA",
    dapi_channel_name: str = "DAPI",
) -> None:
    target_channel = _normalize_channel_id(target_channel_id)
    dapi_channel = _normalize_channel_id(dapi_channel_id)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ch2_scaled = _scale_for_display(ch2_image)
    ch4_scaled = _scale_for_display(ch4_image)
    labels = np.asarray(cell_labels)
    region_mask = labels > 0
    display_threshold = float(
        metrics.get("excluded_signal_display_threshold", _excluded_signal_display_threshold(ch2_image))
    )
    excluded_signal = (~region_mask) & (np.asarray(ch2_image, dtype=np.float64) > display_threshold)
    cell_boundaries = segmentation.find_boundaries(labels, mode="outer")
    nucleus_boundaries = segmentation.find_boundaries(nuclei_mask, mode="outer")

    dapi_rgb = np.dstack([ch4_scaled, ch4_scaled, ch4_scaled])
    dapi_rgb[nucleus_boundaries] = [0.0, 0.9, 0.1]
    _draw_centroid_crosses(
        dapi_rgb,
        nuclei_mask,
        color=np.array([0.0, 1.0, 0.0]),
        halo_color=np.array([0.0, 0.0, 0.0]),
    )

    retained_overlay = _ch2_rgb(ch2_scaled)
    _blend_color(retained_overlay, region_mask, np.array([0.55, 1.0, 0.55]), alpha=0.45)
    retained_overlay[cell_boundaries] = [0.0, 1.0, 0.0]

    combined = _ch2_rgb(ch2_scaled)
    _blend_color(combined, region_mask, np.array([0.55, 1.0, 0.55]), alpha=0.38)
    combined[cell_boundaries] = [0.0, 1.0, 0.0]
    combined[nucleus_boundaries] = [0.0, 0.75, 1.0]
    _draw_centroid_crosses(
        combined,
        nuclei_mask,
        color=np.array([0.0, 1.0, 0.0]),
        halo_color=np.array([0.0, 0.0, 0.0]),
    )

    excluded_rgb = _ch2_rgb(ch2_scaled)
    _blend_color(excluded_rgb, region_mask, np.array([0.55, 1.0, 0.55]), alpha=0.25)
    _blend_color(excluded_rgb, excluded_signal, np.array([1.0, 0.1, 0.7]), alpha=0.6)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10.8), constrained_layout=True)
    flat_axes = axes.ravel()
    flat_axes[0].imshow(dapi_rgb)
    flat_axes[0].set_title(
        f"{image_id} {dapi_channel}/{dapi_channel_name}\n"
        f"green X count={metrics.get('dapi_positive_nucleus_count', 'NA')}"
    )
    flat_axes[1].imshow(ch2_scaled, cmap="gray")
    flat_axes[1].set_title(f"{target_channel}/{target_channel_name} raw display")
    flat_axes[2].imshow(retained_overlay)
    flat_axes[2].set_title(
        f"Cellpose candidate region over {target_channel}\n"
        "light green fill, green edge"
    )
    flat_axes[3].imshow(combined)
    flat_axes[3].set_title(
        "Combined QC\n"
        "green X nuclei, cyan nuclear edges"
    )
    for ax in flat_axes:
        ax.axis("off")
    qc_flags = textwrap.fill(
        str(metrics.get("qc_flags", "NA")).replace(";", "; "),
        width=120,
    )
    fig.suptitle(
        f"Cellpose {target_channel}+{dapi_channel} candidate {target_channel_name}-associated object QC; "
        "not validated whole-cell segmentation\n"
        f"raw in region={_format_scientific(metrics.get('target_integrated_raw_in_cellpose_region'))}, "
        "region/nucleus="
        f"{_format_scientific(metrics.get('target_integrated_intensity_per_DAPI_positive_nucleus'))}, "
        f"area={metrics.get('candidate_region_fraction', float('nan')):.1%}, "
        f"DAPI coverage={metrics.get('dapi_nuclei_centroid_coverage_fraction', float('nan')):.1%}\n"
        f"QC={metrics.get('qc_status', 'NA')}; "
        "excluded corrected CH2="
        f"{metrics.get('excluded_region_background_corrected_fraction', float('nan')):.1%}; "
        f"background={metrics.get('background_method', 'NA')}\n"
        f"flags={qc_flags}",
        fontsize=9.5,
    )
    fig.savefig(output_path, dpi=160)
    plt.close(fig)

    excluded_panel_path = output_path.with_name(output_path.stem + "_excluded_signal_check.png")
    _write_excluded_signal_check(
        output_path=excluded_panel_path,
        image_id=image_id,
        excluded_rgb=excluded_rgb,
        combined=combined,
        display_threshold=display_threshold,
        target_channel_id=target_channel,
    )


def write_cellpose_cell_region_contact_sheet(qc_paths: list[Path], output_path: Path) -> None:
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
        thumbs.append((path.stem.replace("_cellpose_cell_region_qc", ""), thumb))

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


def _cellpose_cell_segmenter(
    *,
    model_name: str,
    gpu: bool,
    flow_threshold: float,
    cellprob_threshold: float,
    diameter: float | None,
) -> CellRegionSegmenter:
    from cellpose import models

    model = models.CellposeModel(
        gpu=gpu,
        pretrained_model=model_name,
        nchan=2,
        use_bfloat16=False,
    )

    def segment(image: np.ndarray) -> np.ndarray:
        return model.eval(
            image,
            channel_axis=0,
            diameter=diameter,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
        )[0]

    return segment


def _mask_path_for_pair(pair: ImagePair, mask_lookup: dict[str, Path]) -> Path:
    for key in [pair.source_id, pair.location]:
        normalized = _lookup_key(key)
        if normalized in mask_lookup:
            return Path(mask_lookup[normalized])
    raise KeyError(
        f"No DAPI nuclei mask found for {pair.source_id} / {pair.location}. "
        "Run Cellpose nuclei counting first or provide a mask lookup with source_id/location keys."
    )


def _lookup_key(value: str) -> str:
    return value.strip().replace("\\", "/").upper()


def _count_nonzero_labels(mask: np.ndarray) -> int:
    labels = np.unique(mask)
    return int(np.count_nonzero(labels))


def _count_nuclei_with_centroid_inside_region(nuclei_mask: np.ndarray, region_mask: np.ndarray) -> int:
    nuclei = np.asarray(nuclei_mask)
    region = np.asarray(region_mask, dtype=bool)
    count = 0
    for label in np.unique(nuclei):
        if label == 0:
            continue
        ys, xs = np.where(nuclei == label)
        if len(ys) == 0:
            continue
        y = int(round(float(np.mean(ys))))
        x = int(round(float(np.mean(xs))))
        y = min(max(y, 0), region.shape[0] - 1)
        x = min(max(x, 0), region.shape[1] - 1)
        if region[y, x]:
            count += 1
    return count


def _cellpose_object_dapi_centroid_occupancy(
    nuclei_mask: np.ndarray,
    cellpose_labels: np.ndarray,
) -> dict[str, int]:
    nuclei = np.asarray(nuclei_mask)
    objects = np.asarray(cellpose_labels)
    object_ids = np.unique(objects)
    object_ids = object_ids[object_ids > 0]
    centroid_counts_by_object = {int(label): 0 for label in object_ids}
    for nucleus_label in np.unique(nuclei):
        if nucleus_label == 0:
            continue
        ys, xs = np.where(nuclei == nucleus_label)
        if len(ys) == 0:
            continue
        y = min(max(int(round(float(np.mean(ys)))), 0), objects.shape[0] - 1)
        x = min(max(int(round(float(np.mean(xs)))), 0), objects.shape[1] - 1)
        object_label = int(objects[y, x])
        if object_label > 0 and object_label in centroid_counts_by_object:
            centroid_counts_by_object[object_label] += 1

    objects_with_dapi = sum(count > 0 for count in centroid_counts_by_object.values())
    objects_with_multiple = sum(count > 1 for count in centroid_counts_by_object.values())
    return {
        "objects_with_dapi_centroid": int(objects_with_dapi),
        "objects_without_dapi_centroid": int(len(object_ids) - objects_with_dapi),
        "objects_with_multiple_dapi_centroids": int(objects_with_multiple),
    }


def _cellpose_object_ids_with_dapi_centroid(
    nuclei_mask: np.ndarray,
    cellpose_labels: np.ndarray,
) -> set[int]:
    nuclei = np.asarray(nuclei_mask)
    objects = np.asarray(cellpose_labels)
    object_ids: set[int] = set()
    for nucleus_label in np.unique(nuclei):
        if nucleus_label == 0:
            continue
        ys, xs = np.where(nuclei == nucleus_label)
        if len(ys) == 0:
            continue
        y = min(max(int(round(float(np.mean(ys)))), 0), objects.shape[0] - 1)
        x = min(max(int(round(float(np.mean(xs)))), 0), objects.shape[1] - 1)
        object_label = int(objects[y, x])
        if object_label > 0:
            object_ids.add(object_label)
    return object_ids


def _mask_for_object_ids(labels: np.ndarray, object_ids: set[int]) -> np.ndarray:
    if not object_ids:
        return np.zeros_like(labels, dtype=bool)
    return np.isin(labels, list(object_ids))


def _background_method(background_value: float) -> str:
    value = float(background_value)
    if value == 0.0:
        return "none"
    return f"constant_value_{value:g}"


def _excluded_signal_display_threshold(image: np.ndarray) -> float:
    return float(np.percentile(np.asarray(image, dtype=np.float64), 75))


def _draw_centroid_crosses(
    rgb: np.ndarray,
    mask: np.ndarray,
    *,
    color: np.ndarray,
    halo_color: np.ndarray | None = None,
) -> None:
    labels = np.unique(mask)
    labels = labels[labels > 0]
    size = max(3, min(8, round(min(mask.shape) / 140)))
    for label in labels:
        ys, xs = np.where(mask == label)
        if len(ys) == 0:
            continue
        y = int(round(float(np.mean(ys))))
        x = int(round(float(np.mean(xs))))
        if halo_color is not None:
            _draw_centroid_cross(rgb, y=y, x=x, size=size + 1, color=halo_color)
        _draw_centroid_cross(rgb, y=y, x=x, size=size, color=color)


def _draw_centroid_cross(rgb: np.ndarray, *, y: int, x: int, size: int, color: np.ndarray) -> None:
    for offset in range(-size, size + 1):
        y1 = y + offset
        x1 = x + offset
        y2 = y + offset
        x2 = x - offset
        if 0 <= y1 < rgb.shape[0] and 0 <= x1 < rgb.shape[1]:
            rgb[y1, x1] = color
        if 0 <= y2 < rgb.shape[0] and 0 <= x2 < rgb.shape[1]:
            rgb[y2, x2] = color


def _scale_for_display(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image, dtype=np.float32)
    if arr.size == 0:
        return arr
    low, high = np.percentile(arr, [1, 99.8])
    if high <= low:
        high = float(np.max(arr))
        low = float(np.min(arr))
    if high <= low:
        return np.zeros_like(arr, dtype=np.float32)
    return np.clip((arr - low) / (high - low), 0, 1)


def _ch2_rgb(ch2_scaled: np.ndarray) -> np.ndarray:
    return np.dstack([ch2_scaled, ch2_scaled, ch2_scaled])


def _blend_color(rgb: np.ndarray, mask: np.ndarray, color: np.ndarray, *, alpha: float) -> None:
    if not np.any(mask):
        return
    rgb[mask] = (1.0 - alpha) * rgb[mask] + alpha * color


def _write_excluded_signal_check(
    *,
    output_path: Path,
    image_id: str,
    excluded_rgb: np.ndarray,
    combined: np.ndarray,
    display_threshold: float,
    target_channel_id: str = "CH2",
) -> None:
    target_channel = _normalize_channel_id(target_channel_id)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    axes[0].imshow(excluded_rgb)
    axes[0].set_title(
        f"High {target_channel} outside candidate region\n"
        f"magenta = excluded display pixels > p75 ({display_threshold:.0f})"
    )
    axes[1].imshow(combined)
    axes[1].set_title("Candidate region + DAPI centroid context")
    for ax in axes:
        ax.axis("off")
    fig.suptitle(f"{image_id} exclusion check", fontsize=12)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def _write_summary_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CELLPOSE_CELL_REGION_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in CELLPOSE_CELL_REGION_COLUMNS})


def _write_summary_plots(rows: list[dict[str, Any]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    labels = [str(row["source_id"]) for row in rows]
    target_channel = str(rows[0].get("target_channel_id", "CH2"))
    dapi_channel = str(rows[0].get("dapi_channel_id", "CH4"))
    per_nucleus = [
        float(row["target_integrated_intensity_per_DAPI_positive_nucleus"]) for row in rows
    ]
    region_fraction = [float(row["candidate_region_fraction"]) for row in rows]
    coverage = [float(row["dapi_nuclei_centroid_coverage_fraction"]) for row in rows]
    excluded_fraction = [float(row["excluded_region_background_corrected_fraction"]) for row in rows]

    fig, axes = plt.subplots(4, 1, figsize=(max(8, 0.45 * len(rows)), 11), constrained_layout=True)
    x = np.arange(len(rows))
    axes[0].bar(x, per_nucleus, color="#2f9e44")
    axes[0].set_ylabel(f"{target_channel} in Cellpose candidate region / DAPI nuclei")
    axes[1].bar(x, region_fraction, color="#1971c2")
    axes[1].set_ylabel("Candidate region fraction")
    axes[1].set_ylim(0, 1)
    axes[2].bar(x, coverage, color="#862e9c")
    axes[2].set_ylabel("DAPI centroid coverage")
    axes[2].set_ylim(0, 1)
    axes[3].bar(x, excluded_fraction, color="#c2255c")
    axes[3].set_ylabel(f"Excluded {target_channel} corrected fraction")
    axes[3].set_ylim(0, 1)
    axes[3].set_xticks(x)
    axes[3].set_xticklabels(labels, rotation=90)
    fig.suptitle(
        f"Cellpose {target_channel}+{dapi_channel} candidate region summary; exploratory only",
        fontsize=12,
    )
    fig.savefig(output_dir / "cellpose_cell_region_summary.png", dpi=160)
    plt.close(fig)


def _write_cellpose_cell_region_metadata(
    *,
    logs_dir: Path,
    output_dir: Path,
    model_name: str,
    gpu: bool,
    background_value: float,
    flow_threshold: float,
    cellprob_threshold: float,
    diameter: float | None,
    target_channel_id: str,
    dapi_channel_id: str,
    target_channel_name: str,
    dapi_channel_name: str,
    method_name: str,
    rows: list[dict[str, Any]],
    image_records: list[dict[str, Any]],
    internal_qc_generated: bool = True,
) -> None:
    target_channel = _normalize_channel_id(target_channel_id)
    dapi_channel = _normalize_channel_id(dapi_channel_id)
    logs_dir.mkdir(parents=True, exist_ok=True)
    config = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "method": {
            "name": method_name,
            "whole_cell_claim": False,
            "interpretation": (
                f"Candidate {target_channel_name}-associated object regions from pretrained Cellpose "
                f"using {target_channel} as target/signal and {dapi_channel} as DAPI context. Manual "
                "validation is required before treating masks as cell segmentation."
            ),
        },
        "channel_extraction": {
            "source": "primary TIFF series",
            "target_channel_id": target_channel,
            "target_channel_name": target_channel_name,
            "dapi_channel_id": dapi_channel,
            "dapi_channel_name": dapi_channel_name,
            "ch2_channel": f"{target_channel}/{target_channel_name} candidate cytoplasm"
            if target_channel == "CH2"
            else "legacy field; see target_channel_id",
            "ch4_channel": f"{dapi_channel}/{dapi_channel_name} nuclei"
            if dapi_channel == "CH4"
            else "legacy field; see dapi_channel_id",
            "z_projection": "none",
        },
        "model": {
            "name": model_name,
            "local_model_path": _resolved_model_path(model_name),
        },
        "segmentation_parameters": {
            "channel_axis": 0,
            "input_channel_order": [
                f"{target_channel}/{target_channel_name}",
                f"{dapi_channel}/{dapi_channel_name}",
            ],
            "diameter": diameter,
            "flow_threshold": flow_threshold,
            "cellprob_threshold": cellprob_threshold,
        },
        "measurement": {
            "background_value_per_px": float(background_value),
            "background_method": _background_method(background_value),
            "raw_measurement": (
                f"sum of {target_channel} pixel intensities inside Cellpose-labeled regions"
            ),
            "background_corrected_measurement": (
                f"sum of max({target_channel} - background_value_per_px, 0) inside "
                "Cellpose-labeled regions"
            ),
            "primary_endpoint": "target_integrated_intensity_per_DAPI_positive_nucleus",
            "normalization_denominator": (
                "count of nonzero labels in the supplied DAPI nuclei mask, not DAPI brightness"
            ),
            "dapi_anchored_region_definition": (
                "Cellpose objects retained only when at least one DAPI nucleus centroid falls "
                "inside the object"
            ),
            "nuclei_filtering_applied": False,
            "nuclei_filtering_policy": "none_count_nonzero_labels_in_supplied_mask",
        },
        "validation_status": {
            "whole_cell_segmentation_validated": False,
            "visual_qc_required": True,
            "manual_mask_validation_required_for_precision_recall_f1_iou": True,
        },
        "software_versions": _software_versions(),
        "device": _device_metadata(gpu),
        "outputs": {
            "image_metrics_csv": str(
                output_dir / "summaries" / "cellpose_cell_region_image_metrics.csv"
            ),
            "masks_dir": str(output_dir / "masks"),
            "qc_dir": str(output_dir / "qc"),
            "internal_qc_generated": bool(internal_qc_generated),
            "qc_contact_sheet": str(output_dir / "qc_contact_sheet.png")
            if internal_qc_generated
            else "",
            "summary_plot": str(output_dir / "plots" / "cellpose_cell_region_summary.png"),
        },
        "image_inputs": image_records,
    }
    (logs_dir / "config_resolved.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False),
        encoding="utf-8",
    )
    status_counts: dict[str, int] = {}
    for row in rows:
        status = str(row["qc_status"])
        status_counts[status] = status_counts.get(status, 0) + 1
    (logs_dir / "run_log.txt").write_text(
        "\n".join(
            [
                f"generated_at_utc: {config['generated_at_utc']}",
                f"method: {method_name}",
                f"model_name: {model_name}",
                f"requested_gpu: {gpu}",
                f"target_channel_id: {target_channel}",
                f"dapi_channel_id: {dapi_channel}",
                f"images_processed: {len(rows)}",
                f"qc_status_counts: {status_counts}",
                "whole_cell_segmentation_validated: false",
                (
                    f"warnings: exploratory output; {target_channel}/{target_channel_name} is target "
                    "signal, not pan-cell stain"
                ),
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _resolved_model_path(model_name: str) -> str | None:
    direct_path = Path(model_name)
    if direct_path.exists():
        return str(direct_path)
    local_model_root = os.environ.get("CELLPOSE_LOCAL_MODELS_PATH")
    candidates = []
    if local_model_root:
        candidates.append(Path(local_model_root) / model_name)
    candidates.append(Path.home() / ".cellpose" / "models" / model_name)
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return None


def _cell_region_method_name(target_channel_id: str, dapi_channel_id: str) -> str:
    target_channel = _normalize_channel_id(target_channel_id).lower()
    dapi_channel = _normalize_channel_id(dapi_channel_id).lower()
    return f"cellpose_{target_channel}_{dapi_channel}_candidate_asma_associated_region"


def _software_versions() -> dict[str, str]:
    versions = {"python": sys.version.split()[0], "platform": platform.platform()}
    for package_name in ["cellpose", "torch", "numpy", "scikit-image", "tifffile"]:
        versions[package_name] = _package_version(package_name)
    return versions


def _package_version(package_name: str) -> str:
    try:
        return importlib_metadata.version(package_name)
    except importlib_metadata.PackageNotFoundError:
        return "not_installed"


def _device_metadata(gpu: bool) -> dict[str, Any]:
    metadata: dict[str, Any] = {"requested_gpu": gpu}
    try:
        import torch
    except ImportError:
        metadata["torch_importable"] = False
        return metadata
    metadata["torch_importable"] = True
    metadata["mps_available"] = bool(torch.backends.mps.is_available())
    metadata["cuda_available"] = bool(torch.cuda.is_available())
    metadata["cuda_device_count"] = int(torch.cuda.device_count())
    if torch.cuda.is_available():
        metadata["cuda_device_name"] = torch.cuda.get_device_name(0)
    return metadata


def _format_scientific(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(number):
        return "NA"
    return f"{number:.2e}"


def _safe_output_stem(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return token or "image"
