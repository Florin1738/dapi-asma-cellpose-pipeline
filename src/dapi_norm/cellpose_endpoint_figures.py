from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from skimage import segmentation
import tifffile

from dapi_norm.image_arrays import read_primary_intensity_plane


DEFAULT_SELECTED_FIELDS_PER_PLATE = 12
DEFAULT_REPRESENTATIVE_FIELDS_PER_PLATE = 4
OLD_CELLPOSE_RETAINED_INTENSITY = "dapi_anchored_cellpose_ch2_integrated_background_corrected"
OLD_CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS = (
    "dapi_anchored_cellpose_ch2_integrated_background_corrected_per_DAPI_positive_nucleus"
)
CELLPOSE_RETAINED_INTENSITY = "dapi_anchored_cellpose_target_integrated_background_corrected"
CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS = (
    "dapi_anchored_cellpose_target_integrated_background_corrected_per_DAPI_positive_nucleus"
)
CELLPOSE_RETAINED_AREA = "dapi_anchored_cellpose_masked_area_px"
CELLPOSE_RETAINED_AREA_PER_NUCLEUS = "dapi_anchored_cellpose_masked_area_per_DAPI_positive_nucleus"


def render_cellpose_endpoint_figures(
    *,
    summary_csv: Path | str,
    output_dir: Path | str,
    panel_page_size: int = 12,
    selected_fields_per_plate: int = DEFAULT_SELECTED_FIELDS_PER_PLATE,
    max_overlay_images: int | None = None,
) -> dict[str, Any]:
    rows = _read_summary(Path(summary_csv))
    if not rows:
        raise ValueError(f"No rows found in {summary_csv}")
    output_path = Path(output_dir)
    panels_dir = output_path / "cellpose_overlay_pages"
    output_path.mkdir(parents=True, exist_ok=True)
    panels_dir.mkdir(parents=True, exist_ok=True)

    selected_rows = _select_rows_for_bar_figures(
        rows,
        per_plate=selected_fields_per_plate,
        metric=CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS,
    )
    sorted_rows = sorted(
        rows,
        key=lambda row: _number(row, CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS),
        reverse=True,
    )
    if max_overlay_images is not None:
        sorted_rows = sorted_rows[:max_overlay_images]

    metric_contrast_path = output_path / "figure_1_endpoint_metric_contrast.png"
    masking_effect_path = output_path / "figure_2_masking_effect_contrast.png"
    plate_summary_path = output_path / "figure_3_plate_level_endpoint_summary.png"
    representative_segmentation_path = output_path / "figure_4_representative_cell_segmentation_examples.png"
    captions_path = output_path / "FIGURE_CAPTIONS.md"

    _write_metric_contrast_figure(selected_rows, metric_contrast_path)
    _write_masking_effect_figure(selected_rows, masking_effect_path)
    _write_plate_summary_figure(rows, plate_summary_path)
    _write_representative_cell_segmentation_panel(
        rows,
        representative_segmentation_path,
        per_plate=DEFAULT_REPRESENTATIVE_FIELDS_PER_PLATE,
    )
    overlay_pages, overlay_index_rows = _write_overlay_pages(
        sorted_rows,
        panels_dir,
        page_size=panel_page_size,
    )
    overlay_index_path = _write_overlay_index(panels_dir / "overlay_index.csv", overlay_index_rows)
    _write_captions(captions_path, rows=rows, overlay_pages=overlay_pages)

    return {
        "metric_contrast": metric_contrast_path,
        "masking_effect": masking_effect_path,
        "plate_summary": plate_summary_path,
        "representative_cell_segmentation": representative_segmentation_path,
        "overlay_pages": overlay_pages,
        "overlay_index": overlay_index_path,
        "captions_markdown": captions_path,
    }


def render_dapi_nuclei_qc_pages(
    *,
    summary_csv: Path | str,
    output_dir: Path | str,
    page_size: int = 12,
) -> dict[str, Any]:
    rows = _read_summary(Path(summary_csv))
    if not rows:
        raise ValueError(f"No rows found in {summary_csv}")
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    all_pages: list[Path] = []
    all_index_rows: list[dict[str, Any]] = []
    for plate, plate_rows in _rows_by_plate(rows).items():
        plate_dir = output_path / plate.replace(" ", "_")
        sorted_rows = sorted(plate_rows, key=_field_sort_key)
        pages, index_rows = _write_dapi_nuclei_overlay_pages(
            sorted_rows,
            plate_dir,
            page_size=page_size,
        )
        all_pages.extend(pages)
        all_index_rows.extend(index_rows)

    index_path = _write_dapi_overlay_index(output_path / "dapi_nuclei_overlay_index.csv", all_index_rows)
    return {
        "pages": all_pages,
        "index": index_path,
        "field_count": len(rows),
        "plate_count": len(_rows_by_plate(rows)),
    }


def render_representative_cell_segmentation_panel(
    *,
    summary_csv: Path | str,
    output_path: Path | str,
    per_plate: int = DEFAULT_REPRESENTATIVE_FIELDS_PER_PLATE,
) -> Path:
    rows = _read_summary(Path(summary_csv))
    if not rows:
        raise ValueError(f"No rows found in {summary_csv}")
    output = Path(output_path)
    _write_representative_cell_segmentation_panel(rows, output, per_plate=per_plate)
    return output


def _read_summary(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        _ensure_alias(row, "whole_field_target_integrated_raw", "whole_field_ch2_integrated_raw")
        _ensure_alias(
            row,
            "whole_field_target_integrated_raw_per_DAPI_positive_nucleus",
            "whole_field_ch2_integrated_raw_per_DAPI_positive_nucleus",
        )
        _ensure_alias(
            row,
            "cellpose_masked_target_integrated_raw",
            "cellpose_masked_ch2_integrated_raw",
        )
        _ensure_alias(
            row,
            "cellpose_masked_target_integrated_raw_per_DAPI_positive_nucleus",
            "cellpose_masked_ch2_integrated_raw_per_DAPI_positive_nucleus",
        )
        _ensure_alias(
            row,
            "cellpose_masked_target_integrated_background_corrected",
            "cellpose_masked_ch2_integrated_background_corrected",
        )
        _ensure_alias(
            row,
            "cellpose_masked_target_integrated_background_corrected_per_DAPI_positive_nucleus",
            "cellpose_masked_ch2_integrated_background_corrected_per_DAPI_positive_nucleus",
        )
        _ensure_alias(row, CELLPOSE_RETAINED_INTENSITY, OLD_CELLPOSE_RETAINED_INTENSITY)
        _ensure_alias(
            row,
            CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS,
            OLD_CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS,
        )
        _ensure_alias(row, "target_path", "ch2_path")
        _ensure_alias(row, "dapi_path", "ch4_path")
        row.setdefault("target_channel_id", "CH2")
        row.setdefault("dapi_channel_id", "CH4")
    return rows


def _ensure_alias(row: dict[str, str], preferred: str, fallback: str) -> None:
    preferred_value = row.get(preferred, "")
    fallback_value = row.get(fallback, "")
    if not preferred_value and fallback_value:
        row[preferred] = fallback_value
    if not fallback_value and preferred_value:
        row[fallback] = preferred_value


def _rows_by_plate(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    by_plate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_plate[row.get("plate", "Plate")].append(row)
    return dict(sorted(by_plate.items()))


def _select_rows_for_bar_figures(
    rows: list[dict[str, str]],
    *,
    per_plate: int,
    metric: str,
) -> list[dict[str, str]]:
    by_plate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_plate[row.get("plate", "Plate")].append(row)

    selected: list[dict[str, str]] = []
    for plate in sorted(by_plate):
        plate_rows = sorted(by_plate[plate], key=lambda row: _number(row, metric))
        if len(plate_rows) <= per_plate:
            selected.extend(plate_rows)
            continue
        indices = np.linspace(0, len(plate_rows) - 1, per_plate)
        seen: set[str] = set()
        for index in indices:
            row = plate_rows[int(round(index))]
            key = row.get("source_id", row.get("location", ""))
            if key not in seen:
                selected.append(row)
                seen.add(key)
    return selected


def _write_metric_contrast_figure(rows: list[dict[str, str]], path: Path) -> None:
    labels = [_short_label(row) for row in rows]
    target = _target_label(rows)
    x = np.arange(len(rows))
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(max(16, 0.72 * len(rows)), 15),
        sharex=True,
        constrained_layout=True,
    )
    _apply_figure_style(fig, axes)

    axes[0].bar(
        x,
        [_number(row, "whole_field_target_integrated_raw") for row in rows],
        color="#6b7280",
        edgecolor="#111827",
        linewidth=0.7,
    )
    axes[0].set_title(f"1. Whole-field {target} signal before Cellpose masking", fontweight="bold")
    axes[0].set_ylabel(f"Total raw {target} intensity", fontweight="bold")

    axes[1].bar(
        x,
        [_number(row, CELLPOSE_RETAINED_INTENSITY) for row in rows],
        color="#39b54a",
        alpha=0.78,
        edgecolor="#0f5f26",
        linewidth=0.7,
    )
    axes[1].set_title(f"2. Cellpose retained-region {target} signal", fontweight="bold")
    axes[1].set_ylabel(f"Retained {target} intensity", fontweight="bold")

    axes[2].bar(
        x,
        [_number(row, CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS) for row in rows],
        color="#1d4ed8",
        alpha=0.78,
        edgecolor="#111827",
        linewidth=0.7,
    )
    axes[2].set_title(
        f"3. Cellpose retained-region {target} signal per DAPI-positive nucleus",
        fontweight="bold",
    )
    axes[2].set_ylabel(f"Retained {target} intensity per nucleus", fontweight="bold")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, rotation=60, ha="right", fontweight="bold")

    for ax in axes:
        ax.yaxis.set_major_formatter(FuncFormatter(_scientific_tick))
    fig.suptitle(
        "How Cellpose retained-region masking and DAPI-positive nucleus counting change the reported aSMA signal",
        fontsize=22,
        fontweight="bold",
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _write_masking_effect_figure(rows: list[dict[str, str]], path: Path) -> None:
    labels = [_short_label(row) for row in rows]
    target = _target_label(rows)
    x = np.arange(len(rows))
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(max(16, 0.72 * len(rows)), 14),
        sharex=True,
        constrained_layout=True,
    )
    _apply_figure_style(fig, axes)

    whole = np.asarray([_number(row, "whole_field_target_integrated_raw") for row in rows])
    retained = np.asarray([_number(row, CELLPOSE_RETAINED_INTENSITY) for row in rows])
    retained_area = [_number(row, CELLPOSE_RETAINED_AREA) for row in rows]
    nuclei_counts = [_number(row, "dapi_positive_nucleus_count") for row in rows]

    axes[0].bar(x, _safe_fraction(retained, whole) * 100.0, color="#39b54a", edgecolor="#0f5f26")
    axes[0].set_title(f"1. Fraction of whole-field {target} signal retained by Cellpose", fontweight="bold")
    axes[0].set_ylabel(f"Retained {target} (%)", fontweight="bold")

    axes[1].bar(x, retained_area, color="#1d4ed8", edgecolor="#172554")
    axes[1].set_title("2. Cellpose retained positive area", fontweight="bold")
    axes[1].set_ylabel("Retained area (px)", fontweight="bold")

    axes[2].bar(x, nuclei_counts, color="#f59e0b", edgecolor="#78350f")
    axes[2].set_title("3. DAPI-positive nuclei used as the denominator", fontweight="bold")
    axes[2].set_ylabel("Nuclei count", fontweight="bold")
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, rotation=60, ha="right", fontweight="bold")

    axes[0].yaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f"{value:.0f}%"))
    axes[1].yaxis.set_major_formatter(FuncFormatter(_scientific_tick))
    fig.suptitle(
        "Cellpose retained-region signal, retained area, and DAPI-positive nucleus denominator",
        fontsize=22,
        fontweight="bold",
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _write_plate_summary_figure(rows: list[dict[str, str]], path: Path) -> None:
    plates = sorted({row.get("plate", "Plate") for row in rows})
    target = _target_label(rows)
    intensity_metrics = [
        ("Whole field / nucleus", "whole_field_target_integrated_raw_per_DAPI_positive_nucleus"),
        ("Cellpose retained region / nucleus", CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS),
    ]
    area_metrics = [
        ("Cellpose retained area / nucleus", CELLPOSE_RETAINED_AREA_PER_NUCLEUS),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(16, 8), constrained_layout=True)
    _apply_figure_style(fig, axes)
    _plot_grouped_plate_summary(axes[0], rows, plates, intensity_metrics, ylabel="Median intensity per nucleus")
    _plot_grouped_plate_summary(axes[1], rows, plates, area_metrics, ylabel="Median area per nucleus (px)")
    axes[0].set_title(f"{target} intensity endpoints by plate", fontweight="bold")
    axes[1].set_title("Cellpose retained positive area by plate", fontweight="bold")
    axes[0].yaxis.set_major_formatter(FuncFormatter(_scientific_tick))
    fig.suptitle("Plate-level median endpoint summary", fontsize=22, fontweight="bold")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _write_representative_cell_segmentation_panel(
    rows: list[dict[str, str]],
    path: Path,
    *,
    per_plate: int,
) -> None:
    selected_rows = _select_rows_for_bar_figures(
        rows,
        per_plate=per_plate,
        metric=CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS,
    )
    selected_rows = sorted(
        selected_rows,
        key=lambda row: (
            _plate_sort_value(row.get("plate", "")),
            _number(row, CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS),
            row.get("source_id", ""),
        ),
    )
    columns = 2
    target = _target_label(rows)
    grid_rows = int(np.ceil(len(selected_rows) / columns))
    fig, axes = plt.subplots(
        grid_rows,
        columns,
        figsize=(18, max(5.4, 4.9 * grid_rows)),
        constrained_layout=True,
    )
    axes_array = np.atleast_1d(axes).ravel()
    for ax, row in zip(axes_array, selected_rows):
        ax.imshow(_render_paired_review_image(row))
        ax.set_title(_representative_overlay_title(row), fontsize=12.5, fontweight="bold")
        ax.axis("off")
    for ax in axes_array[len(selected_rows) :]:
        ax.axis("off")
    fig.suptitle(
        f"Representative Cellpose retained-region examples across low, middle, and high {target} signal",
        fontsize=21,
        fontweight="bold",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _plot_grouped_plate_summary(
    ax: plt.Axes,
    rows: list[dict[str, str]],
    plates: list[str],
    metrics: list[tuple[str, str]],
    *,
    ylabel: str,
) -> None:
    x = np.arange(len(plates))
    width = 0.78 / len(metrics)
    colors = ["#6b7280", "#39b54a", "#1d4ed8"]
    for index, (label, column) in enumerate(metrics):
        centers = x - 0.39 + width / 2 + index * width
        medians = []
        lows = []
        highs = []
        for plate in plates:
            values = np.asarray([_number(row, column) for row in rows if row.get("plate") == plate])
            values = values[np.isfinite(values)]
            medians.append(float(np.median(values)) if len(values) else np.nan)
            lows.append(float(np.percentile(values, 25)) if len(values) else np.nan)
            highs.append(float(np.percentile(values, 75)) if len(values) else np.nan)
        yerr = np.vstack([
            np.asarray(medians) - np.asarray(lows),
            np.asarray(highs) - np.asarray(medians),
        ])
        ax.bar(
            centers,
            medians,
            width=width,
            color=colors[index % len(colors)],
            alpha=0.75,
            edgecolor="#111827",
            linewidth=0.7,
            yerr=yerr,
            capsize=4,
            label=label,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(plates, fontweight="bold")
    ax.set_ylabel(ylabel, fontweight="bold")
    ax.legend(loc="upper left", frameon=False, prop={"weight": "bold", "size": 12})


def _write_overlay_pages(
    rows: list[dict[str, str]],
    output_dir: Path,
    *,
    page_size: int,
) -> tuple[list[Path], list[dict[str, Any]]]:
    pages: list[Path] = []
    index_rows: list[dict[str, Any]] = []
    target = _target_label(rows)
    dapi = _dapi_label(rows)
    for page_index, start in enumerate(range(0, len(rows), page_size), start=1):
        page_rows = rows[start : start + page_size]
        columns = 2
        grid_rows = int(np.ceil(len(page_rows) / columns))
        fig, axes = plt.subplots(
            grid_rows,
            columns,
            figsize=(18, max(5.2, 4.9 * grid_rows)),
            constrained_layout=True,
        )
        axes_array = np.atleast_1d(axes).ravel()
        for ax, row in zip(axes_array, page_rows):
            ax.imshow(_render_paired_review_image(row))
            ax.set_title(_overlay_title(row), fontsize=13, fontweight="bold")
            ax.axis("off")
            tile_index = len(index_rows)
            index_rows.append(
                {
                    "source_id": row.get("source_id", ""),
                    "plate": row.get("plate", ""),
                    "location": row.get("location", ""),
                    "page": output_path_name(page_index),
                    "page_path": str(output_dir / output_path_name(page_index)),
                    "tile_number_on_page": (tile_index % page_size) + 1,
                    "sort_metric": row.get(CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS, ""),
                    "qc_status": row.get("qc_status", ""),
                    "qc_flags": row.get("qc_flags", ""),
                }
            )
        for ax in axes_array[len(page_rows) :]:
            ax.axis("off")
        fig.suptitle(
            f"Cellpose retained-region QC: raw {target} in red, retained regions in green, "
            f"{dapi} nuclei in blue",
            fontsize=20,
            fontweight="bold",
        )
        output_path = output_dir / output_path_name(page_index)
        fig.savefig(output_path, dpi=180)
        plt.close(fig)
        pages.append(output_path)
    return pages, index_rows


def output_path_name(page_index: int) -> str:
    return f"cellpose_overlay_sorted_by_endpoint_page_{page_index:02d}.png"


def _write_overlay_index(path: Path, rows: list[dict[str, Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "source_id",
        "plate",
        "location",
        "page",
        "page_path",
        "tile_number_on_page",
        "sort_metric",
        "qc_status",
        "qc_flags",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _write_dapi_nuclei_overlay_pages(
    rows: list[dict[str, str]],
    output_dir: Path,
    *,
    page_size: int,
) -> tuple[list[Path], list[dict[str, Any]]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    pages: list[Path] = []
    index_rows: list[dict[str, Any]] = []
    dapi = _dapi_label(rows)
    for page_index, start in enumerate(range(0, len(rows), page_size), start=1):
        page_rows = rows[start : start + page_size]
        columns = 2
        grid_rows = int(np.ceil(len(page_rows) / columns))
        fig, axes = plt.subplots(
            grid_rows,
            columns,
            figsize=(18, max(5.2, 4.9 * grid_rows)),
            constrained_layout=True,
        )
        axes_array = np.atleast_1d(axes).ravel()
        for ax, row in zip(axes_array, page_rows):
            ax.imshow(_render_dapi_paired_review_image(row))
            ax.set_title(_dapi_overlay_title(row), fontsize=13, fontweight="bold")
            ax.axis("off")
            tile_index = len(index_rows)
            page_name = dapi_output_path_name(page_index)
            index_rows.append(
                {
                    "source_id": row.get("source_id", ""),
                    "plate": row.get("plate", ""),
                    "location": row.get("location", ""),
                    "page": page_name,
                    "page_path": str(output_dir / page_name),
                    "tile_number_on_page": (tile_index % page_size) + 1,
                    "dapi_positive_nucleus_count": row.get("dapi_positive_nucleus_count", ""),
                    "dapi_path": row.get("dapi_path", row.get("ch4_path", "")),
                    "ch4_path": row.get("ch4_path", ""),
                    "dapi_nuclei_mask_path": row.get("dapi_nuclei_mask_path", ""),
                }
            )
        for ax in axes_array[len(page_rows) :]:
            ax.axis("off")
        fig.suptitle(
            f"DAPI nuclei QC: raw {dapi} in blue, detected nuclei outlined in yellow",
            fontsize=20,
            fontweight="bold",
        )
        output_path = output_dir / dapi_output_path_name(page_index)
        fig.savefig(output_path, dpi=180)
        plt.close(fig)
        pages.append(output_path)
    return pages, index_rows


def dapi_output_path_name(page_index: int) -> str:
    return f"dapi_nuclei_overlay_page_{page_index:02d}.png"


def _write_dapi_overlay_index(path: Path, rows: list[dict[str, Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "source_id",
        "plate",
        "location",
        "page",
        "page_path",
        "tile_number_on_page",
        "dapi_positive_nucleus_count",
        "dapi_path",
        "ch4_path",
        "dapi_nuclei_mask_path",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _render_overlay(row: dict[str, str]) -> np.ndarray:
    ch2 = _read_plane(Path(row.get("target_path") or row["ch2_path"]))
    dapi_path = row.get("dapi_path") or row.get("ch4_path")
    ch4 = _read_plane(Path(dapi_path)) if dapi_path else np.zeros_like(ch2)
    labels = np.asarray(tifffile.imread(row["cellpose_mask_path"]))
    nuclei = np.asarray(tifffile.imread(row["dapi_nuclei_mask_path"]))
    if ch2.shape != labels.shape or ch2.shape != nuclei.shape:
        raise ValueError(f"Overlay shape mismatch for {row.get('source_id')}")

    rgb = np.dstack([_scale_for_display(ch2)] * 3)
    region_mask = labels > 0
    anchored_ids = _cellpose_object_ids_with_dapi_centroid(nuclei, labels)
    anchored_mask = _mask_for_object_ids(labels, anchored_ids)
    no_dapi_mask = region_mask & ~anchored_mask
    dapi_display = _scale_for_display(ch4)
    dapi_signal = dapi_display > np.percentile(dapi_display, 99.2) if np.any(dapi_display) else False

    _blend_color(rgb, anchored_mask, np.array([0.30, 1.0, 0.38]), alpha=0.18)
    _blend_color(rgb, dapi_signal, np.array([0.10, 0.38, 1.0]), alpha=0.28)
    anchored_labels = np.where(anchored_mask, labels, 0)
    no_dapi_labels = np.where(no_dapi_mask, labels, 0)
    rgb[segmentation.find_boundaries(anchored_labels, mode="outer")] = [0.0, 0.80, 0.12]
    if np.any(no_dapi_mask):
        rgb[segmentation.find_boundaries(no_dapi_labels, mode="outer")] = [
            1.0,
            0.55,
            0.0,
        ]
    _draw_nucleus_marks(rgb, nuclei)
    return np.clip(rgb, 0, 1)


def _render_red_intensity(row: dict[str, str]) -> np.ndarray:
    ch2 = _read_plane(Path(row.get("target_path") or row["ch2_path"]))
    scaled = _scale_for_display(ch2)
    rgb = np.zeros((*scaled.shape, 3), dtype=np.float32)
    rgb[..., 0] = scaled
    rgb[..., 1] = scaled * 0.06
    rgb[..., 2] = scaled * 0.06
    return np.clip(rgb, 0, 1)


def _render_paired_review_image(row: dict[str, str]) -> np.ndarray:
    left = _render_red_intensity(row)
    right = _render_overlay(row)
    if left.shape != right.shape:
        raise ValueError(f"Paired image shape mismatch for {row.get('source_id')}")
    gutter = np.ones((left.shape[0], max(8, left.shape[1] // 45), 3), dtype=np.float32)
    return np.concatenate([left, gutter, right], axis=1)


def _render_blue_dapi(row: dict[str, str]) -> np.ndarray:
    ch4 = _read_plane(Path(row.get("dapi_path") or row["ch4_path"]))
    scaled = _scale_for_display(ch4)
    rgb = np.zeros((*scaled.shape, 3), dtype=np.float32)
    rgb[..., 0] = scaled * 0.05
    rgb[..., 1] = scaled * 0.12
    rgb[..., 2] = scaled
    return np.clip(rgb, 0, 1)


def _render_dapi_nuclei_overlay(row: dict[str, str]) -> np.ndarray:
    rgb = _render_blue_dapi(row)
    nuclei = np.asarray(tifffile.imread(row["dapi_nuclei_mask_path"]))
    if rgb.shape[:2] != nuclei.shape:
        raise ValueError(f"DAPI overlay shape mismatch for {row.get('source_id')}")
    nuclei_mask = nuclei > 0
    _blend_color(rgb, nuclei_mask, np.array([0.0, 1.0, 0.95]), alpha=0.16)
    rgb[segmentation.find_boundaries(nuclei, mode="outer")] = [1.0, 0.86, 0.0]
    _draw_nucleus_marks_colored(
        rgb,
        nuclei,
        color=np.array([1.0, 0.80, 0.0]),
        halo=np.array([0.0, 0.0, 0.0]),
    )
    return np.clip(rgb, 0, 1)


def _render_dapi_paired_review_image(row: dict[str, str]) -> np.ndarray:
    left = _render_blue_dapi(row)
    right = _render_dapi_nuclei_overlay(row)
    if left.shape != right.shape:
        raise ValueError(f"DAPI paired image shape mismatch for {row.get('source_id')}")
    gutter = np.ones((left.shape[0], max(8, left.shape[1] // 45), 3), dtype=np.float32)
    return np.concatenate([left, gutter, right], axis=1)


def _write_captions(path: Path, *, rows: list[dict[str, str]], overlay_pages: list[Path]) -> None:
    target = _target_label(rows)
    dapi = _dapi_label(rows)
    lines = [
        "# Figure Captions",
        "",
        f"**Figure 1. Endpoint metric contrast across representative fields.** Bars use the same field order across rows. The first row shows whole-field raw {target} integrated intensity before masking. The second row shows {target} integrated intensity inside Cellpose retained regions. The third row divides the Cellpose retained-region intensity by the DAPI-positive nucleus count.",
        "",
        f"**Figure 2. Cellpose retained-region signal, retained area, and DAPI-positive nucleus denominator.** The first row shows the percentage of whole-field raw {target} signal retained inside Cellpose regions. The second row shows the retained positive area in pixels. The third row shows the DAPI-positive nucleus count used as the denominator.",
        "",
        "**Figure 3. Plate-level endpoint summary.** Bars show plate-level medians; error bars show the interquartile range. DAPI is used only as a DAPI-positive nucleus count denominator.",
        "",
        f"**Figure 4. Representative Cellpose retained-region examples.** Examples were selected from low-to-high quantiles of the Cellpose retained-region {target} integrated intensity per DAPI-positive nucleus within each plate. Each tile has raw {target} in red on the left and the Cellpose retained-region overlay on the right.",
        "",
        f"**Cellpose overlay pages.** Each tile has two views of the same field: raw {target} intensity in red on the left, and the Cellpose retained-region overlay on the right. The overlay shows retained regions in transparent green, {dapi} signal/centroids in blue, and orange boundaries for Cellpose objects that contain no DAPI nucleus centroid. The overlay pages are sorted from high to low Cellpose retained-region {target} integrated intensity per DAPI-positive nucleus.",
        "",
    ]
    for page in overlay_pages:
        lines.append(f"- `{page.name}`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_plane(path: Path) -> np.ndarray:
    try:
        image, _pages = read_primary_intensity_plane(path)
        return np.asarray(image)
    except Exception:
        return np.asarray(tifffile.imread(path))


def _scale_for_display(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image, dtype=np.float32)
    if arr.size == 0:
        return arr
    low, high = np.percentile(arr, [1, 99.8])
    if high <= low:
        low = float(np.min(arr))
        high = float(np.max(arr))
    if high <= low:
        return np.zeros_like(arr, dtype=np.float32)
    return np.clip((arr - low) / (high - low), 0, 1)


def _blend_color(rgb: np.ndarray, mask: np.ndarray | bool, color: np.ndarray, *, alpha: float) -> None:
    if isinstance(mask, bool):
        if not mask:
            return
    elif not np.any(mask):
        return
    rgb[mask] = (1.0 - alpha) * rgb[mask] + alpha * color


def _draw_nucleus_marks(rgb: np.ndarray, nuclei: np.ndarray) -> None:
    _draw_nucleus_marks_colored(
        rgb,
        nuclei,
        color=np.array([0.0, 0.28, 1.0]),
        halo=np.array([1.0, 1.0, 1.0]),
    )


def _draw_nucleus_marks_colored(
    rgb: np.ndarray,
    nuclei: np.ndarray,
    *,
    color: np.ndarray,
    halo: np.ndarray,
) -> None:
    labels = np.unique(nuclei)
    labels = labels[labels > 0]
    size = max(1, min(3, round(min(nuclei.shape) / 260)))
    radius = size + 2
    for label in labels:
        ys, xs = np.where(nuclei == label)
        if len(ys) == 0:
            continue
        y = int(round(float(np.mean(ys))))
        x = int(round(float(np.mean(xs))))
        _draw_cross(rgb, y, x, size + 1, halo)
        _draw_cross(rgb, y, x, size, color)
        _draw_circle(rgb, y, x, radius + 1, halo)
        _draw_circle(rgb, y, x, radius, color)


def _draw_cross(rgb: np.ndarray, y: int, x: int, size: int, color: np.ndarray) -> None:
    for offset in range(-size, size + 1):
        for yy, xx in [(y + offset, x + offset), (y + offset, x - offset)]:
            if 0 <= yy < rgb.shape[0] and 0 <= xx < rgb.shape[1]:
                rgb[yy, xx] = color


def _draw_circle(rgb: np.ndarray, y: int, x: int, radius: int, color: np.ndarray) -> None:
    for angle in np.linspace(0, 2 * np.pi, max(16, radius * 8), endpoint=False):
        yy = int(round(y + radius * np.sin(angle)))
        xx = int(round(x + radius * np.cos(angle)))
        if 0 <= yy < rgb.shape[0] and 0 <= xx < rgb.shape[1]:
            rgb[yy, xx] = color


def _cellpose_object_ids_with_dapi_centroid(nuclei_mask: np.ndarray, cellpose_labels: np.ndarray) -> set[int]:
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


def _apply_figure_style(fig: plt.Figure, axes: np.ndarray | list[plt.Axes]) -> None:
    for ax in np.atleast_1d(axes).ravel():
        ax.tick_params(axis="both", labelsize=12, width=1.2)
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontweight("bold")
        ax.grid(axis="y", alpha=0.24, linewidth=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.patch.set_facecolor("white")


def _overlay_title(row: dict[str, str]) -> str:
    return (
        f"{_short_label(row)} | nuclei={_number(row, 'dapi_positive_nucleus_count'):.0f}\n"
        "Cellpose intensity="
        f"{_format_scientific(_number(row, CELLPOSE_RETAINED_INTENSITY))}; "
        "per nucleus="
        f"{_format_scientific(_number(row, CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS))}"
    )


def _representative_overlay_title(row: dict[str, str]) -> str:
    return (
        f"{_short_label(row)} | nuclei={_number(row, 'dapi_positive_nucleus_count'):.0f} | "
        f"area={_number(row, CELLPOSE_RETAINED_AREA):.0f}px\n"
        "retained intensity="
        f"{_format_scientific(_number(row, CELLPOSE_RETAINED_INTENSITY))}; "
        "per nucleus="
        f"{_format_scientific(_number(row, CELLPOSE_RETAINED_INTENSITY_PER_NUCLEUS))}"
    )


def _target_label(rows: list[dict[str, str]]) -> str:
    row = rows[0] if rows else {}
    channel = (row.get("target_channel_id") or "CH2").upper()
    return f"{channel}/aSMA"


def _dapi_label(rows: list[dict[str, str]]) -> str:
    row = rows[0] if rows else {}
    channel = (row.get("dapi_channel_id") or "CH4").upper()
    return f"{channel}/DAPI"


def _short_label(row: dict[str, str]) -> str:
    plate = row.get("plate", "").replace("Plate ", "P")
    source = row.get("source_id", row.get("location", ""))
    if "/" in source:
        run, xy = source.rsplit("/", 1)
        return f"{plate} {run[-6:]}/{xy}"
    return f"{plate} {source}"


def _dapi_overlay_title(row: dict[str, str]) -> str:
    dapi = _dapi_label([row])
    return (
        f"{_short_label(row)} | DAPI nuclei="
        f"{_number(row, 'dapi_positive_nucleus_count'):.0f}\n"
        f"left: raw {dapi}; right: detected nuclei overlay"
    )


def _field_sort_key(row: dict[str, str]) -> tuple[int, str]:
    location = row.get("location", "")
    digits = "".join(character for character in location if character.isdigit())
    return (int(digits) if digits else 9999, row.get("source_id", ""))


def _plate_sort_value(value: str) -> tuple[int, str]:
    digits = "".join(character for character in value if character.isdigit())
    return (int(digits) if digits else 9999, value)


def _number(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "nan"))
    except (TypeError, ValueError):
        return float("nan")


def _safe_fraction(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    return np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=np.float64),
        where=np.isfinite(denominator) & (denominator != 0),
    )


def _scientific_tick(value: float, _pos: int) -> str:
    if not np.isfinite(value):
        return ""
    if value == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / (10**exponent)
    if -2 <= exponent <= 3:
        return f"{value:.0f}"
    return f"{mantissa:.1f}x10^{exponent}"


def _format_scientific(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / (10**exponent)
    return f"{mantissa:.2f}x10^{exponent}"
