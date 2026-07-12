from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from skimage import segmentation


METHOD_REVIEW_COLUMNS = [
    "image_id",
    "dapi_nucleus_count",
    "propagation_region_area_px",
    "cellpose_region_area_px",
    "both_region_area_px",
    "propagation_only_area_px",
    "cellpose_only_area_px",
    "union_area_px",
    "method_region_jaccard",
    "crop_box",
    "propagation_per_DAPI_positive_nucleus",
    "cellpose_per_DAPI_positive_nucleus",
    "propagation_qc_status",
    "cellpose_qc_status",
    "propagation_qc_flags",
    "cellpose_qc_flags",
    "interpretation",
]


@dataclass(frozen=True)
class MethodReviewRecord:
    image_id: str
    ch2_image: np.ndarray
    ch4_image: np.ndarray
    nuclei_mask: np.ndarray
    propagation_labels: np.ndarray
    cellpose_labels: np.ndarray
    propagation_metrics: dict[str, Any]
    cellpose_metrics: dict[str, Any]


def write_method_review_package(
    *,
    records: Iterable[MethodReviewRecord],
    output_dir: Path,
    crop_size: int = 280,
) -> dict[str, Path]:
    record_list = list(records)
    if not record_list:
        raise ValueError("Cannot write a method review package without records")
    output_dir.mkdir(parents=True, exist_ok=True)

    summaries: list[dict[str, Any]] = []
    for record in record_list:
        _validate_record(record)
        crop_box = select_review_crop_box(record, crop_size=crop_size)
        summaries.append(
            summarize_method_overlap(
                image_id=record.image_id,
                nuclei_mask=record.nuclei_mask,
                propagation_labels=record.propagation_labels,
                cellpose_labels=record.cellpose_labels,
                propagation_metrics=record.propagation_metrics,
                cellpose_metrics=record.cellpose_metrics,
                crop_box=crop_box,
            )
        )

    full_field_panel = output_dir / "method_comparison_full_field_panel.png"
    crop_panel = output_dir / "method_comparison_crop_panel.png"
    summary_csv = output_dir / "method_comparison_review_summary.csv"
    readme = output_dir / "README.md"
    _write_review_panel(
        records=record_list,
        summaries=summaries,
        output_path=full_field_panel,
        crop=False,
    )
    _write_review_panel(
        records=record_list,
        summaries=summaries,
        output_path=crop_panel,
        crop=True,
    )
    _write_summary_csv(summary_csv, summaries)
    _write_readme(readme, image_count=len(record_list))
    return {
        "full_field_panel": full_field_panel,
        "crop_panel": crop_panel,
        "summary_csv": summary_csv,
        "readme": readme,
    }


def summarize_method_overlap(
    *,
    image_id: str,
    nuclei_mask: np.ndarray,
    propagation_labels: np.ndarray,
    cellpose_labels: np.ndarray,
    propagation_metrics: dict[str, Any],
    cellpose_metrics: dict[str, Any],
    crop_box: tuple[int, int, int, int],
) -> dict[str, Any]:
    nuclei = np.asarray(nuclei_mask)
    propagation = np.asarray(propagation_labels) > 0
    cellpose = np.asarray(cellpose_labels) > 0
    both = propagation & cellpose
    prop_only = propagation & ~cellpose
    cellpose_only = cellpose & ~propagation
    union = propagation | cellpose
    union_area = int(np.count_nonzero(union))
    return {
        "image_id": image_id,
        "dapi_nucleus_count": _count_nonzero_labels(nuclei),
        "propagation_region_area_px": int(np.count_nonzero(propagation)),
        "cellpose_region_area_px": int(np.count_nonzero(cellpose)),
        "both_region_area_px": int(np.count_nonzero(both)),
        "propagation_only_area_px": int(np.count_nonzero(prop_only)),
        "cellpose_only_area_px": int(np.count_nonzero(cellpose_only)),
        "union_area_px": union_area,
        "method_region_jaccard": float(np.count_nonzero(both) / union_area)
        if union_area
        else float("nan"),
        "crop_box": ",".join(str(part) for part in crop_box),
        "propagation_per_DAPI_positive_nucleus": _float_or_nan(
            propagation_metrics.get("seeded_region_intensity_per_DAPI_positive_nucleus")
        ),
        "cellpose_per_DAPI_positive_nucleus": _float_or_nan(
            cellpose_metrics.get("target_integrated_intensity_per_DAPI_positive_nucleus")
        ),
        "propagation_qc_status": propagation_metrics.get("qc_status", ""),
        "cellpose_qc_status": cellpose_metrics.get("qc_status", ""),
        "propagation_qc_flags": propagation_metrics.get("qc_flags", ""),
        "cellpose_qc_flags": cellpose_metrics.get("qc_flags", ""),
        "interpretation": "qualitative_qc_only_not_manual_validation",
    }


def select_review_crop_box(
    record: MethodReviewRecord,
    *,
    crop_size: int,
) -> tuple[int, int, int, int]:
    ch2 = np.asarray(record.ch2_image, dtype=np.float64)
    propagation = np.asarray(record.propagation_labels) > 0
    cellpose = np.asarray(record.cellpose_labels) > 0
    nuclei = np.asarray(record.nuclei_mask) > 0
    disagreement = propagation ^ cellpose
    high_ch2 = ch2 >= np.percentile(ch2, 97) if ch2.size else np.zeros_like(ch2, dtype=bool)
    target = disagreement | (high_ch2 & (propagation | cellpose)) | nuclei
    if not np.any(target):
        target = np.ones_like(ch2, dtype=bool)
    ys, xs = np.nonzero(target)
    weights = ch2[ys, xs] + 1.0
    center_y = int(round(float(np.average(ys, weights=weights))))
    center_x = int(round(float(np.average(xs, weights=weights))))
    return _crop_box_around(center_y, center_x, ch2.shape, crop_size=crop_size)


def _write_review_panel(
    *,
    records: list[MethodReviewRecord],
    summaries: list[dict[str, Any]],
    output_path: Path,
    crop: bool,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = len(records)
    fig, axes = plt.subplots(
        n_rows,
        6,
        figsize=(18, max(3.2, 3.15 * n_rows)),
        squeeze=False,
        constrained_layout=True,
        width_ratios=[1, 1, 1, 1, 1, 0.82],
    )
    for row_index, (record, summary) in enumerate(zip(records, summaries, strict=True)):
        arrays = _arrays_for_display(record, crop_box=_parse_crop_box(summary["crop_box"]) if crop else None)
        ch2_scaled = _scale_for_display(arrays["ch2"])
        ch4_scaled = _scale_for_display(arrays["ch4"])
        nuclei = arrays["nuclei"]
        propagation = arrays["propagation"]
        cellpose = arrays["cellpose"]

        dapi_rgb = np.dstack([ch4_scaled, ch4_scaled, ch4_scaled])
        _draw_nuclei_context(dapi_rgb, nuclei)
        prop_overlay = _region_overlay(ch2_scaled, propagation, fill_color=np.array([0.55, 1.0, 0.55]))
        cellpose_overlay = _region_overlay(
            ch2_scaled,
            cellpose,
            fill_color=np.array([0.40, 0.85, 1.0]),
        )
        disagreement = _disagreement_overlay(ch2_scaled, propagation, cellpose)

        titles = [
            f"{record.image_id} CH4/DAPI\nGreen X nuclei",
            "CH2/aSMA raw",
            "Propagation/Otsu\nlight green",
            "Cellpose CH2+CH4\nlight blue",
            "Method disagreement\nboth green, prop magenta, Cellpose cyan",
            "Metrics",
        ]
        panels = [dapi_rgb, ch2_scaled, prop_overlay, cellpose_overlay, disagreement]
        for col_index, panel in enumerate(panels):
            if col_index == 1:
                axes[row_index, col_index].imshow(panel, cmap="gray")
            else:
                axes[row_index, col_index].imshow(panel)
            axes[row_index, col_index].set_title(titles[col_index], fontsize=8)
            axes[row_index, col_index].axis("off")
        axes[row_index, 5].axis("off")
        axes[row_index, 5].text(
            0,
            1,
            _summary_text(summary),
            ha="left",
            va="top",
            fontsize=7.4,
            transform=axes[row_index, 5].transAxes,
            color="#1f2933",
        )
    title = "Full-field" if not crop else "Matched-crop"
    fig.suptitle(
        f"{title} comparison: candidate aSMA-region methods. Qualitative QC only, not validation.",
        fontsize=12,
    )
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def _arrays_for_display(
    record: MethodReviewRecord,
    *,
    crop_box: tuple[int, int, int, int] | None,
) -> dict[str, np.ndarray]:
    arrays = {
        "ch2": np.asarray(record.ch2_image),
        "ch4": np.asarray(record.ch4_image),
        "nuclei": np.asarray(record.nuclei_mask),
        "propagation": np.asarray(record.propagation_labels) > 0,
        "cellpose": np.asarray(record.cellpose_labels) > 0,
    }
    if crop_box is None:
        return arrays
    y0, x0, y1, x1 = crop_box
    return {name: value[y0:y1, x0:x1] for name, value in arrays.items()}


def _summary_text(summary: dict[str, Any]) -> str:
    return (
        f"DAPI nuclei: {summary['dapi_nucleus_count']}\n"
        f"Propagation / nucleus: {_format_scientific(summary['propagation_per_DAPI_positive_nucleus'])}\n"
        f"Cellpose / nucleus: {_format_scientific(summary['cellpose_per_DAPI_positive_nucleus'])}\n"
        f"Jaccard: {_format_float(summary['method_region_jaccard'])}\n"
        f"Prop only px: {summary['propagation_only_area_px']}\n"
        f"Cellpose only px: {summary['cellpose_only_area_px']}\n"
        f"Prop QC: {summary['propagation_qc_status']}\n"
        f"Cellpose QC: {summary['cellpose_qc_status']}\n"
        "No manual ground truth here."
    )


def _region_overlay(ch2_scaled: np.ndarray, region: np.ndarray, *, fill_color: np.ndarray) -> np.ndarray:
    rgb = np.dstack([ch2_scaled, ch2_scaled, ch2_scaled])
    if np.any(region):
        rgb[region] = 0.55 * rgb[region] + 0.45 * fill_color
        rgb[segmentation.find_boundaries(region, mode="outer")] = fill_color
    return rgb


def _disagreement_overlay(ch2_scaled: np.ndarray, propagation: np.ndarray, cellpose: np.ndarray) -> np.ndarray:
    rgb = np.dstack([ch2_scaled, ch2_scaled, ch2_scaled])
    both = propagation & cellpose
    prop_only = propagation & ~cellpose
    cellpose_only = cellpose & ~propagation
    _blend(rgb, both, np.array([0.30, 1.0, 0.35]), alpha=0.55)
    _blend(rgb, prop_only, np.array([1.0, 0.10, 0.85]), alpha=0.60)
    _blend(rgb, cellpose_only, np.array([0.10, 0.85, 1.0]), alpha=0.60)
    return rgb


def _draw_nuclei_context(rgb: np.ndarray, nuclei_mask: np.ndarray) -> None:
    nuclei = np.asarray(nuclei_mask)
    boundaries = segmentation.find_boundaries(nuclei, mode="outer")
    rgb[boundaries] = [0.0, 0.75, 1.0]
    for label in np.unique(nuclei):
        if label == 0:
            continue
        ys, xs = np.nonzero(nuclei == label)
        if len(ys) == 0:
            continue
        y = int(round(float(np.mean(ys))))
        x = int(round(float(np.mean(xs))))
        _draw_cross(rgb, y=y, x=x, size=max(2, min(7, min(rgb.shape[:2]) // 70)))


def _draw_cross(rgb: np.ndarray, *, y: int, x: int, size: int) -> None:
    for offset in range(-size, size + 1):
        for y_pos, x_pos in [(y + offset, x + offset), (y + offset, x - offset)]:
            if 0 <= y_pos < rgb.shape[0] and 0 <= x_pos < rgb.shape[1]:
                rgb[y_pos, x_pos] = [0.0, 1.0, 0.0]


def _blend(rgb: np.ndarray, mask: np.ndarray, color: np.ndarray, *, alpha: float) -> None:
    if np.any(mask):
        rgb[mask] = (1.0 - alpha) * rgb[mask] + alpha * color


def _validate_record(record: MethodReviewRecord) -> None:
    shapes = {
        "ch2": np.asarray(record.ch2_image).shape,
        "ch4": np.asarray(record.ch4_image).shape,
        "nuclei": np.asarray(record.nuclei_mask).shape,
        "propagation": np.asarray(record.propagation_labels).shape,
        "cellpose": np.asarray(record.cellpose_labels).shape,
    }
    if len(set(shapes.values())) != 1:
        raise ValueError(f"{record.image_id} shape mismatch: {shapes}")
    for name in ["ch2", "ch4", "nuclei", "propagation", "cellpose"]:
        if len(shapes[name]) != 2:
            raise ValueError(f"{record.image_id} expected 2-D arrays, got {name} shape {shapes[name]}")


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


def _crop_box_around(
    center_y: int,
    center_x: int,
    shape: tuple[int, int],
    *,
    crop_size: int,
) -> tuple[int, int, int, int]:
    height, width = shape
    size = max(8, int(crop_size))
    half = size // 2
    y0 = max(0, int(center_y) - half)
    x0 = max(0, int(center_x) - half)
    y1 = min(height, y0 + size)
    x1 = min(width, x0 + size)
    y0 = max(0, y1 - size)
    x0 = max(0, x1 - size)
    return y0, x0, y1, x1


def _parse_crop_box(value: str) -> tuple[int, int, int, int]:
    parts = [int(part) for part in str(value).split(",")]
    if len(parts) != 4:
        raise ValueError(f"Expected crop box with four comma-separated values, got {value!r}")
    return tuple(parts)  # type: ignore[return-value]


def _write_summary_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=METHOD_REVIEW_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in METHOD_REVIEW_COLUMNS})


def _write_readme(path: Path, *, image_count: int) -> None:
    path.write_text(
        "\n".join(
            [
                "# Candidate Method Review Panels",
                "",
                f"Images reviewed: `{image_count}`",
                "",
                "These panels compare exploratory candidate aSMA-region methods on the same fields and crops.",
                "",
                "- Propagation/Otsu is the CellProfiler-style DAPI-seeded aSMA-associated region method.",
                "- Cellpose CH2+CH4 is the pretrained Cellpose candidate-region method.",
                "- The disagreement panel shows overlap and method-only pixels.",
                "",
                "This package is qualitative QC only. It is not manual validation and does not report segmentation accuracy. Precision, recall, F1, and IoU require completed manual/reference masks.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _float_or_nan(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _format_scientific(value: Any) -> str:
    number = _float_or_nan(value)
    if not np.isfinite(number):
        return "NA"
    return f"{number:.3e}"


def _format_float(value: Any) -> str:
    number = _float_or_nan(value)
    if not np.isfinite(number):
        return "NA"
    return f"{number:.3f}"


def _count_nonzero_labels(mask: np.ndarray) -> int:
    labels = np.unique(mask)
    return int(np.count_nonzero(labels))
