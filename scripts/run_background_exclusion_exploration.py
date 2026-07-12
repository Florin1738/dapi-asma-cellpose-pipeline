#!/usr/bin/env python
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw
from scipy.stats import pearsonr, spearmanr
from skimage.filters import threshold_otsu
import typer

from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.pi_simple_summary import (
    DEFAULT_PLATES,
    _count_lookup_key,
    _resolve_plate_roots,
    find_image_pairs,
    load_counts_by_plate,
)

app = typer.Typer(help="Compare raw aSMA intensity with simple background-excluded variants.")


METHODS = [
    ("raw_per_nucleus", "Raw CH2 / nucleus"),
    ("p10_corrected_per_nucleus", "P10-corrected CH2 / nucleus"),
    ("global_foreground_per_nucleus", "Global foreground CH2 / nucleus"),
    ("per_image_otsu_foreground_per_nucleus", "Per-image Otsu foreground CH2 / nucleus"),
]


@app.command()
def main(
    input_root: Path = typer.Option(Path("data/aSMA_DAPI_plates"), "--input"),
    counts_root: Path = typer.Option(Path("output/pi_simple_summary/cellpose_counts"), "--counts"),
    output_dir: Path = typer.Option(Path("output/background_exclusion_exploration"), "--output"),
    background_percentile: float = typer.Option(10.0, "--background-percentile"),
    sample_pixels_per_image: int = typer.Option(50_000, "--sample-pixels-per-image", min=1000),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    pairs_by_plate = _load_pairs_by_plate(input_root)
    counts_by_plate = load_counts_by_plate(counts_root)
    global_threshold, sample_count = _estimate_global_threshold(
        pairs_by_plate=pairs_by_plate,
        background_percentile=background_percentile,
        sample_pixels_per_image=sample_pixels_per_image,
    )
    rows = _measure_rows(
        pairs_by_plate=pairs_by_plate,
        counts_by_plate=counts_by_plate,
        background_percentile=background_percentile,
        global_threshold=global_threshold,
    )
    _write_csv(output_dir / "background_comparison_image_level.csv", rows)
    _write_csv(output_dir / "background_comparison_method_correlations.csv", _method_correlations(rows))
    _write_csv(output_dir / "background_comparison_plate_summary.csv", _plate_summary(rows))
    _write_plots(rows, output_dir / "plots")
    _write_foreground_qc(rows, output_dir / "qc_global_foreground_examples.png")
    _write_report(
        output_dir / "background_exclusion_report.md",
        rows=rows,
        global_threshold=global_threshold,
        sample_count=sample_count,
        background_percentile=background_percentile,
    )
    typer.echo(f"rows={len(rows)}")
    typer.echo(f"global_foreground_threshold={global_threshold:.6g}")
    typer.echo(f"wrote={output_dir}")


def _load_pairs_by_plate(input_root: Path) -> dict[str, list[Any]]:
    plate_roots = _resolve_plate_roots(input_root, DEFAULT_PLATES)
    return {
        plate_name: find_image_pairs(plate_root) if plate_root is not None else []
        for plate_name, plate_root in plate_roots.items()
    }


def _estimate_global_threshold(
    *,
    pairs_by_plate: dict[str, list[Any]],
    background_percentile: float,
    sample_pixels_per_image: int,
) -> tuple[float, int]:
    samples: list[np.ndarray] = []
    for pairs in pairs_by_plate.values():
        for pair in pairs:
            image, _ = read_primary_intensity_plane(pair.ch2_path)
            corrected = _p_background_corrected(image, background_percentile)
            flat = corrected.ravel()
            if flat.size <= sample_pixels_per_image:
                sample = flat
            else:
                indices = np.linspace(0, flat.size - 1, sample_pixels_per_image, dtype=np.int64)
                sample = flat[indices]
            samples.append(sample.astype(np.float32, copy=False))
    pooled = np.concatenate(samples) if samples else np.array([], dtype=np.float32)
    positive = pooled[pooled > 0]
    if positive.size == 0:
        return 0.0, 0
    return float(threshold_otsu(positive)), int(positive.size)


def _measure_rows(
    *,
    pairs_by_plate: dict[str, list[Any]],
    counts_by_plate: dict[str, dict[str, int]],
    background_percentile: float,
    global_threshold: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for plate_name, pairs in pairs_by_plate.items():
        count_lookup = counts_by_plate.get(plate_name, {})
        for pair in pairs:
            count = count_lookup.get(_count_lookup_key(pair.source_id))
            if count is None:
                raise ValueError(f"No nucleus count for {plate_name} {pair.source_id}")
            image, _ = read_primary_intensity_plane(pair.ch2_path)
            target = image.astype(np.float64)
            raw_sum = float(np.sum(target))
            bg_value = float(np.percentile(target, background_percentile))
            p10_corrected = np.clip(target - bg_value, 0, None)
            p10_sum = float(np.sum(p10_corrected))
            global_mask = p10_corrected > global_threshold
            global_sum = float(np.sum(p10_corrected[global_mask]))
            per_image_threshold = _safe_otsu(p10_corrected)
            per_image_mask = p10_corrected > per_image_threshold
            per_image_sum = float(np.sum(p10_corrected[per_image_mask]))
            rows.append(
                {
                    "plate": plate_name,
                    "source_id": pair.source_id,
                    "location": pair.location,
                    "ch2_path": str(pair.ch2_path),
                    "nuclei_count": int(count),
                    "pixel_count": int(target.size),
                    "raw_integrated": raw_sum,
                    "raw_per_nucleus": _divide(raw_sum, count),
                    "background_percentile": background_percentile,
                    "background_value_per_px": bg_value,
                    "p10_corrected_integrated": p10_sum,
                    "p10_corrected_per_nucleus": _divide(p10_sum, count),
                    "p10_removed_fraction_of_raw": _removed_fraction(raw_sum, p10_sum),
                    "global_foreground_threshold_after_p10": global_threshold,
                    "global_foreground_area_px": int(np.count_nonzero(global_mask)),
                    "global_foreground_area_fraction": float(np.count_nonzero(global_mask) / target.size),
                    "global_foreground_integrated": global_sum,
                    "global_foreground_per_nucleus": _divide(global_sum, count),
                    "global_foreground_removed_fraction_of_p10": _removed_fraction(p10_sum, global_sum),
                    "per_image_otsu_threshold_after_p10": per_image_threshold,
                    "per_image_otsu_area_px": int(np.count_nonzero(per_image_mask)),
                    "per_image_otsu_area_fraction": float(np.count_nonzero(per_image_mask) / target.size),
                    "per_image_otsu_foreground_integrated": per_image_sum,
                    "per_image_otsu_foreground_per_nucleus": _divide(per_image_sum, count),
                }
            )
    return rows


def _p_background_corrected(image: np.ndarray, percentile: float) -> np.ndarray:
    target = image.astype(np.float64)
    background = float(np.percentile(target, percentile))
    return np.clip(target - background, 0, None)


def _safe_otsu(image: np.ndarray) -> float:
    positive = image[image > 0]
    if positive.size == 0:
        return 0.0
    if np.min(positive) == np.max(positive):
        return float(np.min(positive))
    return float(threshold_otsu(positive))


def _method_correlations(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for x_key, x_label in METHODS:
        for y_key, y_label in METHODS:
            if x_key >= y_key:
                continue
            x = np.array([float(row[x_key]) for row in rows], dtype=np.float64)
            y = np.array([float(row[y_key]) for row in rows], dtype=np.float64)
            output.append(
                {
                    "x_metric": x_key,
                    "y_metric": y_key,
                    "x_label": x_label,
                    "y_label": y_label,
                    "pearson_r": _safe_corr(pearsonr, x, y),
                    "spearman_r": _safe_corr(spearmanr, x, y),
                }
            )
    for metric, label in METHODS:
        x = np.array([float(row["nuclei_count"]) for row in rows], dtype=np.float64)
        y = np.array([float(row[metric]) for row in rows], dtype=np.float64)
        output.append(
            {
                "x_metric": "nuclei_count",
                "y_metric": metric,
                "x_label": "Nuclei count",
                "y_label": label,
                "pearson_r": _safe_corr(pearsonr, x, y),
                "spearman_r": _safe_corr(spearmanr, x, y),
            }
        )
    return output


def _plate_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for plate in sorted({row["plate"] for row in rows}):
        plate_rows = [row for row in rows if row["plate"] == plate]
        entry: dict[str, Any] = {"plate": plate, "rows": len(plate_rows)}
        for metric, _label in METHODS:
            values = np.array([float(row[metric]) for row in plate_rows], dtype=np.float64)
            entry[f"{metric}_median"] = float(np.median(values))
            entry[f"{metric}_mean"] = float(np.mean(values))
        entry["median_p10_removed_fraction_of_raw"] = float(
            np.median([float(row["p10_removed_fraction_of_raw"]) for row in plate_rows])
        )
        entry["median_global_foreground_area_fraction"] = float(
            np.median([float(row["global_foreground_area_fraction"]) for row in plate_rows])
        )
        output.append(entry)
    return output


def _write_plots(rows: list[dict[str, Any]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _scatter_raw_vs_adjusted(rows, output_dir / "raw_vs_background_adjusted_per_nucleus.png")
    _scatter_nuclei_vs_metrics(rows, output_dir / "nuclei_count_vs_intensity_metrics.png")
    _bar_background_removed(rows, output_dir / "background_removed_fraction_by_image.png")
    _rank_shift_plot(rows, output_dir / "rank_shift_raw_to_global_foreground.png")


def _scatter_raw_vs_adjusted(rows: list[dict[str, Any]], output_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    comparisons = [
        ("p10_corrected_per_nucleus", "P10-corrected / nucleus"),
        ("global_foreground_per_nucleus", "Global foreground / nucleus"),
    ]
    colors = {"Plate 1": "#2563eb", "Plate 2": "#dc2626"}
    for ax, (metric, ylabel) in zip(axes, comparisons, strict=True):
        for plate in sorted(colors):
            plate_rows = [row for row in rows if row["plate"] == plate]
            x = [float(row["raw_per_nucleus"]) for row in plate_rows]
            y = [float(row[metric]) for row in plate_rows]
            ax.scatter(x, y, s=22, alpha=0.78, label=plate, color=colors[plate])
        max_value = max(
            max(float(row["raw_per_nucleus"]) for row in rows),
            max(float(row[metric]) for row in rows),
        )
        ax.plot([0, max_value], [0, max_value], color="#6b7280", linewidth=1, linestyle="--")
        ax.set_xlabel("Raw CH2 / nucleus")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel + " vs raw")
        ax.grid(alpha=0.25)
    axes[0].legend(frameon=False)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _scatter_nuclei_vs_metrics(rows: list[dict[str, Any]], output_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    metrics = [
        ("raw_integrated", "Raw integrated CH2"),
        ("p10_corrected_integrated", "P10-corrected integrated CH2"),
        ("global_foreground_integrated", "Global foreground integrated CH2"),
    ]
    colors = {"Plate 1": "#2563eb", "Plate 2": "#dc2626"}
    for ax, (metric, title) in zip(axes, metrics, strict=True):
        for plate in sorted(colors):
            plate_rows = [row for row in rows if row["plate"] == plate]
            ax.scatter(
                [float(row["nuclei_count"]) for row in plate_rows],
                [float(row[metric]) / 1_000_000_000 for row in plate_rows],
                s=22,
                alpha=0.78,
                label=plate,
                color=colors[plate],
            )
        x = np.array([float(row["nuclei_count"]) for row in rows])
        y = np.array([float(row[metric]) for row in rows])
        ax.set_title(f"{title}\nPearson r={_safe_corr(pearsonr, x, y):.2f}")
        ax.set_xlabel("Cellpose CH4 nuclei count")
        ax.set_ylabel("Integrated CH2, billions")
        ax.grid(alpha=0.25)
    axes[0].legend(frameon=False)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _bar_background_removed(rows: list[dict[str, Any]], output_path: Path) -> None:
    sorted_rows = sorted(rows, key=lambda row: float(row["p10_removed_fraction_of_raw"]))
    labels = [row["source_id"] for row in sorted_rows]
    values = [100 * float(row["p10_removed_fraction_of_raw"]) for row in sorted_rows]
    fig, ax = plt.subplots(figsize=(24, 5), constrained_layout=True)
    ax.bar(labels, values, color="#0f766e")
    ax.set_title("Fraction of raw CH2 removed by P10 background subtraction")
    ax.set_ylabel("Removed from raw integrated intensity (%)")
    ax.set_xlabel("Run / LOCATION")
    ax.tick_params(axis="x", labelrotation=90, labelsize=5)
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _rank_shift_plot(rows: list[dict[str, Any]], output_path: Path) -> None:
    raw_rank = _rank_map(rows, "raw_per_nucleus")
    fg_rank = _rank_map(rows, "global_foreground_per_nucleus")
    plot_rows = sorted(
        rows,
        key=lambda row: abs(raw_rank[row["source_id"]] - fg_rank[row["source_id"]]),
        reverse=True,
    )[:30]
    labels = [row["source_id"] for row in plot_rows]
    shifts = [raw_rank[row["source_id"]] - fg_rank[row["source_id"]] for row in plot_rows]
    colors = ["#dc2626" if shift < 0 else "#2563eb" for shift in shifts]
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    ax.barh(labels, shifts, color=colors)
    ax.axvline(0, color="#111827", linewidth=1)
    ax.set_title("Largest rank shifts after global foreground-only background exclusion")
    ax.set_xlabel("Raw rank - foreground rank; positive means higher after exclusion")
    ax.invert_yaxis()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def _write_foreground_qc(rows: list[dict[str, Any]], output_path: Path) -> None:
    selected = _representative_rows(rows)
    thumbs: list[tuple[str, Image.Image]] = []
    for row in selected:
        image, _ = read_primary_intensity_plane(Path(row["ch2_path"]))
        target = image.astype(np.float64)
        corrected = np.clip(target - float(row["background_value_per_px"]), 0, None)
        mask = corrected > float(row["global_foreground_threshold_after_p10"])
        thumbs.append((row["source_id"], _foreground_overlay(corrected, mask)))
    _write_contact_sheet(thumbs, output_path)


def _representative_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ranked_raw = sorted(rows, key=lambda row: float(row["raw_per_nucleus"]), reverse=True)
    ranked_fg = sorted(rows, key=lambda row: float(row["global_foreground_per_nucleus"]), reverse=True)
    ranked_removed = sorted(rows, key=lambda row: float(row["p10_removed_fraction_of_raw"]), reverse=True)
    picks: list[dict[str, Any]] = []
    for candidate in (
        ranked_raw[:3]
        + ranked_fg[:3]
        + ranked_removed[:3]
        + [ranked_raw[len(ranked_raw) // 2], ranked_raw[-1]]
    ):
        if candidate["source_id"] not in {row["source_id"] for row in picks}:
            picks.append(candidate)
        if len(picks) == 12:
            break
    return picks


def _foreground_overlay(corrected: np.ndarray, mask: np.ndarray) -> Image.Image:
    low, high = np.percentile(corrected, [1, 99.8])
    scaled = np.clip((corrected.astype(np.float32) - low) / max(high - low, 1), 0, 1)
    rgb = np.dstack([scaled, scaled, scaled])
    rgb[mask] = [1.0, 0.85, 0.0]
    return Image.fromarray((rgb * 255).astype(np.uint8), mode="RGB")


def _write_contact_sheet(thumbs: list[tuple[str, Image.Image]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    columns = 3
    thumb_width = 420
    padding = 18
    label_height = 28
    resized = []
    for label, image in thumbs:
        scale = thumb_width / image.width
        resized.append((label, image.resize((thumb_width, int(image.height * scale)))))
    cell_height = max(image.height for _, image in resized) + label_height
    rows = (len(resized) + columns - 1) // columns
    sheet = Image.new(
        "RGB",
        (columns * thumb_width + (columns + 1) * padding, rows * cell_height + (rows + 1) * padding),
        "white",
    )
    draw = ImageDraw.Draw(sheet)
    for index, (label, image) in enumerate(resized):
        row, column = divmod(index, columns)
        x_pos = padding + column * (thumb_width + padding)
        y_pos = padding + row * (cell_height + padding)
        draw.text((x_pos, y_pos), label, fill=(0, 0, 0))
        sheet.paste(image, (x_pos, y_pos + label_height))
    sheet.save(output_path)


def _write_report(
    path: Path,
    *,
    rows: list[dict[str, Any]],
    global_threshold: float,
    sample_count: int,
    background_percentile: float,
) -> None:
    corr_rows = _method_correlations(rows)
    by_pair = {(row["x_metric"], row["y_metric"]): row for row in corr_rows}
    raw_p10 = by_pair[("p10_corrected_per_nucleus", "raw_per_nucleus")]
    raw_global = by_pair[("global_foreground_per_nucleus", "raw_per_nucleus")]
    raw_count = by_pair[("nuclei_count", "raw_per_nucleus")]
    p10_count = by_pair[("nuclei_count", "p10_corrected_per_nucleus")]
    fg_count = by_pair[("nuclei_count", "global_foreground_per_nucleus")]
    top_raw = _top_ids(rows, "raw_per_nucleus")
    top_p10 = _top_ids(rows, "p10_corrected_per_nucleus")
    top_global = _top_ids(rows, "global_foreground_per_nucleus")
    overlap_p10 = len(set(top_raw).intersection(top_p10))
    overlap_global = len(set(top_raw).intersection(top_global))
    median_removed = float(np.median([float(row["p10_removed_fraction_of_raw"]) for row in rows]))
    median_fg_area = float(np.median([float(row["global_foreground_area_fraction"]) for row in rows]))
    lines = [
        "# Background Exclusion Exploration",
        "",
        "This exploratory analysis compares the PI workbook's raw `CH2`/aSMA integrated "
        "intensity with background-adjusted variants. It does not replace the PI workbook.",
        "",
        "## Methods",
        "",
        "- Raw: sum of all `CH2` pixels.",
        (
            f"- P{background_percentile:g} background-corrected: subtract the image-wide "
            f"{background_percentile:g}th percentile CH2 value from every pixel, clip at zero, "
            "then sum the whole field."
        ),
        (
            "- Global foreground-only: after P10 correction, keep only pixels above one global "
            f"Otsu threshold estimated from {sample_count:,} sampled positive corrected pixels "
            f"across the batch. Threshold: {global_threshold:.6g}."
        ),
        "- Per-image Otsu foreground-only is included as a sensitivity check.",
        "",
        "## Main Checks",
        "",
        f"- Median fraction removed by P10 baseline subtraction: {median_removed:.1%}.",
        f"- Median image area kept by global foreground threshold: {median_fg_area:.1%}.",
        (
            "- Raw per-nucleus vs P10-corrected per-nucleus: "
            f"Pearson r={float(raw_p10['pearson_r']):.3f}, "
            f"Spearman r={float(raw_p10['spearman_r']):.3f}."
        ),
        (
            "- Raw per-nucleus vs global-foreground per-nucleus: "
            f"Pearson r={float(raw_global['pearson_r']):.3f}, "
            f"Spearman r={float(raw_global['spearman_r']):.3f}."
        ),
        (
            f"- Top-10 overlap raw vs P10-corrected: {overlap_p10}/10; "
            f"raw vs global foreground: {overlap_global}/10."
        ),
        (
            "- Nuclei-count correlation with per-nucleus metrics: "
            f"raw Pearson r={float(raw_count['pearson_r']):.3f}; "
            f"P10-corrected r={float(p10_count['pearson_r']):.3f}; "
            f"global foreground r={float(fg_count['pearson_r']):.3f}."
        ),
        "",
        "## Interpretation Guardrails",
        "",
        "- P10 correction removes a scalar image baseline; it does not spatially segment cells.",
        "- Global foreground thresholding excludes dim CH2 pixels, so it is stricter and more sensitive to threshold choice.",
        "- These are image-level exploratory metrics; biological interpretation still needs plate map/treatment metadata.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _top_ids(rows: list[dict[str, Any]], metric: str, n: int = 10) -> list[str]:
    return [
        row["source_id"]
        for row in sorted(rows, key=lambda row: float(row[metric]), reverse=True)[:n]
    ]


def _rank_map(rows: list[dict[str, Any]], metric: str) -> dict[str, int]:
    ranked = sorted(rows, key=lambda row: float(row[metric]), reverse=True)
    return {row["source_id"]: index + 1 for index, row in enumerate(ranked)}


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _safe_corr(fn: Any, x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan")
    return float(fn(x, y).statistic)


def _divide(value: float, denominator: int) -> float:
    return float(value / denominator) if denominator > 0 else float("nan")


def _removed_fraction(before: float, after: float) -> float:
    return float(1 - after / before) if before > 0 else float("nan")


if __name__ == "__main__":
    app()
