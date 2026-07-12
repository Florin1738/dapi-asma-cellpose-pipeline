#!/usr/bin/env python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
import typer


DISPLAY_LOW_PERCENTILE = 1.0
DISPLAY_HIGH_PERCENTILE = 99.8
HISTOGRAM_BINS = 256


app = typer.Typer(
    help=(
        "Search CH2/aSMA images for candidate blank-like or low-aSMA fields and "
        "render visual QC panels."
    )
)


@dataclass(frozen=True)
class ImagePair:
    plate: str
    acquisition: str
    location: str
    ch2_path: Path
    ch4_path: Path | None


@dataclass(frozen=True)
class Candidate:
    category: str
    pair: ImagePair
    metrics: pd.Series


@app.command()
def main(
    data_root: Path = typer.Option(
        Path("data/aSMA_DAPI_plates"),
        "--data-root",
        help="Root folder containing plate image folders.",
    ),
    output_dir: Path = typer.Option(
        Path("output/background_candidate_search"),
        "--output-dir",
        help="Directory for candidate metrics and panels.",
    ),
    candidates_per_category: int = typer.Option(4, "--candidates-per-category", min=2, max=8),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    pairs = _find_image_pairs(data_root)
    metrics = _compute_metrics(pairs)
    selected = _select_candidates(metrics, candidates_per_category=candidates_per_category)

    metrics_path = output_dir / "ch2_background_candidate_all_image_metrics.csv"
    candidates_path = output_dir / "ch2_background_candidate_front_runners.csv"
    estimates_path = output_dir / "ch2_background_candidate_estimates.csv"
    panel_path = output_dir / "ch2_background_candidate_front_runners_panel.png"
    overlay_path = output_dir / "ch2_background_candidate_histogram_overlay.png"

    metrics.to_csv(metrics_path, index=False)
    selected.to_csv(candidates_path, index=False)
    _write_estimates(selected, estimates_path)
    _write_front_runner_panel(selected, panel_path)
    _write_overlay_panel(selected, overlay_path)

    typer.echo(f"metrics={metrics_path}")
    typer.echo(f"front_runners={candidates_path}")
    typer.echo(f"estimates={estimates_path}")
    typer.echo(f"panel={panel_path}")
    typer.echo(f"overlay={overlay_path}")


def _find_image_pairs(data_root: Path) -> list[ImagePair]:
    pairs = []
    for ch2_path in sorted(data_root.glob("**/*_CH2.tif")):
        location = ch2_path.parent.name
        acquisition = ch2_path.parent.parent.name
        plate = ch2_path.parent.parent.parent.name
        ch4_path = ch2_path.with_name(ch2_path.name.replace("_CH2.tif", "_CH4.tif"))
        pairs.append(
            ImagePair(
                plate=plate,
                acquisition=acquisition,
                location=location,
                ch2_path=ch2_path,
                ch4_path=ch4_path if ch4_path.exists() else None,
            )
        )
    if not pairs:
        raise FileNotFoundError(f"No *_CH2.tif files found below {data_root}")
    return pairs


def _compute_metrics(pairs: list[ImagePair]) -> pd.DataFrame:
    rows = []
    for pair in pairs:
        ch2 = _read_channel(pair.ch2_path, channel_index=0)
        ch2_values = ch2.astype(np.float64)
        row = {
            "plate": pair.plate,
            "acquisition": pair.acquisition,
            "location": pair.location,
            "image_id": f"{pair.plate}/{pair.acquisition}/{pair.location}",
            "ch2_path": str(pair.ch2_path),
            "ch4_path": "" if pair.ch4_path is None else str(pair.ch4_path),
            "pixel_count": int(ch2_values.size),
            "ch2_sum": float(np.sum(ch2_values, dtype=np.float64)),
            "ch2_mean": float(np.mean(ch2_values)),
            "ch2_p01": float(np.percentile(ch2_values, 1)),
            "ch2_p05": float(np.percentile(ch2_values, 5)),
            "ch2_p10": float(np.percentile(ch2_values, 10)),
            "ch2_median": float(np.percentile(ch2_values, 50)),
            "ch2_p90": float(np.percentile(ch2_values, 90)),
            "ch2_p95": float(np.percentile(ch2_values, 95)),
            "ch2_p99": float(np.percentile(ch2_values, 99)),
            "ch2_p995": float(np.percentile(ch2_values, 99.5)),
            "ch2_max": float(np.max(ch2_values)),
            "ch2_saturated_fraction": float(np.mean(ch2_values >= _dtype_max(ch2))),
        }
        if pair.ch4_path is not None:
            ch4 = _read_channel(pair.ch4_path, channel_index=2)
            ch4_values = ch4.astype(np.float64)
            row.update(
                {
                    "ch4_mean": float(np.mean(ch4_values)),
                    "ch4_median": float(np.percentile(ch4_values, 50)),
                    "ch4_p90": float(np.percentile(ch4_values, 90)),
                    "ch4_p95": float(np.percentile(ch4_values, 95)),
                    "ch4_p99": float(np.percentile(ch4_values, 99)),
                    "ch4_saturated_fraction": float(np.mean(ch4_values >= _dtype_max(ch4))),
                }
            )
        else:
            row.update(
                {
                    "ch4_mean": np.nan,
                    "ch4_median": np.nan,
                    "ch4_p90": np.nan,
                    "ch4_p95": np.nan,
                    "ch4_p99": np.nan,
                    "ch4_saturated_fraction": np.nan,
                }
            )
        rows.append(row)
    metrics = pd.DataFrame(rows)
    metrics["low_ch2_score"] = _mean_rank(
        metrics,
        ["ch2_mean", "ch2_p90", "ch2_p95", "ch2_p99", "ch2_saturated_fraction"],
        ascending=True,
    )
    metrics["low_dapi_score"] = _mean_rank(
        metrics,
        ["ch4_mean", "ch4_p90", "ch4_p95"],
        ascending=True,
    )
    metrics["dapi_present_score"] = _mean_rank(
        metrics,
        ["ch4_mean", "ch4_p90", "ch4_p95"],
        ascending=False,
    )
    metrics["blank_like_score"] = metrics["low_ch2_score"] + metrics["low_dapi_score"]
    metrics["low_asma_with_dapi_score"] = metrics["low_ch2_score"] + metrics["dapi_present_score"]
    return metrics.sort_values(["plate", "acquisition", "location"]).reset_index(drop=True)


def _mean_rank(metrics: pd.DataFrame, columns: list[str], *, ascending: bool) -> pd.Series:
    ranks = []
    for column in columns:
        ranks.append(metrics[column].rank(method="average", ascending=ascending, na_option="bottom", pct=True))
    return pd.concat(ranks, axis=1).mean(axis=1)


def _select_candidates(metrics: pd.DataFrame, *, candidates_per_category: int) -> pd.DataFrame:
    selected_frames = []
    blank = metrics.sort_values(["blank_like_score", "ch2_p95", "ch4_p95"], ascending=True).head(
        candidates_per_category
    )
    blank = blank.copy()
    blank.insert(0, "candidate_category", "blank_like_low_CH2_low_DAPI")
    selected_frames.append(blank)

    dapi_threshold = metrics["ch4_p95"].median(skipna=True)
    with_cells_pool = metrics[
        (metrics["ch4_p95"] >= dapi_threshold) & (metrics["ch4_saturated_fraction"] <= 0.01)
    ].copy()
    if with_cells_pool.empty:
        with_cells_pool = metrics.copy()
    low_asma = with_cells_pool.sort_values(
        ["ch2_p95", "ch2_p99", "ch2_mean"], ascending=True
    ).head(candidates_per_category)
    low_asma = low_asma.copy()
    low_asma.insert(0, "candidate_category", "low_aSMA_with_DAPI_present")
    selected_frames.append(low_asma)

    selected = pd.concat(selected_frames, ignore_index=True)
    selected.insert(1, "candidate_rank", selected.groupby("candidate_category").cumcount() + 1)
    return selected


def _read_channel(path: Path, *, channel_index: int) -> np.ndarray:
    image = np.squeeze(np.asarray(tifffile.imread(path)))
    if image.ndim == 2:
        return image
    if image.ndim == 3 and image.shape[-1] > channel_index:
        return image[..., channel_index]
    raise ValueError(f"Expected grayscale or RGB TIFF for {path}; got {image.shape}")


def _dtype_max(image: np.ndarray) -> float:
    if np.issubdtype(image.dtype, np.integer):
        return float(np.iinfo(image.dtype).max)
    return max(float(np.nanmax(image)), 1.0)


def _write_estimates(selected: pd.DataFrame, output_path: Path) -> None:
    rows = []
    for category, group in selected.groupby("candidate_category", sort=False):
        arrays = [_read_channel(Path(path), channel_index=0).astype(np.float64).ravel() for path in group["ch2_path"]]
        pooled = np.concatenate(arrays)
        hist_counts, edges = np.histogram(pooled, bins=HISTOGRAM_BINS, range=(0, 65535))
        mode_index = int(np.argmax(hist_counts))
        mode_center = float((edges[mode_index] + edges[mode_index + 1]) / 2)
        rows.append(
            {
                "candidate_category": category,
                "n_images": int(len(group)),
                "pooled_pixels": int(pooled.size),
                "pooled_ch2_mode_bin_center": mode_center,
                "pooled_ch2_p50": float(np.percentile(pooled, 50)),
                "pooled_ch2_p90": float(np.percentile(pooled, 90)),
                "pooled_ch2_p95": float(np.percentile(pooled, 95)),
                "pooled_ch2_p99": float(np.percentile(pooled, 99)),
                "pooled_ch2_p995": float(np.percentile(pooled, 99.5)),
            }
        )
    pd.DataFrame(rows).to_csv(output_path, index=False)


def _write_front_runner_panel(selected: pd.DataFrame, output_path: Path) -> None:
    records = selected.to_dict("records")
    ch2_images = [_read_channel(Path(record["ch2_path"]), channel_index=0) for record in records]
    ch4_images = [
        None if not record["ch4_path"] else _read_channel(Path(record["ch4_path"]), channel_index=2)
        for record in records
    ]
    ch2_lo, ch2_hi = _shared_display_limits([image for image in ch2_images])
    ch4_lo, ch4_hi = _shared_display_limits([image for image in ch4_images if image is not None])
    x_max = max(25000.0, min(65535.0, selected["ch2_p99"].max() * 1.15))

    n_cols = len(records)
    fig, axes = plt.subplots(
        3,
        n_cols,
        figsize=(3.25 * n_cols, 9.1),
        dpi=180,
        constrained_layout=False,
    )
    if n_cols == 1:
        axes = axes.reshape(3, 1)
    fig.subplots_adjust(top=0.84, bottom=0.08, left=0.04, right=0.99, hspace=0.3, wspace=0.22)
    fig.suptitle(
        "Candidate low-expression / background-reference fields",
        fontsize=16,
        fontweight="bold",
        y=0.965,
    )
    fig.text(
        0.5,
        0.925,
        (
            "Top row: CH2/aSMA. Middle: CH4/DAPI. Bottom: raw CH2 histogram as fraction of pixels per bin. "
            "These are candidates, not confirmed controls."
        ),
        ha="center",
        fontsize=9.5,
    )

    for idx, (record, ch2, ch4) in enumerate(zip(records, ch2_images, ch4_images, strict=True)):
        _show_red(axes[0, idx], ch2, lo=ch2_lo, hi=ch2_hi)
        axes[0, idx].set_title(
            f"{record['candidate_category'].replace('_', ' ')}\n#{record['candidate_rank']}: {record['image_id']}",
            fontsize=8.5,
        )
        axes[0, idx].text(
            0.02,
            0.98,
            "\n".join(
                [
                    f"CH2 mean {record['ch2_mean']:,.0f}",
                    f"p95 {record['ch2_p95']:,.0f}",
                    f"p99 {record['ch2_p99']:,.0f}",
                ]
            ),
            transform=axes[0, idx].transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color="white",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "black", "alpha": 0.6, "edgecolor": "none"},
        )
        if ch4 is None:
            axes[1, idx].text(0.5, 0.5, "missing CH4", ha="center", va="center")
            axes[1, idx].axis("off")
        else:
            _show_blue(axes[1, idx], ch4, lo=ch4_lo, hi=ch4_hi)
        axes[1, idx].text(
            0.02,
            0.98,
            f"DAPI p95 {record['ch4_p95']:,.0f}",
            transform=axes[1, idx].transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color="white",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "black", "alpha": 0.6, "edgecolor": "none"},
        )
        _plot_fraction_histogram(axes[2, idx], ch2, x_max=x_max)
        axes[2, idx].axvline(record["ch2_p95"], color="#22c55e", linewidth=1.0, label="p95")
        axes[2, idx].axvline(record["ch2_p99"], color="#38bdf8", linewidth=1.0, label="p99")
        axes[2, idx].set_title("raw CH2 distribution", fontsize=8.5)
        axes[2, idx].legend(fontsize=7, frameon=False)

    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def _write_overlay_panel(selected: pd.DataFrame, output_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), dpi=180, constrained_layout=True)
    colors = {
        "blank_like_low_CH2_low_DAPI": "#64748b",
        "low_aSMA_with_DAPI_present": "#be123c",
    }
    x_max = max(25000.0, min(65535.0, selected["ch2_p99"].max() * 1.15))
    bins = np.linspace(0, x_max, HISTOGRAM_BINS + 1)
    for category, group in selected.groupby("candidate_category", sort=False):
        pooled_arrays = []
        for _idx, row in group.iterrows():
            ch2 = _read_channel(Path(row["ch2_path"]), channel_index=0).astype(np.float64).ravel()
            counts, edges = np.histogram(ch2, bins=bins)
            fraction = counts / counts.sum()
            centers = 0.5 * (edges[:-1] + edges[1:])
            axes[0].plot(centers, fraction, color=colors[category], alpha=0.35, linewidth=1)
            pooled_arrays.append(ch2)
        pooled = np.concatenate(pooled_arrays)
        counts, edges = np.histogram(pooled, bins=bins)
        fraction = counts / counts.sum()
        centers = 0.5 * (edges[:-1] + edges[1:])
        axes[0].plot(centers, fraction, color=colors[category], linewidth=2.4, label=category)
        axes[1].plot(centers, fraction, color=colors[category], linewidth=2.4, label=category)
        for percentile, line_style in [(95, "--"), (99, ":")]:
            axes[1].axvline(
                np.percentile(pooled, percentile),
                color=colors[category],
                linestyle=line_style,
                linewidth=1.4,
                label=f"{category} p{percentile}",
            )

    axes[0].set_title("Candidate normalized histograms")
    axes[0].set_xlabel("raw CH2 pixel value")
    axes[0].set_ylabel("fraction of pixels per bin")
    axes[0].set_xlim(0, x_max)
    axes[0].grid(axis="y", color="#e5e7eb")
    axes[0].legend(fontsize=8, frameon=False)

    axes[1].set_title("Pooled candidate background-range markers")
    axes[1].set_xlabel("raw CH2 pixel value")
    axes[1].set_ylabel("fraction of pixels per bin")
    axes[1].set_xlim(0, x_max)
    axes[1].grid(axis="y", color="#e5e7eb")
    axes[1].legend(fontsize=7.5, frameon=False)

    fig.suptitle(
        "Raw CH2 low-expression candidate distributions",
        fontsize=15,
        fontweight="bold",
    )
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def _shared_display_limits(images: list[np.ndarray]) -> tuple[float, float]:
    values = np.concatenate([image.ravel().astype(np.float64) for image in images])
    finite = values[np.isfinite(values)]
    lo, hi = np.percentile(finite, [DISPLAY_LOW_PERCENTILE, DISPLAY_HIGH_PERCENTILE])
    if hi <= lo:
        hi = max(float(np.nanmax(finite)), lo + 1.0)
    return float(lo), float(hi)


def _show_red(ax: plt.Axes, image: np.ndarray, *, lo: float, hi: float) -> None:
    scaled = np.clip((image.astype(np.float64) - lo) / (hi - lo), 0, 1)
    rgb = np.zeros((*scaled.shape, 3), dtype=np.float64)
    rgb[..., 0] = scaled
    ax.imshow(rgb)
    ax.axis("off")


def _show_blue(ax: plt.Axes, image: np.ndarray, *, lo: float, hi: float) -> None:
    scaled = np.clip((image.astype(np.float64) - lo) / (hi - lo), 0, 1)
    rgb = np.zeros((*scaled.shape, 3), dtype=np.float64)
    rgb[..., 2] = scaled
    ax.imshow(rgb)
    ax.axis("off")


def _plot_fraction_histogram(ax: plt.Axes, image: np.ndarray, *, x_max: float) -> None:
    counts, edges = np.histogram(image.ravel(), bins=HISTOGRAM_BINS, range=(0, x_max))
    fraction = counts / counts.sum()
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]
    ax.bar(centers, fraction, width=width, color="#7f1d1d", alpha=0.85)
    ax.set_xlim(0, x_max)
    ax.set_xlabel("raw CH2")
    ax.set_ylabel("fraction")
    ax.grid(axis="y", color="#e5e7eb", linewidth=0.7)


if __name__ == "__main__":
    app()
