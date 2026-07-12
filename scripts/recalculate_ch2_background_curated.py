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


HISTOGRAM_BINS = 256
DISPLAY_LOW_PERCENTILE = 1.0
DISPLAY_HIGH_PERCENTILE = 99.8


app = typer.Typer(
    help="Recalculate CH2 background candidates from manually curated images."
)


@dataclass(frozen=True)
class CuratedCandidate:
    user_rank: int
    original_panel_index: int
    image_id: str
    ch2_path: Path
    ch4_path: Path
    note: str


DEFAULT_CANDIDATES = (
    CuratedCandidate(
        user_rank=1,
        original_panel_index=2,
        image_id="plate 1/ApYYM20AGGSMA_02/XY12",
        ch2_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_02/XY12/"
            "ApYYM20AGGSMA_XY12_CH2.tif"
        ),
        ch4_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_02/XY12/"
            "ApYYM20AGGSMA_XY12_CH4.tif"
        ),
        note="best overall by visual review",
    ),
    CuratedCandidate(
        user_rank=2,
        original_panel_index=5,
        image_id="plate 1/ApYYM20AGGSMA_01/XY41",
        ch2_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY41/"
            "ApYYM20AGGSMA_XY41_CH2.tif"
        ),
        ch4_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY41/"
            "ApYYM20AGGSMA_XY41_CH4.tif"
        ),
        note="accepted low CH2 with DAPI present",
    ),
    CuratedCandidate(
        user_rank=3,
        original_panel_index=6,
        image_id="plate 1/ApYYM20AGGSMA_01/XY40",
        ch2_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY40/"
            "ApYYM20AGGSMA_XY40_CH2.tif"
        ),
        ch4_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY40/"
            "ApYYM20AGGSMA_XY40_CH4.tif"
        ),
        note="accepted low CH2 with DAPI present",
    ),
    CuratedCandidate(
        user_rank=4,
        original_panel_index=1,
        image_id="plate 1/ApYYM20AGGSMA_01/XY01",
        ch2_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY01/"
            "ApYYM20AGGSMA_XY01_CH2.tif"
        ),
        ch4_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY01/"
            "ApYYM20AGGSMA_XY01_CH4.tif"
        ),
        note="accepted but faint CH2 structure visible",
    ),
)


@app.command()
def main(
    output_dir: Path = typer.Option(
        Path("output/background_candidate_search/curated_user_selected"),
        "--output-dir",
        help="Directory for curated recalculation outputs.",
    ),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    records = _candidate_records(list(DEFAULT_CANDIDATES))
    records_path = output_dir / "curated_ch2_background_selected_image_metrics.csv"
    estimates_path = output_dir / "curated_ch2_background_estimates.csv"
    panel_path = output_dir / "curated_ch2_background_candidates_panel.png"
    overlay_path = output_dir / "curated_ch2_background_histogram_overlay.png"

    records.to_csv(records_path, index=False)
    estimates = _estimate_sets(records)
    estimates.to_csv(estimates_path, index=False)
    _write_panel(records, panel_path)
    _write_overlay(records, estimates, overlay_path)

    typer.echo(f"selected_metrics={records_path}")
    typer.echo(f"estimates={estimates_path}")
    typer.echo(f"panel={panel_path}")
    typer.echo(f"overlay={overlay_path}")


def _candidate_records(candidates: list[CuratedCandidate]) -> pd.DataFrame:
    rows = []
    for candidate in candidates:
        ch2 = _read_channel(candidate.ch2_path, channel_index=0)
        ch4 = _read_channel(candidate.ch4_path, channel_index=2)
        values = ch2.astype(np.float64)
        ch4_values = ch4.astype(np.float64)
        row = {
            "user_rank": candidate.user_rank,
            "original_panel_index": candidate.original_panel_index,
            "image_id": candidate.image_id,
            "note": candidate.note,
            "ch2_path": str(candidate.ch2_path),
            "ch4_path": str(candidate.ch4_path),
            "pixel_count": int(values.size),
            "ch2_mode_bin_center": _mode_bin_center(values),
            "ch2_mean": float(np.mean(values)),
            "ch2_p01": float(np.percentile(values, 1)),
            "ch2_p05": float(np.percentile(values, 5)),
            "ch2_p10": float(np.percentile(values, 10)),
            "ch2_p50": float(np.percentile(values, 50)),
            "ch2_p90": float(np.percentile(values, 90)),
            "ch2_p95": float(np.percentile(values, 95)),
            "ch2_p99": float(np.percentile(values, 99)),
            "ch2_p995": float(np.percentile(values, 99.5)),
            "ch2_max": float(np.max(values)),
            "ch2_saturated_fraction": float(np.mean(values >= _dtype_max(ch2))),
            "ch4_p95": float(np.percentile(ch4_values, 95)),
            "ch4_saturated_fraction": float(np.mean(ch4_values >= _dtype_max(ch4))),
        }
        rows.append(row)
    return pd.DataFrame(rows).sort_values("user_rank").reset_index(drop=True)


def _estimate_sets(records: pd.DataFrame) -> pd.DataFrame:
    sets = [
        ("best_single_XY12", records[records["user_rank"] == 1]),
        ("top3_XY12_XY41_XY40", records[records["user_rank"].isin([1, 2, 3])]),
        ("accepted4_XY12_XY41_XY40_XY01", records),
    ]
    rows = []
    for name, group in sets:
        arrays = [_read_channel(Path(path), channel_index=0).astype(np.float64).ravel() for path in group["ch2_path"]]
        pooled = np.concatenate(arrays)
        rows.append(
            {
                "estimate_set": name,
                "n_images": int(len(group)),
                "image_ids": ";".join(group["image_id"].tolist()),
                "pooled_pixels": int(pooled.size),
                "ch2_mode_bin_center": _mode_bin_center(pooled),
                "ch2_mean": float(np.mean(pooled)),
                "ch2_p01": float(np.percentile(pooled, 1)),
                "ch2_p05": float(np.percentile(pooled, 5)),
                "ch2_p10": float(np.percentile(pooled, 10)),
                "ch2_p50": float(np.percentile(pooled, 50)),
                "ch2_p90": float(np.percentile(pooled, 90)),
                "ch2_p95": float(np.percentile(pooled, 95)),
                "ch2_p99": float(np.percentile(pooled, 99)),
                "ch2_p995": float(np.percentile(pooled, 99.5)),
                "recommended_use": _recommended_use(name),
            }
        )
    return pd.DataFrame(rows)


def _recommended_use(name: str) -> str:
    if name == "best_single_XY12":
        return "strictest single-image candidate; best visual blank"
    if name == "top3_XY12_XY41_XY40":
        return "preferred pooled baseline candidate from best visual images"
    return "accepted sensitivity set including original image 1"


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


def _mode_bin_center(values: np.ndarray) -> float:
    counts, edges = np.histogram(values, bins=HISTOGRAM_BINS, range=(0, 65535))
    idx = int(np.argmax(counts))
    return float((edges[idx] + edges[idx + 1]) / 2)


def _write_panel(records: pd.DataFrame, output_path: Path) -> None:
    rows = records.to_dict("records")
    ch2_images = [_read_channel(Path(row["ch2_path"]), channel_index=0) for row in rows]
    ch4_images = [_read_channel(Path(row["ch4_path"]), channel_index=2) for row in rows]
    ch2_lo, ch2_hi = _display_limits(ch2_images)
    ch4_lo, ch4_hi = _display_limits(ch4_images)
    x_max = max(18000.0, records["ch2_p99"].max() * 1.15)

    fig, axes = plt.subplots(
        3,
        len(rows),
        figsize=(3.8 * len(rows), 9.5),
        dpi=180,
        constrained_layout=False,
    )
    fig.subplots_adjust(top=0.84, bottom=0.08, left=0.05, right=0.99, hspace=0.3, wspace=0.24)
    fig.suptitle("Curated CH2 background candidates", fontsize=16, fontweight="bold", y=0.965)
    fig.text(
        0.5,
        0.925,
        (
            "Manual selection from prior panel: original columns 2, 5, 6, and 1. "
            "Original columns 3 and 4 were excluded as CH2-positive."
        ),
        ha="center",
        fontsize=9.5,
    )
    for idx, (row, ch2, ch4) in enumerate(zip(rows, ch2_images, ch4_images, strict=True)):
        _show_channel(axes[0, idx], ch2, color="red", lo=ch2_lo, hi=ch2_hi)
        axes[0, idx].set_title(
            f"rank {row['user_rank']} | original #{row['original_panel_index']}\n{row['image_id']}",
            fontsize=9,
        )
        axes[0, idx].text(
            0.02,
            0.98,
            f"CH2 p50 {row['ch2_p50']:,.0f}\np95 {row['ch2_p95']:,.0f}\np99 {row['ch2_p99']:,.0f}",
            transform=axes[0, idx].transAxes,
            ha="left",
            va="top",
            fontsize=8,
            color="white",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "black", "alpha": 0.62, "edgecolor": "none"},
        )
        _show_channel(axes[1, idx], ch4, color="blue", lo=ch4_lo, hi=ch4_hi)
        axes[1, idx].text(
            0.02,
            0.98,
            f"DAPI p95 {row['ch4_p95']:,.0f}",
            transform=axes[1, idx].transAxes,
            ha="left",
            va="top",
            fontsize=8,
            color="white",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "black", "alpha": 0.62, "edgecolor": "none"},
        )
        _plot_fraction_hist(axes[2, idx], ch2, x_max=x_max)
        axes[2, idx].axvline(row["ch2_p50"], color="#fbbf24", linewidth=1.1, label="p50")
        axes[2, idx].axvline(row["ch2_p95"], color="#22c55e", linewidth=1.1, label="p95")
        axes[2, idx].axvline(row["ch2_p99"], color="#38bdf8", linewidth=1.1, label="p99")
        axes[2, idx].legend(fontsize=7, frameon=False)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def _write_overlay(records: pd.DataFrame, estimates: pd.DataFrame, output_path: Path) -> None:
    rows = records.to_dict("records")
    x_max = max(18000.0, records["ch2_p99"].max() * 1.15)
    bins = np.linspace(0, x_max, HISTOGRAM_BINS + 1)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), dpi=180, constrained_layout=True)
    colors = ["#0f766e", "#2563eb", "#9333ea", "#78716c"]
    for row, color in zip(rows, colors, strict=True):
        ch2 = _read_channel(Path(row["ch2_path"]), channel_index=0).astype(np.float64).ravel()
        counts, edges = np.histogram(ch2, bins=bins)
        fraction = counts / counts.sum()
        centers = 0.5 * (edges[:-1] + edges[1:])
        axes[0].plot(
            centers,
            fraction,
            color=color,
            linewidth=1.8,
            label=f"rank {row['user_rank']} original #{row['original_panel_index']}",
        )
    axes[0].set_title("Selected image histograms")
    axes[0].set_xlabel("raw CH2 pixel value")
    axes[0].set_ylabel("fraction of pixels per bin")
    axes[0].set_xlim(0, x_max)
    axes[0].grid(axis="y", color="#e5e7eb")
    axes[0].legend(fontsize=8, frameon=False)

    set_colors = {
        "best_single_XY12": "#0f766e",
        "top3_XY12_XY41_XY40": "#be123c",
        "accepted4_XY12_XY41_XY40_XY01": "#475569",
    }
    for _idx, row in estimates.iterrows():
        arrays = []
        for image_id in row["image_ids"].split(";"):
            path = records.loc[records["image_id"] == image_id, "ch2_path"].iloc[0]
            arrays.append(_read_channel(Path(path), channel_index=0).astype(np.float64).ravel())
        pooled = np.concatenate(arrays)
        counts, edges = np.histogram(pooled, bins=bins)
        fraction = counts / counts.sum()
        centers = 0.5 * (edges[:-1] + edges[1:])
        color = set_colors[row["estimate_set"]]
        axes[1].plot(centers, fraction, color=color, linewidth=2.2, label=row["estimate_set"])
        axes[1].axvline(row["ch2_p50"], color=color, linestyle="-", linewidth=1.2)
        axes[1].axvline(row["ch2_p95"], color=color, linestyle="--", linewidth=1.2)
    axes[1].set_title("Pooled recalculated estimates: solid p50, dashed p95")
    axes[1].set_xlabel("raw CH2 pixel value")
    axes[1].set_ylabel("fraction of pixels per bin")
    axes[1].set_xlim(0, x_max)
    axes[1].grid(axis="y", color="#e5e7eb")
    axes[1].legend(fontsize=7.5, frameon=False)
    fig.suptitle("Curated CH2 background recalculation", fontsize=15, fontweight="bold")
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def _display_limits(images: list[np.ndarray]) -> tuple[float, float]:
    values = np.concatenate([image.ravel().astype(np.float64) for image in images])
    lo, hi = np.percentile(values[np.isfinite(values)], [DISPLAY_LOW_PERCENTILE, DISPLAY_HIGH_PERCENTILE])
    if hi <= lo:
        hi = max(float(np.nanmax(values)), lo + 1)
    return float(lo), float(hi)


def _show_channel(ax: plt.Axes, image: np.ndarray, *, color: str, lo: float, hi: float) -> None:
    scaled = np.clip((image.astype(np.float64) - lo) / (hi - lo), 0, 1)
    rgb = np.zeros((*scaled.shape, 3), dtype=np.float64)
    if color == "red":
        rgb[..., 0] = scaled
    elif color == "blue":
        rgb[..., 2] = scaled
    else:
        raise ValueError(color)
    ax.imshow(rgb)
    ax.axis("off")


def _plot_fraction_hist(ax: plt.Axes, image: np.ndarray, *, x_max: float) -> None:
    counts, edges = np.histogram(image.ravel(), bins=HISTOGRAM_BINS, range=(0, x_max))
    fraction = counts / counts.sum()
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]
    ax.bar(centers, fraction, width=width, color="#7f1d1d", alpha=0.82)
    ax.set_xlabel("raw CH2")
    ax.set_ylabel("fraction")
    ax.set_xlim(0, x_max)
    ax.grid(axis="y", color="#e5e7eb")


if __name__ == "__main__":
    app()
