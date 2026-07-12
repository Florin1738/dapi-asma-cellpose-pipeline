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


THRESHOLD = 8100.0
HISTOGRAM_BINS = 256
DISPLAY_LOW_PERCENTILE = 1.0
DISPLAY_HIGH_PERCENTILE = 99.8


app = typer.Typer(
    help="Render a CH2 threshold-retained-area QC panel for low, medium, and high expression images."
)


@dataclass(frozen=True)
class Example:
    category: str
    image_id: str
    ch2_path: Path


EXAMPLES = (
    Example(
        category="low",
        image_id="plate 1/ApYYM20AGGSMA_02/XY12",
        ch2_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_02/XY12/"
            "ApYYM20AGGSMA_XY12_CH2.tif"
        ),
    ),
    Example(
        category="low",
        image_id="plate 1/ApYYM20AGGSMA_01/XY41",
        ch2_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY41/"
            "ApYYM20AGGSMA_XY41_CH2.tif"
        ),
    ),
    Example(
        category="medium",
        image_id="plate 1/ApYYM20AGGSMA_01/XY08",
        ch2_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY08/"
            "ApYYM20AGGSMA_XY08_CH2.tif"
        ),
    ),
    Example(
        category="medium",
        image_id="plate 1/ApYYM20AGGSMA_01/XY95",
        ch2_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY95/"
            "ApYYM20AGGSMA_XY95_CH2.tif"
        ),
    ),
    Example(
        category="high",
        image_id="plate 1/ApYYM20AGGSMA_02/XY04",
        ch2_path=Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_02/XY04/"
            "ApYYM20AGGSMA_XY04_CH2.tif"
        ),
    ),
    Example(
        category="high",
        image_id="plate 2/APIPIKEALDSMA/XY73",
        ch2_path=Path(
            "data/aSMA_DAPI_plates/plate 2/APIPIKEALDSMA/XY73/"
            "APIPIKEALDSMA_XY73_CH2.tif"
        ),
    ),
)


@dataclass(frozen=True)
class Measurement:
    example: Example
    ch2: np.ndarray
    retained_mask: np.ndarray
    retained_fraction: float
    excluded_fraction: float
    raw_sum: float
    retained_raw_sum: float
    retained_raw_sum_fraction: float
    background_subtracted_sum: float
    mean: float
    median: float
    p95: float
    p99: float


@app.command()
def main(
    output_dir: Path = typer.Option(
        Path("output/ch2_threshold_retained_area/threshold_8100_low_medium_high"),
        "--output-dir",
        help="Directory for threshold panel and summary CSV.",
    ),
    threshold: float = typer.Option(THRESHOLD, "--threshold", help="Raw CH2 threshold to retain."),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    measurements = [_measure(example, threshold=threshold) for example in EXAMPLES]

    panel_path = output_dir / "ch2_threshold_8100_retained_area_low_medium_high_panel.png"
    summary_path = output_dir / "ch2_threshold_8100_retained_area_low_medium_high_summary.csv"
    _write_panel(measurements, threshold=threshold, output_path=panel_path)
    _write_summary(measurements, threshold=threshold, output_path=summary_path)

    typer.echo(f"panel={panel_path}")
    typer.echo(f"summary={summary_path}")


def _measure(example: Example, *, threshold: float) -> Measurement:
    ch2 = _read_ch2(example.ch2_path)
    values = ch2.astype(np.float64)
    retained = values >= threshold
    raw_sum = float(np.sum(values, dtype=np.float64))
    retained_sum = float(np.sum(values[retained], dtype=np.float64))
    return Measurement(
        example=example,
        ch2=ch2,
        retained_mask=retained,
        retained_fraction=float(np.mean(retained)),
        excluded_fraction=float(1.0 - np.mean(retained)),
        raw_sum=raw_sum,
        retained_raw_sum=retained_sum,
        retained_raw_sum_fraction=retained_sum / raw_sum if raw_sum > 0 else float("nan"),
        background_subtracted_sum=float(np.sum(np.maximum(values - threshold, 0.0), dtype=np.float64)),
        mean=float(np.mean(values)),
        median=float(np.percentile(values, 50)),
        p95=float(np.percentile(values, 95)),
        p99=float(np.percentile(values, 99)),
    )


def _read_ch2(path: Path) -> np.ndarray:
    image = np.squeeze(np.asarray(tifffile.imread(path)))
    if image.ndim == 2:
        return image
    if image.ndim == 3 and image.shape[-1] >= 1:
        return image[..., 0]
    raise ValueError(f"Expected grayscale or RGB CH2 TIFF for {path}; got {image.shape}")


def _write_panel(
    measurements: list[Measurement],
    *,
    threshold: float,
    output_path: Path,
) -> None:
    ch2_lo, ch2_hi = _display_limits([measurement.ch2 for measurement in measurements])
    x_max = 65535.0
    y_max = _histogram_ymax(measurements, x_max=x_max)

    fig, axes = plt.subplots(
        4,
        len(measurements),
        figsize=(3.65 * len(measurements), 13.2),
        dpi=180,
        constrained_layout=False,
    )
    fig.subplots_adjust(top=0.88, bottom=0.06, left=0.045, right=0.99, hspace=0.36, wspace=0.22)
    fig.suptitle(
        "CH2/aSMA pixels retained after excluding raw intensity < 8,100",
        fontsize=17,
        fontweight="bold",
        y=0.975,
    )
    fig.text(
        0.5,
        0.935,
        (
            "Rows: original CH2, green retained overlay, binary retained mask, raw histogram. "
            "Keep rule: CH2 >= 8,100. All images use the same display and histogram axes."
        ),
        ha="center",
        fontsize=10,
    )

    for col, measurement in enumerate(measurements):
        title = (
            f"{measurement.example.category.upper()}\n"
            f"{measurement.example.image_id}\n"
            f"retained area {measurement.retained_fraction:.1%}"
        )
        _show_ch2(axes[0, col], measurement.ch2, lo=ch2_lo, hi=ch2_hi)
        axes[0, col].set_title(title, fontsize=8.7)
        axes[0, col].text(
            0.02,
            0.98,
            (
                f"mean {measurement.mean:,.0f}\n"
                f"median {measurement.median:,.0f}\n"
                f"p95 {measurement.p95:,.0f}"
            ),
            transform=axes[0, col].transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color="white",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "black", "alpha": 0.62, "edgecolor": "none"},
        )

        _show_retained_overlay(
            axes[1, col],
            measurement.ch2,
            measurement.retained_mask,
            lo=ch2_lo,
            hi=ch2_hi,
        )
        axes[1, col].text(
            0.02,
            0.98,
            "green = kept",
            transform=axes[1, col].transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color="white",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "black", "alpha": 0.62, "edgecolor": "none"},
        )

        axes[2, col].imshow(measurement.retained_mask, cmap="gray", vmin=0, vmax=1)
        axes[2, col].axis("off")
        axes[2, col].text(
            0.02,
            0.98,
            (
                f"kept pixels {measurement.retained_fraction:.1%}\n"
                f"kept raw sum {measurement.retained_raw_sum_fraction:.1%}\n"
                f"bg-sub sum {_sci(measurement.background_subtracted_sum)}"
            ),
            transform=axes[2, col].transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            color="white",
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "black", "alpha": 0.62, "edgecolor": "none"},
        )

        _plot_histogram(axes[3, col], measurement, threshold=threshold, x_max=x_max, y_max=y_max)

    fig.text(
        0.5,
        0.02,
        (
            "This is a threshold-QC visualization. It shows what remains after an empirical "
            "background floor, not a validated cell mask."
        ),
        ha="center",
        fontsize=9,
    )
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def _display_limits(images: list[np.ndarray]) -> tuple[float, float]:
    values = np.concatenate([image.ravel().astype(np.float64) for image in images])
    lo, hi = np.percentile(values[np.isfinite(values)], [DISPLAY_LOW_PERCENTILE, DISPLAY_HIGH_PERCENTILE])
    if hi <= lo:
        hi = max(float(np.nanmax(values)), lo + 1.0)
    return float(lo), float(hi)


def _show_ch2(ax: plt.Axes, image: np.ndarray, *, lo: float, hi: float) -> None:
    scaled = np.clip((image.astype(np.float64) - lo) / (hi - lo), 0, 1)
    rgb = np.zeros((*scaled.shape, 3), dtype=np.float64)
    rgb[..., 0] = scaled
    ax.imshow(rgb)
    ax.axis("off")


def _show_retained_overlay(
    ax: plt.Axes,
    image: np.ndarray,
    retained: np.ndarray,
    *,
    lo: float,
    hi: float,
) -> None:
    base = np.clip((image.astype(np.float64) - lo) / (hi - lo), 0, 1)
    rgb = np.zeros((*base.shape, 3), dtype=np.float64)
    rgb[..., 0] = base * 0.45
    rgb[..., 1] = base * 0.45
    rgb[..., 2] = base * 0.45
    rgb[retained, 0] = np.maximum(rgb[retained, 0], 0.05)
    rgb[retained, 1] = 1.0
    rgb[retained, 2] = 0.05
    ax.imshow(rgb)
    ax.axis("off")


def _histogram_ymax(measurements: list[Measurement], *, x_max: float) -> float:
    max_count = 1
    for measurement in measurements:
        counts, _edges = np.histogram(
            measurement.ch2.ravel(),
            bins=HISTOGRAM_BINS,
            range=(0, x_max),
        )
        max_count = max(max_count, int(counts.max()))
    return max_count * 1.08


def _plot_histogram(
    ax: plt.Axes,
    measurement: Measurement,
    *,
    threshold: float,
    x_max: float,
    y_max: float,
) -> None:
    ax.hist(
        measurement.ch2.ravel(),
        bins=HISTOGRAM_BINS,
        range=(0, x_max),
        color="#7f1d1d",
        alpha=0.86,
    )
    ax.axvline(threshold, color="#22c55e", linewidth=1.5, label="8,100 cutoff")
    ax.fill_betweenx([0, y_max], 0, threshold, color="#64748b", alpha=0.15)
    ax.set_xlim(0, x_max)
    ax.set_ylim(0, y_max)
    ax.set_xlabel("raw CH2")
    ax.set_ylabel("pixels")
    ax.grid(axis="y", color="#e5e7eb", linewidth=0.7)
    ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    ax.legend(fontsize=7, frameon=False, loc="upper right")


def _write_summary(
    measurements: list[Measurement],
    *,
    threshold: float,
    output_path: Path,
) -> None:
    rows = []
    for measurement in measurements:
        rows.append(
            {
                "category": measurement.example.category,
                "image_id": measurement.example.image_id,
                "ch2_path": str(measurement.example.ch2_path),
                "threshold_raw_ch2": threshold,
                "retained_rule": "CH2 >= threshold",
                "retained_area_fraction": measurement.retained_fraction,
                "excluded_area_fraction": measurement.excluded_fraction,
                "raw_sum": measurement.raw_sum,
                "retained_raw_sum": measurement.retained_raw_sum,
                "retained_raw_sum_fraction": measurement.retained_raw_sum_fraction,
                "background_subtracted_sum": measurement.background_subtracted_sum,
                "mean": measurement.mean,
                "median": measurement.median,
                "p95": measurement.p95,
                "p99": measurement.p99,
            }
        )
    pd.DataFrame(rows).to_csv(output_path, index=False)


def _sci(value: float) -> str:
    mantissa, exponent = f"{value:.2e}".split("e")
    return f"{mantissa}e{int(exponent)}"


if __name__ == "__main__":
    app()
