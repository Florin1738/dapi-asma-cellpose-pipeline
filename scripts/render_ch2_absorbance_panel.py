#!/usr/bin/env python
from __future__ import annotations

from dataclasses import dataclass, replace
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
OD_ZERO_EPS = 1e-12
OD_MAX_8BIT = -np.log(1.0 / 256.0)


app = typer.Typer(
    help="Render CH2/aSMA-only absorbance-style density panels for selected wells."
)


@dataclass(frozen=True)
class WellInput:
    location: str
    path: Path


@dataclass(frozen=True)
class WellMeasurement:
    location: str
    path: Path
    ch2: np.ndarray
    display: np.ndarray
    absorbance: np.ndarray
    raw_integrated_intensity: float
    absorbance_sum: float
    absorbance_mean: float
    absorbance_median: float
    absorbance_p95: float
    absorbance_p99: float
    saturated_pixel_fraction: float


@app.command()
def main(
    output_dir: Path = typer.Option(
        Path("output/ch2_absorbance/plate1_xy22_xy23_xy24"),
        "--output-dir",
        help="Directory for the panel and summary CSV.",
    ),
    wells: list[str] = typer.Option(
        ["XY22", "XY23", "XY24"],
        "--well",
        help="Well/location identifiers to include.",
    ),
    plate_dir: Path = typer.Option(
        Path("data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01"),
        "--plate-dir",
        help="Plate directory containing XY folders.",
    ),
    bins: int = typer.Option(256, "--bins", min=32, max=512),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    inputs = [_resolve_well_input(plate_dir, well) for well in wells]
    measurements = _apply_shared_display_scale([_measure_well(input_) for input_ in inputs])

    panel_path = output_dir / "plate1_XY22_XY23_XY24_CH2_absorbance_panel.png"
    csv_path = output_dir / "plate1_XY22_XY23_XY24_CH2_absorbance_summary.csv"
    _write_panel(measurements=measurements, output_path=panel_path, bins=bins)
    _write_summary(measurements=measurements, output_path=csv_path, bins=bins)

    typer.echo(f"panel={panel_path}")
    typer.echo(f"summary={csv_path}")


def _resolve_well_input(plate_dir: Path, well: str) -> WellInput:
    matches = sorted((plate_dir / well).glob(f"*_{well}_CH2.tif"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one CH2 TIFF for {well}, found {len(matches)}")
    return WellInput(location=well, path=matches[0])


def _measure_well(input_: WellInput) -> WellMeasurement:
    ch2 = _read_ch2(input_.path)
    fraction = _to_dtype_fraction(ch2)
    display = _red_display(ch2)
    absorbance = _fluorescence_to_od_compatible_density(fraction)
    max_value = _dtype_max(ch2)
    return WellMeasurement(
        location=input_.location,
        path=input_.path,
        ch2=ch2,
        display=display,
        absorbance=absorbance,
        raw_integrated_intensity=float(np.sum(ch2, dtype=np.float64)),
        absorbance_sum=float(np.sum(absorbance, dtype=np.float64)),
        absorbance_mean=float(np.mean(absorbance)),
        absorbance_median=float(np.median(absorbance)),
        absorbance_p95=float(np.percentile(absorbance, 95)),
        absorbance_p99=float(np.percentile(absorbance, 99)),
        saturated_pixel_fraction=float(np.mean(ch2 >= max_value)),
    )


def _read_ch2(path: Path) -> np.ndarray:
    image = np.squeeze(np.asarray(tifffile.imread(path)))
    if image.ndim == 2:
        return image
    if image.ndim == 3 and image.shape[-1] >= 1:
        return image[..., 0]
    raise ValueError(f"Expected grayscale or RGB CH2 TIFF for {path}; got shape {image.shape}")


def _dtype_max(image: np.ndarray) -> float:
    if np.issubdtype(image.dtype, np.integer):
        return float(np.iinfo(image.dtype).max)
    return max(float(np.nanmax(image)), 1.0)


def _to_dtype_fraction(image: np.ndarray) -> np.ndarray:
    return np.clip(image.astype(np.float64) / _dtype_max(image), 0.0, 1.0)


def _red_display(image: np.ndarray) -> np.ndarray:
    scaled = _percentile_display_scale(image)
    return _red_rgb(scaled)


def _red_rgb(scaled: np.ndarray) -> np.ndarray:
    rgb = np.zeros((*scaled.shape, 3), dtype=np.float64)
    rgb[..., 0] = scaled
    return rgb


def _apply_shared_display_scale(measurements: list[WellMeasurement]) -> list[WellMeasurement]:
    values = np.concatenate([measurement.ch2.ravel().astype(np.float64) for measurement in measurements])
    finite = values[np.isfinite(values)]
    lo, hi = np.percentile(finite, [DISPLAY_LOW_PERCENTILE, DISPLAY_HIGH_PERCENTILE])
    if hi <= lo:
        hi = max(float(np.nanmax(finite)), lo + 1.0)
    scaled_measurements = []
    for measurement in measurements:
        scaled = np.clip((measurement.ch2.astype(np.float64) - lo) / (hi - lo), 0.0, 1.0)
        scaled_measurements.append(replace(measurement, display=_red_rgb(scaled)))
    return scaled_measurements


def _percentile_display_scale(
    image: np.ndarray,
    low: float = DISPLAY_LOW_PERCENTILE,
    high: float = DISPLAY_HIGH_PERCENTILE,
) -> np.ndarray:
    values = image.astype(np.float64)
    lo, hi = np.percentile(values[np.isfinite(values)], [low, high])
    if hi <= lo:
        hi = max(float(np.nanmax(values)), lo + 1.0)
    return np.clip((values - lo) / (hi - lo), 0.0, 1.0)


def _fluorescence_to_od_compatible_density(fraction: np.ndarray) -> np.ndarray:
    fixed_red_8bit = np.round(np.clip(1.0 - fraction, 0.0, 1.0) * 255.0)
    density = -np.log((fixed_red_8bit + 1.0) / 256.0)
    density[np.abs(density) < OD_ZERO_EPS] = 0.0
    return density


def _write_panel(
    *,
    measurements: list[WellMeasurement],
    output_path: Path,
    bins: int,
) -> None:
    x_max = OD_MAX_8BIT
    y_max = _shared_histogram_ymax(measurements, bins=bins, x_max=x_max)
    fig, axes = plt.subplots(
        2,
        len(measurements),
        figsize=(5.1 * len(measurements), 7.8),
        dpi=180,
        constrained_layout=False,
    )
    if len(measurements) == 1:
        axes = axes.reshape(2, 1)
    fig.subplots_adjust(top=0.83, bottom=0.12, left=0.07, right=0.98, hspace=0.36, wspace=0.24)
    fig.suptitle(
        "Plate 1 CH2 / alpha-smooth muscle actin: absorbance-style density",
        fontsize=15,
        fontweight="bold",
        y=0.965,
    )
    fig.text(
        0.5,
        0.915,
        (
            "All pixels included. Higher density means brighter CH2/aSMA after the "
            "white-background OD-compatible transform. CH2 thumbnails use one shared display scale."
        ),
        ha="center",
        fontsize=9.5,
    )

    for column, measurement in enumerate(measurements):
        image_ax = axes[0, column]
        image_ax.imshow(measurement.display)
        image_ax.set_title(measurement.location, fontsize=13, fontweight="bold")
        image_ax.axis("off")
        image_ax.text(
            0.02,
            0.98,
            "\n".join(
                [
                    f"raw sum: {_sci(measurement.raw_integrated_intensity)}",
                    f"density sum: {_sci(measurement.absorbance_sum)}",
                    f"mean density: {measurement.absorbance_mean:.3f}",
                    f"sat px: {measurement.saturated_pixel_fraction:.1%}",
                ]
            ),
            transform=image_ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color="white",
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "black", "alpha": 0.62, "edgecolor": "none"},
        )

        hist_ax = axes[1, column]
        _plot_histogram(hist_ax, measurement, bins=bins, x_max=x_max, y_max=y_max)

    fig.text(
        0.5,
        0.03,
        (
            "Read this as a CH2 intensity distribution, not cell segmentation. "
            "The spike at the far right reflects saturated/near-saturated CH2 pixels."
        ),
        ha="center",
        fontsize=9,
    )
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def _shared_histogram_ymax(
    measurements: list[WellMeasurement],
    *,
    bins: int,
    x_max: float,
) -> float:
    max_count = 1
    for measurement in measurements:
        counts, _edges = np.histogram(measurement.absorbance.ravel(), bins=bins, range=(0.0, x_max))
        max_count = max(max_count, int(counts.max()))
    return 10 ** np.ceil(np.log10(max_count))


def _plot_histogram(
    ax: plt.Axes,
    measurement: WellMeasurement,
    *,
    bins: int,
    x_max: float,
    y_max: float,
) -> None:
    ax.hist(
        measurement.absorbance.ravel(),
        bins=bins,
        range=(0.0, x_max),
        color="#7f1d1d",
        alpha=0.88,
    )
    ax.axvline(measurement.absorbance_median, color="#fbbf24", linewidth=1.4, label="median")
    ax.axvline(measurement.absorbance_p95, color="#22c55e", linewidth=1.4, label="p95")
    ax.axvline(measurement.absorbance_p99, color="#38bdf8", linewidth=1.4, label="p99")
    ax.set_title(
        (
            f"median {measurement.absorbance_median:.2f} | "
            f"p95 {measurement.absorbance_p95:.2f} | p99 {measurement.absorbance_p99:.2f}"
        ),
        fontsize=10,
    )
    ax.set_xlabel("CH2 density / absorbance-style value")
    ax.set_ylabel("pixels")
    ax.set_xlim(0.0, x_max)
    ax.set_ylim(1.0, y_max)
    ax.set_yscale("log")
    ax.grid(axis="y", color="#e5e7eb", linewidth=0.7)
    ax.legend(loc="upper right", fontsize=8, frameon=False)


def _write_summary(
    *,
    measurements: list[WellMeasurement],
    output_path: Path,
    bins: int,
) -> None:
    rows = []
    for measurement in measurements:
        rows.append(
            {
                "location": measurement.location,
                "path": str(measurement.path),
                "bins": bins,
                "display_low_percentile": DISPLAY_LOW_PERCENTILE,
                "display_high_percentile": DISPLAY_HIGH_PERCENTILE,
                "od_zero_eps": OD_ZERO_EPS,
                "density_formula": "-ln((round((1 - CH2_fraction) * 255) + 1) / 256)",
                "raw_integrated_intensity": measurement.raw_integrated_intensity,
                "absorbance_sum": measurement.absorbance_sum,
                "absorbance_mean": measurement.absorbance_mean,
                "absorbance_median": measurement.absorbance_median,
                "absorbance_p95": measurement.absorbance_p95,
                "absorbance_p99": measurement.absorbance_p99,
                "saturated_pixel_fraction": measurement.saturated_pixel_fraction,
            }
        )
    pd.DataFrame(rows).to_csv(output_path, index=False)


def _sci(value: float) -> str:
    mantissa, exponent = f"{value:.2e}".split("e")
    return f"{mantissa} x 10^{int(exponent)}"


if __name__ == "__main__":
    app()
