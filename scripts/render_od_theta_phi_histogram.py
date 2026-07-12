#!/usr/bin/env python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import tifffile
import typer


app = typer.Typer(
    help=(
        "Render NLTD-style theta-phi histograms from RGB optical-density or "
        "fluorescence-adapted pseudo-OD space."
    )
)


@dataclass(frozen=True)
class RhoRange:
    lower: float
    upper: float | None

    @property
    def label(self) -> str:
        if self.upper is None:
            return f"rho >= {self.lower:g}"
        return f"{self.lower:g} <= rho < {self.upper:g}"


DEFAULT_RHO_RANGES = (
    RhoRange(0.0, 0.5),
    RhoRange(0.5, 1.0),
    RhoRange(1.0, 2.0),
    RhoRange(2.0, 4.0),
    RhoRange(4.0, 8.0),
    RhoRange(8.0, None),
)


@app.command()
def main(
    input_path: Path = typer.Option(
        Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY22/"
            "ApYYM20AGGSMA_XY22_Overlay.tif"
        ),
        "--input",
        help="RGB image path. For this project, use the microscope Overlay TIFF.",
    ),
    output_dir: Path = typer.Option(
        Path("output/od_spherical/plate1_xy22"),
        "--output-dir",
        help="Directory for panel and summary CSV.",
    ),
    title: str = typer.Option(
        "Plate 1 XY22 overlay: theta-phi histograms in fluorescence pseudo-OD space",
        "--title",
        help="Figure title.",
    ),
    mode: str = typer.Option(
        "fluorescence-pseudo-od",
        "--mode",
        help="'fluorescence-pseudo-od' maps bright fluorescence to high rho; "
        "'brightfield-od' uses standard OD = -ln(I/I0).",
    ),
    bins: int = typer.Option(256, "--bins", min=16, max=1024),
    theta_max: float = typer.Option(
        90.0,
        "--theta-max",
        min=0.1,
        max=360.0,
        help="Upper theta display/histogram limit in degrees. Use a small value when theta collapses.",
    ),
    eps: float = typer.Option(1e-6, "--eps", min=1e-12),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    image = _read_rgb(input_path)
    density = _density_image(image, mode=mode, eps=eps)
    theta, phi, rho = _spherical_coordinates(density)
    panel_path = output_dir / "XY22_overlay_theta_phi_by_rho_panel.png"
    csv_path = output_dir / "XY22_overlay_theta_phi_by_rho_summary.csv"
    _write_panel(
        image=image,
        theta=theta,
        phi=phi,
        rho=rho,
        rho_ranges=DEFAULT_RHO_RANGES,
        bins=bins,
        theta_max=theta_max,
        title=title,
        mode=mode,
        output_path=panel_path,
    )
    _write_summary(
        path=csv_path,
        input_path=input_path,
        mode=mode,
        bins=bins,
        theta_max=theta_max,
        rho_ranges=DEFAULT_RHO_RANGES,
        theta=theta,
        phi=phi,
        rho=rho,
    )
    typer.echo(f"panel={panel_path}")
    typer.echo(f"summary={csv_path}")


def _read_rgb(path: Path) -> np.ndarray:
    arr = np.asarray(tifffile.imread(path))
    arr = np.squeeze(arr)
    if arr.ndim != 3 or arr.shape[-1] < 3:
        raise ValueError(f"Expected RGB/YXS image with at least 3 samples: {path}")
    return arr[..., :3]


def _density_image(image: np.ndarray, *, mode: str, eps: float) -> np.ndarray:
    arr = image.astype(np.float64)
    if np.issubdtype(image.dtype, np.integer):
        max_value = float(np.iinfo(image.dtype).max)
    else:
        max_value = float(np.nanmax(arr))
    if max_value <= 0:
        raise ValueError("Image maximum is zero; cannot normalize RGB values")
    normalized = np.clip(arr / max_value, 0.0, 1.0)
    if mode == "fluorescence-pseudo-od":
        return -np.log(np.clip(1.0 - normalized, eps, 1.0))
    if mode == "brightfield-od":
        return -np.log(np.clip(normalized, eps, 1.0))
    raise ValueError(f"Unsupported mode: {mode}")


def _spherical_coordinates(density: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    red = density[..., 0]
    green = density[..., 1]
    blue = density[..., 2]
    rho = np.sqrt(red * red + green * green + blue * blue)
    theta = np.degrees(np.arctan2(green, red))
    theta = np.where(theta < 0, theta + 360.0, theta)
    with np.errstate(invalid="ignore", divide="ignore"):
        phi = np.degrees(np.arccos(np.clip(blue / np.maximum(rho, 1e-15), -1.0, 1.0)))
    phi = np.nan_to_num(phi, nan=0.0)
    return theta, phi, rho


def _write_panel(
    *,
    image: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    rho: np.ndarray,
    rho_ranges: tuple[RhoRange, ...],
    bins: int,
    theta_max: float,
    title: str,
    mode: str,
    output_path: Path,
) -> None:
    histograms = []
    for rho_range in rho_ranges:
        mask = _rho_mask(rho, rho_range)
        histograms.append(_theta_phi_hist(theta[mask], phi[mask], bins=bins, theta_max=theta_max))
    max_log = max(float(np.log10(hist + 1).max()) for hist in histograms)

    fig = plt.figure(figsize=(16, 12), dpi=180, constrained_layout=True)
    subfigs = fig.subfigures(3, 1, height_ratios=[0.12, 0.28, 0.60])
    subfigs[0].suptitle(title, fontsize=16, fontweight="bold")
    subfigs[0].text(
        0.5,
        0.18,
        (
            "Coordinates: theta = atan2(G_density, R_density); "
            "phi = arccos(B_density / rho); rho = ||RGB_density||. "
            "Histograms show log10(pixel count + 1)."
        ),
        ha="center",
        va="center",
        fontsize=9,
    )

    top_axes = subfigs[1].subplots(1, 3)
    _show_rgb_preview(top_axes[0], image)
    top_axes[0].set_title("XY22 microscope RGB overlay")
    _plot_rho_distribution(top_axes[1], rho, rho_ranges)
    top_axes[1].set_title("rho / absorbance distribution")
    _plot_theta_phi_all(top_axes[2], theta, phi, bins=bins, theta_max=theta_max)
    top_axes[2].set_title("all rho values")

    grid_axes = subfigs[2].subplots(2, 3, sharex=True, sharey=True)
    mappable = None
    for ax, rho_range, hist in zip(grid_axes.ravel(), rho_ranges, histograms, strict=True):
        log_hist = np.log10(hist + 1).T
        mappable = ax.imshow(
            log_hist,
            origin="lower",
            extent=[0, theta_max, 0, 90],
            aspect="auto",
            cmap="magma",
            vmin=0,
            vmax=max_log,
        )
        count = int(_rho_mask(rho, rho_range).sum())
        fraction = count / rho.size
        ax.set_title(f"{rho_range.label}\n{count:,} px ({fraction:.1%})", fontsize=9)
        ax.set_xlabel("theta degrees")
        ax.set_ylabel("phi degrees")
        ax.set_xlim(0, theta_max)
        ax.set_ylim(0, 90)
    if mappable is not None:
        subfigs[2].colorbar(mappable, ax=grid_axes.ravel().tolist(), label="log10(pixel count + 1)")

    fig.text(
        0.5,
        0.01,
        (
            f"Mode: {mode}. This is an exploratory NLTD-style visualization. "
            "Because this fluorescence overlay has no green signal, theta is expected to collapse near 0; "
            "phi carries most red-vs-blue separation."
        ),
        ha="center",
        va="bottom",
        fontsize=9,
    )
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def _show_rgb_preview(ax: plt.Axes, image: np.ndarray) -> None:
    arr = image.astype(np.float64)
    if np.issubdtype(image.dtype, np.integer):
        arr /= float(np.iinfo(image.dtype).max)
    else:
        arr /= max(float(np.nanmax(arr)), 1.0)
    ax.imshow(np.clip(arr, 0, 1))
    ax.axis("off")


def _plot_rho_distribution(ax: plt.Axes, rho: np.ndarray, rho_ranges: tuple[RhoRange, ...]) -> None:
    ax.hist(rho.ravel(), bins=160, color="#4b5563", alpha=0.85)
    colors = ["#dbeafe", "#bfdbfe", "#93c5fd", "#60a5fa", "#3b82f6", "#1d4ed8"]
    ymax = ax.get_ylim()[1]
    for rho_range, color in zip(rho_ranges, colors, strict=True):
        lower = rho_range.lower
        upper = rho_range.upper if rho_range.upper is not None else float(np.nanmax(rho))
        ax.axvspan(lower, upper, color=color, alpha=0.25)
    ax.set_xlabel("rho / absorbance")
    ax.set_ylabel("pixels")
    ax.set_yscale("log")
    ax.set_ylim(1, ymax)


def _plot_theta_phi_all(
    ax: plt.Axes, theta: np.ndarray, phi: np.ndarray, *, bins: int, theta_max: float
) -> None:
    hist = _theta_phi_hist(theta.ravel(), phi.ravel(), bins=bins, theta_max=theta_max)
    ax.imshow(
        np.log10(hist + 1).T,
        origin="lower",
        extent=[0, theta_max, 0, 90],
        aspect="auto",
        cmap="magma",
    )
    ax.set_xlabel("theta degrees")
    ax.set_ylabel("phi degrees")
    ax.set_xlim(0, theta_max)
    ax.set_ylim(0, 90)


def _theta_phi_hist(
    theta: np.ndarray, phi: np.ndarray, *, bins: int, theta_max: float
) -> np.ndarray:
    theta = np.asarray(theta)
    phi = np.asarray(phi)
    in_range = (theta >= 0) & (theta <= theta_max) & (phi >= 0) & (phi <= 90)
    hist, _theta_edges, _phi_edges = np.histogram2d(
        theta[in_range],
        phi[in_range],
        bins=[bins, bins],
        range=[[0, theta_max], [0, 90]],
    )
    return hist


def _rho_mask(rho: np.ndarray, rho_range: RhoRange) -> np.ndarray:
    if rho_range.upper is None:
        return rho >= rho_range.lower
    return (rho >= rho_range.lower) & (rho < rho_range.upper)


def _write_summary(
    *,
    path: Path,
    input_path: Path,
    mode: str,
    bins: int,
    theta_max: float,
    rho_ranges: tuple[RhoRange, ...],
    theta: np.ndarray,
    phi: np.ndarray,
    rho: np.ndarray,
) -> None:
    lines = [
        "input_path,mode,bins,theta_max,rho_range,pixel_count,pixel_fraction,"
        "rho_min,rho_max,rho_median,theta_median,phi_median\n"
    ]
    for rho_range in rho_ranges:
        mask = _rho_mask(rho, rho_range)
        count = int(mask.sum())
        if count == 0:
            values = ["nan"] * 5
        else:
            values = [
                f"{float(np.nanmin(rho[mask])):.8g}",
                f"{float(np.nanmax(rho[mask])):.8g}",
                f"{float(np.nanmedian(rho[mask])):.8g}",
                f"{float(np.nanmedian(theta[mask])):.8g}",
                f"{float(np.nanmedian(phi[mask])):.8g}",
            ]
        lines.append(
            ",".join(
                [
                    str(input_path),
                    mode,
                    str(bins),
                    f"{theta_max:.8g}",
                    f'"{rho_range.label}"',
                    str(count),
                    f"{count / rho.size:.8g}",
                    *values,
                ]
            )
            + "\n"
        )
    path.write_text("".join(lines), encoding="utf-8")


if __name__ == "__main__":
    app()
