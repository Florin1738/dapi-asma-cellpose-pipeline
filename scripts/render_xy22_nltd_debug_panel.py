#!/usr/bin/env python
from __future__ import annotations

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
THETA_ZOOM_MAX_DEGREES = 5.0
PHI_MAX_DEGREES = 90.0


app = typer.Typer(
    help=(
        "Render a Plate 1 XY22 diagnostic panel comparing CH2, CH4, and a "
        "white-background overlay prepared for NLTD-style optical-density analysis."
    )
)


@app.command()
def main(
    ch2_path: Path = typer.Option(
        Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY22/"
            "ApYYM20AGGSMA_XY22_CH2.tif"
        ),
        "--ch2",
        help="Plate 1 XY22 CH2/aSMA TIFF.",
    ),
    ch4_path: Path = typer.Option(
        Path(
            "data/aSMA_DAPI_plates/plate 1/ApYYM20AGGSMA_01/XY22/"
            "ApYYM20AGGSMA_XY22_CH4.tif"
        ),
        "--ch4",
        help="Plate 1 XY22 CH4/DAPI TIFF.",
    ),
    output_dir: Path = typer.Option(
        Path("output/od_spherical/plate1_xy22_fixed_review"),
        "--output-dir",
        help="Directory for the corrected diagnostic panel and summary CSV.",
    ),
    eps: float = typer.Option(1e-6, "--eps", min=1e-12),
    bins: int = typer.Option(256, "--bins", min=32, max=512),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    ch2_raw = _read_rgb_channel(ch2_path, 0)
    ch4_raw = _read_rgb_channel(ch4_path, 2)
    _assert_same_shape(ch2_raw, ch4_raw)

    ch2_fraction = _to_dtype_fraction(ch2_raw)
    ch4_fraction = _to_dtype_fraction(ch4_raw)

    columns = [
        _PanelColumn(
            name="CH2 / aSMA only",
            source_paths=[ch2_path],
            source_rgb=_display_rgb(red=ch2_raw),
            fixed_rgb=_fixed_od_rgb(red_fraction=ch2_fraction),
            fluorescence_density=_fluorescence_density(red=ch2_fraction, eps=eps),
        ),
        _PanelColumn(
            name="CH4 / DAPI only",
            source_paths=[ch4_path],
            source_rgb=_display_rgb(blue=ch4_raw),
            fixed_rgb=_fixed_od_rgb(blue_fraction=ch4_fraction),
            fluorescence_density=_fluorescence_density(blue=ch4_fraction, eps=eps),
        ),
        _PanelColumn(
            name="Fixed overlay: CH2 + CH4",
            source_paths=[ch2_path, ch4_path],
            source_rgb=_display_rgb(red=ch2_raw, blue=ch4_raw),
            fixed_rgb=_fixed_od_rgb(red_fraction=ch2_fraction, blue_fraction=ch4_fraction),
            fluorescence_density=_fluorescence_density(
                red=ch2_fraction,
                blue=ch4_fraction,
                eps=eps,
            ),
        ),
    ]
    for column in columns:
        column.od_density = _nltd_brightfield_od(column.fixed_rgb)
        column.theta, column.phi, column.rho = _spherical_coordinates(column.od_density)
        column.fluorescence_rho = _rho(column.fluorescence_density)

    panel_path = output_dir / "XY22_CH2_CH4_fixed_overlay_NLTD_diagnostic.png"
    csv_path = output_dir / "XY22_CH2_CH4_fixed_overlay_NLTD_summary.csv"
    fixed_overlay_path = output_dir / "XY22_fixed_white_background_overlay_for_OD.png"

    plt.imsave(fixed_overlay_path, columns[-1].fixed_rgb)
    _write_panel(columns=columns, output_path=panel_path, bins=bins)
    _write_summary(
        columns=columns,
        output_path=csv_path,
        eps=eps,
        bins=bins,
        ch2_raw=ch2_raw,
        ch4_raw=ch4_raw,
    )

    typer.echo(f"panel={panel_path}")
    typer.echo(f"fixed_overlay={fixed_overlay_path}")
    typer.echo(f"summary={csv_path}")


class _PanelColumn:
    def __init__(
        self,
        *,
        name: str,
        source_paths: list[Path],
        source_rgb: np.ndarray,
        fixed_rgb: np.ndarray,
        fluorescence_density: np.ndarray,
    ) -> None:
        self.name = name
        self.source_paths = source_paths
        self.source_rgb = source_rgb
        self.fixed_rgb = fixed_rgb
        self.fluorescence_density = fluorescence_density
        self.fluorescence_rho: np.ndarray | None = None
        self.od_density: np.ndarray | None = None
        self.theta: np.ndarray | None = None
        self.phi: np.ndarray | None = None
        self.rho: np.ndarray | None = None


def _read_rgb_channel(path: Path, channel_index: int) -> np.ndarray:
    image = np.squeeze(np.asarray(tifffile.imread(path)))
    if image.ndim == 2:
        return image
    if image.ndim == 3 and image.shape[-1] >= channel_index + 1:
        return image[..., channel_index]
    raise ValueError(f"Expected a grayscale or RGB TIFF for {path}; got shape {image.shape}")


def _assert_same_shape(ch2: np.ndarray, ch4: np.ndarray) -> None:
    if ch2.shape != ch4.shape:
        raise ValueError(f"CH2 and CH4 shapes differ: {ch2.shape} versus {ch4.shape}")


def _to_dtype_fraction(image: np.ndarray) -> np.ndarray:
    if np.issubdtype(image.dtype, np.integer):
        max_value = float(np.iinfo(image.dtype).max)
    else:
        max_value = max(float(np.nanmax(image)), 1.0)
    return np.clip(image.astype(np.float64) / max_value, 0.0, 1.0)


def _display_rgb(
    *,
    red: np.ndarray | None = None,
    blue: np.ndarray | None = None,
) -> np.ndarray:
    shape_source = red if red is not None else blue
    if shape_source is None:
        raise ValueError("At least one source channel is required")
    rgb = np.zeros((*shape_source.shape, 3), dtype=np.float64)
    if red is not None:
        rgb[..., 0] = _percentile_display_scale(red)
    if blue is not None:
        rgb[..., 2] = _percentile_display_scale(blue)
    return rgb


def _percentile_display_scale(
    image: np.ndarray,
    low: float = DISPLAY_LOW_PERCENTILE,
    high: float = DISPLAY_HIGH_PERCENTILE,
) -> np.ndarray:
    values = image.astype(np.float64)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.zeros_like(values)
    lo, hi = np.percentile(finite, [low, high])
    if hi <= lo:
        hi = max(float(np.nanmax(finite)), lo + 1.0)
    return np.clip((values - lo) / (hi - lo), 0.0, 1.0)


def _fixed_od_rgb(
    *,
    red_fraction: np.ndarray | None = None,
    blue_fraction: np.ndarray | None = None,
) -> np.ndarray:
    shape_source = red_fraction if red_fraction is not None else blue_fraction
    if shape_source is None:
        raise ValueError("At least one source channel is required")
    rgb = np.ones((*shape_source.shape, 3), dtype=np.float64)
    if red_fraction is not None:
        rgb[..., 0] = 1.0 - red_fraction
    if blue_fraction is not None:
        rgb[..., 2] = 1.0 - blue_fraction
    return np.round(np.clip(rgb, 0.0, 1.0) * 255.0).astype(np.uint8)


def _fluorescence_density(
    *,
    red: np.ndarray | None = None,
    blue: np.ndarray | None = None,
    eps: float,
) -> np.ndarray:
    shape_source = red if red is not None else blue
    if shape_source is None:
        raise ValueError("At least one source channel is required")
    density = np.zeros((*shape_source.shape, 3), dtype=np.float64)
    if red is not None:
        density[..., 0] = -np.log(np.clip(1.0 - red, eps, 1.0))
    if blue is not None:
        density[..., 2] = -np.log(np.clip(1.0 - blue, eps, 1.0))
    return density


def _nltd_brightfield_od(rgb8: np.ndarray) -> np.ndarray:
    od = -np.log((rgb8.astype(np.float64) + 1.0) / 256.0)
    od[np.abs(od) < OD_ZERO_EPS] = 0.0
    return od


def _spherical_coordinates(density: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    red = density[..., 0]
    green = density[..., 1]
    blue = density[..., 2]
    rho = _rho(density)
    valid = _valid_od_mask(rho)
    theta = np.full_like(rho, np.nan, dtype=np.float64)
    phi = np.full_like(rho, np.nan, dtype=np.float64)
    theta[valid] = np.degrees(np.arctan2(green[valid], red[valid]))
    theta = np.where(theta < 0, theta + 360.0, theta)
    with np.errstate(divide="ignore", invalid="ignore"):
        phi[valid] = np.degrees(np.arccos(np.clip(blue[valid] / rho[valid], -1.0, 1.0)))
    return theta, phi, rho


def _rho(density: np.ndarray) -> np.ndarray:
    return np.sqrt(np.sum(density * density, axis=-1))


def _valid_od_mask(rho: np.ndarray) -> np.ndarray:
    return np.isfinite(rho) & (rho > OD_ZERO_EPS)


def _valid_angle_mask(theta: np.ndarray, phi: np.ndarray, rho: np.ndarray) -> np.ndarray:
    return _valid_od_mask(rho) & np.isfinite(theta) & np.isfinite(phi)


def _write_panel(*, columns: list[_PanelColumn], output_path: Path, bins: int) -> None:
    fig, axes = plt.subplots(
        5,
        len(columns),
        figsize=(14, 18),
        dpi=180,
        constrained_layout=False,
    )
    fig.subplots_adjust(top=0.9, bottom=0.045, left=0.06, right=0.96, hspace=0.62, wspace=0.34)
    fig.suptitle(
        "Plate 1 XY22: CH2, CH4, and corrected NLTD-style OD diagnostic",
        fontsize=16,
        fontweight="bold",
        y=0.985,
    )
    fig.text(
        0.5,
        0.953,
        (
            "Source fluorescence is red+blue only. The fixed overlay maps signal to "
            "absorbance on a white background before applying NLTD's OD transform."
        ),
        ha="center",
        va="top",
        fontsize=10,
    )

    max_fluor_rho = max(float(np.nanpercentile(col.fluorescence_rho, 99.9)) for col in columns)
    max_od_rho = max(_valid_percentile(col.rho, _valid_od_mask(col.rho), 99.9) for col in columns)
    max_fluor_rho = max(max_fluor_rho, 0.01)
    max_od_rho = max(max_od_rho, 0.01)

    for index, column in enumerate(columns):
        axes[0, index].imshow(column.source_rgb)
        axes[0, index].set_title(f"{column.name}\nsource display")
        axes[0, index].axis("off")

        axes[1, index].imshow(column.fixed_rgb)
        axes[1, index].set_title("fixed input for OD\nwhite background")
        axes[1, index].axis("off")

        _plot_rho_histogram(
            axes[2, index],
            column.fluorescence_rho,
            x_max=max_fluor_rho,
            title="1D rho, fluorescence density",
        )
        _plot_theta_phi(
            axes[3, index],
            column.theta,
            column.phi,
            column.rho,
            bins=bins,
            title="NLTD theta-phi\nzoomed theta 0-5 deg",
        )
        _plot_phi_rho(
            axes[4, index],
            column.phi,
            column.rho,
            bins=bins,
            rho_max=max_od_rho,
            title="More useful here:\nphi versus OD rho",
        )

    fig.text(
        0.5,
        0.008,
        (
            "Key check: theta is expected to collapse because the fixed OD input has no "
            "green absorbance axis. For these two-channel fluorescence images, rho and "
            "phi/rho views are interpretable; a theta-phi TPOM is mostly a line."
        ),
        ha="center",
        va="bottom",
        fontsize=9,
    )
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def _plot_rho_histogram(ax: plt.Axes, rho: np.ndarray, *, x_max: float, title: str) -> None:
    values = np.ravel(rho)
    ax.hist(values, bins=180, color="#475569", alpha=0.9)
    ax.set_title(title)
    ax.set_xlabel("rho")
    ax.set_ylabel("pixels")
    ax.set_yscale("log")
    ax.set_xlim(0, x_max)


def _plot_theta_phi(
    ax: plt.Axes,
    theta: np.ndarray,
    phi: np.ndarray,
    rho: np.ndarray,
    *,
    bins: int,
    title: str,
) -> None:
    theta_max = THETA_ZOOM_MAX_DEGREES
    valid = _valid_angle_mask(theta, phi, rho)
    theta_values = np.ravel(theta[valid])
    phi_values = np.ravel(phi[valid])
    if theta_values.size == 0:
        ax.text(0.5, 0.5, "no nonzero OD pixels", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        ax.set_xlabel("theta degrees")
        ax.set_ylabel("phi degrees")
        ax.set_xlim(0.0, theta_max)
        ax.set_ylim(0.0, PHI_MAX_DEGREES)
        return
    hist, _x_edges, _y_edges = np.histogram2d(
        theta_values,
        phi_values,
        bins=[bins, bins],
        range=[[0.0, theta_max], [0.0, PHI_MAX_DEGREES]],
    )
    if np.nanpercentile(theta_values, 99) - np.nanpercentile(theta_values, 1) < 1e-9:
        counts, edges = np.histogram(phi_values, bins=bins, range=(0.0, PHI_MAX_DEGREES))
        centers = 0.5 * (edges[:-1] + edges[1:])
        if counts.max() > 0:
            widths = theta_max * 0.9 * np.log10(counts + 1.0) / np.log10(counts.max() + 1.0)
            ax.fill_betweenx(centers, 0.0, widths, step="mid", color="#f97316", alpha=0.9)
        ax.text(
            0.97,
            0.04,
            "theta collapsed;\nwidth shows phi counts",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            color="#111827",
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
        )
    else:
        image = ax.imshow(
            np.log10(hist.T + 1.0),
            origin="lower",
            extent=[0.0, theta_max, 0.0, PHI_MAX_DEGREES],
            aspect="auto",
            cmap="magma",
        )
        plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="log10(count + 1)")
    ax.axvline(0.0, color="cyan", linewidth=0.8)
    ax.set_title(title)
    ax.set_xlabel("theta degrees")
    ax.set_ylabel("phi degrees")
    ax.set_xlim(0.0, theta_max)
    ax.set_ylim(0.0, PHI_MAX_DEGREES)


def _plot_phi_rho(
    ax: plt.Axes,
    phi: np.ndarray,
    rho: np.ndarray,
    *,
    bins: int,
    rho_max: float,
    title: str,
) -> None:
    valid = _valid_angle_mask(np.zeros_like(phi), phi, rho)
    phi_values = np.ravel(phi[valid])
    rho_values = np.ravel(rho[valid])
    if phi_values.size == 0:
        ax.text(0.5, 0.5, "no nonzero OD pixels", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        ax.set_xlabel("phi degrees")
        ax.set_ylabel("OD rho")
        ax.set_xlim(0.0, PHI_MAX_DEGREES)
        ax.set_ylim(0.0, rho_max)
        return
    hist, _x_edges, _y_edges = np.histogram2d(
        phi_values,
        rho_values,
        bins=[bins, bins],
        range=[[0.0, PHI_MAX_DEGREES], [0.0, rho_max]],
    )
    image = ax.imshow(
        np.log10(hist.T + 1.0),
        origin="lower",
        extent=[0.0, PHI_MAX_DEGREES, 0.0, rho_max],
        aspect="auto",
        cmap="viridis",
    )
    ax.set_title(title)
    ax.set_xlabel("phi degrees")
    ax.set_ylabel("OD rho")
    ax.set_xlim(0.0, PHI_MAX_DEGREES)
    ax.set_ylim(0.0, rho_max)
    plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="log10(count + 1)")


def _write_summary(
    *,
    columns: list[_PanelColumn],
    output_path: Path,
    eps: float,
    bins: int,
    ch2_raw: np.ndarray,
    ch4_raw: np.ndarray,
) -> None:
    rows = []
    for column in columns:
        valid = _valid_angle_mask(column.theta, column.phi, column.rho)
        tp_hist, _x_edges, _y_edges = np.histogram2d(
            np.ravel(column.theta[valid]),
            np.ravel(column.phi[valid]),
            bins=[bins, bins],
            range=[[0.0, THETA_ZOOM_MAX_DEGREES], [0.0, PHI_MAX_DEGREES]],
        )
        rows.append(
            {
                "column": column.name,
                "source_paths": ";".join(str(path) for path in column.source_paths),
                "eps": eps,
                "bins": bins,
                "display_low_percentile": DISPLAY_LOW_PERCENTILE,
                "display_high_percentile": DISPLAY_HIGH_PERCENTILE,
                "od_zero_eps": OD_ZERO_EPS,
                "theta_zoom_max_degrees": THETA_ZOOM_MAX_DEGREES,
                "phi_max_degrees": PHI_MAX_DEGREES,
                "source_red_nonzero": int(np.count_nonzero(column.source_rgb[..., 0])),
                "source_green_nonzero": int(np.count_nonzero(column.source_rgb[..., 1])),
                "source_blue_nonzero": int(np.count_nonzero(column.source_rgb[..., 2])),
                "fluorescence_rho_median": float(np.nanmedian(column.fluorescence_rho)),
                "fluorescence_rho_p99": float(np.nanpercentile(column.fluorescence_rho, 99)),
                "nltd_od_valid_pixels": int(np.count_nonzero(valid)),
                "nltd_od_theta_median": _valid_percentile(column.theta, valid, 50),
                "nltd_od_theta_p99": _valid_percentile(column.theta, valid, 99),
                "nltd_od_phi_median": _valid_percentile(column.phi, valid, 50),
                "nltd_od_phi_p01": _valid_percentile(column.phi, valid, 1),
                "nltd_od_phi_p99": _valid_percentile(column.phi, valid, 99),
                "nltd_od_rho_median": _valid_percentile(column.rho, valid, 50),
                "nltd_od_rho_p99": _valid_percentile(column.rho, valid, 99),
                "theta_phi_nonzero_bins": int(np.count_nonzero(tp_hist)),
            }
        )
    rows.append(
        {
            "column": "raw_channel_stats",
            "source_paths": "raw CH2 and CH4 arrays before display scaling",
            "eps": eps,
            "bins": bins,
            "display_low_percentile": DISPLAY_LOW_PERCENTILE,
            "display_high_percentile": DISPLAY_HIGH_PERCENTILE,
            "od_zero_eps": OD_ZERO_EPS,
            "theta_zoom_max_degrees": THETA_ZOOM_MAX_DEGREES,
            "phi_max_degrees": PHI_MAX_DEGREES,
            "source_red_nonzero": int(np.count_nonzero(ch2_raw)),
            "source_green_nonzero": 0,
            "source_blue_nonzero": int(np.count_nonzero(ch4_raw)),
            "fluorescence_rho_median": float("nan"),
            "fluorescence_rho_p99": float("nan"),
            "nltd_od_valid_pixels": 0,
            "nltd_od_theta_median": float("nan"),
            "nltd_od_theta_p99": float("nan"),
            "nltd_od_phi_median": float("nan"),
            "nltd_od_phi_p01": float("nan"),
            "nltd_od_phi_p99": float("nan"),
            "nltd_od_rho_median": float("nan"),
            "nltd_od_rho_p99": float("nan"),
            "theta_phi_nonzero_bins": 0,
        }
    )
    pd.DataFrame(rows).to_csv(output_path, index=False)


def _valid_percentile(values: np.ndarray, valid: np.ndarray, percentile: float) -> float:
    selected = values[valid]
    if selected.size == 0:
        return float("nan")
    return float(np.nanpercentile(selected, percentile))


if __name__ == "__main__":
    app()
