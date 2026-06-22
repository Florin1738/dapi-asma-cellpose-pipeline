from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from dapi_norm.image_arrays import read_primary_intensity_plane


def read_2d(path: Path) -> np.ndarray:
    image, _page_count = read_primary_intensity_plane(path)
    return image


def robust_scale(image: np.ndarray) -> np.ndarray:
    arr = image.astype(np.float32, copy=False)
    low, high = np.percentile(arr, [1, 99.8])
    if high <= low:
        high = float(np.max(arr))
        low = float(np.min(arr))
    if high <= low:
        return np.zeros_like(arr, dtype=np.float32)
    return np.clip((arr - low) / (high - low), 0, 1)


def compose_channel_rgb(scaled_channels: dict[str, np.ndarray]) -> np.ndarray:
    if not scaled_channels:
        raise ValueError("At least one scaled channel is required")

    first_shape = next(iter(scaled_channels.values())).shape
    composite = np.zeros((*first_shape, 3), dtype=np.float32)

    for channel, image in scaled_channels.items():
        channel_upper = channel.upper()
        if image.shape != first_shape:
            raise ValueError("All channels must have the same shape")
        if channel_upper == "CH2":
            composite[..., 0] = np.maximum(composite[..., 0], image)
        elif channel_upper == "CH4":
            composite[..., 2] = np.maximum(composite[..., 2], image)
        else:
            composite[..., 1] = np.maximum(composite[..., 1], image)

    return np.clip(composite, 0, 1)


def generate_position_preview(
    *,
    position_id: str,
    channel_paths: dict[str, Path],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ordered_channels = sorted(channel_paths)
    if not ordered_channels:
        raise ValueError("At least one channel path is required")

    images = {channel: read_2d(path) for channel, path in channel_paths.items()}
    scaled = {channel: robust_scale(image) for channel, image in images.items()}

    n_panels = len(ordered_channels) + 1
    fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, 4), constrained_layout=True)
    if n_panels == 1:
        axes = [axes]

    for ax, channel in zip(axes, ordered_channels, strict=False):
        ax.imshow(scaled[channel], cmap="gray")
        raw = images[channel]
        ax.set_title(
            f"{position_id} {channel}\n"
            f"min={raw.min():.0f} p50={np.percentile(raw, 50):.0f} "
            f"p99={np.percentile(raw, 99):.0f} max={raw.max():.0f}",
            fontsize=9,
        )
        ax.axis("off")

    composite_ax = axes[-1]
    composite = compose_channel_rgb(scaled)
    composite_ax.imshow(composite)
    composite_ax.set_title("Composite\nCH2 red, CH4 blue", fontsize=9)
    composite_ax.axis("off")

    fig.savefig(output_path, dpi=150)
    plt.close(fig)
