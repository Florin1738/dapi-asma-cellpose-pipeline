from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile


def read_primary_intensity_plane(
    path: Path, *, allow_multi_active_rgb: bool = False
) -> tuple[np.ndarray, int]:
    with tifffile.TiffFile(path) as tif:
        page_count = len(tif.pages)
        image = np.asarray(tif.series[0].asarray())

    return to_intensity_plane(
        image, path=path, allow_multi_active_rgb=allow_multi_active_rgb
    ), page_count


def to_intensity_plane(
    image: np.ndarray,
    *,
    path: Path | None = None,
    allow_multi_active_rgb: bool = False,
) -> np.ndarray:
    arr = np.squeeze(np.asarray(image))
    if arr.ndim == 2:
        return arr

    if arr.ndim == 3 and arr.shape[-1] in {3, 4}:
        samples = arr[..., :3]
        sample_nonzero = np.count_nonzero(samples.reshape(-1, samples.shape[-1]), axis=0)
        active_samples = np.flatnonzero(sample_nonzero)
        if len(active_samples) == 1:
            return samples[..., int(active_samples[0])]
        if not allow_multi_active_rgb:
            source = f": {path}" if path is not None else ""
            raise ValueError(
                "Refusing to collapse multi-active RGB image into one intensity plane"
                f"{source}. Use an explicit display-only projection for overlays."
            )
        return np.max(samples, axis=-1)

    source = f": {path}" if path is not None else ""
    raise ValueError(f"Expected a 2-D image or YXS RGB image after squeezing singleton axes{source}")
