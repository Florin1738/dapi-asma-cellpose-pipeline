#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import tifffile
import typer
from skimage.segmentation import find_boundaries

from dapi_norm.image_arrays import read_primary_intensity_plane

app = typer.Typer(help="Run a one-image Cellpose smoke test and save visual QC outputs.")


@app.command()
def main(
    input_path: Path = typer.Option(..., "--input", help="Candidate DAPI TIFF image."),
    output_dir: Path = typer.Option(..., "--output", help="Smoke-test output directory."),
    model_name: str = typer.Option("cpsam_v2", "--model", help="Cellpose pretrained model name."),
    gpu: bool = typer.Option(True, "--gpu/--cpu", help="Request GPU/MPS acceleration."),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    image, _page_count = read_primary_intensity_plane(input_path)
    extracted_path = output_dir / f"{input_path.stem}_candidate_dapi_intensity.tif"
    tifffile.imwrite(extracted_path, image.astype(np.uint16), photometric="minisblack")

    from cellpose import models

    model = models.CellposeModel(gpu=gpu, pretrained_model=model_name, use_bfloat16=False)
    masks = model.eval(
        image,
        channels=[0, 0],
        diameter=None,
        flow_threshold=0.4,
        cellprob_threshold=0.0,
    )[0]

    mask_path = output_dir / f"{input_path.stem}_{model_name}_masks.tif"
    tifffile.imwrite(mask_path, masks.astype(np.uint16), photometric="minisblack")
    montage_path = output_dir / f"{input_path.stem}_{model_name}_smoke_montage.png"
    _write_montage(image=image, masks=masks, output_path=montage_path, model_name=model_name)

    typer.echo(f"extracted_intensity={extracted_path}")
    typer.echo(f"mask={mask_path}")
    typer.echo(f"montage={montage_path}")
    typer.echo(f"label_count={int(masks.max())}")


def _write_montage(*, image: np.ndarray, masks: np.ndarray, output_path: Path, model_name: str) -> None:
    low, high = np.percentile(image, [1, 99.8])
    scaled = np.clip((image.astype(np.float32) - low) / max(high - low, 1), 0, 1)
    overlay = np.dstack([scaled, scaled, scaled])
    boundaries = find_boundaries(masks, mode="outer")
    overlay[boundaries] = [1.0, 0.0, 0.0]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    axes[0].imshow(scaled, cmap="gray")
    axes[0].set_title("Candidate DAPI intensity plane")
    axes[1].imshow(masks, cmap="nipy_spectral")
    axes[1].set_title(f"Cellpose {model_name} labels\ncount={int(masks.max())}")
    axes[2].imshow(overlay)
    axes[2].set_title("Cellpose outlines")
    for ax in axes:
        ax.axis("off")
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    app()

