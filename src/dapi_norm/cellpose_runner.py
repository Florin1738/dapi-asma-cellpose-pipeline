from __future__ import annotations

from collections.abc import Callable
from datetime import UTC, datetime
import hashlib
import importlib.metadata as importlib_metadata
import os
from pathlib import Path
import platform
import re
import sys
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw
import tifffile
from skimage.segmentation import find_boundaries
import yaml

from dapi_norm.data_inventory import inventory_dataset
from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.nuclei_outputs import summarize_labeled_mask, write_nuclei_count_tables

Segmenter = Callable[[np.ndarray], np.ndarray]
CELLPOSE_CHANNELS = [0, 0]
CELLPOSE_DIAMETER = None
CELLPOSE_FLOW_THRESHOLD = 0.4
CELLPOSE_CELLPROB_THRESHOLD = 0.0


def run_nuclei_count_batch(
    *,
    input_root: Path | str,
    output_dir: Path | str,
    channel_id: str = "CH4",
    model_name: str = "cpsam_v2",
    gpu: bool = True,
    max_images: int | None = None,
    segmenter: Segmenter | None = None,
    channel_identity_confirmed: bool = False,
) -> list[dict[str, Any]]:
    input_path = Path(input_root)
    output_path = Path(output_dir)
    inventory = inventory_dataset(input_path)
    channel_key = channel_id.upper()

    image_items = [
        (position_id, position.channels[channel_key].path)
        for position_id, position in inventory.positions.items()
        if channel_key in position.channels
    ]
    if max_images is not None:
        image_items = image_items[:max_images]
    if not image_items:
        raise ValueError(f"No {channel_key} images found under {input_path}")

    masks_dir = output_path / "masks"
    qc_dir = output_path / "qc"
    summaries_dir = output_path / "summaries"
    logs_dir = output_path / "logs"
    masks_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    active_segmenter = segmenter
    if active_segmenter is None:
        active_segmenter = _cellpose_segmenter(model_name=model_name, gpu=gpu)

    summary_rows: list[dict[str, Any]] = []
    nucleus_rows: list[dict[str, Any]] = []
    image_records: list[dict[str, Any]] = []
    montage_paths: list[Path] = []

    for position_id, image_file in image_items:
        image, page_count = read_primary_intensity_plane(image_file)
        mask = validate_label_mask(active_segmenter(image), image_shape=image.shape)
        safe_model = _safe_token(model_name)
        mask_path = masks_dir / f"{position_id}_{channel_key}_{safe_model}_labels.tif"
        montage_path = qc_dir / f"{position_id}_{channel_key}_{safe_model}_montage.png"
        tifffile.imwrite(mask_path, mask, photometric="minisblack")
        montage_paths.append(montage_path)
        write_segmentation_montage(
            image=image,
            mask=mask,
            output_path=montage_path,
            title=f"{position_id} {channel_key} {model_name}",
        )
        summary, rows = summarize_labeled_mask(
            image_id=position_id,
            input_path=image_file,
            mask=mask,
            backend="cellpose",
            model_name=model_name,
            channel_id=channel_key,
            candidate_stain="candidate_DAPI"
            if channel_identity_confirmed or channel_key == "CH4"
            else "unknown",
            channel_identity_confirmed=channel_identity_confirmed,
            mask_path=mask_path,
            qc_montage_path=montage_path,
        )
        summary_rows.append(summary)
        nucleus_rows.extend(rows)
        image_records.append(
            {
                "image_id": position_id,
                "input_path": str(image_file),
                "page_count": int(page_count),
                "extracted_plane_shape": [int(image.shape[0]), int(image.shape[1])],
                "extracted_plane_dtype": str(image.dtype),
                "channel_id": channel_key,
                "candidate_stain": summary["candidate_stain"],
                "channel_identity_confirmed": channel_identity_confirmed,
                "nucleus_count": summary["nucleus_count"],
                "mask_path": str(mask_path),
                "qc_montage_path": str(montage_path),
            }
        )

    write_nuclei_count_tables(
        output_dir=summaries_dir,
        summary_rows=summary_rows,
        nucleus_rows=nucleus_rows,
    )
    contact_sheet_path = output_path / "qc_contact_sheet.png"
    write_qc_contact_sheet(montage_paths, contact_sheet_path)
    write_run_metadata(
        logs_dir=logs_dir,
        input_root=input_path,
        output_dir=output_path,
        channel_id=channel_key,
        model_name=model_name,
        gpu=gpu,
        channel_identity_confirmed=channel_identity_confirmed,
        image_records=image_records,
        contact_sheet_path=contact_sheet_path,
    )
    return summary_rows


def _cellpose_segmenter(*, model_name: str, gpu: bool) -> Segmenter:
    from cellpose import models

    model = models.CellposeModel(gpu=gpu, pretrained_model=model_name, use_bfloat16=False)

    def segment(image: np.ndarray) -> np.ndarray:
        return model.eval(
            image,
            channels=CELLPOSE_CHANNELS,
            diameter=CELLPOSE_DIAMETER,
            flow_threshold=CELLPOSE_FLOW_THRESHOLD,
            cellprob_threshold=CELLPOSE_CELLPROB_THRESHOLD,
        )[0]

    return segment


def write_segmentation_montage(
    *,
    image: np.ndarray,
    mask: np.ndarray,
    output_path: Path,
    title: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    low, high = np.percentile(image, [1, 99.8])
    scaled = np.clip((image.astype(np.float32) - low) / max(high - low, 1), 0, 1)
    overlay = np.dstack([scaled, scaled, scaled])
    boundaries = find_boundaries(mask, mode="outer")
    overlay[boundaries] = [1.0, 0.0, 0.0]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    axes[0].imshow(scaled, cmap="gray")
    axes[0].set_title(f"{title}\nintensity")
    axes[1].imshow(mask, cmap="nipy_spectral")
    axes[1].set_title(f"labels\ncount={count_labels(mask)}")
    axes[2].imshow(overlay)
    axes[2].set_title("outlines")
    for ax in axes:
        ax.axis("off")
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def _safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return token or "model"


def count_labels(mask: np.ndarray) -> int:
    labels = np.unique(mask)
    return int(np.count_nonzero(labels))


def validate_label_mask(mask: np.ndarray, *, image_shape: tuple[int, int]) -> np.ndarray:
    label_mask = np.asarray(mask)
    if label_mask.ndim != 2:
        raise ValueError("Cellpose segmenter returned a non-2-D mask")
    if label_mask.shape != image_shape:
        raise ValueError(
            "Cellpose segmenter returned a mask that is not the same shape as the input image: "
            f"mask={label_mask.shape}, image={image_shape}"
        )
    if np.any(label_mask < 0):
        raise ValueError("Cellpose segmenter returned negative label IDs")
    if label_mask.max(initial=0) > np.iinfo(np.uint32).max:
        raise ValueError("Cellpose segmenter returned label IDs larger than uint32")
    return label_mask.astype(np.uint32, copy=False)


def write_qc_contact_sheet(montage_paths: list[Path], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not montage_paths:
        raise ValueError("Cannot write a contact sheet without QC montage paths")

    thumb_width = 720
    padding = 24
    label_height = 34
    columns = 2
    thumbs: list[tuple[str, Image.Image]] = []
    for path in montage_paths:
        with Image.open(path) as image:
            rgb_image = image.convert("RGB")
            scale = thumb_width / rgb_image.width
            thumb = rgb_image.resize((thumb_width, int(rgb_image.height * scale)))
        thumbs.append((path.stem.replace("_cpsam_v2_montage", ""), thumb))

    cell_height = max(thumb.height for _label, thumb in thumbs) + label_height
    rows = (len(thumbs) + columns - 1) // columns
    sheet = Image.new(
        "RGB",
        (columns * thumb_width + (columns + 1) * padding, rows * cell_height + (rows + 1) * padding),
        "white",
    )
    draw = ImageDraw.Draw(sheet)
    for index, (label, thumb) in enumerate(thumbs):
        row, column = divmod(index, columns)
        x_pos = padding + column * (thumb_width + padding)
        y_pos = padding + row * (cell_height + padding)
        draw.text((x_pos, y_pos), label, fill=(0, 0, 0))
        sheet.paste(thumb, (x_pos, y_pos + label_height))

    sheet.save(output_path)


def write_run_metadata(
    *,
    logs_dir: Path,
    input_root: Path,
    output_dir: Path,
    channel_id: str,
    model_name: str,
    gpu: bool,
    channel_identity_confirmed: bool,
    image_records: list[dict[str, Any]],
    contact_sheet_path: Path,
) -> None:
    logs_dir.mkdir(parents=True, exist_ok=True)
    config = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "input_root": str(input_root),
        "output_dir": str(output_dir),
        "channel_id": channel_id,
        "channel_identity_confirmed": channel_identity_confirmed,
        "channel_extraction": {
            "source": "primary TIFF series",
            "axes_policy": (
                "singleton axes are squeezed; 2-D images are used directly; RGB/YXS exports with "
                "exactly one active sample use that active sample"
            ),
            "z_projection": "none",
            "current_channel_note": (
                f"{channel_id} is treated as DAPI/nuclei for this run when channel identity is "
                "confirmed; otherwise CH4 is treated as candidate DAPI based on blue, punctate "
                "morphology and non-CH4 channels are logged as unknown."
            ),
        },
        "model": _model_metadata(model_name),
        "software_versions": _software_versions(),
        "device": _device_metadata(gpu),
        "segmentation_parameters": {
            "channels": CELLPOSE_CHANNELS,
            "diameter": CELLPOSE_DIAMETER,
            "flow_threshold": CELLPOSE_FLOW_THRESHOLD,
            "cellprob_threshold": CELLPOSE_CELLPROB_THRESHOLD,
        },
        "filtering": {
            "implemented": False,
            "kept_after_filtering_default": True,
            "notes": "No area, border, or shape filtering is applied in the current count-only run.",
        },
        "background_settings": "not_applicable_count_only_run",
        "outputs": {
            "nucleus_counts_csv": str(output_dir / "summaries" / "nucleus_counts.csv"),
            "per_nucleus_locations_csv": str(
                output_dir / "summaries" / "per_nucleus_locations.csv"
            ),
            "masks_dir": str(output_dir / "masks"),
            "qc_dir": str(output_dir / "qc"),
            "qc_contact_sheet": str(contact_sheet_path),
        },
        "image_inputs": image_records,
    }
    (logs_dir / "config_resolved.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False), encoding="utf-8"
    )
    total_count = sum(int(record["nucleus_count"]) for record in image_records)
    (logs_dir / "run_log.txt").write_text(
        "\n".join(
            [
                f"generated_at_utc: {config['generated_at_utc']}",
                f"input_root: {input_root}",
                f"output_dir: {output_dir}",
                f"channel_id: {channel_id}",
                f"model_name: {model_name}",
                f"requested_gpu: {gpu}",
                f"images_processed: {len(image_records)}",
                f"total_nucleus_count: {total_count}",
                f"channel_identity_confirmed: {str(channel_identity_confirmed).lower()}",
                f"warnings: {'' if channel_identity_confirmed else 'channel_identity_unconfirmed'}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _software_versions() -> dict[str, str]:
    versions = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
    }
    for package_name in ["cellpose", "torch", "numpy", "scikit-image", "tifffile"]:
        versions[package_name] = _package_version(package_name)
    return versions


def _package_version(package_name: str) -> str:
    try:
        return importlib_metadata.version(package_name)
    except importlib_metadata.PackageNotFoundError:
        return "not_installed"


def _device_metadata(gpu: bool) -> dict[str, Any]:
    metadata: dict[str, Any] = {"requested_gpu": gpu}
    try:
        import torch
    except ImportError:
        metadata["torch_importable"] = False
        return metadata

    metadata["torch_importable"] = True
    metadata["mps_available"] = bool(torch.backends.mps.is_available())
    metadata["cuda_available"] = bool(torch.cuda.is_available())
    metadata["cuda_device_count"] = int(torch.cuda.device_count())
    if torch.cuda.is_available():
        metadata["cuda_device_name"] = torch.cuda.get_device_name(0)
    return metadata


def _model_metadata(model_name: str) -> dict[str, str | None]:
    model_path = _resolved_model_path(model_name)
    return {
        "name": model_name,
        "path": None if model_path is None else str(model_path),
        "sha256": None if model_path is None else _sha256(model_path),
    }


def _resolved_model_path(model_name: str) -> Path | None:
    direct_path = Path(model_name)
    if direct_path.exists():
        return direct_path
    model_root = os.environ.get("CELLPOSE_LOCAL_MODELS_PATH")
    if model_root:
        candidate = Path(model_root) / model_name
        if candidate.exists():
            return candidate
    return None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
