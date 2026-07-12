from __future__ import annotations

from collections import defaultdict
import csv
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from dapi_norm.image_arrays import read_primary_intensity_plane


def render_centroid_overlays(
    *,
    counts_root: Path,
    output_dir: Path,
    contact_sheet_limit: int = 12,
) -> dict[str, int]:
    output_dir.mkdir(parents=True, exist_ok=True)
    overlay_count = 0
    contact_sheet_count = 0
    for counts_csv in sorted(counts_root.glob("*/*/summaries/nucleus_counts.csv")):
        run_output = counts_csv.parents[1]
        run_id = run_output.name
        plate_id = run_output.parent.name
        grouped = _read_centroids_by_image(run_output / "summaries" / "per_nucleus_locations.csv")
        created_for_run: list[Path] = []
        for count_row in _read_count_rows(counts_csv):
            image_id = count_row["image_id"].upper()
            input_path = Path(count_row["input_path"])
            image_rows = grouped.get(image_id, [])
            centroids = [
                (float(row["x_centroid"]), float(row["y_centroid"])) for row in image_rows
            ]
            output_path = output_dir / plate_id / run_id / f"{image_id}_CH4_green_centroids.png"
            _write_overlay(input_path, centroids, output_path)
            created_for_run.append(output_path)
            overlay_count += 1
        if created_for_run:
            _write_contact_sheet(
                created_for_run[:contact_sheet_limit],
                output_dir / plate_id / run_id / "contact_sheet_green_centroids.png",
            )
            contact_sheet_count += 1
    return {"overlay_count": overlay_count, "contact_sheet_count": contact_sheet_count}


def _read_count_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"image_id", "input_path", "nucleus_count"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
        return sorted(list(reader), key=lambda row: _image_sort_key(row["image_id"]))


def _read_centroids_by_image(path: Path) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"image_id", "input_path", "x_centroid", "y_centroid"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            grouped[row["image_id"].upper()].append(row)
    return dict(grouped)


def _write_overlay(input_path: Path, centroids: list[tuple[float, float]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    image, _ = read_primary_intensity_plane(input_path)
    low, high = np.percentile(image, [1, 99.8])
    scaled = np.clip((image.astype(np.float32) - low) / max(high - low, 1), 0, 1)
    gray = (scaled * 255).astype(np.uint8)
    rgb = Image.fromarray(gray, mode="L").convert("RGB")
    draw = ImageDraw.Draw(rgb)
    size = max(4, min(8, round(min(rgb.size) / 120)))
    for x, y in centroids:
        draw.line((x - size, y - size, x + size, y + size), fill=(0, 255, 0), width=2)
        draw.line((x - size, y + size, x + size, y - size), fill=(0, 255, 0), width=2)
    rgb.save(output_path)


def _write_contact_sheet(paths: list[Path], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    thumb_width = 420
    padding = 18
    label_height = 28
    columns = 3
    thumbs: list[tuple[str, Image.Image]] = []
    for path in paths:
        with Image.open(path) as image:
            rgb = image.convert("RGB")
            scale = thumb_width / rgb.width
            thumb = rgb.resize((thumb_width, int(rgb.height * scale)))
        thumbs.append((path.stem.replace("_CH4_green_centroids", ""), thumb))
    cell_height = max(thumb.height for _, thumb in thumbs) + label_height
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


def _image_sort_key(image_id: str) -> tuple[int, str]:
    digits = "".join(char for char in image_id if char.isdigit())
    return (int(digits) if digits else 10**9, image_id)
