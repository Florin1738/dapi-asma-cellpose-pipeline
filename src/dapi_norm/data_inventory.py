from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import csv
import re

import numpy as np

from dapi_norm.image_arrays import read_primary_intensity_plane

IMAGE_NAME_RE = re.compile(
    r"^(?P<sample_id>.+)_(?P<position_id>XY\d+)_(?P<token>CH\d+|Overlay)\.tiff?$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class ParsedImageName:
    sample_id: str
    position_id: str
    channel_id: str | None
    kind: str


@dataclass(frozen=True)
class ImageFileInfo:
    path: Path
    width: int
    height: int
    dtype: str
    page_count: int
    min_intensity: float
    max_intensity: float
    mean_intensity: float
    p01_intensity: float
    p50_intensity: float
    p99_intensity: float
    nonzero_fraction: float


@dataclass
class PositionInventory:
    position_id: str
    sample_id: str
    channels: dict[str, ImageFileInfo] = field(default_factory=dict)
    overlay: ImageFileInfo | None = None


@dataclass
class DatasetInventory:
    root: Path
    file_counts_by_suffix: dict[str, int]
    positions: dict[str, PositionInventory]
    sidecar_count: int
    unparsed_tiff_paths: list[Path]

    @property
    def tiff_count(self) -> int:
        return sum(
            count
            for suffix, count in self.file_counts_by_suffix.items()
            if suffix in {".tif", ".tiff"}
        )


def parse_image_filename(filename: str) -> ParsedImageName | None:
    match = IMAGE_NAME_RE.match(Path(filename).name)
    if match is None:
        return None

    token = match.group("token")
    if token.lower() == "overlay":
        return ParsedImageName(
            sample_id=match.group("sample_id"),
            position_id=match.group("position_id").upper(),
            channel_id=None,
            kind="overlay",
        )

    return ParsedImageName(
        sample_id=match.group("sample_id"),
        position_id=match.group("position_id").upper(),
        channel_id=token.upper(),
        kind="channel",
    )


def image_file_info(path: Path, *, allow_display_projection: bool = False) -> ImageFileInfo:
    image_2d, page_count = read_primary_intensity_plane(
        path, allow_multi_active_rgb=allow_display_projection
    )

    height, width = image_2d.shape
    return ImageFileInfo(
        path=path,
        width=int(width),
        height=int(height),
        dtype=str(image_2d.dtype),
        page_count=int(page_count),
        min_intensity=float(np.min(image_2d)),
        max_intensity=float(np.max(image_2d)),
        mean_intensity=float(np.mean(image_2d)),
        p01_intensity=float(np.percentile(image_2d, 1)),
        p50_intensity=float(np.percentile(image_2d, 50)),
        p99_intensity=float(np.percentile(image_2d, 99)),
        nonzero_fraction=float(np.count_nonzero(image_2d) / image_2d.size),
    )


def inventory_dataset(root: Path | str) -> DatasetInventory:
    root_path = Path(root)
    if not root_path.exists():
        raise FileNotFoundError(f"Dataset root does not exist: {root_path}")
    if not root_path.is_dir():
        raise NotADirectoryError(f"Dataset root is not a directory: {root_path}")

    file_counts: dict[str, int] = {}
    positions: dict[str, PositionInventory] = {}
    sidecar_count = 0
    unparsed_tiffs: list[Path] = []

    for path in sorted(p for p in root_path.rglob("*") if p.is_file()):
        suffix = path.suffix.lower()
        file_counts[suffix] = file_counts.get(suffix, 0) + 1

        if suffix not in {".tif", ".tiff"}:
            sidecar_count += 1
            continue

        parsed = parse_image_filename(path.name)
        if parsed is None:
            unparsed_tiffs.append(path)
            continue

        position = positions.setdefault(
            parsed.position_id,
            PositionInventory(position_id=parsed.position_id, sample_id=parsed.sample_id),
        )
        if parsed.kind == "overlay":
            info = image_file_info(path, allow_display_projection=True)
            position.overlay = info
        elif parsed.channel_id is not None:
            info = image_file_info(path, allow_display_projection=False)
            position.channels[parsed.channel_id] = info

    return DatasetInventory(
        root=root_path,
        file_counts_by_suffix=dict(sorted(file_counts.items())),
        positions=dict(sorted(positions.items())),
        sidecar_count=sidecar_count,
        unparsed_tiff_paths=unparsed_tiffs,
    )


def write_inventory_reports(inventory: DatasetInventory, output_dir: Path | str) -> None:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    _write_manifest_csv(inventory, output_path / "image_manifest.csv")
    _write_interpretation_manifest_csv(
        inventory, output_path / "channel_interpretation_manifest.csv"
    )
    _write_summary_markdown(inventory, output_path / "dataset_summary.md")


def _write_manifest_csv(inventory: DatasetInventory, path: Path) -> None:
    fieldnames = [
        "position_id",
        "kind",
        "channel_id",
        "path",
        "width",
        "height",
        "dtype",
        "page_count",
        "min_intensity",
        "p01_intensity",
        "p50_intensity",
        "mean_intensity",
        "p99_intensity",
        "max_intensity",
        "nonzero_fraction",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for position in inventory.positions.values():
            for channel_id, info in sorted(position.channels.items()):
                writer.writerow(_manifest_row(inventory.root, position.position_id, "channel", channel_id, info))
            if position.overlay is not None:
                writer.writerow(
                    _manifest_row(inventory.root, position.position_id, "overlay", "", position.overlay)
                )


def _manifest_row(
    root: Path,
    position_id: str,
    kind: str,
    channel_id: str,
    info: ImageFileInfo,
) -> dict[str, str | int | float]:
    return {
        "position_id": position_id,
        "kind": kind,
        "channel_id": channel_id,
        "path": str(info.path.relative_to(root)),
        "width": info.width,
        "height": info.height,
        "dtype": info.dtype,
        "page_count": info.page_count,
        "min_intensity": info.min_intensity,
        "p01_intensity": info.p01_intensity,
        "p50_intensity": info.p50_intensity,
        "mean_intensity": info.mean_intensity,
        "p99_intensity": info.p99_intensity,
        "max_intensity": info.max_intensity,
        "nonzero_fraction": info.nonzero_fraction,
    }


def _write_interpretation_manifest_csv(inventory: DatasetInventory, path: Path) -> None:
    fieldnames = [
        "position_id",
        "source_file",
        "filename_channel",
        "rgb_component",
        "candidate_stain",
        "confirmed_stain",
        "use_for_segmentation",
        "use_for_measurement",
        "requires_confirmation",
        "notes",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for position in inventory.positions.values():
            for channel_id, info in sorted(position.channels.items()):
                writer.writerow(_interpretation_row(inventory.root, position.position_id, channel_id, info))


def _interpretation_row(
    root: Path,
    position_id: str,
    channel_id: str,
    info: ImageFileInfo,
) -> dict[str, str]:
    if channel_id == "CH4":
        rgb_component = "blue"
        candidate_stain = "candidate_DAPI"
        use_for_segmentation = "unconfirmed"
        use_for_measurement = "unconfirmed"
        notes = "Blue punctate/nuclear-looking signal; candidate DAPI channel, not metadata-confirmed."
    elif channel_id == "CH2":
        rgb_component = "red"
        candidate_stain = "candidate_target"
        use_for_segmentation = "unconfirmed"
        use_for_measurement = "unconfirmed"
        notes = "Red broad/fibrous-looking signal; candidate target channel, not metadata-confirmed."
    else:
        rgb_component = "unknown"
        candidate_stain = "unknown"
        use_for_segmentation = "no"
        use_for_measurement = "no"
        notes = "No channel interpretation available."

    return {
        "position_id": position_id,
        "source_file": str(info.path.relative_to(root)),
        "filename_channel": channel_id,
        "rgb_component": rgb_component,
        "candidate_stain": candidate_stain,
        "confirmed_stain": "",
        "use_for_segmentation": use_for_segmentation,
        "use_for_measurement": use_for_measurement,
        "requires_confirmation": "true",
        "notes": notes,
    }


def _write_summary_markdown(inventory: DatasetInventory, path: Path) -> None:
    channel_ids = sorted(
        {channel_id for position in inventory.positions.values() for channel_id in position.channels}
    )
    lines = [
        "# Dataset Summary",
        "",
        f"Root: `{inventory.root}`",
        f"Positions: {len(inventory.positions)}",
        f"TIFF files: {inventory.tiff_count}",
        f"Sidecar/non-image files: {inventory.sidecar_count}",
        f"Channel IDs: {', '.join(channel_ids) if channel_ids else 'none detected'}",
        "",
        "Important: this inventory treats CH4 as candidate DAPI and CH2 as candidate target "
        "based on RGB export colors and visual morphology. The channel map is not metadata-confirmed.",
        "",
        "## File Counts",
        "",
        "| Suffix | Count |",
        "|---|---:|",
    ]
    for suffix, count in inventory.file_counts_by_suffix.items():
        lines.append(f"| `{suffix or '[no suffix]'}` | {count} |")

    lines.extend(
        [
            "",
            "## Positions",
            "",
            "| Position | Channels | Overlay | Dimensions | Dtype |",
            "|---|---|---|---|---|",
        ]
    )
    for position in inventory.positions.values():
        channels = ", ".join(sorted(position.channels)) or "none"
        overlay = "yes" if position.overlay is not None else "no"
        first_info = next(iter(position.channels.values()), position.overlay)
        dimensions = f"{first_info.width} x {first_info.height}" if first_info is not None else "unknown"
        dtype = first_info.dtype if first_info is not None else "unknown"
        lines.append(
            f"| {position.position_id} | {channels} | {overlay} | {dimensions} | {dtype} |"
        )

    if inventory.unparsed_tiff_paths:
        lines.extend(["", "## Unparsed TIFF Files", ""])
        for unparsed in inventory.unparsed_tiff_paths:
            lines.append(f"- `{unparsed.relative_to(inventory.root)}`")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
