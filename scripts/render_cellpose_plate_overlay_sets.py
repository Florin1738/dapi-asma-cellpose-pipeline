#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Any

from dapi_norm.cellpose_endpoint_figures import _write_overlay_pages


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render plate-specific Cellpose overlay QC pages from the merged full-plate summary."
    )
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--panel-page-size", type=int, default=12)
    args = parser.parse_args()

    rows = _read_rows(args.summary)
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row.get("plate", "Plate")].append(row)

    args.output.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, Any]] = []
    for plate in sorted(grouped, key=_plate_sort_key):
        plate_rows = sorted(
            grouped[plate],
            key=lambda row: _number(row, "cellpose_masked_ch2_integrated_raw_per_DAPI_positive_nucleus"),
            reverse=True,
        )
        plate_dir = args.output / plate.replace(" ", "_")
        plate_dir.mkdir(parents=True, exist_ok=True)
        pages, index_rows = _write_overlay_pages(
            plate_rows,
            plate_dir,
            page_size=args.panel_page_size,
        )
        for row in index_rows:
            row["plate_panel_folder"] = str(plate_dir)
            manifest_rows.append(row)
        print(f"{plate}: fields={len(plate_rows)} pages={len(pages)} folder={plate_dir}")

    manifest_path = args.output / "plate_specific_overlay_index.csv"
    _write_manifest(manifest_path, manifest_rows)
    print(f"manifest={manifest_path}")


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_manifest(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = [
        "plate",
        "source_id",
        "location",
        "page",
        "page_path",
        "tile_number_on_page",
        "sort_metric",
        "qc_status",
        "qc_flags",
        "plate_panel_folder",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _number(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "nan"))
    except (TypeError, ValueError):
        return float("nan")


def _plate_sort_key(value: str) -> tuple[int, str]:
    digits = "".join(character for character in value if character.isdigit())
    return (int(digits) if digits else 999, value)


if __name__ == "__main__":
    main()
