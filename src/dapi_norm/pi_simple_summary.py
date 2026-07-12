from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/dapi_norm_matplotlib")

from dapi_norm.image_arrays import read_primary_intensity_plane

HEADERS = ["LOCATION", "aSMA intensity", "Nuclei Count", "Ratio"]
DEFAULT_PLATES = ("Plate 1", "Plate 2")
_SKIP_DIR_NAMES = {
    ".git",
    ".models",
    ".pytest_cache",
    ".ruff_cache",
    ".venv",
    "__pycache__",
    "output",
}


@dataclass(frozen=True)
class ImagePair:
    location: str
    source_id: str
    ch2_path: Path
    ch4_path: Path
    target_channel_id: str = "CH2"
    dapi_channel_id: str = "CH4"

    @property
    def target_path(self) -> Path:
        return self.ch2_path

    @property
    def dapi_path(self) -> Path:
        return self.ch4_path


@dataclass(frozen=True)
class PiSummaryRow:
    location: str
    asma_intensity: float
    nuclei_count: int | None
    source_id: str | None = None

    @property
    def ratio(self) -> float | None:
        if self.nuclei_count is None or self.nuclei_count <= 0:
            return None
        return self.asma_intensity / self.nuclei_count

    @property
    def plot_label(self) -> str:
        return self.source_id or self.location


def build_pi_summary(
    *,
    input_root: Path,
    counts_by_plate: dict[str, dict[str, int]] | None = None,
    plate_names: Iterable[str] = DEFAULT_PLATES,
    target_channel: str = "CH2",
    dapi_channel: str = "CH4",
) -> dict[str, list[PiSummaryRow]]:
    counts_by_plate = counts_by_plate or {}
    plate_names = tuple(plate_names)
    plate_roots = _resolve_plate_roots(input_root, plate_names)
    summary: dict[str, list[PiSummaryRow]] = {}

    for plate_name in plate_names:
        root = plate_roots.get(plate_name)
        if root is None:
            summary[plate_name] = []
            continue
        count_lookup = counts_by_plate.get(plate_name, {})
        pairs = find_image_pairs(
            root,
            target_channel=target_channel,
            dapi_channel=dapi_channel,
        )
        ambiguous_locations = _duplicated_locations(pairs)
        rows = [
            PiSummaryRow(
                location=pair.location,
                asma_intensity=_raw_integrated_intensity(pair.ch2_path),
                nuclei_count=_count_for_pair(count_lookup, pair, ambiguous_locations),
                source_id=pair.source_id if pair.source_id != pair.location else None,
            )
            for pair in pairs
        ]
        summary[plate_name] = sorted(rows, key=_row_sort_key)

    return summary


def find_image_pairs(
    root: Path,
    *,
    target_channel: str = "CH2",
    dapi_channel: str = "CH4",
) -> list[ImagePair]:
    target_channel_id = _normalize_channel_id(target_channel)
    dapi_channel_id = _normalize_channel_id(dapi_channel)
    if target_channel_id == dapi_channel_id:
        raise ValueError("Target and DAPI channels must be different.")
    pairs: list[ImagePair] = []
    for directory in _walk_dirs(root):
        match = re.fullmatch(r"XY\d+", directory.name, flags=re.IGNORECASE)
        if match is None:
            continue
        target_path = _find_channel_file(directory, target_channel_id)
        dapi_path = _find_channel_file(directory, dapi_channel_id)
        if target_path is not None and dapi_path is not None:
            location = directory.name.upper()
            pairs.append(
                ImagePair(
                    location=location,
                    source_id=_source_id_for_position(root, directory, location),
                    ch2_path=target_path,
                    ch4_path=dapi_path,
                    target_channel_id=target_channel_id,
                    dapi_channel_id=dapi_channel_id,
                )
            )
    return sorted(pairs, key=_pair_sort_key)


def load_count_csv(path: Path, *, key_prefix: str | None = None) -> dict[str, int]:
    counts: dict[str, int] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "image_id" not in (reader.fieldnames or []) or "nucleus_count" not in (
            reader.fieldnames or []
        ):
            raise ValueError(f"Expected image_id and nucleus_count columns in {path}")
        for row in reader:
            image_id = row["image_id"].strip()
            count_key = f"{key_prefix}/{image_id}" if key_prefix else image_id
            normalized_key = _count_lookup_key(count_key)
            if normalized_key in counts:
                raise ValueError(f"Duplicate image_id after count-key normalization: {count_key}")
            counts[normalized_key] = int(row["nucleus_count"])
    return counts


def load_counts_by_plate(
    counts_root: Path | None, plate_names: Iterable[str] = DEFAULT_PLATES
) -> dict[str, dict[str, int]]:
    if counts_root is None:
        return {}
    plate_names = tuple(plate_names)
    counts_by_plate: dict[str, dict[str, int]] = {}
    if counts_root.is_file():
        counts_by_plate[plate_names[0]] = load_count_csv(counts_root)
        return counts_by_plate

    single_counts = counts_root / "summaries" / "nucleus_counts.csv"
    if single_counts.exists():
        counts_by_plate[plate_names[0]] = load_count_csv(single_counts)
        return counts_by_plate

    for plate_name in plate_names:
        for plate_root in _count_plate_root_candidates(counts_root, plate_name):
            direct_counts = plate_root / "summaries" / "nucleus_counts.csv"
            if direct_counts.exists():
                counts_by_plate[plate_name] = load_count_csv(direct_counts)
                break
            nested_counts = sorted(plate_root.glob("*/summaries/nucleus_counts.csv"))
            if nested_counts:
                merged: dict[str, int] = {}
                for nested_count in nested_counts:
                    run_id = nested_count.parent.parent.name
                    for key, count in load_count_csv(nested_count, key_prefix=run_id).items():
                        if key in merged:
                            raise ValueError(
                                f"Duplicate count key {key} while loading counts under {plate_root}"
                            )
                        merged[key] = count
                counts_by_plate[plate_name] = merged
                break
    return counts_by_plate


def write_pi_workbook(
    output_path: Path,
    summary: dict[str, list[PiSummaryRow]],
    *,
    metadata: dict[str, str] | None = None,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sheet_names = [name for name in DEFAULT_PLATES if name in summary]
    sheet_names.extend(name for name in summary if name not in sheet_names)
    if not sheet_names:
        sheet_names = list(DEFAULT_PLATES)
        summary = {name: [] for name in sheet_names}
    metadata_sheet_name = "Channel Mapping" if metadata else None
    workbook_sheet_names = [*sheet_names, *([metadata_sheet_name] if metadata_sheet_name else [])]

    with ZipFile(output_path, "w", compression=ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", _content_types_xml(len(workbook_sheet_names)))
        zf.writestr("_rels/.rels", _root_rels_xml())
        zf.writestr("docProps/core.xml", _core_xml())
        zf.writestr("docProps/app.xml", _app_xml(workbook_sheet_names))
        zf.writestr("xl/workbook.xml", _workbook_xml(workbook_sheet_names))
        zf.writestr("xl/_rels/workbook.xml.rels", _workbook_rels_xml(len(workbook_sheet_names)))
        zf.writestr("xl/styles.xml", _styles_xml())
        for index, sheet_name in enumerate(sheet_names, start=1):
            zf.writestr(
                f"xl/worksheets/sheet{index}.xml",
                _sheet_xml(summary.get(sheet_name, [])),
            )
        if metadata is not None:
            zf.writestr(
                f"xl/worksheets/sheet{len(sheet_names) + 1}.xml",
                _metadata_sheet_xml(metadata),
            )


def write_pi_plots(output_dir: Path, summary: dict[str, list[PiSummaryRow]]) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    created: list[Path] = []
    for plate_name, rows in summary.items():
        if not rows:
            continue
        safe_plate_name = plate_name.replace(" ", "_")
        created.append(
            _write_bar_plot(
                output_dir / f"{safe_plate_name}_aSMA_intensity_by_location.png",
                rows,
                title=f"{plate_name}: raw aSMA intensity by location",
                ylabel="Raw integrated CH2/aSMA intensity",
                value_getter=lambda row: row.asma_intensity,
                color="#7c3aed",
            )
        )
        created.append(
            _write_bar_plot(
                output_dir / f"{safe_plate_name}_ratio_by_location.png",
                rows,
                title=f"{plate_name}: aSMA intensity / nuclei count",
                ylabel="Raw aSMA intensity per CH4 nucleus",
                value_getter=lambda row: row.ratio,
                color="#2563eb",
            )
        )
    return created


def _resolve_plate_roots(input_root: Path, plate_names: tuple[str, ...]) -> dict[str, Path | None]:
    direct_children = {child.name.lower(): child for child in input_root.iterdir() if child.is_dir()}
    has_plate_folder = any(name.lower() in direct_children for name in plate_names)
    if not has_plate_folder:
        return {plate_names[0]: input_root, **{name: None for name in plate_names[1:]}}
    return {name: direct_children.get(name.lower()) for name in plate_names}


def _count_plate_root_candidates(counts_root: Path, plate_name: str) -> list[Path]:
    return [
        counts_root / plate_name,
        counts_root / plate_name.replace(" ", "_"),
        counts_root / plate_name.lower().replace(" ", "_"),
        counts_root / plate_name.lower(),
    ]


def _raw_integrated_intensity(path: Path) -> float:
    image, _ = read_primary_intensity_plane(path)
    return float(np.sum(image.astype(np.float64)))


def _walk_dirs(root: Path) -> Iterable[Path]:
    stack = [root]
    while stack:
        current = stack.pop()
        if current.name in _SKIP_DIR_NAMES:
            continue
        yield current
        try:
            children = [child for child in current.iterdir() if child.is_dir()]
        except OSError:
            continue
        stack.extend(reversed(sorted(children)))


def _find_channel_file(directory: Path, channel: str) -> Path | None:
    channel_id = _normalize_channel_id(channel)
    candidates = [
        path
        for path in directory.iterdir()
        if path.is_file()
        and path.suffix.lower() in {".tif", ".tiff"}
        and re.search(rf"(?:^|_){re.escape(channel_id)}(?:_|$)", path.stem, flags=re.IGNORECASE)
    ]
    return sorted(candidates)[0] if candidates else None


def _normalize_channel_id(channel: str) -> str:
    normalized = channel.strip().upper()
    if not re.fullmatch(r"CH\d+", normalized):
        raise ValueError(f"Channel must look like CH1, CH2, CH4, etc.; got {channel!r}.")
    return normalized


def _source_id_for_position(root: Path, directory: Path, location: str) -> str:
    try:
        relative_parts = directory.relative_to(root).parts
    except ValueError:
        return location
    if len(relative_parts) >= 2:
        return f"{relative_parts[-2]}/{location}"
    return location


def _count_for_pair(
    count_lookup: dict[str, int], pair: ImagePair, ambiguous_locations: set[str]
) -> int | None:
    source_key = _count_lookup_key(pair.source_id)
    if source_key in count_lookup:
        return count_lookup[source_key]

    location_key = _count_lookup_key(pair.location)
    if pair.location in ambiguous_locations and location_key in count_lookup:
        raise ValueError(
            f"Unprefixed count key {pair.location} is ambiguous because that LOCATION appears "
            "in multiple run folders. Use per-run Cellpose count outputs so counts can be matched "
            "as run_folder/XY##."
        )
    if pair.location not in ambiguous_locations and location_key in count_lookup:
        return count_lookup[location_key]
    return None


def _duplicated_locations(pairs: list[ImagePair]) -> set[str]:
    seen: set[str] = set()
    duplicated: set[str] = set()
    for pair in pairs:
        if pair.location in seen:
            duplicated.add(pair.location)
        seen.add(pair.location)
    return duplicated


def _count_lookup_key(value: str) -> str:
    return value.strip().replace("\\", "/").upper()


def _pair_sort_key(pair: ImagePair) -> tuple[str, int, str]:
    return _source_sort_key(pair.source_id)


def _row_sort_key(row: PiSummaryRow) -> tuple[str, int, str]:
    return _source_sort_key(row.plot_label)


def _source_sort_key(source_id: str) -> tuple[str, int, str]:
    prefix, _, location = source_id.rpartition("/")
    if not location:
        location = source_id
    return (prefix.lower(), *_location_sort_key(location))


def _location_sort_key(location: str) -> tuple[int, str]:
    match = re.search(r"(\d+)$", location)
    return (int(match.group(1)) if match else 10**9, location)


def _sheet_xml(rows: list[PiSummaryRow]) -> str:
    row_xml = [_row_xml(1, HEADERS, style_id=1)]
    for row_index, row in enumerate(rows, start=2):
        values = [row.location, row.asma_intensity, row.nuclei_count, None]
        formula = f'IF(C{row_index}>0,B{row_index}/C{row_index},"")'
        row_xml.append(_data_row_xml(row_index, values, formula=formula, cached_ratio=row.ratio))

    dimension_last_row = max(len(rows) + 1, 1)
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<dimension ref="A1:D{dimension_last_row}"/>'
        '<sheetViews><sheetView showGridLines="0" workbookViewId="0">'
        '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
        "</sheetView></sheetViews>"
        '<sheetFormatPr defaultRowHeight="15"/>'
        '<cols><col min="1" max="1" width="14" customWidth="1"/>'
        '<col min="2" max="2" width="20" customWidth="1"/>'
        '<col min="3" max="3" width="14" customWidth="1"/>'
        '<col min="4" max="4" width="18" customWidth="1"/></cols>'
        f"<sheetData>{''.join(row_xml)}</sheetData>"
        '<autoFilter ref="A1:D1"/>'
        "</worksheet>"
    )


def _metadata_sheet_xml(metadata: dict[str, str]) -> str:
    rows = [_row_xml(1, ["Field", "Value"], style_id=1)]
    for row_index, (field, value) in enumerate(metadata.items(), start=2):
        rows.append(_row_xml(row_index, [field, value], style_id=0))
    last_row = max(len(metadata) + 1, 1)
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f'<dimension ref="A1:B{last_row}"/>'
        '<sheetViews><sheetView showGridLines="0" workbookViewId="0">'
        '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
        "</sheetView></sheetViews>"
        '<sheetFormatPr defaultRowHeight="15"/>'
        '<cols><col min="1" max="1" width="34" customWidth="1"/>'
        '<col min="2" max="2" width="72" customWidth="1"/></cols>'
        f"<sheetData>{''.join(rows)}</sheetData>"
        '<autoFilter ref="A1:B1"/>'
        "</worksheet>"
    )


def _row_xml(row_index: int, values: list[str], *, style_id: int) -> str:
    cells = [
        _cell_xml(row_index, col_index, value, style_id=style_id)
        for col_index, value in enumerate(values, start=1)
    ]
    return f'<row r="{row_index}" ht="20" customHeight="1">{"".join(cells)}</row>'


def _data_row_xml(
    row_index: int, values: list[str | float | int | None], *, formula: str, cached_ratio: float | None
) -> str:
    cells = [
        _cell_xml(row_index, 1, values[0], style_id=0),
        _cell_xml(row_index, 2, values[1], style_id=2),
        _cell_xml(row_index, 3, values[2], style_id=2),
        _formula_cell_xml(row_index, 4, formula, cached_ratio, style_id=3),
    ]
    return f'<row r="{row_index}">{"".join(cells)}</row>'


def _cell_xml(row: int, col: int, value: str | float | int | None, *, style_id: int) -> str:
    if value is None:
        return ""
    ref = f"{_column_name(col)}{row}"
    style = f' s="{style_id}"' if style_id else ""
    if isinstance(value, str):
        return (
            f'<c r="{ref}" t="inlineStr"{style}><is><t>'
            f"{escape(value)}"
            "</t></is></c>"
        )
    return f'<c r="{ref}"{style}><v>{_number_text(float(value))}</v></c>'


def _formula_cell_xml(
    row: int, col: int, formula: str, cached_value: float | None, *, style_id: int
) -> str:
    ref = f"{_column_name(col)}{row}"
    value_xml = "" if cached_value is None else f"<v>{_number_text(cached_value)}</v>"
    return f'<c r="{ref}" s="{style_id}"><f>{escape(formula)}</f>{value_xml}</c>'


def _column_name(index: int) -> str:
    letters = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(65 + remainder) + letters
    return letters


def _number_text(value: float) -> str:
    if value.is_integer():
        return str(int(value))
    return f"{value:.12g}"


def _write_bar_plot(
    output_path: Path,
    rows: list[PiSummaryRow],
    *,
    title: str,
    ylabel: str,
    value_getter,
    color: str,
) -> Path:
    import matplotlib.pyplot as plt

    labels = [row.plot_label for row in rows]
    values = [value_getter(row) for row in rows]
    numeric_values = [float(value) if value is not None else np.nan for value in values]

    row_count = len(rows)
    width = max(8.0, min(30.0, 0.24 * row_count + 5.0))
    height = 5.4 if row_count > 40 else 4.8
    label_size = 5 if row_count > 80 else 6 if row_count > 40 else 8
    rotation = 90 if row_count > 12 else 0
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
    bars = ax.bar(labels, numeric_values, color=color)
    ax.set_title(title)
    ax.set_xlabel("LOCATION")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", linewidth=0.4, alpha=0.35)
    ax.tick_params(axis="x", labelrotation=rotation, labelsize=label_size)
    if row_count <= 40:
        for bar, value in zip(bars, numeric_values, strict=True):
            if np.isfinite(value):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height(),
                    _plot_label(value),
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    rotation=0,
                )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def _plot_label(value: float) -> str:
    if abs(value) >= 1_000_000_000:
        return f"{value / 1_000_000_000:.1f}B"
    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    return f"{value:.1f}"


def _content_types_xml(sheet_count: int) -> str:
    sheets = "".join(
        f'<Override PartName="/xl/worksheets/sheet{i}.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        for i in range(1, sheet_count + 1)
    )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/docProps/core.xml" '
        'ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
        '<Override PartName="/docProps/app.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>'
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/styles.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        f"{sheets}</Types>"
    )


def _root_rels_xml() -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        '<Relationship Id="rId2" '
        'Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" '
        'Target="docProps/core.xml"/>'
        '<Relationship Id="rId3" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" '
        'Target="docProps/app.xml"/>'
        "</Relationships>"
    )


def _workbook_xml(sheet_names: list[str]) -> str:
    sheets = "".join(
        f'<sheet name="{escape(name)}" sheetId="{index}" r:id="rId{index}"/>'
        for index, name in enumerate(sheet_names, start=1)
    )
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        "<sheets>"
        f"{sheets}"
        "</sheets>"
        '<calcPr calcMode="auto" fullCalcOnLoad="1"/>'
        "</workbook>"
    )


def _workbook_rels_xml(sheet_count: int) -> str:
    sheets = "".join(
        f'<Relationship Id="rId{i}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        f'Target="worksheets/sheet{i}.xml"/>'
        for i in range(1, sheet_count + 1)
    )
    styles_id = sheet_count + 1
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f"{sheets}"
        f'<Relationship Id="rId{styles_id}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
        "</Relationships>"
    )


def _styles_xml() -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<numFmts count="2">'
        '<numFmt numFmtId="164" formatCode="#,##0"/>'
        '<numFmt numFmtId="165" formatCode="#,##0.00"/>'
        "</numFmts>"
        '<fonts count="2">'
        '<font><sz val="11"/><color theme="1"/><name val="Calibri"/><family val="2"/></font>'
        '<font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/><family val="2"/></font>'
        "</fonts>"
        '<fills count="3">'
        '<fill><patternFill patternType="none"/></fill>'
        '<fill><patternFill patternType="gray125"/></fill>'
        '<fill><patternFill patternType="solid"><fgColor rgb="FF1F4E78"/><bgColor indexed="64"/></patternFill></fill>'
        "</fills>"
        '<borders count="2">'
        "<border><left/><right/><top/><bottom/><diagonal/></border>"
        '<border><left style="thin"><color rgb="FFD9E2F3"/></left>'
        '<right style="thin"><color rgb="FFD9E2F3"/></right>'
        '<top style="thin"><color rgb="FFD9E2F3"/></top>'
        '<bottom style="thin"><color rgb="FFD9E2F3"/></bottom><diagonal/></border>'
        "</borders>"
        '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
        '<cellXfs count="4">'
        '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
        '<xf numFmtId="0" fontId="1" fillId="2" borderId="1" xfId="0" applyFont="1" applyFill="1" applyBorder="1" applyAlignment="1">'
        '<alignment horizontal="center"/></xf>'
        '<xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
        '<xf numFmtId="165" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/>'
        "</cellXfs>"
        '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
        "</styleSheet>"
    )


def _core_xml() -> str:
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/" '
        'xmlns:dcterms="http://purl.org/dc/terms/" '
        'xmlns:dcmitype="http://purl.org/dc/dcmitype/" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        "<dc:creator>Codex</dc:creator>"
        "<cp:lastModifiedBy>Codex</cp:lastModifiedBy>"
        "</cp:coreProperties>"
    )


def _app_xml(sheet_names: list[str]) -> str:
    titles = "".join(f"<vt:lpstr>{escape(name)}</vt:lpstr>" for name in sheet_names)
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" '
        'xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">'
        "<Application>Codex</Application>"
        "<HeadingPairs><vt:vector size=\"2\" baseType=\"variant\">"
        "<vt:variant><vt:lpstr>Worksheets</vt:lpstr></vt:variant>"
        f"<vt:variant><vt:i4>{len(sheet_names)}</vt:i4></vt:variant>"
        "</vt:vector></HeadingPairs>"
        f"<TitlesOfParts><vt:vector size=\"{len(sheet_names)}\" baseType=\"lpstr\">{titles}</vt:vector></TitlesOfParts>"
        "</Properties>"
    )
