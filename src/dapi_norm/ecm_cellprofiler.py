from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import tifffile
from skimage import segmentation

from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.user_cellpose_batch import discover_acquisitions
from dapi_norm.pi_simple_summary import find_image_pairs


ECM_ENDPOINT = "ecm_positive_integrated_background_corrected"


@dataclass(frozen=True)
class EcmImageRecord:
    acquisition_rel: str
    acquisition_name: str
    acquisition_slug: str
    location: str
    source_id: str
    ecm_path: Path
    dapi_path: Path

    @property
    def staged_name(self) -> str:
        return f"{self.acquisition_slug}__{self.location}_ECM.tif"

    @property
    def staged_stem(self) -> str:
        return self.staged_name.removesuffix(".tif")


@dataclass(frozen=True)
class EcmMeasurement:
    acquisition: str
    acquisition_rel: str
    location: str
    source_id: str
    dapi_positive_nucleus_count: int | None
    ecm_channel_id: str
    dapi_channel_id: str
    threshold_method: str
    threshold_deviations: float
    image_area_px: int
    whole_field_ecm_integrated_raw: float
    whole_field_ecm_mean_raw: float
    whole_field_ecm_median_raw: float
    ecm_background_value_per_px: float
    ecm_background_source: str
    ecm_positive_area_px: int
    ecm_positive_area_fraction: float
    ecm_positive_integrated_raw: float
    ecm_positive_integrated_background_corrected: float
    ecm_positive_mean_raw: float | None
    ecm_positive_mean_background_corrected: float | None
    ecm_mask_object_count: int
    qc_status: str
    qc_flags: str
    ecm_path: str
    dapi_path: str
    cellprofiler_mask_path: str
    cellprofiler_overlay_path: str
    qc_panel_path: str

    def as_row(self) -> dict[str, object]:
        return {
            "acquisition": self.acquisition,
            "acquisition_rel": self.acquisition_rel,
            "location": self.location,
            "source_id": self.source_id,
            "dapi_positive_nucleus_count": self.dapi_positive_nucleus_count,
            "ecm_channel_id": self.ecm_channel_id,
            "dapi_channel_id": self.dapi_channel_id,
            "threshold_method": self.threshold_method,
            "threshold_deviations": self.threshold_deviations,
            "image_area_px": self.image_area_px,
            "whole_field_ecm_integrated_raw": self.whole_field_ecm_integrated_raw,
            "whole_field_ecm_mean_raw": self.whole_field_ecm_mean_raw,
            "whole_field_ecm_median_raw": self.whole_field_ecm_median_raw,
            "ecm_background_value_per_px": self.ecm_background_value_per_px,
            "ecm_background_source": self.ecm_background_source,
            "ecm_positive_area_px": self.ecm_positive_area_px,
            "ecm_positive_area_fraction": self.ecm_positive_area_fraction,
            "ecm_positive_integrated_raw": self.ecm_positive_integrated_raw,
            ECM_ENDPOINT: self.ecm_positive_integrated_background_corrected,
            "ecm_positive_mean_raw": self.ecm_positive_mean_raw,
            "ecm_positive_mean_background_corrected": self.ecm_positive_mean_background_corrected,
            "ecm_mask_object_count": self.ecm_mask_object_count,
            "qc_status": self.qc_status,
            "qc_flags": self.qc_flags,
            "ecm_path": self.ecm_path,
            "dapi_path": self.dapi_path,
            "cellprofiler_mask_path": self.cellprofiler_mask_path,
            "cellprofiler_overlay_path": self.cellprofiler_overlay_path,
            "qc_panel_path": self.qc_panel_path,
        }


def discover_ecm_records(
    input_root: Path,
    *,
    ecm_channel_id: str = "CH1",
    dapi_channel_id: str = "CH4",
) -> list[EcmImageRecord]:
    input_root = input_root.expanduser().resolve()
    acquisitions = discover_acquisitions(
        input_root,
        target_channel_id=ecm_channel_id,
        dapi_channel_id=dapi_channel_id,
    )
    records: list[EcmImageRecord] = []
    used_slugs: set[str] = set()
    for acquisition in acquisitions:
        try:
            acquisition_rel = acquisition.input_root.relative_to(input_root).as_posix()
        except ValueError:
            acquisition_rel = acquisition.input_root.name
        if acquisition_rel == ".":
            acquisition_rel = acquisition.input_root.name
        slug = _unique_slug(acquisition.display_name, used_slugs)
        pairs = find_image_pairs(
            acquisition.input_root,
            target_channel=ecm_channel_id,
            dapi_channel=dapi_channel_id,
        )
        for pair in pairs:
            records.append(
                EcmImageRecord(
                    acquisition_rel=acquisition_rel,
                    acquisition_name=acquisition.display_name,
                    acquisition_slug=slug,
                    location=pair.location,
                    source_id=f"{acquisition.display_name}/{pair.location}",
                    ecm_path=pair.target_path.resolve(),
                    dapi_path=pair.dapi_path.resolve(),
                )
            )
    return sorted(records, key=lambda row: (_natural_key(row.acquisition_rel), _natural_key(row.location)))


def stage_ecm_images(records: list[EcmImageRecord], staging_dir: Path, manifest_path: Path) -> None:
    staging_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    for record in records:
        ecm, _ = read_primary_intensity_plane(record.ecm_path)
        out_path = staging_dir / record.staged_name
        tifffile.imwrite(out_path, ecm.astype(np.uint16), photometric="minisblack")
        rows.append(
            {
                "acquisition": record.acquisition_name,
                "acquisition_rel": record.acquisition_rel,
                "location": record.location,
                "source_id": record.source_id,
                "source_ecm_path": str(record.ecm_path),
                "source_dapi_path": str(record.dapi_path),
                "staged_ecm_path": str(out_path),
                "staging_note": "RGB pseudocolor TIFF collapsed to active single channel before CellProfiler",
            }
        )
    write_csv(manifest_path, rows)


def select_representative_records(
    records: list[EcmImageRecord],
    *,
    representative_count: int,
) -> list[EcmImageRecord]:
    if representative_count >= len(records):
        return records
    scored: list[tuple[float, EcmImageRecord]] = []
    for record in records:
        image, _ = read_primary_intensity_plane(record.ecm_path)
        scored.append((float(np.percentile(image.astype(float), 99.0)), record))
    scored.sort(key=lambda item: item[0])
    indices = np.linspace(0, len(scored) - 1, representative_count)
    selected: list[EcmImageRecord] = []
    seen: set[tuple[str, str]] = set()
    for idx in indices:
        record = scored[int(round(idx))][1]
        key = (record.acquisition_rel, record.location)
        if key not in seen:
            selected.append(record)
            seen.add(key)
    for _, record in scored:
        if len(selected) >= representative_count:
            break
        key = (record.acquisition_rel, record.location)
        if key not in seen:
            selected.append(record)
            seen.add(key)
    return sorted(selected, key=lambda row: (_natural_key(row.acquisition_rel), _natural_key(row.location)))


def load_dapi_counts(cellpose_report_root: Path) -> dict[tuple[str, str], int]:
    lookup: dict[tuple[str, str], int] = {}
    for csv_path in cellpose_report_root.glob("**/cellpose_counts/*/*/summaries/nucleus_counts.csv"):
        acquisition = csv_path.parent.parent.name
        with csv_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                lookup[(acquisition, row["image_id"].strip().upper())] = int(row["nucleus_count"])
    return lookup


def measure_ecm_from_mask(
    *,
    record: EcmImageRecord,
    image: np.ndarray,
    mask: np.ndarray,
    dapi_count: int | None,
    threshold_deviations: float,
    ecm_channel_id: str,
    dapi_channel_id: str,
    mask_path: Path,
    overlay_path: Path,
    qc_panel_path: Path,
    root_for_paths: Path,
) -> EcmMeasurement:
    image_float = image.astype(float)
    mask_bool = np.asarray(mask) > 0
    background_pixels = image_float[~mask_bool]
    if background_pixels.size >= max(100, int(image_float.size * 0.01)):
        background_value = float(np.median(background_pixels))
        background_source = "median_cellprofiler_mask_negative_pixels"
    else:
        background_value = float(np.percentile(image_float, 10))
        background_source = "fallback_p10_whole_image_mask_too_large"
    corrected = np.clip(image_float - background_value, 0, None)
    positive_area = int(np.count_nonzero(mask_bool))
    positive_raw = float(np.sum(image_float[mask_bool]))
    positive_corrected = float(np.sum(corrected[mask_bool]))
    image_area = int(image_float.size)
    flags = ecm_qc_flags(
        positive_area_fraction=positive_area / image_area,
        background_source=background_source,
        positive_corrected=positive_corrected,
    )
    status = "review" if flags else "pass"
    return EcmMeasurement(
        acquisition=record.acquisition_name,
        acquisition_rel=record.acquisition_rel,
        location=record.location,
        source_id=record.source_id,
        dapi_positive_nucleus_count=dapi_count,
        ecm_channel_id=ecm_channel_id,
        dapi_channel_id=dapi_channel_id,
        threshold_method="cellprofiler_robust_background",
        threshold_deviations=threshold_deviations,
        image_area_px=image_area,
        whole_field_ecm_integrated_raw=float(np.sum(image_float)),
        whole_field_ecm_mean_raw=float(np.mean(image_float)),
        whole_field_ecm_median_raw=float(np.median(image_float)),
        ecm_background_value_per_px=background_value,
        ecm_background_source=background_source,
        ecm_positive_area_px=positive_area,
        ecm_positive_area_fraction=positive_area / image_area,
        ecm_positive_integrated_raw=positive_raw,
        ecm_positive_integrated_background_corrected=positive_corrected,
        ecm_positive_mean_raw=(positive_raw / positive_area if positive_area else None),
        ecm_positive_mean_background_corrected=(
            positive_corrected / positive_area if positive_area else None
        ),
        ecm_mask_object_count=count_labels(mask),
        qc_status=status,
        qc_flags=";".join(flags),
        ecm_path=str(record.ecm_path),
        dapi_path=str(record.dapi_path),
        cellprofiler_mask_path=_relpath(mask_path, root_for_paths),
        cellprofiler_overlay_path=_relpath(overlay_path, root_for_paths) if overlay_path.exists() else "",
        qc_panel_path=_relpath(qc_panel_path, root_for_paths) if qc_panel_path else "",
    )


def ecm_qc_flags(
    *,
    positive_area_fraction: float,
    background_source: str,
    positive_corrected: float,
) -> list[str]:
    flags: list[str] = []
    if positive_area_fraction < 0.001:
        flags.append("near_empty_ecm_positive_mask")
    if positive_area_fraction > 0.60:
        flags.append("near_full_field_ecm_positive_mask")
    if "fallback" in background_source:
        flags.append("background_estimate_fallback_used")
    if positive_corrected <= 0:
        flags.append("nonpositive_background_corrected_ecm_signal")
    return flags


def count_labels(labels: np.ndarray) -> int:
    values = np.unique(labels)
    return int(np.count_nonzero(values))


def find_cellprofiler_mask(mask_dir: Path, record: EcmImageRecord) -> Path:
    return _find_one(mask_dir, f"{record.staged_stem}_CellProfilerECMLabels*.tif*")


def find_cellprofiler_overlay(qc_dir: Path, record: EcmImageRecord) -> Path:
    matches = sorted(qc_dir.glob(f"{record.staged_stem}_CellProfilerECMOverlay*.png"))
    return matches[0] if matches else Path("")


def robust_background_threshold_preview(
    image: np.ndarray,
    *,
    deviations: float,
    lower_outlier_fraction: float = 0.05,
    upper_outlier_fraction: float = 0.05,
) -> float:
    arr = image.astype(float).ravel()
    lo = float(np.quantile(arr, lower_outlier_fraction))
    hi = float(np.quantile(arr, 1.0 - upper_outlier_fraction))
    trimmed = arr[(arr >= lo) & (arr <= hi)]
    if trimmed.size == 0:
        trimmed = arr
    return float(np.mean(trimmed) + deviations * np.std(trimmed))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_ecm_workbook(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["acquisition"]), []).append(row)
    acquisition_names = list(grouped.keys())
    sheet_names = _unique_sheet_names([*acquisition_names, "Method"])
    with ZipFile(path, "w", compression=ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", _content_types_xml(len(sheet_names)))
        zf.writestr("_rels/.rels", _root_rels_xml())
        zf.writestr("docProps/core.xml", _core_xml())
        zf.writestr("docProps/app.xml", _app_xml(sheet_names))
        zf.writestr("xl/workbook.xml", _workbook_xml(sheet_names))
        zf.writestr("xl/_rels/workbook.xml.rels", _workbook_rels_xml(len(sheet_names)))
        zf.writestr("xl/styles.xml", _styles_xml())
        for idx, sheet_name in enumerate(sheet_names, start=1):
            if idx > len(acquisition_names):
                zf.writestr(f"xl/worksheets/sheet{idx}.xml", _method_sheet_xml(rows))
            else:
                acquisition_name = acquisition_names[idx - 1]
                zf.writestr(
                    f"xl/worksheets/sheet{idx}.xml",
                    _summary_sheet_xml(grouped[acquisition_name]),
                )


def normalize_for_display(image: np.ndarray, *, lo_pct: float = 1, hi_pct: float = 99.5) -> np.ndarray:
    arr = image.astype(float)
    lo, hi = np.percentile(arr, [lo_pct, hi_pct])
    if hi <= lo:
        return np.zeros_like(arr, dtype=float)
    return np.clip((arr - lo) / (hi - lo), 0, 1)


def mask_boundary_rgba(mask: np.ndarray, *, color: tuple[float, float, float]) -> np.ndarray:
    boundaries = segmentation.find_boundaries(mask > 0, mode="outer")
    rgba = np.zeros((*boundaries.shape, 4), dtype=float)
    rgba[..., 0] = color[0]
    rgba[..., 1] = color[1]
    rgba[..., 2] = color[2]
    rgba[..., 3] = boundaries.astype(float) * 0.95
    return rgba


def _unique_slug(name: str, used: set[str]) -> str:
    base = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_") or "acquisition"
    slug = base
    idx = 2
    while slug in used:
        slug = f"{base}_{idx}"
        idx += 1
    used.add(slug)
    return slug


def _natural_key(text: str) -> list[object]:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def _find_one(root: Path, pattern: str) -> Path:
    matches = sorted(root.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one match for {pattern} in {root}; found {len(matches)}")
    return matches[0]


def _relpath(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)


def _content_types_xml(sheet_count: int) -> str:
    sheet_overrides = "\n".join(
        f'<Override PartName="/xl/worksheets/sheet{i}.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        for i in range(1, sheet_count + 1)
    )
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
{sheet_overrides}
</Types>'''


def _root_rels_xml() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>
<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>
</Relationships>'''


def _core_xml() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<dc:title>CellProfiler ECM CH1 Summary</dc:title>
<dc:creator>DAPI intensity quantification pipeline</dc:creator>
</cp:coreProperties>'''


def _app_xml(sheet_names: list[str]) -> str:
    titles = "".join(f"<vt:lpstr>{_xml_escape(name)}</vt:lpstr>" for name in sheet_names)
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
<Application>DAPI intensity quantification pipeline</Application>
<TitlesOfParts><vt:vector size="{len(sheet_names)}" baseType="lpstr">{titles}</vt:vector></TitlesOfParts>
</Properties>'''


def _workbook_xml(sheet_names: list[str]) -> str:
    sheets = "\n".join(
        f'<sheet name="{_xml_escape(_sheet_name(name))}" sheetId="{idx}" r:id="rId{idx}"/>'
        for idx, name in enumerate(sheet_names, start=1)
    )
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
<sheets>{sheets}</sheets>
</workbook>'''


def _workbook_rels_xml(sheet_count: int) -> str:
    rels = "\n".join(
        f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>'
        for i in range(1, sheet_count + 1)
    )
    rels += (
        f'\n<Relationship Id="rId{sheet_count + 1}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
        'Target="styles.xml"/>'
    )
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
{rels}
</Relationships>'''


def _styles_xml() -> str:
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts>
<fills count="1"><fill><patternFill patternType="none"/></fill></fills>
<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0"/></cellXfs>
</styleSheet>'''


def _summary_sheet_xml(rows: list[dict[str, object]]) -> str:
    headers = [
        "LOCATION",
        "ECM integrated intensity",
        "DAPI count (context only)",
        "ECM positive area fraction",
        "ECM background / px",
        "QC Status",
        "QC Flags",
    ]
    body = [headers]
    for row in rows:
        body.append(
            [
                row["location"],
                row["ecm_positive_integrated_background_corrected"],
                row["dapi_positive_nucleus_count"],
                row["ecm_positive_area_fraction"],
                row["ecm_background_value_per_px"],
                row["qc_status"],
                row["qc_flags"],
            ]
        )
    return _worksheet_xml(body)


def _method_sheet_xml(rows: list[dict[str, object]]) -> str:
    ecm_channel_id = _unique_row_value(rows, "ecm_channel_id", default="CH1")
    dapi_channel_id = _unique_row_value(rows, "dapi_channel_id", default="CH4")
    return _worksheet_xml(
        [
            ["field", "value"],
            ["ecm_channel_id", ecm_channel_id],
            ["dapi_channel_id", dapi_channel_id],
            ["mask_method", "Actual CellProfiler IdentifyPrimaryObjects Robust Background"],
            ["background_value", "Median of CellProfiler mask-negative pixels per image"],
            ["endpoint", ECM_ENDPOINT],
            ["normalization_denominator", "none"],
            ["dapi_count_role", "context/QC only; not used for ECM normalization"],
        ]
    )


def _unique_row_value(rows: list[dict[str, object]], key: str, *, default: str) -> str:
    values = sorted({str(row[key]) for row in rows if row.get(key) not in (None, "")})
    if not values:
        return default
    return values[0] if len(values) == 1 else ";".join(values)


def _worksheet_xml(rows: list[list[object]]) -> str:
    row_xml = []
    for row_idx, row in enumerate(rows, start=1):
        cells = []
        for col_idx, value in enumerate(row, start=1):
            ref = f"{_column_letter(col_idx)}{row_idx}"
            style = ' s="1"' if row_idx == 1 else ""
            if value is None or (isinstance(value, float) and math.isnan(value)):
                cells.append(f'<c r="{ref}"{style}/>')
            elif isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool):
                cells.append(f'<c r="{ref}"{style}><v>{float(value):.12g}</v></c>')
            else:
                cells.append(
                    f'<c r="{ref}" t="inlineStr"{style}><is><t>{_xml_escape(str(value))}</t></is></c>'
                )
        row_xml.append(f'<row r="{row_idx}">{"".join(cells)}</row>')
    return f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<sheetData>{"".join(row_xml)}</sheetData>
</worksheet>'''


def _column_letter(index: int) -> str:
    letters = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(65 + remainder) + letters
    return letters


def _sheet_name(name: str) -> str:
    clean = re.sub(r"[\[\]:*?/\\]", "_", name)
    return clean[:31] or "Sheet"


def _unique_sheet_names(names: Sequence[str]) -> list[str]:
    used: set[str] = set()
    unique: list[str] = []
    for name in names:
        base = _sheet_name(name)
        candidate = base
        suffix = 1
        while candidate.casefold() in used:
            suffix_text = f"_{suffix}"
            candidate = f"{base[: 31 - len(suffix_text)]}{suffix_text}"
            suffix += 1
        used.add(candidate.casefold())
        unique.append(candidate)
    return unique


def _xml_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )
