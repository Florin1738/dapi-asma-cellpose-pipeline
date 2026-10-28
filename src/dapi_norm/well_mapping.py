"""Microscope-marker identity; never infer wells from acquisition order."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

WELL_PATTERN = re.compile(r"^_([A-H])(0?[1-9]|1[0-2])$", re.I)
FIELD_PATTERN = re.compile(r"^XY\d+$", re.I)
TREATMENT_COLUMNS = ("compound", "batch", "dose", "condition", "tgf_beta", "serum", "cells_seeded")
MAPPING_COLUMNS = ("assay_folder", "field_id", "well_id", "field_in_well", "fields_in_well",
                   "well_marker", "well_marker_sha256", "treatment_mapping_status",
                   "plate_map_source", "plate_map_sha256", *TREATMENT_COLUMNS)


def load_acquisition_mapping(acquisition: Path, fields, plate_map: Path | None = None):
    acquisition = Path(acquisition).resolve()
    field_ids = list(fields)
    if len(field_ids) != len(set(field_ids)):
        raise ValueError(f"Duplicate image field IDs in {acquisition}")
    records = {}
    for field in field_ids:
        if not FIELD_PATTERN.fullmatch(field):
            raise ValueError(f"Invalid image field ID: {field}")
        folder = acquisition / field
        markers = [p for p in folder.iterdir() if p.is_file() and WELL_PATTERN.fullmatch(p.name)]
        if len(markers) != 1:
            raise ValueError(f"Expected exactly one well marker in {folder}; found {len(markers)}. "
                             "Restore microscope _A01-style markers; do not guess from XY order.")
        marker = markers[0]
        match = WELL_PATTERN.fullmatch(marker.name)
        records[field] = dict(assay_folder=acquisition.name, field_id=field,
            well_id=f"{match[1].upper()}{int(match[2]):02}", well_marker=f"{field}/{marker.name}",
            well_marker_sha256=hashlib.sha256(marker.read_bytes()).hexdigest(),
            treatment_mapping_status="not_provided", plate_map_source="", plate_map_sha256="",
            **{name: "" for name in TREATMENT_COLUMNS})
    counts = Counter(record["well_id"] for record in records.values())
    indices = Counter()
    for field in sorted(records, key=lambda value: int(value[2:])):
        row = records[field]
        indices[row["well_id"]] += 1
        row["field_in_well"] = indices[row["well_id"]]
        row["fields_in_well"] = counts[row["well_id"]]
    map_path = Path(plate_map).resolve() if plate_map else acquisition / "plate_map.csv"
    if plate_map or map_path.exists():
        treatments = {}
        with map_path.open(newline="", encoding="utf-8-sig") as file:
            reader = csv.DictReader(file)
            required = {"assay_folder", "well_id", *TREATMENT_COLUMNS}
            if not required <= set(reader.fieldnames or []):
                raise ValueError(f"Plate map requires columns {sorted(required)}: {map_path}")
            for treatment in reader:
                if treatment["assay_folder"] != acquisition.name:
                    continue
                well = treatment["well_id"]
                if not re.fullmatch(r"[A-H](0[1-9]|1[0-2])", well):
                    raise ValueError(f"Invalid plate-map well: {well}")
                if well in treatments:
                    raise ValueError(f"Duplicate plate-map well {acquisition.name}/{well}")
                if any(not (treatment[name] or '').strip() for name in TREATMENT_COLUMNS):
                    raise ValueError(f"Incomplete treatment metadata for {acquisition.name}/{well}; "
                                     "use an explicit n/a when appropriate")
                treatments[well] = treatment
        digest = hashlib.sha256(map_path.read_bytes()).hexdigest()
        for record in records.values():
            if record["well_id"] not in treatments:
                raise ValueError(f"Plate map missing {acquisition.name}/{record['well_id']}")
            record.update({name: treatments[record["well_id"]][name] for name in TREATMENT_COLUMNS})
            record.update(treatment_mapping_status="provided", plate_map_source=map_path.name,
                          plate_map_sha256=digest)
    return records


def validate_mapping_rows(rows):
    """Reject partial or internally inconsistent provenance at export boundaries."""
    for row in rows:
        if any(name not in row for name in MAPPING_COLUMNS):
            raise ValueError('Verified well mapping with complete provenance is required for workbook export')
        field, well = str(row['field_id']), str(row['well_id'])
        marker = Path(str(row['well_marker']))
        match = WELL_PATTERN.fullmatch(marker.name)
        if (not FIELD_PATTERN.fullmatch(field) or marker.parent.as_posix() != field
                or not match or well != f'{match[1].upper()}{int(match[2]):02}'
                or str(row.get('source_id', '')).rsplit('/', 1)[-1].upper() != field.upper()
                or not re.fullmatch(r'[0-9a-f]{64}', str(row['well_marker_sha256']))
                or not str(row['assay_folder']).strip()):
            raise ValueError('Invalid or inconsistent verified well mapping provenance')
        try:
            index, count = int(row['field_in_well']), int(row['fields_in_well'])
        except (ValueError, TypeError) as exc:
            raise ValueError('Invalid repeated-field mapping counts') from exc
        if not 1 <= index <= count:
            raise ValueError('Invalid repeated-field mapping counts')
        if row['treatment_mapping_status'] == 'provided':
            if (any(not str(row[c] or '').strip() for c in TREATMENT_COLUMNS)
                    or not row['plate_map_source']
                    or not re.fullmatch(r'[0-9a-f]{64}', str(row['plate_map_sha256']))):
                raise ValueError('Provided treatment mapping requires complete plate-map provenance')
        elif row['treatment_mapping_status'] != 'not_provided' or any(
                row[c] for c in (*TREATMENT_COLUMNS, 'plate_map_source', 'plate_map_sha256')):
            raise ValueError('Unverified treatment labels are forbidden')


def mapping_for_acquisitions(acquisitions, plate_map=None):
    acquisitions = list(acquisitions)
    if plate_map and len({Path(a.input_root).name for a in acquisitions}) != len(acquisitions):
        raise ValueError("Shared plate map has ambiguous duplicate assay folder names; use "
                         "a separate plate_map.csv inside each acquisition")
    result = {}
    for acquisition in acquisitions:
        path = Path(acquisition.input_root)
        fields = [p.name for p in path.iterdir() if p.is_dir() and FIELD_PATTERN.fullmatch(p.name)
                  and any(f.suffix.lower() in (".tif", ".tiff") for f in p.iterdir() if f.is_file())]
        result[str(path.resolve())] = load_acquisition_mapping(path, fields, plate_map)
    return result


def attach_mapping(rows, mappings):
    for row in rows:
        target = Path(row.get("target_path", row.get("ch2_path", ""))).resolve()
        record = mappings[str(target.parent.parent)][target.parent.name]
        row.update(record)
    return rows


def write_mapping_artifacts(directory, mappings):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    rows = [row for mapping in mappings.values() for row in mapping.values()]
    with (directory / "image_to_well.csv").open("w", newline="", encoding="utf-8-sig") as file:
        writer = csv.DictWriter(file, fieldnames=MAPPING_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    by_assay = defaultdict(list)
    for row in rows:
        by_assay[row["assay_folder"]].append(row)
    all_wells = {f"{r}{c:02}" for r in "ABCDEFGH" for c in range(1, 13)}
    summary = {}
    for assay, records in by_assay.items():
        counts = Counter(row["well_id"] for row in records)
        summary[assay] = dict(fields=len(records), unique_wells=len(counts),
            repeated_wells={w: n for w, n in counts.items() if n > 1},
            unrepresented_wells=sorted(all_wells - set(counts)),
            treatment_mapping_status=sorted({row["treatment_mapping_status"] for row in records}))
    (directory / "mapping_validation.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
