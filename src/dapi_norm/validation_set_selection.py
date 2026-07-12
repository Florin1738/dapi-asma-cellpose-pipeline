from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from skimage import segmentation


VALIDATION_FEATURE_COLUMNS = [
    "image_id",
    "source_id",
    "ch2_path",
    "ch4_path",
    "nuclei_mask_path",
    "dapi_positive_nucleus_count",
    "target_integrated_raw",
    "target_integrated_raw_per_DAPI_positive_nucleus",
    "target_saturation_fraction",
    "dapi_saturation_fraction",
    "method_region_jaccard",
    "selection_reasons",
]


@dataclass
class ValidationFeatureRecord:
    image_id: str
    source_id: str
    ch2_path: Path
    ch4_path: Path
    nuclei_mask_path: Path
    dapi_positive_nucleus_count: int
    target_integrated_raw: float
    target_integrated_raw_per_DAPI_positive_nucleus: float
    target_saturation_fraction: float
    dapi_saturation_fraction: float
    method_region_jaccard: float | None = None
    selection_reasons: str = ""


def select_validation_set(
    records: Iterable[ValidationFeatureRecord],
    *,
    max_images: int = 16,
    per_bucket: int = 2,
    must_include: Iterable[str] = (),
) -> list[ValidationFeatureRecord]:
    record_list = sorted(list(records), key=_image_sort_key)
    if not record_list:
        raise ValueError("Cannot select a validation set from zero records")
    if max_images <= 0:
        raise ValueError(f"max_images must be > 0, got {max_images}")
    if per_bucket <= 0:
        raise ValueError(f"per_bucket must be > 0, got {per_bucket}")

    by_id = {record.image_id.upper(): record for record in record_list}
    selected: dict[str, ValidationFeatureRecord] = {}
    reason_map: dict[str, list[str]] = {}

    def add(record: ValidationFeatureRecord, reason: str) -> None:
        key = record.image_id.upper()
        if key not in selected and len(selected) >= max_images:
            return
        selected.setdefault(key, _copy_record(record))
        reasons = reason_map.setdefault(key, [])
        if reason not in reasons:
            reasons.append(reason)

    for image_id in must_include:
        record = by_id.get(image_id.strip().upper())
        if record is not None:
            add(record, "must_include")

    buckets = [
        ("low_raw_target_integrated", _sorted_by(record_list, "target_integrated_raw", reverse=False)[:per_bucket]),
        ("high_raw_target_integrated", _sorted_by(record_list, "target_integrated_raw", reverse=True)[:per_bucket]),
        ("mid_raw_target_integrated", _middle_by(record_list, "target_integrated_raw", per_bucket)),
        (
            "low_target_integrated_per_nucleus",
            _sorted_by(
                record_list,
                "target_integrated_raw_per_DAPI_positive_nucleus",
                reverse=False,
            )[:per_bucket],
        ),
        (
            "high_target_integrated_per_nucleus",
            _sorted_by(
                record_list,
                "target_integrated_raw_per_DAPI_positive_nucleus",
                reverse=True,
            )[:per_bucket],
        ),
        (
            "low_dapi_positive_nucleus_count",
            _sorted_by(record_list, "dapi_positive_nucleus_count", reverse=False)[:per_bucket],
        ),
        (
            "high_dapi_positive_nucleus_count",
            _sorted_by(record_list, "dapi_positive_nucleus_count", reverse=True)[:per_bucket],
        ),
        (
            "high_target_saturation_fraction",
            _sorted_by(record_list, "target_saturation_fraction", reverse=True)[:per_bucket],
        ),
        (
            "high_method_disagreement",
            _method_disagreement_records(record_list)[:per_bucket],
        ),
    ]
    for reason, bucket_records in buckets:
        for record in bucket_records:
            add(record, reason)

    selected_rows = list(selected.values())
    for row in selected_rows:
        row.selection_reasons = ";".join(reason_map[row.image_id.upper()])
    return selected_rows


def write_selection_csvs(
    *,
    all_records: Iterable[ValidationFeatureRecord],
    selected_records: Iterable[ValidationFeatureRecord],
    output_dir: Path,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    all_path = output_dir / "validation_candidate_field_features.csv"
    selected_path = output_dir / "selected_manual_validation_fields.csv"
    _write_csv(all_path, list(all_records))
    _write_csv(selected_path, list(selected_records))
    return {"all_features": all_path, "selected": selected_path}


def render_validation_selection_panel(
    *,
    selected_records: Iterable[ValidationFeatureRecord],
    output_path: Path,
    image_loader: Callable[[Path], np.ndarray],
    nuclei_loader: Callable[[Path], np.ndarray],
) -> None:
    records = list(selected_records)
    if not records:
        raise ValueError("Cannot render validation selection panel without records")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(records),
        3,
        figsize=(12.2, max(3.0, 2.75 * len(records))),
        squeeze=False,
        constrained_layout=True,
        width_ratios=[1.0, 1.0, 0.9],
    )
    for row_index, record in enumerate(records):
        ch2 = image_loader(record.ch2_path)
        ch4 = image_loader(record.ch4_path)
        nuclei = nuclei_loader(record.nuclei_mask_path)
        if ch2.shape != ch4.shape or ch2.shape != nuclei.shape:
            raise ValueError(
                f"{record.image_id} shape mismatch for panel: "
                f"ch2={ch2.shape}, ch4={ch4.shape}, nuclei={nuclei.shape}"
            )
        ch2_display = _scale_for_display(ch2)
        ch4_display = _scale_for_display(ch4)
        ch4_rgb = np.dstack([ch4_display, ch4_display, ch4_display])
        _draw_nucleus_crosses(ch4_rgb, nuclei)
        axes[row_index, 0].imshow(ch2_display, cmap="gray")
        axes[row_index, 0].set_title(f"{record.image_id} CH2/aSMA", fontsize=8)
        axes[row_index, 1].imshow(ch4_rgb)
        axes[row_index, 1].set_title("CH4/DAPI nuclei", fontsize=8)
        axes[row_index, 2].text(
            0.0,
            1.0,
            _record_summary(record),
            ha="left",
            va="top",
            fontsize=7.4,
            transform=axes[row_index, 2].transAxes,
            color="#1f2933",
        )
        for col in range(3):
            axes[row_index, col].axis("off")
    fig.suptitle(
        "Recommended manual-validation fields: stratified by intensity, nuclei, and QC risk",
        fontsize=12,
    )
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def _write_csv(path: Path, rows: list[ValidationFeatureRecord]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=VALIDATION_FEATURE_COLUMNS)
        writer.writeheader()
        for row in rows:
            serialized = asdict(row)
            serialized["ch2_path"] = str(row.ch2_path)
            serialized["ch4_path"] = str(row.ch4_path)
            serialized["nuclei_mask_path"] = str(row.nuclei_mask_path)
            writer.writerow(
                {column: serialized.get(column, "") for column in VALIDATION_FEATURE_COLUMNS}
            )


def _copy_record(record: ValidationFeatureRecord) -> ValidationFeatureRecord:
    copied = ValidationFeatureRecord(**asdict(record))
    copied.ch2_path = Path(copied.ch2_path)
    copied.ch4_path = Path(copied.ch4_path)
    copied.nuclei_mask_path = Path(copied.nuclei_mask_path)
    return copied


def _sorted_by(
    records: list[ValidationFeatureRecord],
    field_name: str,
    *,
    reverse: bool,
) -> list[ValidationFeatureRecord]:
    values = [record for record in records if np.isfinite(float(getattr(record, field_name)))]
    return sorted(
        values,
        key=lambda record: (float(getattr(record, field_name)), _image_sort_key(record)),
        reverse=reverse,
    )


def _middle_by(
    records: list[ValidationFeatureRecord],
    field_name: str,
    count: int,
) -> list[ValidationFeatureRecord]:
    values = _sorted_by(records, field_name, reverse=False)
    if not values:
        return []
    numeric = np.array([float(getattr(record, field_name)) for record in values], dtype=np.float64)
    median = float(np.median(numeric))
    return sorted(values, key=lambda record: (abs(float(getattr(record, field_name)) - median), _image_sort_key(record)))[
        :count
    ]


def _method_disagreement_records(records: list[ValidationFeatureRecord]) -> list[ValidationFeatureRecord]:
    values = [
        record
        for record in records
        if record.method_region_jaccard is not None and np.isfinite(float(record.method_region_jaccard))
    ]
    return sorted(values, key=lambda record: (float(record.method_region_jaccard), _image_sort_key(record)))


def _record_summary(record: ValidationFeatureRecord) -> str:
    return (
        f"raw CH2: {_format_scientific(record.target_integrated_raw)}\n"
        f"raw CH2 / DAPI-positive nucleus: "
        f"{_format_scientific(record.target_integrated_raw_per_DAPI_positive_nucleus)}\n"
        f"DAPI-positive nuclei: {record.dapi_positive_nucleus_count}\n"
        f"CH2 saturation: {record.target_saturation_fraction:.3%}\n"
        f"method Jaccard: {_format_optional(record.method_region_jaccard)}\n"
        f"reasons: {record.selection_reasons or 'not selected'}"
    )


def _draw_nucleus_crosses(rgb: np.ndarray, nuclei_mask: np.ndarray) -> None:
    nuclei = np.asarray(nuclei_mask)
    boundaries = segmentation.find_boundaries(nuclei, mode="outer")
    rgb[boundaries] = [0.0, 0.75, 1.0]
    size = max(2, min(7, min(rgb.shape[:2]) // 120))
    for label in np.unique(nuclei):
        if label == 0:
            continue
        ys, xs = np.nonzero(nuclei == label)
        if len(ys) == 0:
            continue
        y = int(round(float(np.mean(ys))))
        x = int(round(float(np.mean(xs))))
        for offset in range(-size, size + 1):
            for y_pos, x_pos in [(y + offset, x + offset), (y + offset, x - offset)]:
                if 0 <= y_pos < rgb.shape[0] and 0 <= x_pos < rgb.shape[1]:
                    rgb[y_pos, x_pos] = [0.0, 1.0, 0.0]


def _scale_for_display(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image, dtype=np.float32)
    if arr.size == 0:
        return arr
    low, high = np.percentile(arr, [1, 99.8])
    if high <= low:
        low = float(np.min(arr))
        high = float(np.max(arr))
    if high <= low:
        return np.zeros_like(arr, dtype=np.float32)
    return np.clip((arr - low) / (high - low), 0, 1)


def _format_scientific(value: float) -> str:
    if not np.isfinite(float(value)):
        return "NA"
    return f"{float(value):.3e}"


def _format_optional(value: float | None) -> str:
    if value is None or not np.isfinite(float(value)):
        return "NA"
    return f"{float(value):.3f}"


def _image_sort_key(record: ValidationFeatureRecord) -> tuple[str, int]:
    image_id = record.image_id.upper()
    digits = "".join(char for char in image_id if char.isdigit())
    return (image_id.rstrip("0123456789"), int(digits) if digits else -1)
