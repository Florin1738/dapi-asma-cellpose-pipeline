from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from scipy import ndimage
from scipy.optimize import linear_sum_assignment
import tifffile
import yaml

from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.seeded_regions import score_seeded_region_qc


SEEDED_REQUIRED_PLOTS = [
    "seeded_region_area_fraction.png",
    "seeded_region_intensity_per_nucleus.png",
]

SEEDED_REQUIRED_WARNING_TOKENS = {
    "exploratory_asma_associated_region_not_validated_cell_mask",
    "ch2_endpoint_signal_used_to_define_region",
}

MASK_VALIDATION_SUMMARY_COLUMNS = [
    "image_id",
    "roi_id",
    "iou_threshold",
    "n_manual",
    "n_predicted",
    "candidate_count",
    "reference_count",
    "true_positives",
    "matched_count",
    "false_positives",
    "false_positive_count",
    "false_negatives",
    "false_negative_count",
    "precision",
    "recall",
    "f1",
    "mean_iou_matched",
    "count_error",
    "count_error_percent",
]

MASK_VALIDATION_MATCH_COLUMNS = [
    "image_id",
    "candidate_label",
    "reference_label",
    "iou",
    "match_status",
    "matched",
]


def validate_seeded_region_outputs(output_dir: Path | str) -> dict[str, Any]:
    output_path = Path(output_dir)
    summary_path = output_path / "summaries" / "seeded_region_image_metrics.csv"
    config_path = output_path / "logs" / "config_resolved.yaml"
    run_log_path = output_path / "logs" / "run_log.txt"
    contact_sheet_path = output_path / "qc_contact_sheet.png"
    for required_path in [summary_path, config_path, run_log_path, contact_sheet_path]:
        if not required_path.exists():
            raise FileNotFoundError(f"Required seeded-region artifact missing: {required_path}")
    for plot_name in SEEDED_REQUIRED_PLOTS:
        plot_path = output_path / "plots" / plot_name
        if not plot_path.exists():
            raise FileNotFoundError(f"Required seeded-region plot missing: {plot_path}")

    rows = _read_csv(summary_path)
    if not rows:
        raise ValueError(f"No rows found in {summary_path}")
    config = _read_config(config_path)
    expected_source_ids = _expected_source_ids(config)
    observed_source_ids = _observed_source_ids(rows)
    missing_rows = sorted(set(expected_source_ids) - set(observed_source_ids))
    extra_rows = sorted(set(observed_source_ids) - set(expected_source_ids))
    if missing_rows:
        raise ValueError(
            "Seeded-region summary has missing rows for configured image IDs: "
            + ", ".join(missing_rows)
        )
    if extra_rows:
        raise ValueError(
            "Seeded-region summary has rows not present in config_resolved.yaml: "
            + ", ".join(extra_rows)
        )
    image_records = _image_records_by_source_id(config)

    for row in rows:
        _validate_seeded_row(row, output_path, image_records)

    validation_status = config.get("validation_status", {})
    return {
        "summary_rows": len(rows),
        "formulas_match": True,
        "masks_exist": True,
        "qc_panels_exist": True,
        "manual_ground_truth_available": bool(
            validation_status.get("manual_ground_truth_available", False)
        ),
        "whole_cell_segmentation_validated": bool(
            validation_status.get("whole_cell_segmentation_validated", False)
        ),
    }


def evaluate_instance_mask_iou(
    *,
    image_id: str,
    candidate_mask: np.ndarray,
    reference_mask: np.ndarray,
    iou_threshold: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if not 0 < iou_threshold <= 1:
        raise ValueError(f"iou_threshold must be greater than 0 and <= 1, got {iou_threshold}")
    candidate = _validate_instance_label_mask(candidate_mask, image_id=image_id, role="candidate")
    reference = _validate_instance_label_mask(reference_mask, image_id=image_id, role="reference")
    if candidate.shape != reference.shape:
        raise ValueError(
            f"Candidate and reference masks must have the same shape: "
            f"candidate={candidate.shape}, reference={reference.shape}"
        )
    candidate_labels = _positive_labels(candidate)
    reference_labels = _positive_labels(reference)
    iou_matrix = _iou_matrix(candidate, reference, candidate_labels, reference_labels)

    matches: list[dict[str, Any]] = []
    used_candidate_indices: set[int] = set()
    used_reference_indices: set[int] = set()
    if iou_matrix.size:
        candidate_indices, reference_indices = linear_sum_assignment(-iou_matrix)
        for candidate_index, reference_index in zip(candidate_indices, reference_indices, strict=True):
            iou = float(iou_matrix[candidate_index, reference_index])
            if iou < iou_threshold:
                continue
            used_candidate_indices.add(int(candidate_index))
            used_reference_indices.add(int(reference_index))
            matches.append(
                {
                    "image_id": image_id,
                    "candidate_label": int(candidate_labels[candidate_index]),
                    "reference_label": int(reference_labels[reference_index]),
                    "iou": iou,
                    "match_status": "true_positive",
                    "matched": True,
                }
            )
    for candidate_index, candidate_label in enumerate(candidate_labels):
        if candidate_index in used_candidate_indices:
            continue
        matches.append(
            {
                "image_id": image_id,
                "candidate_label": int(candidate_label),
                "reference_label": "",
                "iou": 0.0,
                "match_status": "false_positive",
                "matched": False,
            }
        )
    for reference_index, reference_label in enumerate(reference_labels):
        if reference_index in used_reference_indices:
            continue
        matches.append(
            {
                "image_id": image_id,
                "candidate_label": "",
                "reference_label": int(reference_label),
                "iou": 0.0,
                "match_status": "false_negative",
                "matched": False,
            }
        )

    matched_count = len(used_candidate_indices)
    candidate_count = len(candidate_labels)
    reference_count = len(reference_labels)
    false_positive_count = candidate_count - len(used_candidate_indices)
    false_negative_count = reference_count - len(used_reference_indices)
    precision = _safe_divide(matched_count, candidate_count)
    recall = _safe_divide(matched_count, reference_count)
    f1 = _f1(precision, recall)
    matched_ious = [match["iou"] for match in matches if match["matched"]]
    mean_iou = float(np.mean(matched_ious)) if matched_ious else float("nan")
    count_error = candidate_count - reference_count
    count_error_percent = _safe_divide(count_error, reference_count) * 100

    return (
        {
            "image_id": image_id,
            "roi_id": "full_image",
            "iou_threshold": float(iou_threshold),
            "n_manual": reference_count,
            "n_predicted": candidate_count,
            "candidate_count": candidate_count,
            "reference_count": reference_count,
            "true_positives": matched_count,
            "matched_count": matched_count,
            "false_positives": false_positive_count,
            "false_positive_count": false_positive_count,
            "false_negatives": false_negative_count,
            "false_negative_count": false_negative_count,
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "mean_iou_matched": mean_iou,
            "count_error": count_error,
            "count_error_percent": count_error_percent,
        },
        matches,
    )


def write_instance_mask_validation_tables(
    *,
    output_dir: Path | str,
    summaries: list[dict[str, Any]],
    matches: list[dict[str, Any]],
) -> None:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    _write_csv(output_path / "manual_mask_validation_summary.csv", summaries, MASK_VALIDATION_SUMMARY_COLUMNS)
    _write_csv(output_path / "manual_mask_validation_matches.csv", matches, MASK_VALIDATION_MATCH_COLUMNS)


def run_manual_mask_validation(
    *,
    candidate_mask_paths: dict[str, Path],
    reference_mask_paths: dict[str, Path],
    output_dir: Path | str,
    iou_threshold: float,
    reference_completion_status: dict[str, str] | None = None,
    reference_completion_mask_paths: dict[str, Path] | None = None,
    run_metadata: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    if not reference_mask_paths:
        raise ValueError("At least one reference/manual mask is required for validation")
    candidate_without_reference = sorted(set(candidate_mask_paths) - set(reference_mask_paths))
    reference_without_candidate = sorted(set(reference_mask_paths) - set(candidate_mask_paths))
    reference_masks = {
        image_id: _read_instance_label_mask(path, image_id=image_id, role="reference")
        for image_id, path in sorted(reference_mask_paths.items())
    }
    if reference_completion_status is not None:
        _validate_reference_completion_status(
            status_by_image=reference_completion_status,
            reference_masks=reference_masks,
            reference_mask_paths=reference_mask_paths,
            status_mask_paths=reference_completion_mask_paths,
        )
    if not any(np.any(reference > 0) for reference in reference_masks.values()):
        raise ValueError(
            "all reference masks are empty; fill manual/reference labels before validation"
        )
    summaries: list[dict[str, Any]] = []
    all_matches: list[dict[str, Any]] = []
    for image_id, reference in reference_masks.items():
        candidate_path = candidate_mask_paths.get(image_id)
        if candidate_path is None:
            candidate = np.zeros(reference.shape, dtype=np.uint32)
        else:
            candidate = _read_instance_label_mask(candidate_path, image_id=image_id, role="candidate")
        summary, matches = evaluate_instance_mask_iou(
            image_id=image_id,
            candidate_mask=candidate,
            reference_mask=reference,
            iou_threshold=iou_threshold,
        )
        summaries.append(summary)
        all_matches.extend(matches)
    write_instance_mask_validation_tables(
        output_dir=output_dir,
        summaries=summaries,
        matches=all_matches,
    )
    _write_manual_validation_provenance(
        output_dir=output_dir,
        candidate_mask_paths=candidate_mask_paths,
        reference_mask_paths=reference_mask_paths,
        evaluated_image_ids=[summary["image_id"] for summary in summaries],
        candidate_without_reference=candidate_without_reference,
        reference_without_candidate=reference_without_candidate,
        iou_threshold=iou_threshold,
        run_metadata=run_metadata or {},
    )
    return summaries


def _validate_reference_completion_status(
    *,
    status_by_image: dict[str, str],
    reference_masks: dict[str, np.ndarray],
    reference_mask_paths: dict[str, Path],
    status_mask_paths: dict[str, Path] | None,
) -> None:
    required_statuses = {"complete_non_empty", "confirmed_empty"}
    missing = sorted(set(reference_masks) - set(status_by_image))
    if missing:
        raise ValueError(
            "completion status missing for reference masks: " + ", ".join(missing)
        )
    if status_mask_paths is not None:
        missing_paths = sorted(set(reference_masks) - set(status_mask_paths))
        if missing_paths:
            raise ValueError(
                "completion status mask path missing for reference masks: "
                + ", ".join(missing_paths)
            )
        mismatched_paths = []
        for image_id, reference_path in sorted(reference_mask_paths.items()):
            status_path = status_mask_paths[image_id]
            if _normalized_path(status_path) != _normalized_path(reference_path):
                mismatched_paths.append(
                    f"{image_id}: status={status_path} reference={reference_path}"
                )
        if mismatched_paths:
            raise ValueError(
                "completion status mask path does not match reference mask: "
                + "; ".join(mismatched_paths)
            )

    incomplete = []
    empty_not_confirmed = []
    non_empty_not_complete = []
    for image_id, reference in sorted(reference_masks.items()):
        status = status_by_image[image_id].strip().lower()
        if status not in required_statuses:
            incomplete.append(f"{image_id}={status or 'blank'}")
            continue
        has_reference_objects = bool(np.any(reference > 0))
        if has_reference_objects and status != "complete_non_empty":
            non_empty_not_complete.append(image_id)
        if not has_reference_objects and status != "confirmed_empty":
            empty_not_confirmed.append(image_id)

    if incomplete:
        raise ValueError(
            "reference masks are not marked complete: " + ", ".join(incomplete)
        )
    if non_empty_not_complete:
        raise ValueError(
            "non-empty reference masks must be marked complete_non_empty: "
            + ", ".join(non_empty_not_complete)
        )
    if empty_not_confirmed:
        raise ValueError(
            "empty reference masks must be marked confirmed_empty: "
            + ", ".join(empty_not_confirmed)
        )


def _normalized_path(path: Path) -> Path:
    return path.expanduser().resolve()


def _validate_seeded_row(
    row: dict[str, str],
    output_path: Path,
    image_records: dict[str, dict[str, Any]],
) -> None:
    source_id = row["source_id"]
    record = image_records.get(source_id.upper())
    if record is None:
        raise ValueError(f"No image record for {source_id} in config_resolved.yaml")

    ch2_path = _resolve_reference(record["ch2_path"], output_path)
    nuclei_mask_path = _resolve_reference(record["nuclei_mask_path"], output_path)
    mask_path = _resolve_reference(row["mask_path"], output_path)
    qc_path = _resolve_reference(row["qc_panel_path"], output_path)
    for path in [ch2_path, nuclei_mask_path, mask_path, qc_path]:
        if not path.exists():
            raise FileNotFoundError(f"{source_id} referenced artifact missing: {path}")

    ch2, _page_count = read_primary_intensity_plane(ch2_path)
    nuclei_mask = np.asarray(tifffile.imread(nuclei_mask_path))
    seeded_mask = np.asarray(tifffile.imread(mask_path))
    if ch2.shape != nuclei_mask.shape or ch2.shape != seeded_mask.shape:
        raise ValueError(
            f"{source_id} shape mismatch: ch2={ch2.shape}, nuclei={nuclei_mask.shape}, "
            f"seeded={seeded_mask.shape}"
        )
    _validate_seeded_label_integrity(source_id, seeded_mask, nuclei_mask)

    ch2_float = ch2.astype(np.float64)
    region_mask = seeded_mask > 0
    image_area = int(ch2.size)
    seeded_area = int(np.count_nonzero(region_mask))
    raw_integrated = float(np.sum(ch2_float[region_mask]))
    mean_raw = float(np.mean(ch2_float[region_mask])) if seeded_area else float("nan")
    background_value = float(row["background_value_per_px"])
    corrected = np.clip(ch2_float[region_mask] - background_value, 0, None)
    corrected_integrated = float(np.sum(corrected))
    nucleus_count = _count_nonzero_labels(nuclei_mask)
    per_nucleus = corrected_integrated / nucleus_count if nucleus_count else float("nan")

    _assert_int(source_id, "image_area_px", image_area, row["image_area_px"])
    _assert_int(source_id, "seeded_region_area_px", seeded_area, row["seeded_region_area_px"])
    _assert_int(
        source_id,
        "non_seeded_area_px",
        image_area - seeded_area,
        row["non_seeded_area_px"],
    )
    _assert_int(
        source_id,
        "dapi_positive_nucleus_count",
        nucleus_count,
        row["dapi_positive_nucleus_count"],
    )
    _assert_close(source_id, "seeded_region_integrated_raw", raw_integrated, row)
    _assert_close(source_id, "seeded_region_mean_raw", mean_raw, row)
    _assert_close(
        source_id,
        "seeded_region_integrated_background_corrected",
        corrected_integrated,
        row,
    )
    _assert_close(
        source_id,
        "seeded_region_intensity_per_DAPI_positive_nucleus",
        per_nucleus,
        row,
    )

    expected_status, expected_flags = score_seeded_region_qc(row)
    if row["qc_status"] != expected_status:
        raise ValueError(
            f"{source_id} qc_status mismatch: expected {expected_status}, observed {row['qc_status']}"
        )
    if row["qc_flags"] != ";".join(expected_flags):
        raise ValueError(
            f"{source_id} qc_flags mismatch: expected {';'.join(expected_flags)}, "
            f"observed {row['qc_flags']}"
        )
    _validate_seeded_warning_contract(source_id, row)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_config(config_path: Path) -> dict[str, Any]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ValueError(f"Expected mapping config in {config_path}")
    return config


def _validate_seeded_label_integrity(
    source_id: str,
    seeded_mask: np.ndarray,
    nuclei_mask: np.ndarray,
) -> None:
    if seeded_mask.ndim != 2:
        raise ValueError(f"{source_id} seeded mask must be 2-D; got shape {seeded_mask.shape}")
    if not np.issubdtype(seeded_mask.dtype, np.integer):
        raise ValueError(f"{source_id} seeded mask must use integer labels; got {seeded_mask.dtype}")
    if np.any(seeded_mask < 0):
        raise ValueError(f"{source_id} seeded mask contains negative labels")
    seeded_labels = set(int(label) for label in _positive_labels(seeded_mask))
    nuclei_labels = set(int(label) for label in _positive_labels(nuclei_mask))
    missing_labels = sorted(seeded_labels - nuclei_labels)
    if missing_labels:
        raise ValueError(
            f"{source_id} seeded mask labels not present in DAPI nuclei mask: "
            + ", ".join(str(label) for label in missing_labels)
        )


def _validate_seeded_warning_contract(source_id: str, row: dict[str, str]) -> None:
    warning_tokens = set(token for token in row.get("warnings", "").split(";") if token)
    missing = sorted(SEEDED_REQUIRED_WARNING_TOKENS - warning_tokens)
    if missing:
        raise ValueError(
            f"{source_id} missing required warning tokens: " + ", ".join(missing)
        )


def _expected_source_ids(config: dict[str, Any]) -> list[str]:
    records = config.get("image_records")
    if not isinstance(records, list) or not records:
        raise ValueError("No image_records found in config_resolved.yaml")
    expected: list[str] = []
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("Each image_records entry must be a mapping")
        key = record.get("source_id") or record.get("location")
        if not isinstance(key, str) or not key:
            raise ValueError("Each image_records entry must include source_id or location")
        expected.append(key.upper())
    duplicate_ids = sorted({image_id for image_id in expected if expected.count(image_id) > 1})
    if duplicate_ids:
        raise ValueError("Duplicate image IDs in config_resolved.yaml: " + ", ".join(duplicate_ids))
    images_processed = config.get("images_processed")
    if images_processed is not None and int(images_processed) != len(expected):
        raise ValueError(
            f"images_processed mismatch: expected {len(expected)} image_records, "
            f"observed images_processed={images_processed}"
        )
    return expected


def _observed_source_ids(rows: list[dict[str, str]]) -> list[str]:
    observed = [row["source_id"].upper() for row in rows]
    duplicate_ids = sorted({image_id for image_id in observed if observed.count(image_id) > 1})
    if duplicate_ids:
        raise ValueError("Duplicate rows in seeded-region summary: " + ", ".join(duplicate_ids))
    return observed


def _image_records_by_source_id(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    records = config.get("image_records")
    if not isinstance(records, list) or not records:
        raise ValueError("No image_records found in config_resolved.yaml")
    output: dict[str, dict[str, Any]] = {}
    for record in records:
        if not isinstance(record, dict):
            continue
        for key in [record.get("source_id"), record.get("location")]:
            if isinstance(key, str) and key:
                output[key.upper()] = record
    if not output:
        raise ValueError("No usable image_records found in config_resolved.yaml")
    return output


def _read_instance_label_mask(path: Path, *, image_id: str, role: str) -> np.ndarray:
    return _validate_instance_label_mask(
        np.asarray(tifffile.imread(path)),
        image_id=image_id,
        role=role,
    )


def _validate_instance_label_mask(mask: np.ndarray, *, image_id: str, role: str) -> np.ndarray:
    array = np.asarray(mask)
    if array.ndim != 2:
        raise ValueError(f"{image_id} {role} mask must be a 2-D instance label image; got shape {array.shape}")
    if not (np.issubdtype(array.dtype, np.integer) or np.issubdtype(array.dtype, np.bool_)):
        raise ValueError(f"{image_id} {role} mask must use an integer instance-label dtype; got {array.dtype}")
    if np.any(array < 0):
        raise ValueError(f"{image_id} {role} mask contains negative labels")
    labels = _positive_labels(array)
    if len(labels) == 1 and _connected_component_count(array > 0) > 1:
        raise ValueError(
            f"{image_id} {role} mask appears to be a binary or single-label mask with "
            "multiple disconnected objects; use separate instance labels"
        )
    for label in labels:
        if _connected_component_count(array == label) > 1:
            raise ValueError(
                f"{image_id} {role} mask label {int(label)} has multiple disconnected components"
            )
    return array


def _connected_component_count(binary_mask: np.ndarray) -> int:
    return int(ndimage.label(binary_mask)[1])


def _positive_labels(mask: np.ndarray) -> np.ndarray:
    labels = np.unique(mask)
    return labels[labels > 0].astype(np.int64, copy=False)


def _iou_matrix(
    candidate: np.ndarray,
    reference: np.ndarray,
    candidate_labels: np.ndarray,
    reference_labels: np.ndarray,
) -> np.ndarray:
    if len(candidate_labels) == 0 or len(reference_labels) == 0:
        return np.zeros((len(candidate_labels), len(reference_labels)), dtype=np.float64)
    matrix = np.zeros((len(candidate_labels), len(reference_labels)), dtype=np.float64)
    candidate_masks = [candidate == label for label in candidate_labels]
    reference_masks = [reference == label for label in reference_labels]
    for candidate_index, candidate_binary in enumerate(candidate_masks):
        candidate_area = int(np.count_nonzero(candidate_binary))
        for reference_index, reference_binary in enumerate(reference_masks):
            intersection = int(np.count_nonzero(candidate_binary & reference_binary))
            if intersection == 0:
                continue
            reference_area = int(np.count_nonzero(reference_binary))
            union = candidate_area + reference_area - intersection
            matrix[candidate_index, reference_index] = intersection / union
    return matrix


def _assert_int(image_id: str, label: str, expected: int, observed_text: str) -> None:
    observed = int(observed_text)
    if observed != expected:
        raise ValueError(f"{image_id} {label} mismatch: expected {expected}, observed {observed}")


def _assert_close(image_id: str, label: str, expected: float, row: dict[str, str]) -> None:
    observed = float(row[label])
    if np.isnan(expected) and np.isnan(observed):
        return
    if not np.isclose(observed, expected, rtol=1e-9, atol=1e-6):
        raise ValueError(f"{image_id} {label} mismatch: expected {expected}, observed {observed}")


def _count_nonzero_labels(mask: np.ndarray) -> int:
    return int(np.count_nonzero(np.unique(mask)))


def _safe_divide(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return float("nan")
    return float(numerator / denominator)


def _f1(precision: float, recall: float) -> float:
    if not np.isfinite(precision) or not np.isfinite(recall):
        return float("nan")
    if precision + recall == 0:
        return 0.0
    return float(2 * precision * recall / (precision + recall))


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _resolve_reference(value: str, output_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    candidates = [
        path,
        output_path / path,
        output_path.parent / path,
        output_path.parent.parent / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[-1]


def _write_manual_validation_provenance(
    *,
    output_dir: Path | str,
    candidate_mask_paths: dict[str, Path],
    reference_mask_paths: dict[str, Path],
    evaluated_image_ids: list[str],
    candidate_without_reference: list[str],
    reference_without_candidate: list[str],
    iou_threshold: float,
    run_metadata: dict[str, Any],
) -> None:
    output_path = Path(output_dir)
    logs_path = output_path / "logs"
    logs_path.mkdir(parents=True, exist_ok=True)
    config = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "output_dir": str(output_path),
        "iou_threshold": float(iou_threshold),
        "matching_rule": "candidate and reference objects match when IoU >= iou_threshold",
        "evaluated_image_ids": evaluated_image_ids,
        "reference_without_candidate_counted_as_false_negative": reference_without_candidate,
        "candidate_without_reference_not_evaluated": candidate_without_reference,
        "candidate_masks": {key: str(path) for key, path in sorted(candidate_mask_paths.items())},
        "reference_masks": {key: str(path) for key, path in sorted(reference_mask_paths.items())},
        "run_metadata": run_metadata,
    }
    (logs_path / "config_resolved.yaml").write_text(
        yaml.safe_dump(config, sort_keys=False),
        encoding="utf-8",
    )
    run_log_lines = [
        "Manual instance-mask validation",
        f"generated_at_utc: {config['generated_at_utc']}",
        f"iou_threshold: {iou_threshold}",
        f"evaluated_image_count: {len(evaluated_image_ids)}",
        "reference_without_candidate_counted_as_false_negative: "
        + ", ".join(reference_without_candidate),
        "candidate_without_reference_not_evaluated: " + ", ".join(candidate_without_reference),
    ]
    (logs_path / "run_log.txt").write_text("\n".join(run_log_lines) + "\n", encoding="utf-8")
