from __future__ import annotations

import csv
from pathlib import Path
import re
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from skimage import segmentation
import tifffile

from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.segmentation_validation import run_manual_mask_validation


METHOD_SUMMARY_COLUMNS = [
    "candidate_method",
    "n_images",
    "n_manual",
    "n_predicted",
    "true_positives",
    "false_positives",
    "false_negatives",
    "micro_precision",
    "micro_recall",
    "micro_f1",
    "mean_iou_matched",
    "mean_abs_count_error",
    "max_abs_count_error",
    "passes_acceptance_criteria",
    "acceptance_criteria",
    "validation_scope",
]

IMAGE_SUMMARY_COLUMNS = [
    "candidate_method",
    "image_id",
    "iou_threshold",
    "n_manual",
    "n_predicted",
    "true_positives",
    "false_positives",
    "false_negatives",
    "precision",
    "recall",
    "f1",
    "mean_iou_matched",
    "count_error",
    "count_error_percent",
    "overlay_path",
]


def run_manual_validation_report(
    *,
    candidate_dirs: dict[str, Path],
    reference_dir: Path,
    completion_status_path: Path,
    manifest_path: Path | None,
    output_dir: Path,
    iou_threshold: float,
    min_precision: float = 0.8,
    min_recall: float = 0.8,
    min_f1: float = 0.8,
    min_mean_iou: float = 0.5,
) -> dict[str, Path]:
    if not candidate_dirs:
        raise ValueError("At least one candidate method directory is required")
    reference_mask_paths = _mask_lookup(reference_dir)
    completion_status, completion_mask_paths = _read_completion_status(completion_status_path)
    _validate_all_status_rows_complete(completion_status)
    manifest = _read_manifest(manifest_path) if manifest_path is not None else {}
    acceptance_criteria = (
        f"precision>={min_precision:g}; recall>={min_recall:g}; "
        f"f1>={min_f1:g}; mean_iou_matched>={min_mean_iou:g}"
    )

    method_rows: list[dict[str, Any]] = []
    image_rows: list[dict[str, Any]] = []
    for method_name, candidate_dir in sorted(candidate_dirs.items()):
        method_slug = _slug(method_name)
        candidate_mask_paths = _mask_lookup(candidate_dir)
        candidate_output = output_dir / "per_candidate" / method_slug
        summaries = run_manual_mask_validation(
            candidate_mask_paths=candidate_mask_paths,
            reference_mask_paths=reference_mask_paths,
            output_dir=candidate_output,
            iou_threshold=iou_threshold,
            reference_completion_status=completion_status,
            reference_completion_mask_paths=completion_mask_paths,
            run_metadata={
                "candidate_method": method_name,
                "candidate_dir": str(candidate_dir),
                "reference_dir": str(reference_dir),
                "completion_status_path": str(completion_status_path),
                "manifest_path": str(manifest_path) if manifest_path else "",
            },
        )
        matches = _read_csv(candidate_output / "manual_mask_validation_matches.csv")
        overlay_paths = _write_overlays(
            method_slug=method_slug,
            method_name=method_name,
            output_dir=output_dir,
            summaries=summaries,
            matches=matches,
            reference_mask_paths=reference_mask_paths,
            candidate_mask_paths=candidate_mask_paths,
            manifest=manifest,
        )
        method_rows.append(
            _build_method_row(
                method_name=method_name,
                summaries=summaries,
                matches=matches,
                acceptance_criteria=acceptance_criteria,
                min_precision=min_precision,
                min_recall=min_recall,
                min_f1=min_f1,
                min_mean_iou=min_mean_iou,
            )
        )
        for summary in summaries:
            image_rows.append(
                _build_image_row(
                    method_name=method_name,
                    summary=summary,
                    overlay_path=overlay_paths.get(summary["image_id"], ""),
                )
            )

    method_summary_path = output_dir / "method_validation_summary.csv"
    image_summary_path = output_dir / "validation_image_summary.csv"
    report_path = output_dir / "manual_validation_report.md"
    _write_csv(method_summary_path, method_rows, METHOD_SUMMARY_COLUMNS)
    _write_csv(image_summary_path, image_rows, IMAGE_SUMMARY_COLUMNS)
    _write_report(
        report_path,
        method_rows=method_rows,
        image_rows=image_rows,
        iou_threshold=iou_threshold,
        acceptance_criteria=acceptance_criteria,
    )
    return {
        "method_summary": method_summary_path,
        "image_summary": image_summary_path,
        "report": report_path,
        "overlay_dir": output_dir / "overlays",
    }


def parse_candidate_specs(specs: list[str]) -> dict[str, Path]:
    if not specs:
        raise ValueError("At least one candidate must be supplied as name=directory")
    parsed: dict[str, Path] = {}
    for spec in specs:
        if "=" not in spec:
            raise ValueError(
                "Candidate methods must be supplied as name=directory, got "
                f"{spec!r}"
            )
        name, directory = spec.split("=", 1)
        name = name.strip()
        if not name or not directory.strip():
            raise ValueError(
                "Candidate methods must be supplied as name=directory, got "
                f"{spec!r}"
            )
        if name in parsed:
            raise ValueError(f"Duplicate candidate method {name}")
        parsed[name] = Path(directory.strip())
    return parsed


def _build_method_row(
    *,
    method_name: str,
    summaries: list[dict[str, Any]],
    matches: list[dict[str, str]],
    acceptance_criteria: str,
    min_precision: float,
    min_recall: float,
    min_f1: float,
    min_mean_iou: float,
) -> dict[str, Any]:
    n_manual = sum(int(row["n_manual"]) for row in summaries)
    n_predicted = sum(int(row["n_predicted"]) for row in summaries)
    true_positives = sum(int(row["true_positives"]) for row in summaries)
    false_positives = sum(int(row["false_positive_count"]) for row in summaries)
    false_negatives = sum(int(row["false_negative_count"]) for row in summaries)
    precision = _safe_divide(true_positives, n_predicted)
    recall = _safe_divide(true_positives, n_manual)
    f1 = _f1(precision, recall)
    matched_ious = [
        float(row["iou"])
        for row in matches
        if str(row.get("matched", "")).lower() == "true"
    ]
    mean_iou = float(np.mean(matched_ious)) if matched_ious else float("nan")
    count_errors = [abs(int(row["count_error"])) for row in summaries]
    passes = (
        _meets(precision, min_precision)
        and _meets(recall, min_recall)
        and _meets(f1, min_f1)
        and _meets(mean_iou, min_mean_iou)
    )
    return {
        "candidate_method": method_name,
        "n_images": len(summaries),
        "n_manual": n_manual,
        "n_predicted": n_predicted,
        "true_positives": true_positives,
        "false_positives": false_positives,
        "false_negatives": false_negatives,
        "micro_precision": precision,
        "micro_recall": recall,
        "micro_f1": f1,
        "mean_iou_matched": mean_iou,
        "mean_abs_count_error": float(np.mean(count_errors)) if count_errors else float("nan"),
        "max_abs_count_error": max(count_errors) if count_errors else "",
        "passes_acceptance_criteria": passes,
        "acceptance_criteria": acceptance_criteria,
        "validation_scope": "manual/reference instance masks for the configured task",
    }


def _build_image_row(
    *,
    method_name: str,
    summary: dict[str, Any],
    overlay_path: Path | str,
) -> dict[str, Any]:
    return {
        "candidate_method": method_name,
        "image_id": summary["image_id"],
        "iou_threshold": summary["iou_threshold"],
        "n_manual": summary["n_manual"],
        "n_predicted": summary["n_predicted"],
        "true_positives": summary["true_positives"],
        "false_positives": summary["false_positive_count"],
        "false_negatives": summary["false_negative_count"],
        "precision": summary["precision"],
        "recall": summary["recall"],
        "f1": summary["f1"],
        "mean_iou_matched": summary["mean_iou_matched"],
        "count_error": summary["count_error"],
        "count_error_percent": summary["count_error_percent"],
        "overlay_path": str(overlay_path),
    }


def _write_overlays(
    *,
    method_slug: str,
    method_name: str,
    output_dir: Path,
    summaries: list[dict[str, Any]],
    matches: list[dict[str, str]],
    reference_mask_paths: dict[str, Path],
    candidate_mask_paths: dict[str, Path],
    manifest: dict[str, dict[str, str]],
) -> dict[str, Path]:
    overlay_dir = output_dir / "overlays" / method_slug
    overlay_dir.mkdir(parents=True, exist_ok=True)
    overlay_paths: dict[str, Path] = {}
    matches_by_image = _matches_by_image(matches)
    for summary in summaries:
        image_id = str(summary["image_id"])
        manifest_row = manifest.get(image_id)
        if manifest_row is None or not manifest_row.get("ch2_path"):
            raise ValueError(f"{image_id} manifest is missing CH2 path required for overlay")
        candidate_path = candidate_mask_paths.get(image_id)
        if candidate_path is None:
            candidate_path = _write_empty_candidate_mask_for_overlay(
                output_dir=output_dir,
                method_slug=method_slug,
                image_id=image_id,
                reference_mask_path=reference_mask_paths[image_id],
            )
        overlay_path = overlay_dir / f"{image_id}_candidate_vs_reference_overlay.png"
        write_candidate_reference_overlay(
            ch2_path=Path(manifest_row["ch2_path"]),
            candidate_mask_path=candidate_path,
            reference_mask_path=reference_mask_paths[image_id],
            output_path=overlay_path,
            title=f"{image_id} {method_name}",
            metrics=summary,
            match_rows=matches_by_image.get(image_id, []),
        )
        overlay_paths[image_id] = overlay_path
    missing_overlays = [
        str(summary["image_id"])
        for summary in summaries
        if str(summary["image_id"]) not in overlay_paths
    ]
    if missing_overlays:
        raise ValueError(
            "Missing validation overlays for evaluated images: "
            + ", ".join(missing_overlays)
        )
    if overlay_paths:
        _write_contact_sheet(
            image_paths=list(overlay_paths.values()),
            output_path=output_dir / "overlays" / f"{method_slug}_contact_sheet.png",
        )
    return overlay_paths


def _write_empty_candidate_mask_for_overlay(
    *,
    output_dir: Path,
    method_slug: str,
    image_id: str,
    reference_mask_path: Path,
) -> Path:
    reference = np.asarray(tifffile.imread(reference_mask_path))
    empty_candidate = np.zeros(reference.shape, dtype=np.uint32)
    missing_dir = output_dir / "debug" / "missing_candidate_masks" / method_slug
    missing_dir.mkdir(parents=True, exist_ok=True)
    path = missing_dir / f"{image_id}_missing_candidate_as_empty_labels.tif"
    tifffile.imwrite(path, empty_candidate, photometric="minisblack")
    return path


def write_candidate_reference_overlay(
    *,
    ch2_path: Path,
    candidate_mask_path: Path,
    reference_mask_path: Path,
    output_path: Path,
    title: str,
    metrics: dict[str, Any],
    match_rows: list[dict[str, str]] | None = None,
) -> None:
    ch2, _ = read_primary_intensity_plane(ch2_path)
    candidate_labels = np.asarray(tifffile.imread(candidate_mask_path))
    reference_labels = np.asarray(tifffile.imread(reference_mask_path))
    if ch2.shape != candidate_labels.shape or ch2.shape != reference_labels.shape:
        raise ValueError(
            f"Overlay shape mismatch: ch2={ch2.shape}, candidate={candidate_labels.shape}, "
            f"reference={reference_labels.shape}"
        )
    ch2_scaled = _scale_for_display(ch2)
    overlay = np.dstack([ch2_scaled, ch2_scaled, ch2_scaled])
    status_masks = _object_match_status_masks(
        candidate_labels=candidate_labels,
        reference_labels=reference_labels,
        match_rows=match_rows or [],
    )
    overlay_masks = _status_overlay_masks(status_masks)
    _blend(overlay, overlay_masks["false_negative_fill"], np.array([0.0, 1.0, 0.2]), alpha=0.38)
    _blend(overlay, overlay_masks["false_positive_fill"], np.array([1.0, 0.0, 0.9]), alpha=0.38)
    _blend(overlay, overlay_masks["true_positive_fill"], np.array([1.0, 1.0, 0.0]), alpha=0.48)
    _blend(
        overlay,
        overlay_masks["below_threshold_overlap_fill"],
        np.array([1.0, 0.55, 0.0]),
        alpha=0.52,
    )
    overlay[overlay_masks["false_negative_boundary"]] = [0.0, 1.0, 0.0]
    overlay[overlay_masks["false_positive_boundary"]] = [1.0, 0.0, 1.0]
    overlay[overlay_masks["below_threshold_overlap_boundary"]] = [1.0, 0.55, 0.0]
    overlay[overlay_masks["true_positive_boundary"]] = [1.0, 1.0, 0.0]

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(12.2, 4.1),
        constrained_layout=True,
        width_ratios=[1.0, 1.0, 0.64],
    )
    axes[0].imshow(ch2_scaled, cmap="gray")
    axes[0].set_title("CH2/aSMA")
    axes[1].imshow(overlay)
    axes[1].set_title("Object match status\nYellow=TP, magenta=FP, green=FN, orange=FP/FN overlap")
    axes[2].text(
        0.0,
        0.98,
        _metrics_text(metrics),
        va="top",
        ha="left",
        fontsize=8,
        transform=axes[2].transAxes,
    )
    for ax in axes:
        ax.axis("off")
    fig.suptitle(title, fontsize=12)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def _object_match_status_masks(
    *,
    candidate_labels: np.ndarray,
    reference_labels: np.ndarray,
    match_rows: list[dict[str, Any]],
) -> dict[str, np.ndarray]:
    candidate = np.asarray(candidate_labels)
    reference = np.asarray(reference_labels)
    if candidate.shape != reference.shape:
        raise ValueError(
            "candidate and reference labels must have matching shapes for status masks: "
            f"candidate={candidate.shape}, reference={reference.shape}"
        )
    masks = {
        "true_positive": np.zeros(candidate.shape, dtype=bool),
        "false_positive": np.zeros(candidate.shape, dtype=bool),
        "false_negative": np.zeros(candidate.shape, dtype=bool),
    }
    for row in match_rows:
        status = str(row.get("match_status", "")).strip().lower()
        matched = str(row.get("matched", "")).strip().lower() == "true"
        candidate_label = _optional_positive_label(row.get("candidate_label", ""))
        reference_label = _optional_positive_label(row.get("reference_label", ""))
        if matched or status == "true_positive":
            if candidate_label is not None:
                masks["true_positive"] |= candidate == candidate_label
            if reference_label is not None:
                masks["true_positive"] |= reference == reference_label
        elif status == "false_positive" and candidate_label is not None:
            masks["false_positive"] |= candidate == candidate_label
        elif status == "false_negative" and reference_label is not None:
            masks["false_negative"] |= reference == reference_label
    return masks


def _status_overlay_masks(status_masks: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    true_positive = np.asarray(status_masks["true_positive"], dtype=bool)
    false_positive = np.asarray(status_masks["false_positive"], dtype=bool)
    false_negative = np.asarray(status_masks["false_negative"], dtype=bool)
    below_threshold_overlap = false_positive & false_negative
    false_positive_only = false_positive & ~below_threshold_overlap
    false_negative_only = false_negative & ~below_threshold_overlap
    return {
        "true_positive_fill": true_positive,
        "false_positive_fill": false_positive_only,
        "false_negative_fill": false_negative_only,
        "below_threshold_overlap_fill": below_threshold_overlap,
        "true_positive_boundary": segmentation.find_boundaries(true_positive, mode="inner"),
        "false_positive_boundary": segmentation.find_boundaries(false_positive_only, mode="inner"),
        "false_negative_boundary": segmentation.find_boundaries(false_negative_only, mode="inner"),
        "below_threshold_overlap_boundary": segmentation.find_boundaries(
            below_threshold_overlap,
            mode="inner",
        ),
    }


def _matches_by_image(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    by_image: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        image_id = _normalize_image_id(row.get("image_id", ""))
        if not image_id:
            continue
        by_image.setdefault(image_id, []).append(row)
    return by_image


def _optional_positive_label(value: Any) -> int | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    label = int(text)
    if label <= 0:
        return None
    return label


def _metrics_text(metrics: dict[str, Any]) -> str:
    return "\n".join(
        [
            f"n_manual: {metrics['n_manual']}",
            f"n_predicted: {metrics['n_predicted']}",
            f"TP: {metrics['true_positives']}",
            f"FP: {metrics['false_positive_count']}",
            f"FN: {metrics['false_negative_count']}",
            f"precision: {_format_float(metrics['precision'])}",
            f"recall: {_format_float(metrics['recall'])}",
            f"F1: {_format_float(metrics['f1'])}",
            f"mean IoU: {_format_float(metrics['mean_iou_matched'])}",
        ]
    )


def _write_report(
    path: Path,
    *,
    method_rows: list[dict[str, Any]],
    image_rows: list[dict[str, Any]],
    iou_threshold: float,
    acceptance_criteria: str,
) -> None:
    lines = [
        "# Manual Mask Validation Report",
        "",
        "Status: Validation run completed against manual/reference masks for the configured task.",
        "",
        "This validates only the supplied manual-reference task and mask set. It does not prove biological efficacy, does not validate whole-cell segmentation unless the manual task was whole-cell segmentation, and does not replace visual QC.",
        "",
        f"IoU threshold: `{iou_threshold:g}`",
        f"Acceptance criteria: `{acceptance_criteria}`",
        "",
        "## Method Summary",
        "",
        "| Method | Images | Precision | Recall | F1 | Mean IoU | Pass |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in method_rows:
        lines.append(
            "| {candidate_method} | {n_images} | {micro_precision} | {micro_recall} | "
            "{micro_f1} | {mean_iou_matched} | {passes_acceptance_criteria} |".format(
                **_format_report_row(row)
            )
        )
    lines.extend(
        [
            "",
            "## Image-Level Summary",
            "",
            "| Method | Image | Manual | Predicted | TP | FP | FN | Precision | Recall | F1 | Mean IoU |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in image_rows:
        lines.append(
            "| {candidate_method} | {image_id} | {n_manual} | {n_predicted} | "
            "{true_positives} | {false_positives} | {false_negatives} | "
            "{precision} | {recall} | {f1} | {mean_iou_matched} |".format(
                **_format_report_row(row)
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_completion_status(path: Path) -> tuple[dict[str, str], dict[str, Path]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"image_id", "manual_reference_mask_path", "status"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
        statuses: dict[str, str] = {}
        mask_paths: dict[str, Path] = {}
        for row in reader:
            image_id = _normalize_image_id(row["image_id"])
            if not image_id:
                continue
            if image_id in statuses:
                raise ValueError(f"Duplicate completion status for {image_id} in {path}")
            statuses[image_id] = row["status"].strip().lower()
            mask_paths[image_id] = _resolve_status_path(
                row["manual_reference_mask_path"],
                status_path=path,
            )
    if not statuses:
        raise ValueError(f"No completion status rows found in {path}")
    return statuses, mask_paths


def _validate_all_status_rows_complete(statuses: dict[str, str]) -> None:
    allowed = {"complete_non_empty", "confirmed_empty"}
    incomplete = [
        f"{image_id}={status or 'blank'}"
        for image_id, status in sorted(statuses.items())
        if status not in allowed
    ]
    if incomplete:
        raise ValueError("reference masks are not marked complete: " + ", ".join(incomplete))


def _read_manifest(path: Path) -> dict[str, dict[str, str]]:
    rows = _read_csv(path)
    return {_normalize_image_id(row["image_id"]): row for row in rows if row.get("image_id")}


def _mask_lookup(root: Path) -> dict[str, Path]:
    if not root.exists():
        raise FileNotFoundError(f"Mask directory does not exist: {root}")
    lookup: dict[str, Path] = {}
    for path in sorted(root.rglob("*.tif")) + sorted(root.rglob("*.tiff")):
        match = re.search(r"(XY\d+)", path.name, flags=re.IGNORECASE)
        if match is None:
            continue
        image_id = match.group(1).upper()
        if image_id in lookup:
            raise ValueError(f"Multiple masks found for {image_id} under {root}")
        lookup[image_id] = path
    if not lookup:
        raise ValueError(f"No XY## TIFF masks found under {root}")
    return lookup


def _write_contact_sheet(*, image_paths: list[Path], output_path: Path) -> None:
    thumbs: list[Image.Image] = []
    for image_path in image_paths:
        image = Image.open(image_path).convert("RGB")
        image.thumbnail((520, 220))
        thumbs.append(image.copy())
    if not thumbs:
        return
    width = max(image.width for image in thumbs)
    height = sum(image.height for image in thumbs)
    sheet = Image.new("RGB", (width, height), "white")
    y = 0
    for image in thumbs:
        sheet.paste(image, (0, y))
        y += image.height
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(output_path)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _resolve_status_path(value: str, *, status_path: Path) -> Path:
    path = Path(value.strip()).expanduser()
    if path.is_absolute():
        return path
    if path.exists():
        return path
    return status_path.parent / path


def _normalize_image_id(value: str) -> str:
    return value.strip().upper().replace(" ", "")


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("_")
    return slug or "candidate"


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


def _meets(value: float, threshold: float) -> bool:
    return bool(np.isfinite(value) and value >= threshold)


def _scale_for_display(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image, dtype=np.float32)
    low, high = np.percentile(arr, [1, 99.8])
    if high <= low:
        high = float(np.max(arr))
        low = float(np.min(arr))
    if high <= low:
        return np.zeros_like(arr, dtype=np.float32)
    return np.clip((arr - low) / (high - low), 0, 1)


def _blend(rgb: np.ndarray, mask: np.ndarray, color: np.ndarray, *, alpha: float) -> None:
    rgb[mask] = (1 - alpha) * rgb[mask] + alpha * color


def _format_float(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(number):
        return "NA"
    return f"{number:.3f}"


def _format_report_row(row: dict[str, Any]) -> dict[str, Any]:
    formatted = dict(row)
    for key, value in row.items():
        if isinstance(value, float):
            formatted[key] = _format_float(value)
    return formatted
