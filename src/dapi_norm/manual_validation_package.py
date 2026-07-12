from __future__ import annotations

import csv
from pathlib import Path
import textwrap
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from skimage import segmentation
import tifffile
import yaml

from dapi_norm.image_arrays import read_primary_intensity_plane
from dapi_norm.pi_simple_summary import ImagePair, find_image_pairs
from dapi_norm.seeded_regions import _draw_centroid_crosses, _retained_region_overlay_rgb


MANIFEST_COLUMNS = [
    "image_id",
    "source_id",
    "validation_task",
    "ch2_path",
    "ch4_path",
    "candidate_mask_path",
    "nuclei_mask_path",
    "manual_reference_mask_path",
    "annotation_panel_path",
    "guide_panel_path",
    "method",
    "foreground_method",
    "dapi_positive_nucleus_count",
    "candidate_integrated_raw",
    "candidate_intensity_per_DAPI_positive_nucleus",
    "qc_status",
    "qc_flags",
]

STATUS_COLUMNS = [
    "image_id",
    "manual_reference_mask_path",
    "annotation_panel_path",
    "status",
    "labeler",
    "completed_date",
    "notes",
]


def prepare_manual_validation_package(
    *,
    input_root: Path,
    seeded_run_dir: Path,
    output_dir: Path,
    positions: list[str] | None = None,
    iou_threshold: float = 0.5,
    task: str = "asma_associated_region",
    force_overwrite_reference_masks: bool = False,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    reference_mask_dir = output_dir / "reference_masks_to_fill"
    annotation_dir = output_dir / "annotation_panels_raw_only"
    guide_dir = output_dir / "guide_panels"
    reference_mask_dir.mkdir(parents=True, exist_ok=True)
    annotation_dir.mkdir(parents=True, exist_ok=True)
    guide_dir.mkdir(parents=True, exist_ok=True)

    config_path = seeded_run_dir / "logs" / "config_resolved.yaml"
    config = _load_yaml(config_path)
    pair_lookup = _pair_lookup(input_root)
    metric_rows = _read_metrics(seeded_run_dir / "summaries" / "seeded_region_image_metrics.csv")
    selected_positions = _select_positions(metric_rows, positions)
    manifest_rows: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []
    for image_id in selected_positions:
        pair = pair_lookup[image_id]
        metrics = metric_rows[image_id]
        candidate_mask_path = _resolve_recorded_path(
            Path(metrics["mask_path"]),
            seeded_run_dir=seeded_run_dir,
            config=config,
            config_path=config_path,
        )
        nuclei_mask_path = _nuclei_mask_path_for_image(
            image_id=image_id,
            seeded_run_dir=seeded_run_dir,
            config=config,
            config_path=config_path,
        )
        ch2_image, _ = read_primary_intensity_plane(pair.ch2_path)
        ch4_image, _ = read_primary_intensity_plane(pair.ch4_path)
        candidate_mask = np.asarray(tifffile.imread(candidate_mask_path), dtype=np.uint32)
        nuclei_mask = np.asarray(tifffile.imread(nuclei_mask_path), dtype=np.uint32)
        if candidate_mask.shape != ch2_image.shape:
            raise ValueError(
                f"{image_id} candidate mask shape {candidate_mask.shape} does not match "
                f"CH2 image shape {ch2_image.shape}"
        )
        manual_mask_path = reference_mask_dir / f"{image_id}_manual_reference_labels.tif"
        annotation_path = annotation_dir / f"{image_id}_manual_annotation_panel.png"
        guide_path = guide_dir / f"{image_id}_manual_validation_guide.png"
        _write_blank_reference_mask(
            manual_mask_path,
            shape=candidate_mask.shape,
            force=force_overwrite_reference_masks,
        )
        write_manual_annotation_panel(
            image_id=image_id,
            ch2_image=ch2_image,
            ch4_image=ch4_image,
            nuclei_mask=nuclei_mask,
            output_path=annotation_path,
            task=task,
        )
        write_manual_validation_guide_panel(
            image_id=image_id,
            ch2_image=ch2_image,
            ch4_image=ch4_image,
            candidate_mask=candidate_mask,
            nuclei_mask=nuclei_mask,
            output_path=guide_path,
            metrics=metrics,
            task=task,
        )
        manifest_rows.append(
            {
                "image_id": image_id,
                "source_id": metrics.get("source_id", image_id),
                "validation_task": task,
                "ch2_path": str(pair.ch2_path),
                "ch4_path": str(pair.ch4_path),
                "candidate_mask_path": str(candidate_mask_path),
                "nuclei_mask_path": str(nuclei_mask_path),
                "manual_reference_mask_path": str(manual_mask_path),
                "annotation_panel_path": str(annotation_path),
                "guide_panel_path": str(guide_path),
                "method": metrics.get("method", ""),
                "foreground_method": metrics.get("foreground_method", ""),
                "dapi_positive_nucleus_count": metrics.get("dapi_positive_nucleus_count", ""),
                "candidate_integrated_raw": metrics.get("seeded_region_integrated_raw", ""),
                "candidate_intensity_per_DAPI_positive_nucleus": metrics.get(
                    "seeded_region_intensity_per_DAPI_positive_nucleus", ""
                ),
                "qc_status": metrics.get("qc_status", ""),
                "qc_flags": metrics.get("qc_flags", ""),
            }
        )
        status_rows.append(
            {
                "image_id": image_id,
                "manual_reference_mask_path": str(manual_mask_path.resolve()),
                "annotation_panel_path": str(annotation_path.resolve()),
                "status": "not_started",
                "labeler": "",
                "completed_date": "",
                "notes": "",
            }
        )

    manifest_path = output_dir / "manual_validation_manifest.csv"
    status_path = output_dir / "manual_labeling_status.csv"
    readme_path = output_dir / "README.md"
    _write_csv(manifest_path, manifest_rows, MANIFEST_COLUMNS)
    _write_status_csv(status_path, status_rows, force=force_overwrite_reference_masks)
    _write_readme(
        readme_path,
        output_dir=output_dir,
        seeded_run_dir=seeded_run_dir,
        reference_mask_dir=reference_mask_dir,
        status_path=status_path,
        iou_threshold=iou_threshold,
        task=task,
        image_ids=selected_positions,
    )
    return {
        "manifest": manifest_path,
        "status": status_path,
        "readme": readme_path,
        "reference_mask_dir": reference_mask_dir,
        "annotation_dir": annotation_dir,
        "guide_dir": guide_dir,
    }


def write_manual_annotation_panel(
    *,
    image_id: str,
    ch2_image: np.ndarray,
    ch4_image: np.ndarray,
    nuclei_mask: np.ndarray,
    output_path: Path,
    task: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ch2_scaled = _scale_for_display(ch2_image)
    ch4_scaled = _scale_for_display(ch4_image)
    nuclei = np.asarray(nuclei_mask)
    ch4_rgb = np.dstack([ch4_scaled, ch4_scaled, ch4_scaled])
    nucleus_boundaries = segmentation.find_boundaries(nuclei, mode="outer")
    ch4_rgb[nucleus_boundaries] = [0.0, 1.0, 0.2]
    _draw_centroid_crosses(ch4_rgb, nuclei, color=np.array([0.0, 1.0, 0.2]))

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(12.0, 4.3),
        constrained_layout=True,
        width_ratios=[1.0, 1.0, 0.8],
    )
    axes[0].imshow(ch2_scaled, cmap="gray")
    axes[0].set_title(f"{image_id} CH2/aSMA")
    axes[1].imshow(ch4_rgb)
    axes[1].set_title("CH4/DAPI nuclei\nGreen X = detected nucleus")
    axes[2].text(
        0.0,
        0.98,
        _annotation_protocol_text(task),
        ha="left",
        va="top",
        fontsize=8,
        transform=axes[2].transAxes,
        color="#212529",
    )
    for ax in axes:
        ax.axis("off")
    fig.suptitle("Raw-only annotation panel: create manual reference labels from these channels", fontsize=12)
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def write_manual_validation_guide_panel(
    *,
    image_id: str,
    ch2_image: np.ndarray,
    ch4_image: np.ndarray,
    candidate_mask: np.ndarray,
    nuclei_mask: np.ndarray,
    output_path: Path,
    metrics: dict[str, Any],
    task: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ch2_scaled = _scale_for_display(ch2_image)
    ch4_scaled = _scale_for_display(ch4_image)
    candidate = np.asarray(candidate_mask) > 0
    nuclei = np.asarray(nuclei_mask)

    ch4_rgb = np.dstack([ch4_scaled, ch4_scaled, ch4_scaled])
    nucleus_boundaries = segmentation.find_boundaries(nuclei, mode="outer")
    ch4_rgb[nucleus_boundaries] = [0.0, 1.0, 0.2]
    _draw_centroid_crosses(ch4_rgb, nuclei, color=np.array([0.0, 1.0, 0.2]))

    overlay = _retained_region_overlay_rgb(ch2_scaled, candidate)
    candidate_boundaries = segmentation.find_boundaries(candidate_mask, mode="outer")
    overlay[candidate_boundaries] = [0.0, 1.0, 0.2]
    overlay[nucleus_boundaries] = [0.0, 0.75, 1.0]

    fig, axes = plt.subplots(
        1,
        4,
        figsize=(15.5, 4.3),
        constrained_layout=True,
        width_ratios=[1.0, 1.0, 1.0, 0.72],
    )
    axes[0].imshow(ch2_scaled, cmap="gray")
    axes[0].set_title(f"{image_id} CH2/aSMA")
    axes[1].imshow(ch4_rgb)
    axes[1].set_title("CH4/DAPI nuclei\nGreen X = detected nucleus")
    axes[2].imshow(overlay)
    axes[2].set_title("Candidate retained region\nLight green over CH2")
    axes[3].text(
        0.0,
        0.98,
        _guide_text(metrics=metrics, task=task),
        ha="left",
        va="top",
        fontsize=8,
        transform=axes[3].transAxes,
        color="#212529",
    )
    for ax in axes:
        ax.axis("off")
    fig.suptitle("Candidate guide panel: post-label QC comparison only", fontsize=12)
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def _write_blank_reference_mask(path: Path, *, shape: tuple[int, int], force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(
            "manual reference mask already exists; refusing to overwrite possible labels: "
            f"{path}. Use --force-overwrite-reference-masks only when intentionally "
            "recreating blank placeholders."
        )
    tifffile.imwrite(
        path,
        np.zeros(shape, dtype=np.uint32),
        photometric="minisblack",
    )


def _guide_text(*, metrics: dict[str, Any], task: str) -> str:
    lines = [
        f"Task: {task}",
        f"Method: {metrics.get('method', 'NA')}",
        f"Foreground: {metrics.get('foreground_method', 'NA')}",
        f"Nuclei: {metrics.get('dapi_positive_nucleus_count', 'NA')}",
        "Intensity/nucleus: "
        f"{_format_scientific(metrics.get('seeded_region_intensity_per_DAPI_positive_nucleus'))}",
        f"QC: {metrics.get('qc_status', 'NA')}",
        "Flags:",
        textwrap.fill(str(metrics.get("qc_flags", "")), width=34),
        "",
        "Post-label QC only.",
        "Do not trace this candidate overlay.",
        "Use raw-only annotation panels for labels.",
        "",
        "Use separate integer labels.",
        "0 = background.",
        "Do not make one binary mask for many regions.",
    ]
    return "\n".join(lines)


def _write_readme(
    path: Path,
    *,
    output_dir: Path,
    seeded_run_dir: Path,
    reference_mask_dir: Path,
    status_path: Path,
    iou_threshold: float,
    task: str,
    image_ids: list[str],
) -> None:
    validation_output = output_dir / "validation_results"
    lines = [
        "# Manual Validation Package",
        "",
        f"Validation task: `{task}`.",
        "",
        "Fill the TIFFs in `reference_masks_to_fill/` with manual/reference instance labels.",
        "",
        "Rules:",
        "",
        "- `0` must remain background.",
        "- Each manual object/region must have its own positive integer label.",
        "- Do not save a single binary foreground mask for multiple objects.",
        "- Primary instance unit for `asma_associated_region`: one contiguous aSMA-positive cellular region associated with one or more nearby DAPI-positive nuclei.",
        "- Split touching regions when there is a clear intensity valley or cell boundary; keep inseparable merged regions as one label and record that ambiguity outside the mask if needed.",
        "- Do not label CH2-negative DAPI cells as aSMA-associated regions.",
        "- Include edge objects only if most of the aSMA-associated region is visible in the field.",
        "- Annotate rejected/challenge fields too when visible aSMA-associated regions exist; leave them empty only when the manual reference truly has no object.",
        "- Use `annotation_panels_raw_only/` for drawing reference labels.",
        "- Do not trace the candidate overlay; use candidate-overlay guide panels only for post-label QC comparison.",
        "- Update `manual_labeling_status.csv` before validation: use `complete_non_empty` for filled masks and `confirmed_empty` for fields intentionally left empty.",
        "- Do not copy `manual_labeling_status.csv` across packages; each row's `manual_reference_mask_path` must match the TIFF being validated.",
        "- For whole-cell validation, trace whole-cell boundaries only if the images visibly support that task.",
        "",
        "Do not report precision, recall, F1, or IoU until manual masks are filled and the validator has been run.",
        "",
        "Images included:",
        "",
        *[f"- `{image_id}`" for image_id in image_ids],
        "",
        "Run validation after manual labeling:",
        "",
        "```bash",
        ".venv/bin/python scripts/validate_manual_instance_masks.py \\",
        f"  --candidate-dir {seeded_run_dir / 'masks'} \\",
        f"  --reference-dir {reference_mask_dir} \\",
        f"  --completion-status {status_path} \\",
        f"  --output {validation_output} \\",
        f"  --iou-threshold {iou_threshold:g}",
        "```",
        "",
        "The `guide_panels/` folder shows candidate overlays for post-label QC. It is not the primary source for manual tracing.",
        "",
        "This package supports quantitative validation; it does not itself prove the method is validated.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _annotation_protocol_text(task: str) -> str:
    if task == "asma_associated_region":
        lines = [
            "Manual annotation protocol",
            "",
            "Primary instance unit:",
            "one contiguous aSMA-positive cellular region associated with nearby DAPI-positive nucleus/nuclei.",
            "",
            "Do not trace the candidate overlay.",
            "Use this raw-only panel for reference labels.",
            "",
            "0 = background.",
            "Each object gets a separate positive integer label.",
            "Leave CH2-negative DAPI cells unlabeled for this task.",
        ]
    else:
        lines = [
            "Manual annotation protocol",
            "",
            f"Task: {task}",
            "Define the instance unit before labeling.",
            "Do not trace the candidate overlay.",
            "0 = background.",
            "Each object gets a separate positive integer label.",
        ]
    return "\n".join(textwrap.fill(line, width=34) if len(line) > 34 else line for line in lines)


def _pair_lookup(input_root: Path) -> dict[str, ImagePair]:
    lookup: dict[str, ImagePair] = {}
    for pair in find_image_pairs(input_root):
        for key in {pair.location.upper(), pair.source_id.upper()}:
            if key in lookup:
                raise ValueError(f"Duplicate image pair key {key} under {input_root}")
            lookup[key] = pair
    if not lookup:
        raise ValueError(f"No CH2/CH4 image pairs found under {input_root}")
    return lookup


def _read_metrics(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"image_id", "mask_path", "qc_status", "qc_flags"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
        rows: dict[str, dict[str, str]] = {}
        for row in reader:
            image_id = row["image_id"].strip().upper()
            rows[image_id] = row
        return rows


def _select_positions(metrics: dict[str, dict[str, str]], positions: list[str] | None) -> list[str]:
    if positions is None:
        return sorted(metrics, key=_image_sort_key)
    selected = [position.strip().upper().replace(" ", "") for position in positions]
    missing = [position for position in selected if position not in metrics]
    if missing:
        raise KeyError(f"Requested positions missing from seeded metrics: {', '.join(missing)}")
    return selected


def _nuclei_mask_path_for_image(
    *,
    image_id: str,
    seeded_run_dir: Path,
    config: dict,
    config_path: Path,
) -> Path:
    for record in config.get("image_records", []):
        keys = {str(record.get("source_id", "")).upper(), str(record.get("location", "")).upper()}
        if image_id.upper() not in keys:
            continue
        return _resolve_recorded_path(
            Path(record["nuclei_mask_path"]),
            seeded_run_dir=seeded_run_dir,
            config=config,
            config_path=config_path,
        )
    raise KeyError(f"No image record for {image_id} in {config_path}")


def _resolve_recorded_path(
    path: Path,
    *,
    seeded_run_dir: Path,
    config: dict,
    config_path: Path,
) -> Path:
    if path.is_absolute():
        return path
    for base in _relative_path_bases(
        seeded_run_dir=seeded_run_dir,
        config=config,
        config_path=config_path,
    ):
        candidate = base / path
        if candidate.exists():
            return candidate
    return Path.cwd() / path


def _relative_path_bases(*, seeded_run_dir: Path, config: dict, config_path: Path) -> list[Path]:
    bases: list[Path] = []
    project_root = _infer_project_root_from_output_dir(seeded_run_dir, config)
    if project_root is not None:
        bases.append(project_root)
    bases.extend([seeded_run_dir, config_path.parent, config_path.parent.parent, Path.cwd()])
    deduped: list[Path] = []
    seen: set[Path] = set()
    for base in bases:
        resolved = base.resolve()
        if resolved not in seen:
            deduped.append(resolved)
            seen.add(resolved)
    return deduped


def _infer_project_root_from_output_dir(seeded_run_dir: Path, config: dict) -> Path | None:
    recorded_output_dir = Path(str(config.get("output_dir", "")))
    if recorded_output_dir.is_absolute() or not recorded_output_dir.parts:
        return None
    run_dir = seeded_run_dir.resolve()
    recorded_parts = recorded_output_dir.parts
    if len(recorded_parts) > len(run_dir.parts):
        return None
    if run_dir.parts[-len(recorded_parts) :] != recorded_parts:
        return None
    return run_dir.parents[len(recorded_parts) - 1]


def _load_yaml(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Required config file missing: {path}")
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _write_status_csv(path: Path, rows: list[dict[str, Any]], *, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(
            "manual labeling status file already exists; refusing to overwrite possible "
            f"completion annotations: {path}. Use --force-overwrite-reference-masks only "
            "when intentionally recreating blank placeholders."
        )
    _write_csv(path, rows, STATUS_COLUMNS)


def _scale_for_display(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image, dtype=np.float32)
    low, high = np.percentile(arr, [1, 99.8])
    if high <= low:
        high = float(np.max(arr))
        low = float(np.min(arr))
    if high <= low:
        return np.zeros_like(arr, dtype=np.float32)
    return np.clip((arr - low) / (high - low), 0, 1)


def _format_scientific(value: Any) -> str:
    try:
        return f"{float(value):.3e}"
    except (TypeError, ValueError):
        return "NA"


def _image_sort_key(image_id: str) -> tuple[int, str]:
    digits = "".join(char for char in image_id if char.isdigit())
    return (int(digits) if digits else 10**9, image_id)
